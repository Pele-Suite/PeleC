#include <limits>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Reduce.H>
#include <AMReX_Loop.H>
#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#endif

#include "PeleC.H"
#include "prob.H"
#include "NSCBC.H"

namespace {
// Device-resident fallback counters for the characteristic boundary
// treatment.  Created once, host-side, by PeleC::nscbc_init_diagnostics()
// from read_params (before any boundary fill can run); reduced and reported
// by PeleC::nscbc_report_diagnostics().  Not created lazily from the fill:
// the fill runs inside the OpenMP-parallel MFIter loop of
// StateDataPhysBCFunct, where a first-use allocation would race.
//
// Held as a heap pointer released through amrex::ExecOnFinalize rather than as
// a function-local static object.  A static Gpu::DeviceVector destructs at
// program exit, which is AFTER amrex::Finalize() has torn down the arena it
// allocated from -- harmless on a CPU build and a use-after-free of the device
// allocator on a GPU one.
amrex::Gpu::DeviceVector<amrex::Long>* nscbc_diag_p = nullptr;

// Counter layout: the pc_nscbc::Diag fallbacks, then one count of
// characteristic fills per domain face (2*dir + side).
constexpr int nscbc_ncounters = pc_nscbc::Diag::count + 2 * AMREX_SPACEDIM;

// Host bookkeeping for the report (host-only, touched from the report alone).
amrex::Long nscbc_face_fills_total[2 * AMREX_SPACEDIM] = {0};
bool nscbc_face_announced[2 * AMREX_SPACEDIM] = {false};
bool nscbc_none_warned = false;
int nscbc_last_report_step = 0;

// EB-covered masks of the boundary bands, one set per AMR level, keyed by
// the level's domain box (each level has its own).  Geometry-only, so built
// once per level (PeleC::nscbc_init_eb_mask) and kept for the run; released
// through ExecOnFinalize like the counters.
struct NSCBCCoveredMask
{
  amrex::Box domain;
  amrex::Box band[2 * AMREX_SPACEDIM];
  amrex::Gpu::DeviceVector<char> mask[2 * AMREX_SPACEDIM];
};
amrex::Vector<std::unique_ptr<NSCBCCoveredMask>>* nscbc_masks_p = nullptr;

NSCBCCoveredMask const*
nscbc_find_mask(const amrex::Box& domain)
{
  if (nscbc_masks_p == nullptr) {
    return nullptr;
  }
  for (const auto& m : *nscbc_masks_p) {
    if (m->domain == domain) {
      return m.get();
    }
  }
  return nullptr;
}

amrex::Gpu::DeviceVector<amrex::Long>&
nscbc_diag()
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    nscbc_diag_p != nullptr,
    "NSCBC counters used before PeleC::nscbc_init_diagnostics()");
  return *nscbc_diag_p;
}

} // namespace

struct PCHypFillExtDir
{
  ProbParmDevice const* lprobparm;
  bool m_do_turb_inflow{false};
  bool m_nscbc{false};
  amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM> m_nscbc_prm;
  amrex::Long* m_nscbc_diag{nullptr};
  // TurbInflow fluctuations for the characteristic fill: a copy of the
  // layer-1 ghost velocities add_turb wrote, taken before the launch (the
  // fill overwrites those cells while other threads still need them).  One
  // band per domain face, indexed 2*dir + side (0 = lo, 1 = hi); empty on a
  // face that is not ext_dir or that this FAB does not touch.
  amrex::GpuArray<amrex::Array4<const amrex::Real>, 2 * AMREX_SPACEDIM>
    m_turb{};
  bool m_has_turb{false};
  // EB-covered mask of the depth-3 boundary band of each domain face at this
  // level (nonzero = covered), indexed 2*dir + side; empty without EB or on a
  // face that is not ext_dir.  See PeleC::nscbc_init_eb_mask().
  amrex::GpuArray<amrex::Array4<const char>, 2 * AMREX_SPACEDIM> m_covered{};

  AMREX_GPU_HOST
  explicit PCHypFillExtDir(
    const ProbParmDevice* d_prob_parm,
    const bool do_turb_inflow,
    const bool nscbc,
    const amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM>& nscbc_prm,
    amrex::Long* nscbc_diag,
    amrex::GpuArray<amrex::Array4<const amrex::Real>, 2 * AMREX_SPACEDIM> const&
      turb = {},
    const bool has_turb = false,
    amrex::GpuArray<amrex::Array4<const char>, 2 * AMREX_SPACEDIM> const&
      covered = {})
    : lprobparm(d_prob_parm),
      m_do_turb_inflow(do_turb_inflow),
      m_nscbc(nscbc),
      m_nscbc_prm(nscbc_prm),
      m_nscbc_diag(nscbc_diag),
      m_turb(turb),
      m_has_turb(has_turb),
      m_covered(covered)
  {
  }

  // -------------------------------------------------------------------------
  //  Characteristic (NSCBC) fill for one ghost cell.
  //
  //  Returns true if this ghost cell was filled here, in which case the
  //  ordinary bcnormal() path below is skipped for it entirely.
  //
  //  Corner ownership.  A ghost cell may lie outside the domain in more than
  //  one direction.  Such a cell is owned by the LOWEST idir in which it is
  //  outside on an ext_dir face whose problem hook returns a live target
  //  and whose boundary cell holds fluid (not EB-covered); the state is
  //  written exactly once.  With the tangential stencil row chosen below,
  //  the fill is a pure function of data present before the launch, hence
  //  independent of the order in which ghost cells are visited and identical
  //  on CPU and GPU.
  //
  //  (The legacy bcnormal() path further down does NOT have this property:
  //  at a corner it reads dest() at a tangential index that is itself a
  //  ghost cell, filled by the face pass that AMReX runs before its edge and
  //  corner passes, so its corner value depends on that pass order rather
  //  than on interior data alone.  There is no race; that pre-existing
  //  behaviour is left untouched so that no existing result moves.)
  // -------------------------------------------------------------------------
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  bool nscbc_fill(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr) const
  {
    const int* bc = bcr[URHO].data();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* dx = geom.CellSize();

    const amrex::Dim3 lo3 = amrex::lbound(dest);
    const amrex::Dim3 hi3 = amrex::ubound(dest);
    const int fab_lo[3] = {lo3.x, lo3.y, lo3.z};
    const int fab_hi[3] = {hi3.x, hi3.y, hi3.z};

    // Index of a cell's periodic image inside the domain (identity in a
    // non-periodic direction).
    auto wrap = [&](const int d, const int i) {
      if (geom.isPeriodic(d) == 0) {
        return i;
      }
      const int n_d = domhi[d] - domlo[d] + 1;
      return domlo[d] + (((i - domlo[d]) % n_d) + n_d) % n_d;
    };

    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
      int sgn = 0;
      if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
        sgn = +1;
      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
        (iv[idir] > domhi[idir])) {
        sgn = -1;
      } else {
        continue;
      }
      const int face = 2 * idir + ((sgn > 0) ? 0 : 1);

      const int N_pos = (sgn > 0) ? domlo[idir] : domhi[idir];
      const int layer = sgn * (N_pos - iv[idir]); // 1 = nearest the domain

      // The tangential row the stencil walks along.  Per tangential
      // direction d:
      //
      //  * PERIODIC: the row is read in place.  FillPatch has already filled
      //    this FAB's cells beyond the seam with their periodic images
      //    (ParallelCopy / FillBoundary run before the physical BC fill), so
      //    in place IS the image, and the result agrees bit for bit with
      //    the one AMReX's corner fix-up computes for the same cell on its
      //    shifted strip FAB.  (A clamp to the boundary row here left the
      //    corners outside in two periodic directions -- which the strip
      //    fix-up never revisits -- built from the wrong row.)
      //
      //  * WALL or SYMMETRY (reflect_even scalars): the row is the mirror
      //    image of iv across that face, and the reflect_odd components are
      //    negated after the fill, so a corner ghost is the exact wall image
      //    of the characteristic ghost at the mirrored row.
      //
      //  * anything else (FOExtrap, another Hard face, which owns the corner
      //    only when it has the lower idir): clamped to the boundary row.
      //
      // Clamps are into the domain AND the FAB, so only valid cells are read.
      // A FAB lying wholly beyond the domain in a non-periodic tangential
      // direction has no valid row at all; it is left to the legacy path.
      amrex::IntVect base(AMREX_D_DECL(iv[0], iv[1], iv[2]));
      int mirror[AMREX_SPACEDIM] = {0}; // 0 none, 1 lo face, 2 hi face
      bool usable = true;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if ((d == idir) || (geom.isPeriodic(d) != 0)) {
          continue;
        }
        const int lo_in = amrex::max<int>(domlo[d], fab_lo[d]);
        const int hi_in = amrex::min<int>(domhi[d], fab_hi[d]);
        if (lo_in > hi_in) {
          usable = false;
          break;
        }
        int b = iv[d];
        if ((iv[d] < domlo[d]) && (bc[d] == amrex::BCType::reflect_even)) {
          b = 2 * domlo[d] - 1 - iv[d];
          mirror[d] = 1;
        } else if (
          (iv[d] > domhi[d]) &&
          (bc[d + AMREX_SPACEDIM] == amrex::BCType::reflect_even)) {
          b = 2 * domhi[d] + 1 - iv[d];
          mirror[d] = 2;
        }
        base[d] = amrex::min<int>(amrex::max<int>(b, lo_in), hi_in);
      }
      if (!usable) {
        continue;
      }

      // How deep a normal stencil is available, in the domain and in the FAB.
      const int depth_domain = domhi[idir] - domlo[idir] + 1;
      const int depth_fab =
        (sgn > 0) ? (fab_hi[idir] - N_pos + 1) : (N_pos - fab_lo[idir] + 1);
      int n_stencil =
        amrex::min<int>(3, amrex::min<int>(depth_domain, depth_fab));
      if (n_stencil < 1) {
        continue;
      }

      auto stencil_iv = [&](const int step) {
        amrex::IntVect r = base;
        r[idir] = N_pos + sgn * amrex::min<int>(step, n_stencil - 1);
        return r;
      };

      // EB: covered stencil cells are identified from the EB geometry (a
      // per-level mask of the boundary band, built at level construction),
      // not from the data -- PeleC's default body state is a sampled fluid
      // state that no test on the data can tell from interior data.
      auto covered = [&](const amrex::IntVect& c) {
        const auto& m = m_covered[face];
        if (!m) {
          return false;
        }
        amrex::IntVect w = c;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          w[d] = wrap(d, c[d]);
        }
        const amrex::Dim3 w3 = w.dim3();
        return m.contains(w3.x, w3.y, w3.z) && (m(w) != 0);
      };

      // A covered boundary cell is wall, not a characteristic boundary
      // point: the legacy path fills it, as it would without bc_nscbc.
      const amrex::IntVect ivN = stencil_iv(0);
      if (covered(ivN)) {
        continue;
      }
      amrex::Real s_N[NVAR], s_Nm1[NVAR], s_Nm2[NVAR];
      for (int n = 0; n < NVAR; n++) {
        s_N[n] = dest(ivN, n);
      }

      // Query the problem for this boundary POINT.  x is the location on the
      // boundary plane, not the ghost cell centre: the target is a property
      // of the boundary point and must not vary with the ghost layer, or the
      // relaxation would be applied to a different target in each layer.
      // Across a periodic direction it is the image's location.
      amrex::Real x[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        x[d] = prob_lo[d] +
               (static_cast<amrex::Real>(wrap(d, base[d])) + 0.5) * dx[d];
      }
      x[idir] = prob_lo[idir] + static_cast<amrex::Real>(
                                  (sgn > 0) ? domlo[idir] : (domhi[idir] + 1)) *
                                  dx[idir];

      // TurbInflow fluctuations at this boundary point, as bcnormal gets
      // them; the hook decides how they enter its Target (added to the mean,
      // or as the total velocity).
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> turb_fluc{0.0};
      const auto& turb = m_turb[face];
      if (m_has_turb && turb) {
        amrex::IntVect ivt = base;
        ivt[idir] = N_pos - sgn; // the layer-1 ghost, where add_turb wrote
        // add_turb wrote only within the domain's tangential extent: across
        // a periodic direction read the image (the turbulence is periodic
        // too), and clamp into the band this FAB holds.
        const amrex::Dim3 tlo = amrex::lbound(turb);
        const amrex::Dim3 thi = amrex::ubound(turb);
        const int t_lo[3] = {tlo.x, tlo.y, tlo.z};
        const int t_hi[3] = {thi.x, thi.y, thi.z};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          if (d != idir) {
            ivt[d] = amrex::min<int>(
              amrex::max<int>(wrap(d, base[d]), t_lo[d]), t_hi[d]);
          }
        }
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
          turb_fluc[n] = turb(ivt, n);
        }
      }

      pc_nscbc::Target tgt = ProblemSpecificFunctions::bcnormal_nscbc(
        x, s_N, idir, sgn, time, geom, *lprobparm, turb_fluc);
      if (tgt.type == pc_nscbc::Type::off) {
        continue; // this face is not characteristic here; try the next
      }

      // Covered cells deeper in the stencil shorten it (counted).
      auto bump = [&](const int which) {
        if (m_nscbc_diag != nullptr) {
          amrex::HostDevice::Atomic::Add(&m_nscbc_diag[which], amrex::Long(1));
        }
      };
      const int n_full = n_stencil;
      if ((n_stencil > 2) && covered(stencil_iv(2))) {
        n_stencil = 2;
        bump(pc_nscbc::Diag::body_state);
      }
      if ((n_stencil > 1) && covered(stencil_iv(1))) {
        n_stencil = 1;
        bump(pc_nscbc::Diag::body_state);
      }

      // The rest of the stencil, now that this point is characteristic.
      const amrex::IntVect ivNm1 = stencil_iv(1);
      const amrex::IntVect ivNm2 = stencil_iv(2);
      for (int n = 0; n < NVAR; n++) {
        s_Nm1[n] = dest(ivNm1, n);
        s_Nm2[n] = dest(ivNm2, n);
      }

      // Tangential neighbours of the boundary cell, for the transverse terms,
      // clamped into the domain and the FAB: at a corner the clamp collapses
      // the stencil and inv_dt falls to a one-sided spacing, or to zero if
      // both neighbours land on the same cell.  Against a tangential wall or
      // symmetry face the missing neighbour is known exactly -- it is the
      // wall image of the boundary cell itself -- so the centred difference
      // is kept there and the fill stays mirror-symmetric about the wall.
      // An EB-covered neighbour is not fluid data: the difference falls back
      // to one-sided against the boundary cell, and both sides covered leaves
      // nothing to difference.  That is geometry, expected wherever a wall
      // meets the face, so it is not counted (as a covered boundary cell is
      // not).  A neighbour the mask calls fluid but whose state is zero or
      // non-finite is bad data: same fallback, counted as body_state.
      pc_nscbc::Transverse tr;
      amrex::Real s_tm[AMREX_SPACEDIM][NVAR], s_tp[AMREX_SPACEDIM][NVAR];
      if (m_nscbc_prm[idir].beta < 1.0) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          if (d == idir) {
            continue;
          }
          const int dlo = amrex::max<int>(domlo[d], fab_lo[d]);
          const int dhi = amrex::min<int>(domhi[d], fab_hi[d]);
          int jm = amrex::max<int>(ivN[d] - 1, dlo);
          int jp = amrex::min<int>(ivN[d] + 1, dhi);
          amrex::IntVect ivm(ivN), ivp(ivN);
          ivm[d] = jm;
          ivp[d] = jp;
          for (int n = 0; n < NVAR; n++) {
            s_tm[d][n] = dest(ivm, n);
            s_tp[d][n] = dest(ivp, n);
          }
          const bool wall_m =
            (ivN[d] == domlo[d]) && (bc[d] == amrex::BCType::reflect_even);
          const bool wall_p =
            (ivN[d] == domhi[d]) &&
            (bc[d + AMREX_SPACEDIM] == amrex::BCType::reflect_even);
          if (wall_m) {
            jm = ivN[d] - 1;
            for (int n = 0; n < NVAR; n++) {
              s_tm[d][n] =
                (bcr[n].lo(d) == amrex::BCType::reflect_odd) ? -s_N[n] : s_N[n];
            }
          }
          if (wall_p) {
            jp = ivN[d] + 1;
            for (int n = 0; n < NVAR; n++) {
              s_tp[d][n] =
                (bcr[n].hi(d) == amrex::BCType::reflect_odd) ? -s_N[n] : s_N[n];
            }
          }
          const bool cm = (jm != ivN[d]) && !wall_m && covered(ivm);
          const bool cp = (jp != ivN[d]) && !wall_p && covered(ivp);
          const bool xm = (jm != ivN[d]) && !wall_m && !cm &&
                          pc_nscbc::is_body_state(s_tm[d]);
          const bool xp = (jp != ivN[d]) && !wall_p && !cp &&
                          pc_nscbc::is_body_state(s_tp[d]);
          const bool bm = cm || xm;
          const bool bp = cp || xp;
          if (xm || xp) {
            bump(pc_nscbc::Diag::body_state);
          }
          if (bm || bp) {
            if (bm) {
              jm = ivN[d];
              for (int n = 0; n < NVAR; n++) {
                s_tm[d][n] = s_N[n];
              }
            }
            if (bp) {
              jp = ivN[d];
              for (int n = 0; n < NVAR; n++) {
                s_tp[d][n] = s_N[n];
              }
            }
          }
          if (jp == jm) {
            continue;
          }
          tr.sm[d] = s_tm[d];
          tr.sp[d] = s_tp[d];
          tr.inv_dt[d] = 1.0 / (static_cast<amrex::Real>(jp - jm) * dx[d]);
          tr.valid = true;
        }
      }

      amrex::Real s_ghost[NVAR];
      if (!pc_nscbc::apply(
            s_N, s_Nm1, s_Nm2, n_stencil, dx[idir], idir, sgn, layer, tgt,
            m_nscbc_prm[idir], s_ghost, m_nscbc_diag, &tr)) {
        // The boundary cell's data is not fluid (zero EB body state): there
        // is nothing to build a characteristic state from, so this face
        // does not own the cell.  Another direction may; otherwise the
        // ordinary bcnormal() path below fills it exactly as it would
        // without the characteristic treatment.
        continue;
      }
      // A FAB holding fewer interior cells along the normal than the
      // extrapolation wants (the coarse patch of a two-level fill under a
      // fine grid on this face, say) dropped the order to what fits.
      // Counted, because the ghost then depends on the box layout; counted
      // only once the fill is known to be this face's.
      if (n_full < amrex::min<int>(3, depth_domain)) {
        bump(pc_nscbc::Diag::stencil_short);
      }
      // Per-face count of characteristic fills, so the report can say which
      // faces the treatment actually reached.
      bump(pc_nscbc::Diag::count + face);

      // Wall/symmetry image across a tangential face (see base above).
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (mirror[d] != 0) {
          for (int n = UMX; n < UMX + AMREX_SPACEDIM; ++n) {
            const int t = (mirror[d] == 1) ? bcr[n].lo(d) : bcr[n].hi(d);
            if (t == amrex::BCType::reflect_odd) {
              s_ghost[n] = -s_ghost[n];
            }
          }
        }
      }

      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ghost[n];
      }
      return true;
    }
    return false;
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int dcomp,
    const int numcomp,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    // We need the whole state to apply physical BCs
    AMREX_ALWAYS_ASSERT(dcomp == 0 && numcomp == NVAR);

    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    // Characteristic boundary treatment, where the problem asks for it.
    if (m_nscbc && nscbc_fill(iv, dest, geom, time, bcr)) {
      return;
    }

    amrex::Real s_int[NVAR] = {0.0};
    amrex::Real s_ext[NVAR] = {0.0};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> turb_fluc{0.0};

    // Fill external boundaries:
    // bcnormal populates s_ext based on s_int (the state in the 1st domain
    // cell) s_ext is initialized to a reflection (odd for velocity, even for
    // all others) of the interior state, such that if bcnormal does nothing the
    // boundary is equivalent to an adiabatic NoSlipWall. The user can provide
    // any arbitrary bcnormal in the prob.H for each case to define custom
    // combinations of inflows, outflows, and walls on the boundary face.

    // boundary conditions in x, y, [z if 3D]
    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
      if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
        // xlo, ylo, [zlo if 3D]

        // interior state at edge of domain
        amrex::IntVect loc_e{iv};
        loc_e[idir] = domlo[idir];
        for (int n = 0; n < NVAR; n++) {
          s_int[n] = dest(loc_e, n);
        }

        // interior reflected position state (odd for velocity to make
        // NoSlipWall)
        amrex::IntVect loc_r{iv};
        loc_r[idir] = domlo[idir] + (domlo[idir] - iv[idir] - 1);
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(loc_r, n);
        }
        for (int n = UMX; n < UMX + AMREX_SPACEDIM; n++) {
          s_ext[n] *= -1.0;
        }

        // turbulent fluctuations
        if (m_do_turb_inflow && (iv[idir] == domlo[idir] - 1)) {
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
            turb_fluc[n] = dest(iv, UMX + n);
          }
        }

        // Compute and populate ghost cells
        bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm, turb_fluc);
        for (int n = 0; n < NVAR; n++) {
          dest(iv, n) = s_ext[n];
        }

      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
        (iv[idir] > domhi[idir])) {
        // xhi, yhi, [zhi if 3D]

        // interior state at edge of domain
        amrex::IntVect loc_e{iv};
        loc_e[idir] = domhi[idir];
        for (int n = 0; n < NVAR; n++) {
          s_int[n] = dest(loc_e, n);
        }

        // interior reflected position state (odd for velocity to make
        // NoSlipWall)
        amrex::IntVect loc_r{iv};
        loc_r[idir] = domhi[idir] - (iv[idir] - domhi[idir] - 1);
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(loc_r, n);
        }
        for (int n = UMX; n < UMX + AMREX_SPACEDIM; n++) {
          s_ext[n] *= -1.0;
        }

        // turbulent fluctuations
        if (m_do_turb_inflow && (iv[idir] == domhi[idir] + 1)) {
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
            turb_fluc[n] = dest(iv, UMX + n);
          }
        }

        // Compute and populate ghost cells
        bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm, turb_fluc);
        for (int n = 0; n < NVAR; n++) {
          dest(iv, n) = s_ext[n];
        }
      }
    }
  }
};

struct PCReactFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dest*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
  }
};

void
pc_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{

  if (PeleC::turb_inflow.is_initialized()) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(), dir) & bx);
      if (bcr[1].lo()[dir] == amrex::BCType::ext_dir && bndryBoxLO.ok()) {
        // Create box with ghost cells and set them to zero
        amrex::IntVect growVect(PeleC::numGrow());
        growVect[dir] = 0;
        const amrex::Box modDom = amrex::grow(geom.Domain(), growVect);
        const auto bndryBoxLO_ghost =
          amrex::Box(amrex::adjCellLo(modDom, dir) & bx);
        data.setVal<amrex::RunOn::Device>(
          0.0, bndryBoxLO_ghost, UMX, AMREX_SPACEDIM);

        PeleC::turb_inflow.add_turb(
          bndryBoxLO, data, UMX, geom, time, dir, amrex::Orientation::low);
      }

      auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(), dir) & bx);
      if (bcr[1].hi()[dir] == amrex::BCType::ext_dir && bndryBoxHI.ok()) {
        // Create box with ghost cells and set them to zero
        amrex::IntVect growVect(PeleC::numGrow());
        growVect[dir] = 0;
        const amrex::Box modDom = amrex::grow(geom.Domain(), growVect);
        const auto bndryBoxHI_ghost =
          amrex::Box(amrex::adjCellHi(modDom, dir) & bx);
        data.setVal<amrex::RunOn::Device>(
          0.0, bndryBoxHI_ghost, UMX, AMREX_SPACEDIM);

        PeleC::turb_inflow.add_turb(
          bndryBoxHI, data, UMX, geom, time, dir, amrex::Orientation::high);
      }
    }
  }

  const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;

  // Capture the NSCBC parameters HOST-side (a device kernel must
  // never touch ParmParse or a class static).
  const bool nscbc = PeleC::nscbc_active();
  amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM> nscbc_prm;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    nscbc_prm[d] = PeleC::nscbc_params(d);
  }
  amrex::Long* diag = nscbc ? nscbc_diag().data() : nullptr;

  amrex::GpuArray<amrex::Array4<const char>, 2 * AMREX_SPACEDIM> covered{};
  if (nscbc) {
    if (NSCBCCoveredMask const* m = nscbc_find_mask(geom.Domain())) {
      for (int f = 0; f < 2 * AMREX_SPACEDIM; ++f) {
        if (!m->mask[f].empty()) {
          covered[f] =
            amrex::makeArray4<const char>(m->mask[f].data(), m->band[f], 1);
        }
      }
    } else {
      // Every level of an EB run registers its mask at construction; a
      // fill at a geometry without one would read covered cells as fluid.
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !PeleC::nscbc_needs_eb_mask(),
        "NSCBC: no EB covered-cell mask for this level's domain");
    }
  }

  // The characteristic fill reads the TurbInflow fluctuations from a copy
  // (see PCHypFillExtDir::m_turb): only the one-cell band add_turb wrote on
  // each ext_dir face, not the whole FAB.  Async arena: freed in stream
  // order.
  amrex::Array<amrex::FArrayBox, 2 * AMREX_SPACEDIM> turb_copy;
  amrex::GpuArray<amrex::Array4<const amrex::Real>, 2 * AMREX_SPACEDIM>
    turb_arr{};
  const bool has_turb = nscbc && PeleC::turb_inflow.is_initialized();
  if (has_turb) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      for (int side = 0; side < 2; ++side) {
        const bool ext_dir = (side == 0)
                               ? (bcr[1].lo()[dir] == amrex::BCType::ext_dir)
                               : (bcr[1].hi()[dir] == amrex::BCType::ext_dir);
        const amrex::Box band =
          ((side == 0) ? amrex::adjCellLo(geom.Domain(), dir)
                       : amrex::adjCellHi(geom.Domain(), dir)) &
          data.box();
        if (ext_dir && band.ok()) {
          auto& fab = turb_copy[2 * dir + side];
          fab.resize(band, AMREX_SPACEDIM, amrex::The_Async_Arena());
          fab.copy<amrex::RunOn::Device>(
            data, band, UMX, band, 0, AMREX_SPACEDIM);
          turb_arr[2 * dir + side] = fab.const_array();
        }
      }
    }
  }

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(
    PCHypFillExtDir{
      lprobparm, PeleC::turb_inflow.is_initialized(), nscbc, nscbc_prm, diag,
      turb_arr, has_turb, covered});
  hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

void
pc_reactfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  amrex::GpuBndryFuncFab<PCReactFillExtDir> react_bndry_func(
    PCReactFillExtDir{});
  react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

void
pc_nullfill(
  amrex::Box const& /*bx*/,
  amrex::FArrayBox& /*data*/,
  const int /*dcomp*/,
  const int /*numcomp*/,
  amrex::Geometry const& /*geom*/,
  const amrex::Real /*time*/,
  const amrex::Vector<amrex::BCRec>& /*bcr*/,
  const int /*bcomp*/,
  const int /*scomp*/)
{
}

void
PeleC::nscbc_check_fine_faces() const
{
  if (!bc_nscbc || level == 0) {
    return;
  }
  const amrex::Box& dom = geom.Domain();
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    for (int side = 0; side < 2; ++side) {
      const int t = (side == 0) ? phys_bc.lo(dir) : phys_bc.hi(dir);
      if ((t != PCPhysBCType::inflow) && (t != PCPhysBCType::user_bc)) {
        continue; // not a face the characteristic treatment can reach
      }
      const int face = (side == 0) ? dom.smallEnd(dir) : dom.bigEnd(dir);
      // A box touches when its FILLED region reaches the face, not only its
      // valid region: a box ending a few cells short still has its ghost
      // band crossing the face, and the characteristic fill then builds
      // those ghosts from coarse-interpolated stencil cells.
      const int ng = nscbc_fill_depth();
      bool touches = false;
      for (int i = 0; i < grids.size(); ++i) {
        const amrex::Box b = amrex::grow(grids[i], dir, ng);
        if (
          ((side == 0) && (b.smallEnd(dir) <= face)) ||
          ((side == 1) && (b.bigEnd(dir) >= face))) {
          touches = true;
          break;
        }
      }
      // Once per face per run: regridding can recur every few steps, and a
      // warning repeated 1600 times is a warning nobody reads.
      static bool warned[AMREX_SPACEDIM][2] = {};
      if (
        touches && !warned[dir][side] &&
        amrex::ParallelDescriptor::IOProcessor()) {
        warned[dir][side] = true;
        // A warning rather than an abort: the problem hook decides per
        // boundary POINT whether a face is characteristic, and the host
        // cannot know what it will return.  But if any point of this face
        // is characteristic, the fill's extrapolation stencil is
        // level-local, so the fine patch imposes a DIFFERENT boundary
        // condition than the coarse level does on the same face -- a
        // level-dependent artefact that refining cannot remove.
        amrex::Warning(
          "NSCBC: level " + std::to_string(level) +
          " grids (with their ghost cells) reach the domain " +
          (side == 0 ? "lo" : "hi") + " face in direction " +
          std::to_string(dir) +
          ", which is a Hard/UserBC face with pelec.bc_nscbc = 1. The "
          "characteristic fill's stencil is level-local, so a refined patch "
          "on a characteristic face makes the boundary condition "
          "level-dependent. Keep refinement away from characteristic faces "
          "(see Source/NSCBC.H). (This warning is printed once per face.)");
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  Periodic-seam gate.  Measures how far the characteristic fill is from
//  periodic where the domain is: the worst relative mismatch between ghost
//  cells and their images one period away.  The normal stencil row is read
//  in place across a periodic direction, but the transverse terms' tangential
//  neighbours (beta < 1) are clamped into the domain and the FAB (see
//  nscbc_fill), which leaves a deliberate residual at the seam under amrex's
//  corner-strip protocol -- measured 2e-4 on the inert vortex and 1.6e-3
//  with a flame front sitting on the seam corner -- so the gate REPORTS the
//  measured value and aborts only above 1e-2, an order of margin above the
//  worst known-good state and an order below a broken stencil (the naive
//  wrap measured 2.2e-2 here).  The report guards its own blind spots: the
//  tangential spread of the boundary row is printed beside the mismatch (a pass
//  on a boundary-uniform row gates nothing), and a decomposition with no image
//  pair in one FAB says NOT CHECKED instead of passing.  Exercised by
//  NSCBC-COVO/nscbc-wrapgate.inp.
// ---------------------------------------------------------------------------
void
PeleC::nscbc_check_periodic_wrap()
{
  if (!bc_nscbc || (level != 0)) {
    return;
  }
  static bool done = false;
  if (done) {
    return;
  }

  const amrex::Box& dom = geom.Domain();
  auto characteristic = [&](const int dir, const int side) {
    const int t = (side == 0) ? phys_bc.lo(dir) : phys_bc.hi(dir);
    return (t == PCPhysBCType::inflow) || (t == PCPhysBCType::user_bc);
  };

  bool relevant = false;
  for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      relevant =
        relevant || ((characteristic(idir, 0) || characteristic(idir, 1)) &&
                     (d != idir) && geom.isPeriodic(d));
    }
  }
  if (!relevant) {
    return;
  }
  done = true;

  const int ng = numGrow();
  amrex::MultiFab S(grids, dmap, NVAR, ng, amrex::MFInfo(), Factory());
  FillPatch(*this, S, ng, state[State_Type].curTime(), State_Type, 0, NVAR);

  const amrex::Real big = std::numeric_limits<amrex::Real>::max();
  const amrex::Real tiny = std::numeric_limits<amrex::Real>::min();

  // One scale per component, taken over the whole field rather than
  // pointwise: a component that is legitimately near zero at the seam (normal
  // momentum where u_n ~ 0) would otherwise turn an O(1e-3) absolute residual
  // into an O(1) relative one.  Momentum components share the momentum-vector
  // scale and species densities share the density scale, so a trace species or
  // a stagnant velocity component is measured against the flow it sits in.
  amrex::Vector<amrex::Real> h_scale(NVAR);
  for (int n = 0; n < NVAR; n++) {
    h_scale[n] = S.norminf(n, ng);
  }
  amrex::Real mom_scale = 0.0;
  for (int n = UMX; n < UMX + AMREX_SPACEDIM; n++) {
    mom_scale = amrex::max<amrex::Real>(mom_scale, h_scale[n]);
  }
  const amrex::Real rho_scale = h_scale[URHO];
  for (int n = UMX; n < UMX + AMREX_SPACEDIM; n++) {
    h_scale[n] = mom_scale;
  }
  for (int n = UFS; n < UFS + NUM_SPECIES; n++) {
    h_scale[n] = rho_scale;
  }
  for (int n = 0; n < NVAR; n++) {
    h_scale[n] = amrex::max<amrex::Real>(h_scale[n], tiny);
  }
  amrex::Gpu::DeviceVector<amrex::Real> d_scale(NVAR);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, h_scale.begin(), h_scale.end(), d_scale.begin());
  const amrex::Real* scale = d_scale.data();

  for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    for (int side = 0; side < 2; ++side) {
      if (!characteristic(idir, side)) {
        continue;
      }
      const int N_pos = (side == 0) ? dom.smallEnd(idir) : dom.bigEnd(idir);

      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if ((d == idir) || (!geom.isPeriodic(d))) {
          continue;
        }
        const int n_d = dom.length(d);
        if (n_d < ng) {
          continue; // the image slab would fold over itself
        }

        // The ghost cells of this face that sit above the domain in d.  The
        // image of each is n_d cells below, and is a ghost of this face too.
        amrex::Box reg = amrex::grow(dom, ng);
        if (side == 0) {
          reg.setBig(idir, dom.smallEnd(idir) - 1);
        } else {
          reg.setSmall(idir, dom.bigEnd(idir) + 1);
        }
        reg.setSmall(d, dom.bigEnd(d) + 1);
        reg.setBig(d, dom.bigEnd(d) + ng);

        const amrex::IntVect img = -n_d * amrex::IntVect::TheDimensionVector(d);

        amrex::ReduceOps<
          amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpMax,
          amrex::ReduceOpMin>
          op;
        amrex::ReduceData<amrex::Real, amrex::Long, amrex::Real, amrex::Real>
          rd(op);
        using RT = typename decltype(rd)::Type;

        for (amrex::MFIter mfi(S); mfi.isValid(); ++mfi) {
          auto const& a = S.const_array(mfi);
          const amrex::Box& fbx = mfi.fabbox();

          // The image pairs this FAB holds both halves of.
          amrex::Box sh(fbx);
          sh.shift(-img);
          const amrex::Box pbx = fbx & reg & sh;
          if (!pbx.isEmpty()) {
            op.eval(pbx, rd, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT {
              amrex::ignore_unused(k); // 2-D: AMREX_D_DECL drops it
              const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
              const amrex::IntVect iw = iv + img;
              amrex::Real e = 0.0;
              for (int n = 0; n < NVAR; n++) {
                const amrex::Real u = a(iv, n);
                const amrex::Real v = a(iw, n);
                // A non-finite ghost must fail the check, not vanish into
                // max(e, NaN) == e.
                if (!amrex::Math::isfinite(u) || !amrex::Math::isfinite(v)) {
                  e = big;
                  break;
                }
                e = amrex::max<amrex::Real>(
                  e, amrex::Math::abs(u - v) / scale[n]);
              }
              return {e, amrex::Long(1), -big, big};
            });
          }

          // The tangential structure of the boundary row itself: valid data,
          // so this measures whether the check above could have failed.
          amrex::Box row = dom;
          row.setSmall(idir, N_pos);
          row.setBig(idir, N_pos);
          row &= mfi.validbox();
          if (!row.isEmpty()) {
            op.eval(row, rd, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT {
              const amrex::Real r = a(i, j, k, URHO);
              return {0.0, amrex::Long(0), r, r};
            });
          }
        }

        auto hv = rd.value(op);
        amrex::Real worst = amrex::get<0>(hv);
        amrex::Long npairs = amrex::get<1>(hv);
        amrex::Real rmax = amrex::get<2>(hv);
        amrex::Real rmin = amrex::get<3>(hv);
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceLongSum(npairs);
        amrex::ParallelDescriptor::ReduceRealMax(rmax);
        amrex::ParallelDescriptor::ReduceRealMin(rmin);
        const amrex::Real spread =
          (rmax > -big)
            ? (rmax - rmin) / amrex::max<amrex::Real>(
                                amrex::Math::abs(rmax),
                                std::numeric_limits<amrex::Real>::min())
            : 0.0;

        // Relative to each component's field scale (see h_scale above).
        constexpr amrex::Real tol = 1.0e-2;
        if (worst > tol) {
          amrex::Abort(
            "NSCBC periodic-seam check FAILED on direction " +
            std::to_string(idir) + " " + (side == 0 ? "lo" : "hi") +
            " with periodic tangential direction " + std::to_string(d) +
            ": worst relative mismatch " + std::to_string(worst) + " over " +
            std::to_string(npairs) +
            " image pairs.  The characteristic fill is not periodic where the "
            "domain is.");
        }
        if (amrex::ParallelDescriptor::IOProcessor() && (verbose > 0)) {
          amrex::Print() << "  NSCBC periodic-seam check: dir " << idir
                         << (side == 0 ? " lo" : " hi")
                         << ", periodic tangential dir " << d << " -- ";
          if (npairs == 0) {
            amrex::Print()
              << "NOT CHECKED: no box holds a ghost cell and its image "
                 "together.  Raise amr.max_grid_size in direction "
              << d << " to span the domain if you want this gated.\n";
          } else {
            amrex::Print() << npairs << " image pairs agree to " << worst
                           << " (boundary-row density spread " << spread
                           << ")\n";
            if (spread == 0.0) {
              amrex::Print()
                << "    (that row is uniform along the boundary, so this pass "
                   "is vacuous -- a clamped stencil would pass it too.)\n";
            }
          }
        }
      }
    }
  }
}

pc_nscbc::Params
PeleC::nscbc_params(const int idir)
{
  pc_nscbc::Params p;
  p.sigma = bc_nscbc_sigma;
  p.relax_u = bc_nscbc_relax_u;
  p.relax_t = bc_nscbc_relax_t;
  p.order = bc_nscbc_order;
  p.beta = bc_nscbc_beta;
  p.beta_s = bc_nscbc_beta_s;
  p.pin_farfield = bc_nscbc_pin_farfield;
  // Only the ratio sigma/L_ref is physical.  L_ref is fixed to the domain
  // extent along the boundary normal so that sigma keeps the meaning it has
  // in the literature, rather than being one of two dials for one degree of
  // freedom.  Uses probhi - problo, not probhi: the legacy Fortran used
  // probhi(idir) and was therefore silently wrong for any domain not
  // anchored at the origin.
  const auto& geom = amrex::DefaultGeometry();
  p.extrap_temperature = bc_nscbc_extrap_temperature;
  p.extrap_material = bc_nscbc_extrap_material;
  p.backflow_material = bc_nscbc_backflow_material;
  p.L_ref = geom.ProbHi(idir) - geom.ProbLo(idir);
  return p;
}

int
PeleC::nscbc_fill_depth() const
{
  // The ghost depth of the level's advance FillPatch: the state's own ghost
  // cells, the hydro + explicit-filter border, and Sborder (which also
  // carries the spray ghosts).
  return amrex::max<int>(
    get_new_data(State_Type).nGrow(),
    amrex::max<int>(numGrow() + nGrowF, Sborder.nGrow()));
}

int
PeleC::nscbc_last_report()
{
  return nscbc_last_report_step;
}

bool
PeleC::nscbc_needs_eb_mask()
{
  return bc_nscbc && eb_in_domain;
}

void
PeleC::nscbc_init_eb_mask()
{
#ifdef AMREX_USE_EB
  if (!nscbc_needs_eb_mask() || (nscbc_find_mask(geom.Domain()) != nullptr)) {
    return;
  }
  if (nscbc_masks_p == nullptr) {
    nscbc_masks_p = new amrex::Vector<std::unique_ptr<NSCBCCoveredMask>>();
    amrex::ExecOnFinalize([]() {
      delete nscbc_masks_p;
      nscbc_masks_p = nullptr;
    });
  }
  auto m = std::make_unique<NSCBCCoveredMask>();
  const amrex::Box& dom = geom.Domain();
  m->domain = dom;
  const auto& eb_level = amrex::EB2::IndexSpace::top().getLevel(geom);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    for (int side = 0; side < 2; ++side) {
      const int t = (side == 0) ? phys_bc.lo(dir) : phys_bc.hi(dir);
      if ((t != PCPhysBCType::inflow) && (t != PCPhysBCType::user_bc)) {
        continue; // not a face the characteristic treatment can reach
      }
      // The stencil's cells: three deep along the normal, the whole domain
      // across it (the fill only ever reads in-domain rows; a periodic
      // index is wrapped into the domain before the lookup).
      amrex::Box band = dom;
      if (side == 0) {
        band.setBig(
          dir, amrex::min<int>(dom.smallEnd(dir) + 2, dom.bigEnd(dir)));
      } else {
        band.setSmall(
          dir, amrex::max<int>(dom.bigEnd(dir) - 2, dom.smallEnd(dir)));
      }
      // Collective: every rank fills its share of the band's volume
      // fraction from the EB index space, and an all-reduce assembles the
      // whole band on every rank (a boundary fill may need any part of it).
      amrex::BoxArray ba(band);
      ba.maxSize(64);
      const amrex::DistributionMapping dm(ba);
      amrex::MultiFab vf(ba, dm, 1, 0);
      eb_level.fillVolFrac(vf, geom);
      const auto npts = static_cast<std::size_t>(band.numPts());
      amrex::Vector<int> h(npts, 0);
      for (amrex::MFIter mfi(vf); mfi.isValid(); ++mfi) {
        const amrex::Box& b = mfi.validbox();
        amrex::FArrayBox hf(b, 1, amrex::The_Pinned_Arena());
        hf.copy<amrex::RunOn::Device>(vf[mfi], b, 0, b, 0, 1);
        amrex::Gpu::streamSynchronize();
        const auto ha = hf.const_array();
        const auto ma = amrex::makeArray4<int>(h.data(), band, 1);
        amrex::LoopOnCpu(b, [&](int i, int j, int k) {
          ma(i, j, k) = (ha(i, j, k) > 0.0) ? 0 : 1;
        });
      }
      amrex::ParallelDescriptor::ReduceIntMax(h.data(), static_cast<int>(npts));
      // Stored one byte per cell: the mask is replicated on every rank.
      amrex::Vector<char> hc(npts);
      for (std::size_t i = 0; i < npts; ++i) {
        hc[i] = static_cast<char>(h[i] != 0);
      }
      const int f = 2 * dir + side;
      m->band[f] = band;
      m->mask[f].resize(npts);
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, hc.begin(), hc.end(), m->mask[f].begin());
    }
  }
  amrex::Gpu::streamSynchronize();
  nscbc_masks_p->push_back(std::move(m));
#endif
}

void
PeleC::nscbc_init_diagnostics()
{
  if (nscbc_diag_p != nullptr) {
    return;
  }
  nscbc_diag_p = new amrex::Gpu::DeviceVector<amrex::Long>(nscbc_ncounters, 0);
  amrex::ExecOnFinalize([]() {
    // Whatever counted since the last report is reported here, before the
    // counters go: a run that ends between two reports, or before the
    // first, must not drop its fallbacks.  The finalize callbacks run before
    // AMReX tears down the arenas and MPI.
    nscbc_report_diagnostics(-1);
    delete nscbc_diag_p;
    nscbc_diag_p = nullptr;
  });
}

void
PeleC::nscbc_report_diagnostics(const int step)
{
  if (!bc_nscbc || (nscbc_diag_p == nullptr)) {
    return;
  }
  if (step >= 0) {
    nscbc_last_report_step = step;
  }
  // Settle any counter atomics still in flight BEFORE the copy: the fills
  // bump these counters from MFIter's rotating streams, and the blocking
  // Gpu::copy below synchronises only its own stream.  A count landing
  // between the copy and the reset would never be reported.
  amrex::Gpu::Device::synchronize();
  std::vector<amrex::Long> h(nscbc_ncounters, 0);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, nscbc_diag().begin(), nscbc_diag().end(),
    h.begin());
  amrex::ParallelDescriptor::ReduceLongSum(h.data(), nscbc_ncounters);

  // Which faces the treatment actually reached.  bc_nscbc = 1 with a hook
  // that is never called (an outlet left FOExtrap, or a bcnormal_nscbc that
  // is not a static member of the problem's ProblemSpecificFunctions) or
  // that returns off everywhere would otherwise look exactly like a healthy
  // boundary: only fallbacks print.
  static const char* face_name[6] = {"x-lo", "x-hi", "y-lo",
                                     "y-hi", "z-lo", "z-hi"};
  amrex::Long all_fills = 0;
  for (int f = 0; f < 2 * AMREX_SPACEDIM; ++f) {
    nscbc_face_fills_total[f] += h[pc_nscbc::Diag::count + f];
    all_fills += nscbc_face_fills_total[f];
    if (
      (nscbc_face_fills_total[f] > 0) && !nscbc_face_announced[f] &&
      amrex::ParallelDescriptor::IOProcessor()) {
      nscbc_face_announced[f] = true;
      amrex::Print() << "  NSCBC: characteristic fills active on the "
                     << face_name[f] << " face\n";
    }
  }
  // Warn once the run has advanced past the first report (the end-of-run
  // report passes -1) without a single characteristic fill.
  static int first_step = -2;
  if (first_step == -2) {
    first_step = step;
  }
  if (
    ((step < 0) || (step > first_step)) && (all_fills == 0) &&
    !nscbc_none_warned) {
    nscbc_none_warned = true;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Warning(
        "NSCBC: pelec.bc_nscbc = 1 but no boundary point has been filled by "
        "the characteristic treatment. It acts only on Hard/UserBC faces, "
        "and only where ProblemSpecificFunctions::bcnormal_nscbc -- a static "
        "member of the struct the problem's prob.H aliases to "
        "ProblemSpecificFunctions -- returns a type other than off.");
    }
  }

  // The supersonic path is exact, not a degradation, so it is reported but is
  // not a warning.  The others mean the boundary is being asked for something
  // it cannot cleanly provide.
  const amrex::Long total =
    h[pc_nscbc::Diag::reversed] + h[pc_nscbc::Diag::body_state] +
    h[pc_nscbc::Diag::eos_failure] + h[pc_nscbc::Diag::floored] +
    h[pc_nscbc::Diag::transverse_drop] + h[pc_nscbc::Diag::source_drop] +
    h[pc_nscbc::Diag::stencil_short] + h[pc_nscbc::Diag::target_incomplete] +
    h[pc_nscbc::Diag::passive_unset];
  const amrex::Long structure = h[pc_nscbc::Diag::structure];
  if (
    amrex::ParallelDescriptor::IOProcessor() &&
    (total > 0 || structure > 0 || verbose > 1)) {
    amrex::Print() << "  NSCBC fallbacks since last report:" << "  supersonic "
                   << h[pc_nscbc::Diag::supersonic] << ",  flow reversal "
                   << h[pc_nscbc::Diag::reversed] << ",  EB covered stencil "
                   << h[pc_nscbc::Diag::body_state] << ",  EOS failure "
                   << h[pc_nscbc::Diag::eos_failure] << ",  floored "
                   << h[pc_nscbc::Diag::floored] << ",  transverse dropped "
                   << h[pc_nscbc::Diag::transverse_drop] << ",  source dropped "
                   << h[pc_nscbc::Diag::source_drop] << ",  short stencil "
                   << h[pc_nscbc::Diag::stencil_short]
                   << ",  target incomplete "
                   << h[pc_nscbc::Diag::target_incomplete]
                   << ",  inflow passives unset "
                   << h[pc_nscbc::Diag::passive_unset] << "\n";
    if (structure > 0) {
      // Advisory, not a fallback: a front is in the outflow boundary cells.
      // A flame there is what the flame closures exist for; a front without
      // heat release is not (NSCBC-FlameOutflow/README.md, the exit test).
      amrex::Print()
        << "  NSCBC: material structure (|dS| > 5% of rho per cell) sat in "
        << structure << " outflow boundary-cell fills since the last report.\n"
        << "         A front is on this outflow.  If it burns, set "
           "bc_nscbc_extrap_temperature = 1 and bc_nscbc_beta_s = 0\n"
        << "         before it crosses; if it does not, raise bc_nscbc_sigma "
           "toward 16 instead (extrap_temperature alone\n"
        << "         over-vents an inert front at low sigma).  Keep "
           "bc_nscbc_extrap_material off during a transit.\n"
        << "         See the BCs chapter and NSCBC-FlameOutflow/README.md.\n";
    }
  }
  nscbc_diag().assign(nscbc_ncounters, 0);
  // The reset itself is a device fill on the current stream; the next
  // advance's fills bump these counters from MFIter's rotating streams,
  // which are not ordered against it.  Settle the reset before returning
  // so a zero can never land on top of a fresh count.
  amrex::Gpu::streamSynchronize();
}
