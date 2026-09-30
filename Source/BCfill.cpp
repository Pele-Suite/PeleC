#include <limits>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Reduce.H>

#include "PeleC.H"
#include "prob.H"
#include "NSCBC.H"

namespace {
// Device-resident fallback counters for the characteristic boundary
// treatment.  Allocated on first use, reduced and reported by
// PeleC::nscbc_report_diagnostics().
//
// Held as a heap pointer released through amrex::ExecOnFinalize rather than as
// a function-local static object.  A static Gpu::DeviceVector destructs at
// program exit, which is AFTER amrex::Finalize() has torn down the arena it
// allocated from -- harmless on a CPU build and a use-after-free of the device
// allocator on a GPU one.
amrex::Gpu::DeviceVector<amrex::Long>* nscbc_diag_p = nullptr;

amrex::Gpu::DeviceVector<amrex::Long>&
nscbc_diag()
{
  if (nscbc_diag_p == nullptr) {
    nscbc_diag_p =
      new amrex::Gpu::DeviceVector<amrex::Long>(pc_nscbc::Diag::count, 0);
    amrex::ExecOnFinalize([]() {
      delete nscbc_diag_p;
      nscbc_diag_p = nullptr;
    });
  }
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

  AMREX_GPU_HOST
  explicit PCHypFillExtDir(
    const ProbParmDevice* d_prob_parm,
    const bool do_turb_inflow,
    const bool nscbc,
    const amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM>& nscbc_prm,
    amrex::Long* nscbc_diag)
    : lprobparm(d_prob_parm),
      m_do_turb_inflow(do_turb_inflow),
      m_nscbc(nscbc),
      m_nscbc_prm(nscbc_prm),
      m_nscbc_diag(nscbc_diag)
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
  //  outside on an ext_dir face whose problem hook returns a live target;
  //  the state is written exactly once.  Combined with the clamped stencil
  //  below this makes the fill a pure function of valid interior data, hence
  //  independent of the order in which ghost cells are visited and identical
  //  on CPU and GPU.
  //
  //  (Note that the legacy bcnormal() path further down does NOT have this
  //  property: at a corner it reads dest() at a tangential index that is
  //  itself a ghost cell, which another thread in the same launch may be
  //  writing.  That is pre-existing behaviour and is left untouched here so
  //  that no existing result moves.)
  // -------------------------------------------------------------------------
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  bool nscbc_fill(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const int* bc) const
  {
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* dx = geom.CellSize();

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

      const int N_pos = (sgn > 0) ? domlo[idir] : domhi[idir];
      const int layer = sgn * (N_pos - iv[idir]); // 1 = nearest the domain

      const amrex::Dim3 lo3 = amrex::lbound(dest);
      const amrex::Dim3 hi3 = amrex::ubound(dest);
      const int fab_lo[3] = {lo3.x, lo3.y, lo3.z};
      const int fab_hi[3] = {hi3.x, hi3.y, hi3.z};

      // Tangential indices are clamped into the domain (and the FAB), in
      // every tangential direction: only valid interior cells are read, and
      // the fill is a pure function of pre-launch data.
      //
      // In a PERIODIC tangential direction this leaves a small, measured
      // seam residual, and that is a deliberate trade.  amrex's corner
      // protocol (StateDataPhysBCFunct) recomputes corner ghosts on a strip
      // FAB holding only their image band, so a seam-adjacent cell is filled
      // twice on two different FABs; exact agreement requires restricting
      // the tangential stencil to the band the strip can see.  Both
      // alternatives were built and measured: a wrap through the resident
      // images leaves the array aperiodic at the 2e-2 level, and the
      // strip-consistent band is bitwise-periodic and measured equally
      // stable (NSCBC-FlameOutflow-DRM, beta = 0.5) -- but costs ~100 lines
      // of decomposition-sensitive index machinery to remove a residual that
      // measures 2e-4 (inert) to 1.6e-3 (flame on the seam corner).  The
      // clamp was kept for simplicity; nscbc_check_periodic_wrap() reports
      // the residual and aborts above 1e-2, which a broken stencil fails.
      amrex::IntVect base(AMREX_D_DECL(iv[0], iv[1], iv[2]));
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != idir) {
          base[d] = amrex::min<int>(
            amrex::max<int>(iv[d], amrex::max<int>(domlo[d], fab_lo[d])),
            amrex::min<int>(domhi[d], fab_hi[d]));
        }
      }

      // How deep a normal stencil is available, in the domain and in the FAB.
      const int depth_domain = domhi[idir] - domlo[idir] + 1;
      const int depth_fab =
        (sgn > 0) ? (fab_hi[idir] - N_pos + 1) : (N_pos - fab_lo[idir] + 1);
      const int n_stencil =
        amrex::min<int>(3, amrex::min<int>(depth_domain, depth_fab));
      if (n_stencil < 1) {
        continue;
      }

      auto stencil_iv = [&](const int step) {
        amrex::IntVect r = base;
        r[idir] = N_pos + sgn * amrex::min<int>(step, n_stencil - 1);
        return r;
      };
      amrex::Real s_N[NVAR], s_Nm1[NVAR], s_Nm2[NVAR];
      const amrex::IntVect ivN = stencil_iv(0);
      const amrex::IntVect ivNm1 = stencil_iv(1);
      const amrex::IntVect ivNm2 = stencil_iv(2);
      for (int n = 0; n < NVAR; n++) {
        s_N[n] = dest(ivN, n);
        s_Nm1[n] = dest(ivNm1, n);
        s_Nm2[n] = dest(ivNm2, n);
      }

      // Query the problem for this boundary POINT.  x is the location on the
      // boundary plane, not the ghost cell centre: the target is a property
      // of the boundary point and must not vary with the ghost layer, or the
      // relaxation would be applied to a different target in each layer.
      amrex::Real x[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        x[d] = prob_lo[d] + (static_cast<amrex::Real>(base[d]) + 0.5) * dx[d];
      }
      x[idir] = prob_lo[idir] + static_cast<amrex::Real>(
                                  (sgn > 0) ? domlo[idir] : (domhi[idir] + 1)) *
                                  dx[idir];

      pc_nscbc::Target tgt = ProblemSpecificFunctions::bcnormal_nscbc(
        x, s_N, idir, sgn, time, geom, *lprobparm);
      if (tgt.type == pc_nscbc::Type::off) {
        continue; // this face is not characteristic here; try the next
      }

      // Boundary-register composition (level 0 only; the registers change
      // once per advance and are read frozen here).  The kernel below sees
      // only the composed Target -- it stays a pure function.
      // Tangential neighbours of the boundary cell, for the transverse terms,
      // clamped into the domain: at a corner the clamp collapses the stencil
      // and inv_dt falls to a one-sided spacing, or to zero if both
      // neighbours land on the same cell.
      pc_nscbc::Transverse tr;
      amrex::Real s_tm[AMREX_SPACEDIM][NVAR], s_tp[AMREX_SPACEDIM][NVAR];
      if (m_nscbc_prm[idir].beta < 1.0) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          if (d == idir) {
            continue;
          }
          const int dlo = amrex::max<int>(domlo[d], fab_lo[d]);
          const int dhi = amrex::min<int>(domhi[d], fab_hi[d]);
          const int jm = amrex::max<int>(ivN[d] - 1, dlo);
          const int jp = amrex::min<int>(ivN[d] + 1, dhi);
          if (jp == jm) {
            continue;
          }
          amrex::IntVect ivm(ivN), ivp(ivN);
          ivm[d] = jm;
          ivp[d] = jp;
          for (int n = 0; n < NVAR; n++) {
            s_tm[d][n] = dest(ivm, n);
            s_tp[d][n] = dest(ivp, n);
          }
          tr.sm[d] = s_tm[d];
          tr.sp[d] = s_tp[d];
          tr.inv_dt[d] = 1.0 / (static_cast<amrex::Real>(jp - jm) * dx[d]);
          tr.valid = true;
        }
      }

      amrex::Real s_ghost[NVAR];
      pc_nscbc::apply(
        s_N, s_Nm1, s_Nm2, n_stencil, dx[idir], idir, sgn, layer, tgt,
        m_nscbc_prm[idir], s_ghost, m_nscbc_diag, &tr);
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
    if (m_nscbc && nscbc_fill(iv, dest, geom, time, bc)) {
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

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(
    PCHypFillExtDir{
      lprobparm, PeleC::turb_inflow.is_initialized(), nscbc, nscbc_prm, diag});
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
      bool touches = false;
      for (int i = 0; i < grids.size(); ++i) {
        const amrex::Box& b = grids[i];
        if (
          ((side == 0) && (b.smallEnd(dir) == face)) ||
          ((side == 1) && (b.bigEnd(dir) == face))) {
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
          "NSCBC: level " + std::to_string(level) + " grids touch the domain " +
          (side == 0 ? "lo" : "hi") + " face in direction " +
          std::to_string(dir) +
          ", which is a Hard/UserBC face with pelec.bc_nscbc = 1. The "
          "characteristic fill's stencil is level-local, so a refined patch "
          "on a characteristic face makes the boundary condition "
          "level-dependent. Keep refinement away from characteristic faces "
          "(see the BCs chapter). (This warning is printed once per face.)");
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  Periodic-seam gate.  Measures how far the characteristic fill is from
//  periodic where the domain is: the worst relative mismatch between ghost
//  cells and their images one period away.  The clamped tangential stencil
//  (see nscbc_fill) leaves a deliberate residual at the seam under amrex's
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
                const amrex::Real s = amrex::max<amrex::Real>(
                  amrex::Math::abs(u),
                  amrex::max<amrex::Real>(
                    amrex::Math::abs(v),
                    std::numeric_limits<amrex::Real>::min()));
                e = amrex::max<amrex::Real>(e, amrex::Math::abs(u - v) / s);
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

void
PeleC::nscbc_report_diagnostics()
{
  if (!bc_nscbc) {
    return;
  }
  std::vector<amrex::Long> h(pc_nscbc::Diag::count, 0);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, nscbc_diag().begin(), nscbc_diag().end(),
    h.begin());
  amrex::ParallelDescriptor::ReduceLongSum(h.data(), pc_nscbc::Diag::count);

  // The supersonic path is exact, not a degradation, so it is reported but is
  // not a warning.  The others mean the boundary is being asked for something
  // it cannot cleanly provide.
  // transverse_drop and source_drop belong here as much as the rest: a
  // beta or beta_s that is silently not being applied looks exactly like a
  // beta or beta_s that does nothing, and the only way to tell the two apart
  // is to count it.
  const amrex::Long total =
    h[pc_nscbc::Diag::reversed] + h[pc_nscbc::Diag::body_state] +
    h[pc_nscbc::Diag::eos_failure] + h[pc_nscbc::Diag::floored] +
    h[pc_nscbc::Diag::transverse_drop] + h[pc_nscbc::Diag::source_drop] +
    h[pc_nscbc::Diag::target_incomplete];
  const amrex::Long structure = h[pc_nscbc::Diag::structure];
  if (
    amrex::ParallelDescriptor::IOProcessor() &&
    (total > 0 || structure > 0 || verbose > 1)) {
    amrex::Print() << "  NSCBC fallbacks since last report:" << "  supersonic "
                   << h[pc_nscbc::Diag::supersonic] << ",  flow reversal "
                   << h[pc_nscbc::Diag::reversed] << ",  EB body state "
                   << h[pc_nscbc::Diag::body_state] << ",  EOS failure "
                   << h[pc_nscbc::Diag::eos_failure] << ",  floored "
                   << h[pc_nscbc::Diag::floored] << ",  transverse dropped "
                   << h[pc_nscbc::Diag::transverse_drop] << ",  source dropped "
                   << h[pc_nscbc::Diag::source_drop] << ",  target incomplete "
                   << h[pc_nscbc::Diag::target_incomplete] << "\n";
    if (structure > 0) {
      // Advisory, not a fallback: a front is in the outflow boundary cells,
      // which is the configuration the flame closures exist for.
      amrex::Print()
        << "  NSCBC: material structure (|dS| > 5% of rho per cell) sat in "
        << structure << " outflow boundary-cell fills since the last report.\n"
        << "         A flame or front is on this outflow: set "
           "bc_nscbc_extrap_temperature = 1 and bc_nscbc_beta_s = 0 (any "
           "sigma then\n"
        << "         survives a crossing), and keep bc_nscbc_extrap_material "
           "off during a transit.  See the BCs chapter and "
           "NSCBC-FlameOutflow/README.md.\n";
    }
  }
  // Settle any counter atomics still in flight on other streams before the
  // reset; the blocking Gpu::copy above synchronised only its own stream.
  amrex::Gpu::Device::synchronize();
  nscbc_diag().assign(pc_nscbc::Diag::count, 0);
  // The reset itself is a device fill on the current stream; the next
  // advance's fills bump these counters from MFIter's rotating streams,
  // which are not ordered against it.  Settle the reset before returning
  // so a zero can never land on top of a fresh count.
  amrex::Gpu::streamSynchronize();
}
