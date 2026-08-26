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

// ---------------------------------------------------------------------------
//  Boundary registers (Docs/NSCBC-boundary-registers-design.md).
//
//  Per-point state on the LEVEL-0 boundary-adjacent band of every domain
//  face: the learned reference means (EMAs of the outgoing acoustic
//  invariant and of the pressure, tau = 3 t_a) and the ready-to-add
//  composition terms derived from them once per advance.  The band of a
//  face is a structured (SPACEDIM-1)-dimensional strip over the domain's
//  tangential extent, so the whole store is one small flat device array
//  indexed by face and flattened tangential index -- no MultiFab, no
//  communication: every register a fill reads belongs to a band cell
//  inside the same FAB (the tangential clamp guarantees it), and every
//  band cell is updated by the one rank that owns it.  The invariant that
//  makes the checkpoint gather trivial: entries a rank does not own are
//  EXACTLY ZERO on that rank, always.
//
//  The fill reads the registers frozen (they change once per advance,
//  outside every fill path), and the kernel in NSCBC.H never sees them at
//  all: BCfill composes the effective Target -- u += u_minus at a
//  registered inflow (NRI), p += dp_ac at a registered outflow (NDNR) --
//  and hands the same pure function the same pure arguments.
// ---------------------------------------------------------------------------
namespace reg {
enum : int {
  ema_Rout = 0, // EMA of the outgoing invariant u_out + p/(rho c)
  ema_p,        // EMA of the boundary-cell pressure
  u_minus,      // ready-to-add LAB-frame velocity correction (NRI inflow)
  dp_ac,        // ready-to-add pressure correction (NDNR outflow)
  trend_dS,     // EMA of |dS|/rho: the slow material-structure trend
                // (phase C, advisory -- reported, consumed by nothing)
  ema_rhoc,     // slowly-EMA'd reference impedance.  The invariant MUST be
                // referenced to a time-frozen rho c: the mean p/(rho c) is
                // ~1e4 cm/s, so a coherent 0.1% oscillation of the local
                // impedance aliases tens of cm/s into R_out -- signal
                // amplitude.  (The kernel's per-fill frozen rho_c is safe
                // because it differences across the stencil at one instant;
                // the register differences across time.)  Measured: local
                // impedance degrades I_in from 1.03 to 1.6 in the driver
                // and identically in PeleC.
  init,         // 1 once this entry has been seeded
  ncomp
};
}

amrex::Gpu::DeviceVector<amrex::Real>* nscbc_reg_p = nullptr;
amrex::Box nscbc_reg_dom; // the domain the store was sized for

long
nscbc_reg_total(const amrex::Box& dom)
{
  long total = 0;
  for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    long ext = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (d != idir) {
        ext *= dom.length(d);
      }
    }
    total += 2 * ext; // lo and hi face
  }
  return total;
}

// Allocation happens ONLY where the domain is authoritative -- the update
// and the restart, which pass the level-0 geometry's Domain().  Everything
// else peeks: a null peek simply means "no update has run yet", and the
// pre-update registers are all-zero classical behaviour anyway.
amrex::Gpu::DeviceVector<amrex::Real>&
nscbc_registers(const amrex::Box& dom)
{
  const long want = nscbc_reg_total(dom) * reg::ncomp;
  if (nscbc_reg_p == nullptr) {
    nscbc_reg_p = new amrex::Gpu::DeviceVector<amrex::Real>(want, 0.0);
    amrex::ExecOnFinalize([]() {
      delete nscbc_reg_p;
      nscbc_reg_p = nullptr;
    });
  } else if (static_cast<long>(nscbc_reg_p->size()) != want) {
    nscbc_reg_p->clear();
    nscbc_reg_p->resize(want, 0.0);
  }
  nscbc_reg_dom = dom;
  return *nscbc_reg_p;
}

// The peek doubles as the level test: only the geometry the store was
// sized for -- level 0, recorded at update/restart time -- gets a pointer.
// (amrex::DefaultGeometry() proved unreliable for this comparison, which
// is how both the zero-size allocation and a silently-classical fill
// happened; matching the recorded domain is self-consistent.)
const amrex::Real*
nscbc_registers_peek(const amrex::Box& dom)
{
  return ((nscbc_reg_p != nullptr) && !nscbc_reg_p->empty() &&
          (dom == nscbc_reg_dom))
           ? nscbc_reg_p->data()
           : nullptr;
}

// Offset of face (idir, side) into the flat register array, and the
// flattened tangential index of a band cell.  Both host- and device-safe,
// domain passed explicitly so the device never touches globals.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE long
nscbc_reg_face_offset(const int idir, const int side, const amrex::Box& dom)
{
  long off = 0;
  for (int fd = 0; fd < AMREX_SPACEDIM; ++fd) {
    long ext = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (d != fd) {
        ext *= dom.length(d);
      }
    }
    for (int s = 0; s < 2; ++s) {
      if ((fd == idir) && (s == side)) {
        return off;
      }
      off += ext;
    }
  }
  return off; // unreachable
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE long
nscbc_reg_index(
  const int idir,
  const int side,
  const amrex::IntVect& iv,
  const amrex::Box& dom)
{
  long idx = nscbc_reg_face_offset(idir, side, dom);
  long stride = 1;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (d != idir) {
      idx += (iv[d] - dom.smallEnd(d)) * stride;
      stride *= dom.length(d);
    }
  }
  return idx * reg::ncomp;
}
} // namespace

struct PCHypFillExtDir
{
  ProbParmDevice const* lprobparm;
  bool m_do_turb_inflow{false};
  bool m_nscbc{false};
  amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM> m_nscbc_prm;
  amrex::Long* m_nscbc_diag{nullptr};
  // Boundary registers: null when neither NRI nor NDNR is enabled; the
  // mode bits say which composition applies (1 = NRI inflow, 2 = NDNR
  // outflow).  Level-0 only -- see the register banner in this file.
  const amrex::Real* m_nscbc_reg{nullptr};
  int m_reg_mode{0};

  AMREX_GPU_HOST
  explicit PCHypFillExtDir(
    const ProbParmDevice* d_prob_parm,
    const bool do_turb_inflow,
    const bool nscbc,
    const amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM>& nscbc_prm,
    amrex::Long* nscbc_diag,
    const amrex::Real* nscbc_reg = nullptr,
    const int reg_mode = 0)
    : lprobparm(d_prob_parm),
      m_do_turb_inflow(do_turb_inflow),
      m_nscbc(nscbc),
      m_nscbc_prm(nscbc_prm),
      m_nscbc_diag(nscbc_diag),
      m_nscbc_reg(nscbc_reg),
      m_reg_mode(reg_mode)
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
      if (m_nscbc_reg != nullptr) {
        const int side = (sgn > 0) ? 0 : 1;
        const long ridx =
          nscbc_reg_index(idir, side, base, geom.Domain());
        if (
          (tgt.type == pc_nscbc::Type::inflow) && ((m_reg_mode & 1) != 0)) {
          tgt.u[idir] += m_nscbc_reg[ridx + reg::u_minus];
        }
        if (
          (tgt.type == pc_nscbc::Type::outflow) && ((m_reg_mode & 2) != 0)) {
          tgt.p += m_nscbc_reg[ridx + reg::dp_ac];
        }
      }

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

  // Boundary registers: the peek returns a pointer only for the geometry
  // the store was sized for (level 0, recorded at update time), and only
  // when a register consumer is enabled.
  const amrex::Real* reg_ptr = nullptr;
  int reg_mode = 0;
  if (nscbc) {
    reg_mode = PeleC::nscbc_register_mode();
    if (reg_mode != 0) {
      reg_ptr = nscbc_registers_peek(geom.Domain());
    }
  }

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(PCHypFillExtDir{
    lprobparm, PeleC::turb_inflow.is_initialized(), nscbc, nscbc_prm, diag,
    reg_ptr, reg_mode});
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
      if (touches && !warned[dir][side] &&
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
          " grids touch the domain " + (side == 0 ? "lo" : "hi") +
          " face in direction " + std::to_string(dir) +
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
//  wrap measured 2.2e-2 here).  The report guards its own blind spots: the tangential spread
//  of the boundary row is printed beside the mismatch (a pass on a
//  boundary-uniform row gates nothing), and a decomposition with no image
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
                     (d != idir) && (geom.isPeriodic(d) != 0));
    }
  }
  if (!relevant) {
    return;
  }
  done = true;

  const int ng = numGrow();
  amrex::MultiFab S(grids, dmap, NVAR, ng, amrex::MFInfo(), Factory());
  FillPatch(
    *this, S, ng, state[State_Type].curTime(), State_Type, 0, NVAR);

  const amrex::Real big = std::numeric_limits<amrex::Real>::max();

  for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    for (int side = 0; side < 2; ++side) {
      if (!characteristic(idir, side)) {
        continue;
      }
      const int N_pos = (side == 0) ? dom.smallEnd(idir) : dom.bigEnd(idir);

      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if ((d == idir) || (geom.isPeriodic(d) == 0)) {
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
            op.eval(
              pbx, rd, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const amrex::IntVect iw = iv + img;
                amrex::Real e = 0.0;
                for (int n = 0; n < NVAR; n++) {
                  const amrex::Real u = a(iv, n);
                  const amrex::Real v = a(iw, n);
                  const amrex::Real s = amrex::max<amrex::Real>(
                    amrex::Math::abs(u), amrex::max<amrex::Real>(amrex::Math::abs(v), 1.0e-300));
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
            op.eval(
              row, rd, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT {
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
          (rmax > -big) ? (rmax - rmin) / amrex::max<amrex::Real>(
                                            amrex::Math::abs(rmax), 1.0e-300)
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
          amrex::Print()
            << "  NSCBC periodic-seam check: dir " << idir
            << (side == 0 ? " lo" : " hi") << ", periodic tangential dir " << d
            << " -- ";
          if (npairs == 0) {
            amrex::Print()
              << "NOT CHECKED: no box holds a ghost cell and its image "
                 "together.  Raise amr.max_grid_size in direction "
              << d << " to span the domain if you want this gated.\n";
          } else {
            amrex::Print()
              << npairs << " image pairs agree to " << worst
              << " (boundary-row density spread " << spread << ")\n";
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

// The once-per-advance register update (design note, phase A/B/C).  Runs on
// level 0 after the new-time state exists; per-point local, so it is
// decomposition-independent, and it writes only the band entries this rank
// owns -- the zero-elsewhere invariant the checkpoint gather relies on.
// Every domain face band is updated regardless of its BC type: registers on
// non-characteristic faces are simply never read.
void
PeleC::nscbc_update_registers(const amrex::Real dt)
{
  if (
    (level != 0) || !bc_nscbc || !(bc_nscbc_nri || bc_nscbc_ndnr) ||
    (dt <= 0.0)) {
    return;
  }
  // One-time announce: silent stateful machinery is unverifiable (the
  // counters lesson) -- say once that the registers are live and for whom.
  static bool announced = false;
  if (!announced && amrex::ParallelDescriptor::IOProcessor()) {
    announced = true;
    amrex::Print() << "NSCBC boundary registers active (NRI=" << bc_nscbc_nri
                   << ", NDNR=" << bc_nscbc_ndnr << ")\n";
  }
  const amrex::Box dom = geom.Domain();
  const amrex::MultiFab& S = get_new_data(State_Type);
  amrex::Real* rp = nscbc_registers(dom).data();

  for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    const amrex::Real L_ref = geom.ProbHi(idir) - geom.ProbLo(idir);
    for (int side = 0; side < 2; ++side) {
      amrex::Box band = dom;
      if (side == 0) {
        band.setBig(idir, dom.smallEnd(idir));
      } else {
        band.setSmall(idir, dom.bigEnd(idir));
      }
      const int sgn = (side == 0) ? +1 : -1;
      const amrex::Real n_sgn = static_cast<amrex::Real>(-sgn);
      const int e_int = sgn; // step toward the interior along idir

      for (amrex::MFIter mfi(S, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox() & band;
        if (!b.ok()) {
          continue;
        }
        auto const& s = S.const_array(mfi);
        const amrex::Box sbox = mfi.validbox();
        amrex::ParallelFor(
          b, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            auto eos = pele::physics::PhysicsType::eos();
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            auto prims = [&](const amrex::IntVect& p_iv, amrex::Real& rho,
                             amrex::Real& u_n, amrex::Real& p,
                             amrex::Real& c) {
              rho = s(p_iv, URHO);
              u_n = s(p_iv, UMX + idir) / rho;
              const amrex::Real T = s(p_iv, UTEMP);
              amrex::Real Y[NUM_SPECIES];
              for (int n = 0; n < NUM_SPECIES; ++n) {
                Y[n] = s(p_iv, UFS + n) / rho;
              }
              eos.RTY2P(rho, T, Y, p);
              eos.RTY2Cs(rho, T, Y, c);
            };
            amrex::Real rho = 0.0, u_n = 0.0, p = 0.0, c = 0.0;
            prims(iv, rho, u_n, p, c);
            if (!(amrex::Math::isfinite(p) && (c > 0.0) && (rho > 0.0))) {
              return; // leave the register untouched on a sick state
            }
            const amrex::Real rhoc = rho * c;
            const amrex::Real rc_ref = (rp == nullptr)
                                         ? rhoc
                                         : ((rp[nscbc_reg_index(
                                               idir, side, iv, dom) +
                                             reg::init] == 0.0)
                                              ? rhoc
                                              : rp[nscbc_reg_index(
                                                     idir, side, iv, dom) +
                                                   reg::ema_rhoc]);
            const amrex::Real R_out = n_sgn * u_n + p / rc_ref;

            const long ridx = nscbc_reg_index(idir, side, iv, dom);
            const amrex::Real tau = 6.0 * L_ref / c; // 3 t_a = 6 L/c
            const amrex::Real w = dt / (tau + dt);

            // Phase C's advisory trend: the entropy-family slope at the
            // band cell, frozen impedance, same limiter as the kernel.
            amrex::Real dS_abs = 0.0;
            {
              amrex::IntVect iv1 = iv, iv2 = iv;
              iv1[idir] += e_int;
              iv2[idir] += 2 * e_int;
              if (sbox.contains(iv1) && sbox.contains(iv2)) {
                const amrex::Real inv_c2 = 1.0 / (c * c);
                amrex::Real r1, u1, p1, c1, r2, u2, p2, c2;
                prims(iv1, r1, u1, p1, c1);
                prims(iv2, r2, u2, p2, c2);
                const amrex::Real S0 = rho - p * inv_c2;
                const amrex::Real S1 = r1 - p1 * inv_c2;
                const amrex::Real S2 = r2 - p2 * inv_c2;
                dS_abs =
                  std::abs(pc_nscbc::minmod(S0 - S1, S1 - S2)) / rho;
              }
            }

            if (rp[ridx + reg::init] == 0.0) {
              rp[ridx + reg::ema_Rout] = R_out;
              rp[ridx + reg::ema_p] = p;
              rp[ridx + reg::trend_dS] = dS_abs;
              rp[ridx + reg::ema_rhoc] = rhoc;
              rp[ridx + reg::init] = 1.0;
            } else {
              rp[ridx + reg::ema_Rout] += w * (R_out - rp[ridx + reg::ema_Rout]);
              rp[ridx + reg::ema_p] += w * (p - rp[ridx + reg::ema_p]);
              rp[ridx + reg::trend_dS] +=
                w * (dS_abs - rp[ridx + reg::trend_dS]);
              rp[ridx + reg::ema_rhoc] += w * (rhoc - rp[ridx + reg::ema_rhoc]);
            }
            rp[ridx + reg::u_minus] =
              n_sgn * 0.5 * (R_out - rp[ridx + reg::ema_Rout]);
            rp[ridx + reg::dp_ac] = p - rp[ridx + reg::ema_p];
          });
      }
    }
  }
  amrex::Gpu::streamSynchronize();
}

// Checkpoint the registers: gather (sum -- the zero-elsewhere invariant
// makes it exact) and write one small file from the IOProcessor.
void
PeleC::nscbc_registers_checkpoint(const std::string& dir)
{
  if (!bc_nscbc || !(bc_nscbc_nri || bc_nscbc_ndnr)) {
    return;
  }
  if ((nscbc_reg_p == nullptr) || nscbc_reg_p->empty()) {
    return; // no update has run yet; nothing to checkpoint
  }
  auto& regs = *nscbc_reg_p;
  std::vector<amrex::Real> h(regs.size());
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, regs.begin(), regs.end(), h.begin());
  amrex::ParallelDescriptor::ReduceRealSum(
    h.data(), static_cast<int>(h.size()));
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream f(dir + "/NSCBCRegisters", std::ios::binary);
    const long n = static_cast<long>(h.size());
    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
    f.write(
      reinterpret_cast<const char*>(h.data()),
      static_cast<std::streamsize>(n * sizeof(amrex::Real)));
  }
}

// Restore the registers: read, broadcast, then zero the entries this rank
// does not own so the invariant (and the next checkpoint's gather) holds.
// A checkpoint without the file (pre-register, or registers then disabled)
// restarts with cold registers, which re-seed on the first advance.
void
PeleC::nscbc_registers_restart(const std::string& dir)
{
  if (!bc_nscbc || !(bc_nscbc_nri || bc_nscbc_ndnr) || (level != 0)) {
    return;
  }
  auto& regs = nscbc_registers(geom.Domain());
  std::vector<amrex::Real> h(regs.size(), 0.0);
  int have = 0;
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ifstream f(dir + "/NSCBCRegisters", std::ios::binary);
    if (f.good()) {
      long n = 0;
      f.read(reinterpret_cast<char*>(&n), sizeof(n));
      if (n == static_cast<long>(h.size())) {
        f.read(
          reinterpret_cast<char*>(h.data()),
          static_cast<std::streamsize>(n * sizeof(amrex::Real)));
        have = 1;
      }
    }
  }
  amrex::ParallelDescriptor::Bcast(
    &have, 1, amrex::ParallelDescriptor::IOProcessorNumber());
  if (have == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "  NSCBC: no register file in the checkpoint; "
                        "registers re-seed on the first advance.\n";
    }
    return;
  }
  amrex::ParallelDescriptor::Bcast(
    h.data(), static_cast<int>(h.size()),
    amrex::ParallelDescriptor::IOProcessorNumber());

  // Zero what this rank does not own.
  const amrex::Box dom = geom.Domain();
  std::vector<char> owned(h.size() / reg::ncomp, 0);
  for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    for (int side = 0; side < 2; ++side) {
      amrex::Box band = dom;
      if (side == 0) {
        band.setBig(idir, dom.smallEnd(idir));
      } else {
        band.setSmall(idir, dom.bigEnd(idir));
      }
      for (int bi = 0; bi < grids.size(); ++bi) {
        if (dmap[bi] != amrex::ParallelDescriptor::MyProc()) {
          continue;
        }
        const amrex::Box b = grids[bi] & band;
        if (!b.ok()) {
          continue;
        }
        amrex::LoopOnCpu(b, [&](int i, int j, int k) {
          const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
          owned[nscbc_reg_index(idir, side, iv, dom) / reg::ncomp] = 1;
        });
      }
    }
  }
  for (size_t e = 0; e < owned.size(); ++e) {
    if (owned[e] == 0) {
      for (int n = 0; n < reg::ncomp; ++n) {
        h[e * reg::ncomp + n] = 0.0;
      }
    }
  }
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, h.begin(), h.end(), regs.begin());
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
    h[pc_nscbc::Diag::reversed] +
    h[pc_nscbc::Diag::body_state] + h[pc_nscbc::Diag::eos_failure] +
    h[pc_nscbc::Diag::floored] + h[pc_nscbc::Diag::transverse_drop] +
    h[pc_nscbc::Diag::source_drop] + h[pc_nscbc::Diag::target_incomplete];
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
                   << h[pc_nscbc::Diag::transverse_drop]
                   << ",  source dropped "
                   << h[pc_nscbc::Diag::source_drop]
                   << ",  target incomplete "
                   << h[pc_nscbc::Diag::target_incomplete] << "\n";
    if (structure > 0) {
      // Advisory, not a fallback: a front is in the outflow boundary cells,
      // which is the configuration the flame closures exist for.
      amrex::Print()
        << "  NSCBC: material structure (|dS| > 5% of rho per cell) sat in "
        << structure
        << " outflow boundary-cell fills since the last report.\n"
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
