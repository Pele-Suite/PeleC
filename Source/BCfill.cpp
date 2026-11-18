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

  AMREX_GPU_HOST
  explicit PCHypFillExtDir(
    const ProbParmDevice* d_prob_parm,
    const bool do_turb_inflow,
    const bool nscbc,
    const amrex::GpuArray<pc_nscbc::Params, AMREX_SPACEDIM>& nscbc_prm,
    amrex::Long* nscbc_diag,
    amrex::GpuArray<amrex::Array4<const amrex::Real>, 2 * AMREX_SPACEDIM> const&
      turb = {},
    const bool has_turb = false)
    : lprobparm(d_prob_parm),
      m_do_turb_inflow(do_turb_inflow),
      m_nscbc(nscbc),
      m_nscbc_prm(nscbc_prm),
      m_nscbc_diag(nscbc_diag),
      m_turb(turb),
      m_has_turb(has_turb)
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
  //  and whose boundary cell holds fluid (not EB body state); the state is
  //  written exactly once.  Combined with the clamped stencil
  //  below this makes the fill a pure function of valid interior data, hence
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
      // the fill is a pure function of pre-launch data.  Next to a WALL in a
      // tangential direction this reads the domain's boundary row, not the
      // wall's reflected ghost (which the legacy path would have used), so a
      // corner ghost owned by this face carries that row's tangential
      // velocity rather than the wall-flipped one.
      //
      // In a PERIODIC tangential direction this leaves a small seam
      // residual, and that is a deliberate trade.  amrex's corner protocol
      // (StateDataPhysBCFunct) recomputes corner ghosts on a strip FAB
      // holding only their image band, so a seam-adjacent cell is filled
      // twice on two different FABs; exact agreement requires restricting
      // the tangential stencil to the band the strip can see.  Both
      // alternatives were built and measured: a wrap through the resident
      // images leaves the array aperiodic at the 2e-2 level, and the
      // strip-consistent band is bitwise-periodic and equally stable -- but
      // costs ~100 lines of decomposition-sensitive index machinery to
      // remove a residual that measures 2e-4 (inert) to 1.6e-3 (a flame on
      // the seam corner).  The clamp was kept for simplicity.  Nothing in
      // this code checks the residual; the verification cases of a later PR
      // in this series measure it.
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
      for (int n = 0; n < NVAR; n++) {
        s_N[n] = dest(ivN, n);
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

      // TurbInflow fluctuations at this boundary point, as bcnormal gets
      // them; the hook decides how they enter its Target (added to the mean,
      // or as the total velocity).
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> turb_fluc{0.0};
      const auto& turb = m_turb[2 * idir + ((sgn > 0) ? 0 : 1)];
      if (m_has_turb && turb) {
        amrex::IntVect ivt = base;
        ivt[idir] = N_pos - sgn; // the layer-1 ghost, where add_turb wrote
        // Across a periodic tangential direction the turbulence is periodic
        // too: a seam ghost reads its image, not the clamped boundary row
        // (which would put an O(u') seam into the inflow), whenever this FAB
        // holds the image.
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          if ((d != idir) && (geom.isPeriodic(d) != 0)) {
            const int n_d = domhi[d] - domlo[d] + 1;
            const int w = domlo[d] + (((iv[d] - domlo[d]) % n_d) + n_d) % n_d;
            if ((w >= fab_lo[d]) && (w <= fab_hi[d])) {
              ivt[d] = w;
            }
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

      // The rest of the stencil, now that this point is characteristic.
      const amrex::IntVect ivNm1 = stencil_iv(1);
      const amrex::IntVect ivNm2 = stencil_iv(2);
      for (int n = 0; n < NVAR; n++) {
        s_Nm1[n] = dest(ivNm1, n);
        s_Nm2[n] = dest(ivNm2, n);
      }
      // A FAB holding fewer interior cells along the normal than the
      // extrapolation wants (the coarse patch of a two-level fill under a
      // fine grid on this face, say) drops the order to what fits.  Counted,
      // because the ghost then depends on the box layout.
      if (
        (m_nscbc_diag != nullptr) &&
        (depth_fab < amrex::min<int>(3, depth_domain))) {
        amrex::HostDevice::Atomic::Add(
          &m_nscbc_diag[pc_nscbc::Diag::stencil_short], amrex::Long(1));
      }

      amrex::Real s_ghost[NVAR];
      if (!pc_nscbc::apply(
            s_N, s_Nm1, s_Nm2, n_stencil, dx[idir], idir, sgn, layer, tgt,
            m_nscbc_prm[idir], s_ghost, m_nscbc_diag)) {
        // The boundary cell of the stencil is EB body state: there is
        // nothing to build a characteristic state from, so this face does
        // not own the cell.  Another direction may; otherwise the ordinary
        // bcnormal() path below fills it exactly as it would without the
        // characteristic treatment.
        continue;
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
      turb_arr, has_turb});
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
          "(see Source/NSCBC.H). (This warning is printed once per face.)");
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
  p.pin_farfield = bc_nscbc_pin_farfield;
  // Only the ratio sigma/L_ref is physical.  L_ref is fixed to the domain
  // extent along the boundary normal so that sigma keeps the meaning it has
  // in the literature, rather than being one of two dials for one degree of
  // freedom.  Uses probhi - problo, not probhi: the legacy Fortran used
  // probhi(idir) and was therefore silently wrong for any domain not
  // anchored at the origin.
  const auto& geom = amrex::DefaultGeometry();
  p.L_ref = geom.ProbHi(idir) - geom.ProbLo(idir);
  return p;
}

void
PeleC::nscbc_init_diagnostics()
{
  if (nscbc_diag_p != nullptr) {
    return;
  }
  nscbc_diag_p =
    new amrex::Gpu::DeviceVector<amrex::Long>(pc_nscbc::Diag::count, 0);
  amrex::ExecOnFinalize([]() {
    delete nscbc_diag_p;
    nscbc_diag_p = nullptr;
  });
}

void
PeleC::nscbc_report_diagnostics()
{
  if (!bc_nscbc) {
    return;
  }
  // Settle any counter atomics still in flight BEFORE the copy: the fills
  // bump these counters from MFIter's rotating streams, and the blocking
  // Gpu::copy below synchronises only its own stream.  A count landing
  // between the copy and the reset would never be reported.
  amrex::Gpu::Device::synchronize();
  std::vector<amrex::Long> h(pc_nscbc::Diag::count, 0);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, nscbc_diag().begin(), nscbc_diag().end(),
    h.begin());
  amrex::ParallelDescriptor::ReduceLongSum(h.data(), pc_nscbc::Diag::count);

  // The supersonic path is exact, not a degradation, so it is reported but is
  // not a warning.  The others mean the boundary is being asked for something
  // it cannot cleanly provide.
  const amrex::Long total =
    h[pc_nscbc::Diag::reversed] + h[pc_nscbc::Diag::body_state] +
    h[pc_nscbc::Diag::eos_failure] + h[pc_nscbc::Diag::floored] +
    h[pc_nscbc::Diag::stencil_short] + h[pc_nscbc::Diag::target_incomplete] +
    h[pc_nscbc::Diag::passive_unset];
  if (amrex::ParallelDescriptor::IOProcessor() && (total > 0 || verbose > 1)) {
    amrex::Print() << "  NSCBC fallbacks since last report:" << "  supersonic "
                   << h[pc_nscbc::Diag::supersonic] << ",  flow reversal "
                   << h[pc_nscbc::Diag::reversed] << ",  EB body state "
                   << h[pc_nscbc::Diag::body_state] << ",  EOS failure "
                   << h[pc_nscbc::Diag::eos_failure] << ",  floored "
                   << h[pc_nscbc::Diag::floored] << ",  short stencil "
                   << h[pc_nscbc::Diag::stencil_short]
                   << ",  target incomplete "
                   << h[pc_nscbc::Diag::target_incomplete]
                   << ",  inflow passives unset "
                   << h[pc_nscbc::Diag::passive_unset] << "\n";
  }
  nscbc_diag().assign(pc_nscbc::Diag::count, 0);
  // The reset itself is a device fill on the current stream; the next
  // advance's fills bump these counters from MFIter's rotating streams,
  // which are not ordered against it.  Settle the reset before returning
  // so a zero can never land on top of a fresh count.
  amrex::Gpu::streamSynchronize();
}
