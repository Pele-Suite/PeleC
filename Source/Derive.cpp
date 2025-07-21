#include "mechanism.H"

#include "PeleCUtilities.H"
#include "Derive.H"
#include "PeleC.H"
#include "IndexDefines.H"
#include "TransCoeff.H"

void
pc_dermagvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto magvel = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real datinv = 1.0 / dat(i, j, k, URHO);
    const amrex::Real dat1 = (dat(i, j, k, UMX) * datinv);
    const amrex::Real dat2 = (dat(i, j, k, UMY) * datinv);
    const amrex::Real dat3 = (dat(i, j, k, UMZ) * datinv);
    magvel(i, j, k) = sqrt((dat1 * dat1) + (dat2 * dat2) + (dat3 * dat3));
  });
}

void
PeleC::pc_dermagvel_if(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
#if AMREX_SPACEDIM > 2
  const amrex::Geometry& geomdata,
#else
  const amrex::Geometry& /*geomdata*/,
#endif
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto magvel = derfab.array();

  if (do_rf) {
#if AMREX_SPACEDIM > 2
    int axis = rf_axis;
    amrex::Real omega = rf_omega;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};
    auto const gdata = geomdata.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::IntVect iv(AMREX_D_DECL(i, j, k));
      amrex::RealVect r = get_rotaxis_vec(iv, axis_loc, gdata);
      amrex::RealVect w(0.0, 0.0, 0.0);
      w[axis] = omega;
      amrex::RealVect w_cross_r = w.crossProduct(r);

      const amrex::Real datinv = 1.0 / dat(i, j, k, URHO);
      const amrex::Real dat1 = dat(i, j, k, UMX) * datinv + w_cross_r[0];
      const amrex::Real dat2 = dat(i, j, k, UMY) * datinv + w_cross_r[1];
      const amrex::Real dat3 = dat(i, j, k, UMZ) * datinv + w_cross_r[2];
      magvel(i, j, k) = sqrt((dat1 * dat1) + (dat2 * dat2) + (dat3 * dat3));
    });
#endif
  } else {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real datinv = 1.0 / dat(i, j, k, URHO);
      const amrex::Real dat1 = dat(i, j, k, UMX) * datinv;
      const amrex::Real dat2 = dat(i, j, k, UMY) * datinv;
      const amrex::Real dat3 = dat(i, j, k, UMZ) * datinv;
      magvel(i, j, k) = sqrt((dat1 * dat1) + (dat2 * dat2) + (dat3 * dat3));
    });
  }
}

void
pc_dermagmom(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto magmom = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    magmom(i, j, k) = sqrt(
      dat(i, j, k, UMX) * dat(i, j, k, UMX) +
      dat(i, j, k, UMY) * dat(i, j, k, UMY) +
      dat(i, j, k, UMZ) * dat(i, j, k, UMZ));
  });
}

void
pc_derkineng(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto kineng = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real datxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX);
    const amrex::Real datysq = dat(i, j, k, UMY) * dat(i, j, k, UMY);
    const amrex::Real datzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ);
    kineng(i, j, k) = 0.5 / dat(i, j, k, URHO) * (datxsq + datysq + datzsq);
  });
}

void
PeleC::pc_dereint1(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  //  Compute internal energy from (rho E).
  auto const dat = datfab.const_array();
  auto e = derfab.array();
  if (do_rf) {
#if AMREX_SPACEDIM > 2
    int axis = rf_axis;
    amrex::Real omega = rf_omega;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};
    auto const gdata = geomdata.data();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::IntVect iv(AMREX_D_DECL(i, j, k));
      amrex::Real roteng = get_rot_energy(iv, omega, axis, axis_loc, gdata);
      const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
      const amrex::Real ux = dat(i, j, k, UMX) * rhoInv;
      const amrex::Real uy = dat(i, j, k, UMY) * rhoInv;
      const amrex::Real uz = dat(i, j, k, UMZ) * rhoInv;
      e(i, j, k) =
        dat(i, j, k, UEDEN) * rhoInv - 0.5 * (ux * ux + uy * uy + uz * uz);
      // see Blazek, Appendix A3, Navier-Stokes in rotating frame of reference
      e(i, j, k) -= roteng;
    });
#else
    amrex::ignore_unused(geomdata);
#endif
  } else {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
      const amrex::Real ux = dat(i, j, k, UMX) * rhoInv;
      const amrex::Real uy = dat(i, j, k, UMY) * rhoInv;
      const amrex::Real uz = dat(i, j, k, UMZ) * rhoInv;
      e(i, j, k) =
        dat(i, j, k, UEDEN) * rhoInv - 0.5 * (ux * ux + uy * uy + uz * uz);
    });
  }
}

void
pc_derlogden(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto logden = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    logden(i, j, k) = log10(dat(i, j, k));
  });
}

void
pc_dermagvort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto vort = derfab.array();

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3, amrex::The_Async_Arena());
  auto larr = local.array();

  const auto& flag_fab = amrex::getEBCellFlagFab(datfab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;

  // Convert momentum to velocity.
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    AMREX_D_TERM(larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
                 , larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
                 , larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoInv;)
  });

  AMREX_D_TERM(
    const amrex::Real dx = geomdata.CellSize(0);
    , const amrex::Real dy = geomdata.CellSize(1);
    , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate vorticity.
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)

    // if fab is all regular -> call regular idx and weights
    // otherwise
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(
      const amrex::Real wi = get_weight(im, ip);
      , const amrex::Real wj = get_weight(jm, jp);
      , const amrex::Real wk = get_weight(km, kp);)

    AMREX_D_TERM(
      vort(i, j, k) = 0.0 * dx;
      ,
      const amrex::Real vx = wi * (larr(ip, j, k, 1) - larr(im, j, k, 1)) / dx;
      const amrex::Real uy = wj * (larr(i, jp, k, 0) - larr(i, jm, k, 0)) / dy;
      const amrex::Real v3 = vx - uy;
      ,
      const amrex::Real wx = wi * (larr(ip, j, k, 2) - larr(im, j, k, 2)) / dx;
      const amrex::Real wy = wj * (larr(i, jp, k, 2) - larr(i, jm, k, 2)) / dy;
      const amrex::Real uz = wk * (larr(i, j, kp, 0) - larr(i, j, km, 0)) / dz;
      const amrex::Real vz = wk * (larr(i, j, kp, 1) - larr(i, j, km, 1)) / dz;
      const amrex::Real v1 = wy - vz; const amrex::Real v2 = uz - wx;);
    vort(i, j, k) = sqrt(AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void
pc_derdivu(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto divu = derfab.array();

  const auto& flag_fab = amrex::getEBCellFlagFab(datfab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;

  AMREX_D_TERM(
    const amrex::Real dx = geomdata.CellSize(0);
    , const amrex::Real dy = geomdata.CellSize(1);
    , const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(
      const amrex::Real wi = get_weight(im, ip);
      , const amrex::Real wj = get_weight(jm, jp);
      , const amrex::Real wk = get_weight(km, kp);)

    AMREX_D_TERM(
      const amrex::Real uhi = dat(ip, j, k, UMX) / dat(ip, j, k, URHO);
      const amrex::Real ulo = dat(im, j, k, UMX) / dat(im, j, k, URHO);
      , const amrex::Real vhi = dat(i, jp, k, UMY) / dat(i, jp, k, URHO);
      const amrex::Real vlo = dat(i, jm, k, UMY) / dat(i, jm, k, URHO);
      , const amrex::Real whi = dat(i, j, kp, UMZ) / dat(i, j, kp, URHO);
      const amrex::Real wlo = dat(i, j, km, UMZ) / dat(i, j, km, URHO););
    divu(i, j, k) = AMREX_D_TERM(
      wi * (uhi - ulo) / dx, +wj * (vhi - vlo) / dy, +wk * (whi - wlo) / dz);
  });
}

void
pc_derenstrophy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  // This routine will derive enstrophy  = 1/2 rho (x_vorticity^2 +
  // y_vorticity^2 + z_vorticity^2)
  auto const dat = datfab.const_array();
  auto enstrophy = derfab.array();

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3, amrex::The_Async_Arena());
  auto larr = local.array();

  const auto& flag_fab = amrex::getEBCellFlagFab(datfab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;

  // Convert momentum to velocity.
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
    larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
    larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoInv;
  });

  AMREX_D_TERM(
    const amrex::Real dx = geomdata.CellSize(0);
    , const amrex::Real dy = geomdata.CellSize(1);
    , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate enstrophy.
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(
      const amrex::Real wi = get_weight(im, ip);
      , const amrex::Real wj = get_weight(jm, jp);
      , const amrex::Real wk = get_weight(km, kp);)

    AMREX_D_TERM(
      enstrophy(i, j, k) = 0.0 * dx;
      ,
      const amrex::Real vx = wi * (larr(ip, j, k, 1) - larr(im, j, k, 1)) / dx;
      const amrex::Real uy = wj * (larr(i, jp, k, 0) - larr(i, jm, k, 0)) / dy;
      const amrex::Real v3 = vx - uy;
      ,
      const amrex::Real wx = wi * (larr(ip, j, k, 2) - larr(im, j, k, 2)) / dx;

      const amrex::Real wy = wj * (larr(i, jp, k, 2) - larr(i, jm, k, 2)) / dy;

      const amrex::Real uz = wk * (larr(i, j, kp, 0) - larr(i, j, km, 0)) / dz;
      const amrex::Real vz = wk * (larr(i, j, kp, 1) - larr(i, j, km, 1)) / dz;

      const amrex::Real v1 = wy - vz; const amrex::Real v2 = uz - wx;);
    enstrophy(i, j, k) = 0.5 * dat(i, j, k, URHO) *
                         (AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void
pc_dernull(
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // This routine is used by particle_count.  Yes it does nothing.
}

void
pc_dermolefrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // Derive the mole fractions of the species
  auto const dat = datfab.const_array();
  auto spec = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real mass[NUM_SPECIES];
    amrex::Real mole[NUM_SPECIES];
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);

    for (int n = 0; n < NUM_SPECIES; n++) {
      mass[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.Y2X(mass, mole);
    for (int n = 0; n < NUM_SPECIES; n++) {
      spec(i, j, k, n) = mole[n];
    }
  });
}

void
pc_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto cfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real c;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2Cs(rho, T, massfrac, c);
    cfab(i, j, k) = c;
  });
}

void
pc_derentropy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto sfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real s;
    auto eos = pele::physics::PhysicsType::eos();
    eos.S(s);
    sfab(i, j, k) = s;
  });
}

void
pc_dermachnumber(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto mach = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real c;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2Cs(rho, T, massfrac, c);
    const amrex::Real datxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX);
    const amrex::Real datysq = dat(i, j, k, UMY) * dat(i, j, k, UMY);
    const amrex::Real datzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ);
    mach(i, j, k) = sqrt(datxsq + datysq + datzsq) / dat(i, j, k, URHO) / c;
  });
}

void
pc_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto pfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real p;
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2P(rho, T, massfrac, p);
    pfab(i, j, k) = p;
  });
}

void
pc_dertemp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto tfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    tfab(i, j, k) = dat(i, j, k, UTEMP);
  });
}

void
pc_derspectrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/,
  const int idx)
{
  auto const dat = datfab.const_array();
  auto spectrac = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    spectrac(i, j, k) = dat(i, j, k, UFS + idx) / dat(i, j, k, URHO);
  });
}

void
pc_derradialvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto rvel = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
    geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();
  AMREX_D_TERM(
    const amrex::Real centerx = 0.5 * (prob_lo[0] + prob_hi[0]);
    , const amrex::Real centery = 0.5 * (prob_lo[1] + prob_hi[1]);
    , const amrex::Real centerz = 0.5 * (prob_lo[2] + prob_hi[2]));

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0] - centerx;
      , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1] - centery;
      , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2] - centerz;)
    const amrex::Real r = sqrt(AMREX_D_TERM(x * x, +y * y, +z * z));
    rvel(i, j, k) = (AMREX_D_TERM(
                      dat(i, j, k, UMX) * x, +dat(i, j, k, UMY) * y,
                      +dat(i, j, k, UMZ) * z)) /
                    (dat(i, j, k, URHO) * r);
  });
}

void
pc_dercp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto cp_arr = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real mass[NUM_SPECIES];
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);

    for (int n = 0; n < NUM_SPECIES; n++) {
      mass[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real cp = 0.0;
    eos.RTY2Cp(dat(i, j, k, URHO), dat(i, j, k, UTEMP), mass, cp);
    cp_arr(i, j, k) = cp;
  });
}

void
pc_dercv(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto cv_arr = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real mass[NUM_SPECIES];
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);

    for (int n = 0; n < NUM_SPECIES; n++) {
      mass[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real cv = 0.0;
    eos.RTY2Cv(dat(i, j, k, URHO), dat(i, j, k, UTEMP), mass, cv);
    cv_arr(i, j, k) = cv;
  });
}

void
pc_dercoord(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto coord_arr = derfab.array();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(coord_arr(i, j, k, 0) = prob_lo[0] + (i + 0.5) * dx[0];
                 , coord_arr(i, j, k, 1) = prob_lo[1] + (j + 0.5) * dx[1];
                 , coord_arr(i, j, k, 2) = prob_lo[2] + (k + 0.5) * dx[2];);
  });
}

void
PeleC::pc_derviscosity(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto mu_arr = derfab.array();
  auto const* ltransparm = trans_parms.device_parm();
  const ProbParmDevice* lprobparm = d_prob_parm_device;
  const auto& gdata = geomdata.data();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real massfrac[NUM_SPECIES];
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    amrex::Real mu = 0.0, dum1 = 0.0, dum2 = 0.0;
    const bool get_xi = false, get_mu = true, get_lam = false,
               get_Ddiag = false, get_chi = false;
    const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, gdata);
    pc_transcoeff(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac, nullptr,
      nullptr, mu, dum1, dum2, ltransparm, *lprobparm, x);
    mu_arr(i, j, k) = mu;
  });
}

void
PeleC::pc_derbulkviscosity(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto xi_arr = derfab.array();
  auto const* ltransparm = trans_parms.device_parm();
  const ProbParmDevice* lprobparm = d_prob_parm_device;
  const auto& gdata = geomdata.data();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real massfrac[NUM_SPECIES];
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    amrex::Real xi = 0.0, dum1 = 0.0, dum2 = 0.0;
    const bool get_xi = true, get_mu = false, get_lam = false,
               get_Ddiag = false, get_chi = false;
    const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, gdata);
    pc_transcoeff(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac, nullptr,
      nullptr, dum1, xi, dum2, ltransparm, *lprobparm, x);
    xi_arr(i, j, k) = xi;
  });
}

void
PeleC::pc_derconductivity(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto lam_arr = derfab.array();
  auto const* ltransparm = trans_parms.device_parm();
  const ProbParmDevice* lprobparm = d_prob_parm_device;
  const auto& gdata = geomdata.data();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real massfrac[NUM_SPECIES];
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    amrex::Real lam = 0.0, dum1 = 0.0, dum2 = 0.0;
    const bool get_xi = false, get_mu = false, get_lam = true,
               get_Ddiag = false, get_chi = false;
    const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, gdata);
    pc_transcoeff(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac, nullptr,
      nullptr, dum1, dum2, lam, ltransparm, *lprobparm, x);
    lam_arr(i, j, k) = lam;
  });
}

void
PeleC::pc_derdiffusivity(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto d_arr = derfab.array();
  auto const* ltransparm = trans_parms.device_parm();
  const ProbParmDevice* lprobparm = d_prob_parm_device;
  const auto& gdata = geomdata.data();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real ddiag[NUM_SPECIES] = {0.0};
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    amrex::Real dum1 = 0.0, dum2 = 0.0, dum3 = 0.0;
    const bool get_xi = false, get_mu = false, get_lam = false,
               get_Ddiag = true, get_chi = false;
    const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, gdata);
    pc_transcoeff(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac, ddiag,
      nullptr, dum1, dum2, dum3, ltransparm, *lprobparm, x);
    for (int n = 0; n < NUM_SPECIES; n++) {
      d_arr(i, j, k, n) = ddiag[n];
    }
  });
}

#ifdef PELE_USE_MASA
void
pc_derrhommserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto rhommserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real rho = masa_eval_3d_exact_rho(x, y, z);
    rhommserror(i, j, k) = dat(i, j, k, URHO) - rho;
  });
}

void
pc_derummserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto ummserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real u = masa_eval_3d_exact_u(x, y, z);
    ummserror(i, j, k) = dat(i, j, k, UMX) / dat(i, j, k, URHO) - u;
  });
}

void
pc_dervmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto vmmserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real v = masa_eval_3d_exact_v(x, y, z);
    vmmserror(i, j, k) = dat(i, j, k, UMY) / dat(i, j, k, URHO) - v;
  });
}

void
pc_derwmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto wmmserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real w = masa_eval_3d_exact_w(x, y, z);
    wmmserror(i, j, k) = dat(i, j, k, UMZ) / dat(i, j, k, URHO) - w;
  });
}

void
pc_derpmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.const_array();
  auto pmmserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    auto eos = pele::physics::PhysicsType::eos();
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    amrex::Real eint = dat(i, j, k, UEINT) * rhoinv;
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }

    amrex::Real pdat;
    eos.RYET2P(rho, massfrac, eint, T, pdat);
    const amrex::Real p = masa_eval_3d_exact_p(x, y, z);
    pmmserror(i, j, k) = pdat - p;
  });
}
#endif
