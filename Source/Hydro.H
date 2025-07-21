#ifndef HYDRO_H
#define HYDRO_H

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "PeleCUtilities.H"
#include "Godunov.H"
#include "AMReX_EB_Redistribution.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_srctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& src,
  amrex::Array4<amrex::Real> const& srcq)
{
  // Note (Hari S):
  // rotational frame not incorporated here, but shouldn't be called because
  // RF required MOL and this is only used for SDC/Godunov

  // Local loading Vars ==> Load once from global memory use multiple times
  const amrex::Real srcrho = src(i, j, k, URHO);
  const amrex::Real rhoinv = 1.0 / q(i, j, k, QRHO);
  const amrex::Real vx = q(i, j, k, QU);
  const amrex::Real vy = q(i, j, k, QV);
  const amrex::Real vz = q(i, j, k, QW);
  const amrex::Real umx = src(i, j, k, UMX);
  const amrex::Real umy = src(i, j, k, UMY);
  const amrex::Real umz = src(i, j, k, UMZ);
  // SrcQ Calc
  srcq(i, j, k, QRHO) = srcrho;
  srcq(i, j, k, QU) = (umx - vx * srcrho) * rhoinv;
  srcq(i, j, k, QV) = (umy - vy * srcrho) * rhoinv;
  srcq(i, j, k, QW) = (umz - vz * srcrho) * rhoinv;
  srcq(i, j, k, QREINT) = src(i, j, k, UEDEN) - vx * umx - vy * umy - vz * umz +
                          srcrho * (vx * vx + vy * vy + vz * vz);
  srcq(i, j, k, QPRES) =
    qa(i, j, k, QDPDE) *
      (srcq(i, j, k, QREINT) - q(i, j, k, QREINT) * srcrho * rhoinv) * rhoinv +
    qa(i, j, k, QDPDR) * srcrho;

#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    srcq(i, j, k, QFA + n) =
      (src(i, j, k, UFA + n) - q(i, j, k, QFA + n) * srcrho) * rhoinv;
  }
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    srcq(i, j, k, QFS + n) =
      (src(i, j, k, UFS + n) - q(i, j, k, QFS + n) * srcrho) * rhoinv;
  }
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    srcq(i, j, k, QFX + n) = src(i, j, k, UFX + n);
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    srcq(i, j, k, QLIN + n) = src(i, j, k, ULIN + n);
  }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pc_compute_hyp_wallflux(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real const> const& q,
  amrex::Array4<amrex::Real const> const& qa,
  const amrex::Real axm,
  const amrex::Real axp,
  const amrex::Real aym,
  const amrex::Real ayp,
#if (AMREX_SPACEDIM == 3)
  const amrex::Real azm,
  const amrex::Real azp,
#endif
  amrex::GpuArray<amrex::Real, NVAR>& fluxw) noexcept
{

  // Here we have 2 options. Option 1: call the Riemann solver with dummy data.
  // Option 2 is set the known flux at the wall, as we do in MOL. We are going
  // with option 2.
  // clang-format off
  // Option 1
  // AMREX_D_TERM(amrex::Real u = q(i, j, k, QU);, amrex::Real v = q(i, j, k, QV);
  //              , amrex::Real w = q(i, j, k, QW););
  // const amrex::Real apnorm = std::sqrt(AMREX_D_TERM(
  //   (axm - axp) * (axm - axp), +(aym - ayp) * (aym - ayp),
  //   +(azm - azp) * (azm - azp)));
  // const amrex::Real apnorminv = 1. / apnorm;
  // amrex::Real un = AMREX_D_TERM(
  //   u * (axm - axp) * apnorminv, +v * (aym - ayp) * apnorminv,
  //   +w * (azm - azp) * apnorminv);

  // amrex::Real momfluxn, momfluxt1, momfluxt2;
  // amrex::Real dummy_gu, dummy_gv, dummy_gv2, dummy_gp, dummy_gg;
  // amrex::Real rho = q(i, j, k, QRHO);
  // amrex::Real p = q(i, j, k, QPRES);
  // amrex::Real cav = qa(i, j, k, QC);
  // amrex::Real ustar;
  // amrex::Real sp[NUM_SPECIES] = {0.0};
  // for (int n = 0; n < NUM_SPECIES; ++n) {
  //   sp[n] = q(i, j, k, QFS + n);
  // }

  // const int bc_test_val = 1;
  // riemann(
  //   rho, un, 0., 0., p, sp, rho, un, 0., 0., p, sp, bc_test_val, cav, ustar,
  //   fluxw[URHO], &fluxw[UFS], momfluxn, momfluxt1, momfluxt2, fluxw[UEDEN],
  //   fluxw[UEINT], dummy_gu, dummy_gv, dummy_gv2, dummy_gp, dummy_gg);

  // // enforce that no spurious flux was created by the Riemann solver
  // for (int n = 0; n < NVAR; n++) {
  //   fluxw[n] = 0.;
  // }
  // AMREX_D_TERM(fluxw[UMX] = (axm - axp) * momfluxn;
  //              , fluxw[UMY] = (aym - ayp) * momfluxn;
  //              , fluxw[UMZ] = (azm - azp) * momfluxn;);
  // clang-format on

  // Option 2
  amrex::ignore_unused(qa);
  amrex::Real p = q(i, j, k, QPRES);
  AMREX_D_TERM(
    fluxw[UMX] = (axm - axp) * p;, fluxw[UMY] = (aym - ayp) * p;
    , fluxw[UMZ] = (azm - azp) * p;);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
eb_compute_div(
  int i,
  int j,
  int k,
  int n,
  amrex::IntVect const& blo,
  amrex::IntVect const& bhi,
  amrex::Array4<amrex::Real const> const& q,
  amrex::Array4<amrex::Real const> const& qaux,
  amrex::Array4<amrex::Real> const& divu,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& u,
    amrex::Array4<amrex::Real const> const& v,
    amrex::Array4<amrex::Real const> const& w),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real> const& fx,
    amrex::Array4<amrex::Real> const& fy,
    amrex::Array4<amrex::Real> const& fz),
  amrex::Array4<amrex::EBCellFlag const> const& flag,
  amrex::Array4<amrex::Real const> const& vfrc,
  amrex::Array4<amrex::Real> const& redistwgt,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& apx,
    amrex::Array4<amrex::Real const> const& apy,
    amrex::Array4<amrex::Real const> const& apz),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& fcx,
    amrex::Array4<amrex::Real const> const& fcy,
    amrex::Array4<amrex::Real const> const& fcz),
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
  const int eb_weights_type)
{
  AMREX_D_TERM(
    bool x_high = (i == bhi[0]);, bool y_high = (j == bhi[1]);
    , bool z_high = (k == bhi[2]));
  bool valid_cell = AMREX_D_TERM(
    (blo[0] <= i) && (i <= bhi[0]), &&(blo[1] <= j) && (j <= bhi[1]),
    &&(blo[2] <= k) && (k <= bhi[2]));

#if (AMREX_SPACEDIM == 2)
  amrex::Real cell_vol_inv = dxinv[0] * dxinv[1];
  if (flag(i, j, k).isCovered()) {
    divu(i, j, k, n) = 0.0;
    if (valid_cell) {
      fx(i, j, k, n) = 0.;
      fy(i, j, k, n) = 0.;
      if (x_high) {
        fx(i + 1, j, k, n) = 0.;
      }
      if (y_high) {
        fy(i, j + 1, k, n) = 0.;
      }
    }
  } else if (flag(i, j, k).isRegular()) {
    divu(i, j, k, n) = cell_vol_inv * (u(i + 1, j, k, n) - u(i, j, k, n)) +
                       cell_vol_inv * (v(i, j + 1, k, n) - v(i, j, k, n));
    if (valid_cell) {
      fx(i, j, k, n) = u(i, j, k, n);
      fy(i, j, k, n) = v(i, j, k, n);
      if (x_high) {
        fx(i + 1, j, k, n) = u(i + 1, j, k, n);
      }
      if (y_high) {
        fy(i, j + 1, k, n) = v(i, j + 1, k, n);
      }
    }
  } else {
    amrex::Real fxm = u(i, j, k, n);
    if (apx(i, j, k) != 0.0 && apx(i, j, k) != 1.0) {
      int jj =
        j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i, j, k, 0)));
      amrex::Real fracy = flag(i, j, k).isConnected(0, jj - j, 0)
                            ? amrex::Math::abs(fcx(i, j, k, 0))
                            : 0.0;
      fxm = (1.0 - fracy) * fxm + fracy * u(i, jj, k, n);
    }
    if (valid_cell) {
      fx(i, j, k, n) = fxm;
    }

    amrex::Real fxp = u(i + 1, j, k, n);
    if (apx(i + 1, j, k) != 0.0 && apx(i + 1, j, k) != 1.0) {
      int jj =
        j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i + 1, j, k, 0)));
      amrex::Real fracy = flag(i + 1, j, k).isConnected(0, jj - j, 0)
                            ? amrex::Math::abs(fcx(i + 1, j, k, 0))
                            : 0.0;
      fxp = (1.0 - fracy) * fxp + fracy * u(i + 1, jj, k, n);
    }
    if (valid_cell && x_high) {
      fx(i + 1, j, k, n) = fxp;
    }

    amrex::Real fym = v(i, j, k, n);
    if (apy(i, j, k) != 0.0 && apy(i, j, k) != 1.0) {
      int ii =
        i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j, k, 0)));
      amrex::Real fracx = flag(i, j, k).isConnected(ii - i, 0, 0)
                            ? amrex::Math::abs(fcy(i, j, k, 0))
                            : 0.0;
      fym = (1.0 - fracx) * fym + fracx * v(ii, j, k, n);
    }
    if (valid_cell) {
      fy(i, j, k, n) = fym;
    }

    amrex::Real fyp = v(i, j + 1, k, n);
    if (apy(i, j + 1, k) != 0.0 && apy(i, j + 1, k) != 1.0) {
      int ii =
        i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j + 1, k, 0)));
      amrex::Real fracx = flag(i, j + 1, k).isConnected(ii - i, 0, 0)
                            ? amrex::Math::abs(fcy(i, j + 1, k, 0))
                            : 0.0;
      fyp = (1.0 - fracx) * fyp + fracx * v(ii, j + 1, k, n);
    }
    if (valid_cell && y_high) {
      fy(i, j + 1, k, n) = fyp;
    }

    divu(i, j, k, n) = (1.0 / vfrc(i, j, k)) * cell_vol_inv *
                       (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm +
                        apy(i, j + 1, k) * fyp - apy(i, j, k) * fym);

    amrex::GpuArray<amrex::Real, NVAR> flux_hyp_wall = {0.0};
    pc_compute_hyp_wallflux(
      i, j, k, q, qaux, apx(i, j, k), apx(i + 1, j, k), apy(i, j, k),
      apy(i, j + 1, k), flux_hyp_wall);

    // Here we assume dx == dy == dz
    divu(i, j, k, n) += flux_hyp_wall[n] * dxinv[0] / vfrc(i, j, k);
  }

#else // 3-d starts here

  amrex::Real cell_vol_inv = dxinv[0] * dxinv[1] * dxinv[2];
  if (flag(i, j, k).isCovered()) {
    divu(i, j, k, n) = 0.0;
    if (valid_cell) {
      fx(i, j, k, n) = 0.;
      fy(i, j, k, n) = 0.;
      fz(i, j, k, n) = 0.;
      if (x_high) {
        fx(i + 1, j, k, n) = 0.;
      }
      if (y_high) {
        fy(i, j + 1, k, n) = 0.;
      }
      if (z_high) {
        fz(i, j, k + 1, n) = 0.;
      }
    }
  } else if (flag(i, j, k).isRegular()) {
    divu(i, j, k, n) = cell_vol_inv * (u(i + 1, j, k, n) - u(i, j, k, n)) +
                       cell_vol_inv * (v(i, j + 1, k, n) - v(i, j, k, n)) +
                       cell_vol_inv * (w(i, j, k + 1, n) - w(i, j, k, n));
    if (valid_cell) {
      fx(i, j, k, n) = u(i, j, k, n);
      fy(i, j, k, n) = v(i, j, k, n);
      fz(i, j, k, n) = w(i, j, k, n);
      if (x_high) {
        fx(i + 1, j, k, n) = u(i + 1, j, k, n);
      }
      if (y_high) {
        fy(i, j + 1, k, n) = v(i, j + 1, k, n);
      }
      if (z_high) {
        fz(i, j, k + 1, n) = w(i, j, k + 1, n);
      }
    }
  } else {
    amrex::Real fxm = u(i, j, k, n);
    if (apx(i, j, k) != 0.0 && apx(i, j, k) != 1.0) {
      int jj =
        j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i, j, k, 0)));
      int kk =
        k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i, j, k, 1)));
      amrex::Real fracy = flag(i, j, k).isConnected(0, jj - j, 0)
                            ? amrex::Math::abs(fcx(i, j, k, 0))
                            : 0.0;
      amrex::Real fracz = flag(i, j, k).isConnected(0, 0, kk - k)
                            ? amrex::Math::abs(fcx(i, j, k, 1))
                            : 0.0;
      fxm = (1.0 - fracy) * (1.0 - fracz) * fxm +
            fracy * (1.0 - fracz) * u(i, jj, k, n) +
            fracz * (1.0 - fracy) * u(i, j, kk, n) +
            fracy * fracz * u(i, jj, kk, n);
    }
    if (valid_cell) {
      fx(i, j, k, n) = fxm;
    }

    amrex::Real fxp = u(i + 1, j, k, n);
    if (apx(i + 1, j, k) != 0.0 && apx(i + 1, j, k) != 1.0) {
      int jj =
        j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i + 1, j, k, 0)));
      int kk =
        k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i + 1, j, k, 1)));
      amrex::Real fracy = flag(i + 1, j, k).isConnected(0, jj - j, 0)
                            ? amrex::Math::abs(fcx(i + 1, j, k, 0))
                            : 0.0;
      amrex::Real fracz = flag(i + 1, j, k).isConnected(0, 0, kk - k)
                            ? amrex::Math::abs(fcx(i + 1, j, k, 1))
                            : 0.0;
      fxp = (1.0 - fracy) * (1.0 - fracz) * fxp +
            fracy * (1.0 - fracz) * u(i + 1, jj, k, n) +
            fracz * (1.0 - fracy) * u(i + 1, j, kk, n) +
            fracy * fracz * u(i + 1, jj, kk, n);
    }
    if (valid_cell && x_high) {
      fx(i + 1, j, k, n) = fxp;
    }

    amrex::Real fym = v(i, j, k, n);
    if (apy(i, j, k) != 0.0 && apy(i, j, k) != 1.0) {
      int ii =
        i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j, k, 0)));
      int kk =
        k + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j, k, 1)));
      amrex::Real fracx = flag(i, j, k).isConnected(ii - i, 0, 0)
                            ? amrex::Math::abs(fcy(i, j, k, 0))
                            : 0.0;
      amrex::Real fracz = flag(i, j, k).isConnected(0, 0, kk - k)
                            ? amrex::Math::abs(fcy(i, j, k, 1))
                            : 0.0;
      fym = (1.0 - fracx) * (1.0 - fracz) * fym +
            fracx * (1.0 - fracz) * v(ii, j, k, n) +
            fracz * (1.0 - fracx) * v(i, j, kk, n) +
            fracx * fracz * v(ii, j, kk, n);
    }
    if (valid_cell) {
      fy(i, j, k, n) = fym;
    }

    amrex::Real fyp = v(i, j + 1, k, n);
    if (apy(i, j + 1, k) != 0.0 && apy(i, j + 1, k) != 1.0) {
      int ii =
        i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j + 1, k, 0)));
      int kk =
        k + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j + 1, k, 1)));
      amrex::Real fracx = flag(i, j + 1, k).isConnected(ii - i, 0, 0)
                            ? amrex::Math::abs(fcy(i, j + 1, k, 0))
                            : 0.0;
      amrex::Real fracz = flag(i, j + 1, k).isConnected(0, 0, kk - k)
                            ? amrex::Math::abs(fcy(i, j + 1, k, 1))
                            : 0.0;
      fyp = (1.0 - fracx) * (1.0 - fracz) * fyp +
            fracx * (1.0 - fracz) * v(ii, j + 1, k, n) +
            fracz * (1.0 - fracx) * v(i, j + 1, kk, n) +
            fracx * fracz * v(ii, j + 1, kk, n);
    }
    if (valid_cell && y_high) {
      fy(i, j + 1, k, n) = fyp;
    }

    amrex::Real fzm = w(i, j, k, n);
    if (apz(i, j, k) != 0.0 && apz(i, j, k) != 1.0) {
      int ii =
        i + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k, 0)));
      int jj =
        j + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k, 1)));
      amrex::Real fracx = flag(i, j, k).isConnected(ii - i, 0, 0)
                            ? amrex::Math::abs(fcz(i, j, k, 0))
                            : 0.0;
      amrex::Real fracy = flag(i, j, k).isConnected(0, jj - j, 0)
                            ? amrex::Math::abs(fcz(i, j, k, 1))
                            : 0.0;

      fzm = (1.0 - fracx) * (1.0 - fracy) * fzm +
            fracx * (1.0 - fracy) * w(ii, j, k, n) +
            fracy * (1.0 - fracx) * w(i, jj, k, n) +
            fracx * fracy * w(ii, jj, k, n);
    }
    if (valid_cell) {
      fz(i, j, k, n) = fzm;
    }

    amrex::Real fzp = w(i, j, k + 1, n);
    if (apz(i, j, k + 1) != 0.0 && apz(i, j, k + 1) != 1.0) {
      int ii =
        i + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k + 1, 0)));
      int jj =
        j + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k + 1, 1)));
      amrex::Real fracx = flag(i, j, k + 1).isConnected(ii - i, 0, 0)
                            ? amrex::Math::abs(fcz(i, j, k + 1, 0))
                            : 0.0;
      amrex::Real fracy = flag(i, j, k + 1).isConnected(0, jj - j, 0)
                            ? amrex::Math::abs(fcz(i, j, k + 1, 1))
                            : 0.0;
      fzp = (1.0 - fracx) * (1.0 - fracy) * fzp +
            fracx * (1.0 - fracy) * w(ii, j, k + 1, n) +
            fracy * (1.0 - fracx) * w(i, jj, k + 1, n) +
            fracx * fracy * w(ii, jj, k + 1, n);
    }
    if (valid_cell && z_high) {
      fz(i, j, k + 1, n) = fzp;
    }

    amrex::GpuArray<amrex::Real, NVAR> flux_hyp_wall = {0.0};
    pc_compute_hyp_wallflux(
      i, j, k, q, qaux, apx(i, j, k), apx(i + 1, j, k), apy(i, j, k),
      apy(i, j + 1, k), apz(i, j, k), apz(i, j, k + 1), flux_hyp_wall);

    // With EB we assume dx == dy == dz
    // NOTE: we have already made the fluxes extensive so we define by dx*dy*dz
    // here ...
    divu(i, j, k, n) =
      cell_vol_inv / vfrc(i, j, k) *
      (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm + apy(i, j + 1, k) * fyp -
       apy(i, j, k) * fym + apz(i, j, k + 1) * fzp - apz(i, j, k) * fzm);
    divu(i, j, k, n) += flux_hyp_wall[n] * dxinv[0] / vfrc(i, j, k);
  }
#endif

  // The operations following this assume we have returned the negative of the
  // divergence of fluxes.
  divu(i, j, k, n) *= -1.0;

  // Go ahead and make the redistwgt array here since we'll need it in
  // flux_redistribute
  if (eb_weights_type == 0) {
    redistwgt(i, j, k) = 1.0;
  } else if (eb_weights_type == 1) {
    redistwgt(i, j, k) =
      q(i, j, k, QRHO) *
      (q(i, j, k, QREINT) + 0.5 * (AMREX_D_TERM(
                                    q(i, j, k, QU) * q(i, j, k, QU),
                                    +q(i, j, k, QV) * q(i, j, k, QV),
                                    +q(i, j, k, QW) * q(i, j, k, QW))));
  } else if (eb_weights_type == 2) {
    redistwgt(i, j, k) = q(i, j, k, QRHO);
  } else if (eb_weights_type == 3) {
    redistwgt(i, j, k) = vfrc(i, j, k);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_divu(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& q,
  AMREX_D_DECL(
    const amrex::Real dx, const amrex::Real dy, const amrex::Real dz),
  amrex::Array4<amrex::Real> const& divu)
{
#if AMREX_SPACEDIM == 1
  amrex::Real ux = (q(i, j, k, QU) - q(i - 1, j, k, QU)) / dx;
#elif AMREX_SPACEDIM == 2
  amrex::Real ux = 0.5 *
                   (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j - 1, k, QU) -
                    q(i - 1, j - 1, k, QU)) /
                   dx;
  amrex::Real vy = 0.5 *
                   (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i - 1, j, k, QV) -
                    q(i - 1, j - 1, k, QV)) /
                   dy;
#elif AMREX_SPACEDIM == 3
  amrex::Real ux =
    0.25 *
    (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j, k - 1, QU) -
     q(i - 1, j, k - 1, QU) + q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU) +
     q(i, j - 1, k - 1, QU) - q(i - 1, j - 1, k - 1, QU)) /
    dx;
  amrex::Real vy =
    0.25 *
    (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i, j, k - 1, QV) -
     q(i, j - 1, k - 1, QV) + q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV) +
     q(i - 1, j, k - 1, QV) - q(i - 1, j - 1, k - 1, QV)) /
    dy;
  amrex::Real wz =
    0.25 *
    (q(i, j, k, QW) - q(i, j, k - 1, QW) + q(i - 1, j, k, QW) -
     q(i - 1, j, k - 1, QW) + q(i, j - 1, k, QW) - q(i, j - 1, k - 1, QW) +
     q(i - 1, j - 1, k, QW) - q(i - 1, j - 1, k - 1, QW)) /
    dz;
#endif
  divu(i, j, k) = AMREX_D_TERM(ux, +vy, +wz);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
eb_divu(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real const> const& q,
  amrex::Array4<amrex::Real> const& divu,
  amrex::Array4<amrex::Real const> const& vfrac,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxInv)
{
#if (AMREX_SPACEDIM == 2)
  const amrex::Real surrounding_vol = vfrac(i, j, k) + vfrac(i - 1, j, k) +
                                      vfrac(i, j - 1, k) +
                                      vfrac(i - 1, j - 1, k);
  const amrex::Real cutoff = 4.0 * (1.0 - 1.e-12);
#elif (AMREX_SPACEDIM == 3)
  const amrex::Real surrounding_vol =
    vfrac(i, j, k) + vfrac(i - 1, j, k) + vfrac(i, j - 1, k) +
    vfrac(i - 1, j - 1, k) + vfrac(i, j, k - 1) + vfrac(i - 1, j, k - 1) +
    vfrac(i, j - 1, k - 1) + vfrac(i - 1, j - 1, k - 1);
  const amrex::Real cutoff = 8.0 * (1.0 - 1.e-12);
#endif

  divu(i, j, k) = 0.0;

  if (surrounding_vol >= cutoff) {

#if (AMREX_SPACEDIM == 2)
    const amrex::Real ux = 0.5 *
                           (q(i, j, k, QU) - q(i - 1, j, k, QU) +
                            q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU)) *
                           dxInv[0];
    const amrex::Real vy = 0.5 *
                           (q(i, j, k, QV) - q(i, j - 1, k, QV) +
                            q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV)) *
                           dxInv[1];
#elif (AMREX_SPACEDIM == 3)
    const amrex::Real ux =
      0.25 *
      (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j, k - 1, QU) -
       q(i - 1, j, k - 1, QU) + q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU) +
       q(i, j - 1, k - 1, QU) - q(i - 1, j - 1, k - 1, QU)) *
      dxInv[0];

    const amrex::Real vy =
      0.25 *
      (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i, j, k - 1, QV) -
       q(i, j - 1, k - 1, QV) + q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV) +
       q(i - 1, j, k - 1, QV) - q(i - 1, j - 1, k - 1, QV)) *
      dxInv[1];

    const amrex::Real wz =
      0.25 *
      (q(i, j, k, QW) - q(i, j, k - 1, QW) + q(i - 1, j, k, QW) -
       q(i - 1, j, k - 1, QW) + q(i, j - 1, k, QW) - q(i, j - 1, k - 1, QW) +
       q(i - 1, j - 1, k, QW) - q(i - 1, j - 1, k - 1, QW)) *
      dxInv[2];
#endif
    divu(i, j, k) = AMREX_D_TERM(ux, +vy, +wz);
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pc_eb_pdivu(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real const> const& q,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& q1,
    amrex::Array4<amrex::Real const> const& q2,
    amrex::Array4<amrex::Real const> const& q3),
  amrex::Array4<amrex::Real> const& divu,
  amrex::Array4<amrex::EBCellFlag const> const& flag,
  amrex::Array4<amrex::Real const> const& vfrc,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& apx,
    amrex::Array4<amrex::Real const> const& apy,
    amrex::Array4<amrex::Real const> const& apz),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& fcx,
    amrex::Array4<amrex::Real const> const& fcy,
    amrex::Array4<amrex::Real const> const& fcz),
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv)
{
  amrex::Real pdivu;

#if (AMREX_SPACEDIM == 2)
  if (flag(i, j, k).isRegular()) {
    pdivu = 0.5 * ((q1(i + 1, j, k, GDPRES) + q1(i, j, k, GDPRES)) *
                     (q1(i + 1, j, k, GDU) - q1(i, j, k, GDU)) * dxinv[0] +
                   (q2(i, j + 1, k, GDPRES) + q2(i, j, k, GDPRES)) *
                     (q2(i, j + 1, k, GDV) - q2(i, j, k, GDV)) * dxinv[1]);
    divu(i, j, k, UEINT) -= pdivu;
  } else if (!flag(i, j, k).isCovered()) {
    amrex::Real fxm, fxp, fym, fyp;
    if (apx(i, j, k) > 0.0) {
      fxm = q1(i, j, k, GDU);
      if (apx(i, j, k) != 0.0 && apx(i, j, k) != 1.0) {
        int jj =
          j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i, j, k, 0)));
        amrex::Real fracy = flag(i, j, k).isConnected(0, jj - j, 0)
                              ? amrex::Math::abs(fcx(i, j, k, 0))
                              : 0.0;
        fxm = (1.0 - fracy) * fxm + fracy * q1(i, jj, k, GDU);
      }
    } else {
      fxm = amrex::Real(0.0);
    }

    if (apx(i + 1, j, k) > 0.0) {
      fxp = q1(i + 1, j, k, GDU);
      if (apx(i + 1, j, k) != 0.0 && apx(i + 1, j, k) != 1.0) {
        int jj =
          j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i + 1, j, k, 0)));
        amrex::Real fracy = flag(i + 1, j, k).isConnected(0, jj - j, 0)
                              ? amrex::Math::abs(fcx(i + 1, j, k, 0))
                              : 0.0;
        fxp = (1.0 - fracy) * fxp + fracy * q1(i + 1, jj, k, GDU);
      }
    } else {
      fxp = amrex::Real(0.0);
    }

    if (apy(i, j, k) > 0.0) {
      fym = q2(i, j, k, GDV);
      if (apy(i, j, k) != 0.0 && apy(i, j, k) != 1.0) {
        int ii =
          i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j, k, 0)));
        amrex::Real fracx = flag(i, j, k).isConnected(ii - i, 0, 0)
                              ? amrex::Math::abs(fcy(i, j, k, 0))
                              : 0.0;
        fym = (1.0 - fracx) * fym + fracx * q2(ii, j, k, GDV);
      }
    } else {
      fym = amrex::Real(0.0);
    }

    if (apy(i, j + 1, k) > 0.0) {
      fyp = q2(i, j + 1, k, GDV);
      if (apy(i, j + 1, k) != 0.0 && apy(i, j + 1, k) != 1.0) {
        int ii =
          i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j + 1, k, 0)));
        amrex::Real fracx = flag(i, j + 1, k).isConnected(ii - i, 0, 0)
                              ? amrex::Math::abs(fcy(i, j + 1, k, 0))
                              : 0.0;
        fyp = (1.0 - fracx) * fyp + fracx * q2(ii, j + 1, k, GDV);
      }
    } else {
      fyp = amrex::Real(0.0);
    }

    divu(i, j, k, UEINT) -=
      (1.0 / vfrc(i, j, k)) * q(i, j, k, QPRES) *
      (dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
       dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym));
  }

#else // 3-d starts here

  if (flag(i, j, k).isRegular()) {
    pdivu = 0.5 * ((q1(i + 1, j, k, GDPRES) + q1(i, j, k, GDPRES)) *
                     (q1(i + 1, j, k, GDU) - q1(i, j, k, GDU)) * dxinv[0] +
                   (q2(i, j + 1, k, GDPRES) + q2(i, j, k, GDPRES)) *
                     (q2(i, j + 1, k, GDV) - q2(i, j, k, GDV)) * dxinv[1] +
                   (q3(i, j, k + 1, GDPRES) + q3(i, j, k, GDPRES)) *
                     (q3(i, j, k + 1, GDW) - q3(i, j, k, GDW)) * dxinv[2]);

    divu(i, j, k, UEINT) -= pdivu;
  } else if (!flag(i, j, k).isCovered()) {
    amrex::Real fxm, fxp, fym, fyp, fzm, fzp;
    if (apx(i, j, k) > 0.0) {
      fxm = q1(i, j, k, GDU);
      if (apx(i, j, k) != 0.0 && apx(i, j, k) != 1.0) {
        int jj =
          j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i, j, k, 0)));
        int kk =
          k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i, j, k, 1)));
        amrex::Real fracy = flag(i, j, k).isConnected(0, jj - j, 0)
                              ? amrex::Math::abs(fcx(i, j, k, 0))
                              : 0.0;
        amrex::Real fracz = flag(i, j, k).isConnected(0, 0, kk - k)
                              ? amrex::Math::abs(fcx(i, j, k, 1))
                              : 0.0;
        fxm = (1.0 - fracy) * (1.0 - fracz) * fxm +
              fracy * (1.0 - fracz) * q1(i, jj, k, GDU) +
              fracz * (1.0 - fracy) * q1(i, j, kk, GDU) +
              fracy * fracz * q1(i, jj, kk, GDU);
      }
    } else {
      fxm = amrex::Real(0.0);
    }

    if (apx(i + 1, j, k) > 0.0) {
      fxp = q1(i + 1, j, k, GDU);
      if (apx(i + 1, j, k) != 0.0 && apx(i + 1, j, k) != 1.0) {
        int jj =
          j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i + 1, j, k, 0)));
        int kk =
          k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i + 1, j, k, 1)));
        amrex::Real fracy = flag(i + 1, j, k).isConnected(0, jj - j, 0)
                              ? amrex::Math::abs(fcx(i + 1, j, k, 0))
                              : 0.0;
        amrex::Real fracz = flag(i + 1, j, k).isConnected(0, 0, kk - k)
                              ? amrex::Math::abs(fcx(i + 1, j, k, 1))
                              : 0.0;
        fxp = (1.0 - fracy) * (1.0 - fracz) * fxp +
              fracy * (1.0 - fracz) * q1(i + 1, jj, k, GDU) +
              fracz * (1.0 - fracy) * q1(i + 1, j, kk, GDU) +
              fracy * fracz * q1(i + 1, jj, kk, GDU);
      }
    } else {
      fxp = amrex::Real(0.0);
    }

    if (apy(i, j, k) > 0.0) {
      fym = q2(i, j, k, GDV);
      if (apy(i, j, k) != 0.0 && apy(i, j, k) != 1.0) {
        int ii =
          i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j, k, 0)));
        int kk =
          k + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j, k, 1)));
        amrex::Real fracx = flag(i, j, k).isConnected(ii - i, 0, 0)
                              ? amrex::Math::abs(fcy(i, j, k, 0))
                              : 0.0;
        amrex::Real fracz = flag(i, j, k).isConnected(0, 0, kk - k)
                              ? amrex::Math::abs(fcy(i, j, k, 1))
                              : 0.0;
        fym = (1.0 - fracx) * (1.0 - fracz) * fym +
              fracx * (1.0 - fracz) * q2(ii, j, k, GDV) +
              fracz * (1.0 - fracx) * q2(i, j, kk, GDV) +
              fracx * fracz * q2(ii, j, kk, GDV);
      }
    } else {
      fym = amrex::Real(0.0);
    }

    if (apy(i, j + 1, k) > 0.0) {
      fyp = q2(i, j + 1, k, GDV);
      if (apy(i, j + 1, k) != 0.0 && apy(i, j + 1, k) != 1.0) {
        int ii =
          i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j + 1, k, 0)));
        int kk =
          k + static_cast<int>(amrex::Math::copysign(1.0, fcy(i, j + 1, k, 1)));
        amrex::Real fracx = flag(i, j + 1, k).isConnected(ii - i, 0, 0)
                              ? amrex::Math::abs(fcy(i, j + 1, k, 0))
                              : 0.0;
        amrex::Real fracz = flag(i, j + 1, k).isConnected(0, 0, kk - k)
                              ? amrex::Math::abs(fcy(i, j + 1, k, 1))
                              : 0.0;
        fyp = (1.0 - fracx) * (1.0 - fracz) * fyp +
              fracx * (1.0 - fracz) * q2(ii, j + 1, k, GDV) +
              fracz * (1.0 - fracx) * q2(i, j + 1, kk, GDV) +
              fracx * fracz * q2(ii, j + 1, kk, GDV);
      }
    } else {
      fyp = amrex::Real(0.0);
    }

    if (apz(i, j, k) > 0.0) {
      fzm = q3(i, j, k, GDW);
      if (apz(i, j, k) != 0.0 && apz(i, j, k) != 1.0) {
        int ii =
          i + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k, 0)));
        int jj =
          j + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k, 1)));
        amrex::Real fracx = flag(i, j, k).isConnected(ii - i, 0, 0)
                              ? amrex::Math::abs(fcz(i, j, k, 0))
                              : 0.0;
        amrex::Real fracy = flag(i, j, k).isConnected(0, jj - j, 0)
                              ? amrex::Math::abs(fcz(i, j, k, 1))
                              : 0.0;

        fzm = (1.0 - fracx) * (1.0 - fracy) * fzm +
              fracx * (1.0 - fracy) * q3(ii, j, k, GDW) +
              fracy * (1.0 - fracx) * q3(i, jj, k, GDW) +
              fracx * fracy * q3(ii, jj, k, GDW);
      }
    } else {
      fzm = amrex::Real(0.0);
    }

    if (apz(i, j, k + 1) > 0.0) {
      fzp = q3(i, j, k + 1, GDW);
      if (apz(i, j, k + 1) != 0.0 && apz(i, j, k + 1) != 1.0) {
        int ii =
          i + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k + 1, 0)));
        int jj =
          j + static_cast<int>(amrex::Math::copysign(1.0, fcz(i, j, k + 1, 1)));
        amrex::Real fracx = flag(i, j, k + 1).isConnected(ii - i, 0, 0)
                              ? amrex::Math::abs(fcz(i, j, k + 1, 0))
                              : 0.0;
        amrex::Real fracy = flag(i, j, k + 1).isConnected(0, jj - j, 0)
                              ? amrex::Math::abs(fcz(i, j, k + 1, 1))
                              : 0.0;
        fzp = (1.0 - fracx) * (1.0 - fracy) * fzp +
              fracx * (1.0 - fracy) * q3(ii, j, k + 1, GDW) +
              fracy * (1.0 - fracx) * q3(i, jj, k + 1, GDW) +
              fracx * fracy * q3(ii, jj, k + 1, GDW);
      }
    } else {
      fzp = amrex::Real(0.0);
    }

    divu(i, j, k, UEINT) -=
      (1.0 / vfrc(i, j, k)) * q(i, j, k, QPRES) *
      (dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
       dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) +
       dxinv[2] * (apz(i, j, k + 1) * fzp - apz(i, j, k) * fzm));
  }
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_norm_spec_flx(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& flx)
{
  const int num = NUM_SPECIES + UFS;
  amrex::Real summ, fac;
  summ = 0.0;
  for (int n = UFS; n < num; ++n) {
    summ += flx(i, j, k, n);
  }

  fac = (summ != 0.0) ? flx(i, j, k, URHO) / summ : 1.0;
  for (int n = UFS; n < num; ++n) {
    flx(i, j, k, n) *= fac;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ext_flx(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<const amrex::Real> const& a)
{
  amrex::Real area = a(i, j, k); // avoid reading from global memory NVAR times
  for (int n = 0; n < NVAR; ++n) {
    flx(i, j, k, n) *= area;
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pc_ext_flx_eb(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Real area,
  amrex::Array4<const amrex::Real> const& aper)
{
  // The fluxes are multiplied by area fractions inside amrex
  //  so here we just weight them by the regular areas
  if (aper(i, j, k) > 0.) {
    for (int n = 0; n < NVAR; ++n) {
      flx(i, j, k, n) *= area;
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_update(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& pdivu)
{
  amrex::Real v = vol(i, j, k);
  for (int n = 0; n < NVAR; ++n) {
    update(i, j, k, n) += (AMREX_D_TERM(
                            flx[0](i, j, k, n) - flx[0](i + 1, j, k, n),
                            +flx[1](i, j, k, n) - flx[1](i, j + 1, k, n),
                            +flx[2](i, j, k, n) - flx[2](i, j, k + 1, n))) /
                          v;
  }

  update(i, j, k, UEINT) -= pdivu(i, j, k);
}

// Host functions
void pc_umdrv(
  const amrex::Real time,
  amrex::Box const& bx,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    src_q, // amrex::IArrayBox const& bcMask,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
  const amrex::Real dt,
  const int ppm_type,
  const int plm_iorder,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme,
  const amrex::Real difmag,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
    a,
  amrex::Array4<amrex::Real> const& vol,
  amrex::Real cflLoc);

void pc_umdrv_eb(
  const amrex::Box& bx,
  const amrex::Box& bxg_i,
  const amrex::MFIter& mfi,
  const amrex::Geometry& geom,
  const amrex::EBFArrayBoxFactory* fact,
  const int* bclo,
  const int* bchi,
  const amrex::Array4<const amrex::Real>& uin,
  const amrex::Array4<amrex::Real>& uout,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::Array4<const amrex::Real>& src_q,
  const amrex::Array4<const amrex::Real>& vf,
  const amrex::Array4<amrex::EBCellFlag const>& flag,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const int as_crse,
  const amrex::Array4<amrex::Real>& drho_as_crse,
  const amrex::Array4<int const>& rrflag_as_crse,
  const int as_fine,
  const amrex::Array4<amrex::Real>& dm_as_fine,
  const amrex::Array4<int const>& lev_mask,
  const amrex::Real dt,
  const int ppm_type,
  const int plm_iorder,
  const bool use_flattening,
  const amrex::Real difmag,
  amrex::BCRec const* bcs_d_ptr,
  const std::string& redistribution_type,
  const int eb_weights_type,
  const int eb_srd_max_order,
  const bool eb_clean_massfrac,
  const amrex::Real eb_clean_massfrac_threshold,
  amrex::Real cflLoc);

void pc_adjust_fluxes(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
    a,
  const amrex::Array4<const amrex::Real>& divu,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& del,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  const amrex::Real difmag);

void pc_adjust_fluxes_eb(
  const amrex::Box& /*bx*/,
  amrex::Array4<const amrex::Real> const& q_arr,
  amrex::Array4<const amrex::Real> const& u_arr,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& apx,
    amrex::Array4<amrex::Real const> const& apy,
    amrex::Array4<amrex::Real const> const& apz),
  amrex::Array4<const amrex::Real> const& vfrac,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flux,
  const amrex::Geometry& geom,
  const int* bclo,
  const int* bchi,
  amrex::Real difmag);

void pc_consup(
  const amrex::Box& bx,
  const amrex::Array4<amrex::Real>& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::Array4<const amrex::Real>& vol,
  const amrex::Array4<const amrex::Real>& pdivu);

void pc_consup_eb(
  const amrex::Box& bx,
  amrex::Array4<const amrex::Real> const& q_arr,
  amrex::Array4<const amrex::Real> const& qaux_arr,
  amrex::Array4<amrex::Real> const& divc_arr,
  amrex::Array4<amrex::Real> const& redistwgt_arr,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real> const& q1,
    amrex::Array4<amrex::Real> const& q2,
    amrex::Array4<amrex::Real> const& q3),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& apx,
    amrex::Array4<amrex::Real const> const& apy,
    amrex::Array4<amrex::Real const> const& apz),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& fcx,
    amrex::Array4<amrex::Real const> const& fcy,
    amrex::Array4<amrex::Real const> const& fcz),
  amrex::Array4<const amrex::Real> const& vfrac,
  amrex::Array4<amrex::EBCellFlag const> const& flag,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>&
    flux_tmp,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flux,
  const int eb_weights_type);

#endif
