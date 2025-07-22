#ifndef LES_H
#define LES_H

#include "PeleC.H"
#include "PeleCUtilities.H"
#include "Diffusion.H"
#include "Diffterm.H"
#include "GradUtil.H"
#include "PelePhysics.H"
#include "Filter.H"

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#if AMREX_SPACEDIM == 3
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_smagorinsky_sfs_stresses(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  const int dir,
  amrex::Real& alpha0,
  amrex::Real& alpha1,
  amrex::Real& alpha2,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  amrex::Real dUdx[3][3] = {{0.0}};

  // Normal velocity gradients
  dUdx[0][dir] = dxinv * (q(iv, QU) - q(ivm, QU));
  dUdx[1][dir] = dxinv * (q(iv, QV) - q(ivm, QV));
  dUdx[2][dir] = dxinv * (q(iv, QW) - q(ivm, QW));
  // Tangential components using precomputed derivatives
  const int c1 = dir == 0 ? 1 : 0;
  const int c2 = dir == 2 ? 1 : 2;

  dUdx[0][c1] = td(iv, 0);
  dUdx[0][c2] = td(iv, 3);
  dUdx[1][c1] = td(iv, 1);
  dUdx[1][c2] = td(iv, 4);
  dUdx[2][c1] = td(iv, 2);
  dUdx[2][c2] = td(iv, 5);

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0.0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (dUdx[m][n] + dUdx[n][m]);
      Sijmag += S[m][n] * S[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut =
    0.5 * (q(iv, QRHO) + q(ivm, QRHO)) * deltabar * deltabar * Sijmag;

  alpha0 = 2.0 * mut * (S[dir][0] - (dir == 0 ? Skk / 3.0 : 0.0));
  alpha1 = 2.0 * mut * (S[dir][1] - (dir == 1 ? Skk / 3.0 : 0.0));
  alpha2 = 2.0 * mut * (S[dir][2] - (dir == 2 ? Skk / 3.0 : 0.0));
  alpha = 2.0 * mut * Sijmag;

  const amrex::Real dTdx = dxinv * (q(iv, QTEMP) - q(ivm, QTEMP));
  flux_T = mut * dTdx;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_smagorinsky_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const int dir,
  const amrex::Real Cs,
  const amrex::Real CI,
  const amrex::Real PrT,
  const amrex::Array4<amrex::Real>& flx)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real dxinv = 1.0 / dx[dir];
  const amrex::Real deltabar =
    std::pow(dx[0] * dx[1] * dx[2], static_cast<amrex::Real>(1.0 / 3.0));
  const amrex::Real Cs2 = Cs * Cs;
  const amrex::Real third = 1.0 / 3.0;
  amrex::Real alphaij[AMREX_SPACEDIM] = {0.0}, alpha, flux_T;
  get_smagorinsky_sfs_stresses(
    i, j, k, q, td, dxinv, deltabar, dir, alphaij[0], alphaij[1], alphaij[2],
    alpha, flux_T);
  const amrex::Real sigmadx =
    Cs2 * alphaij[0] - third * CI * alpha * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    Cs2 * alphaij[1] - third * CI * alpha * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    Cs2 * alphaij[2] - third * CI * alpha * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * Cs2 / PrT * flux_T;

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_sfs_stresses_cc(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Real dxinv[AMREX_SPACEDIM],
  const amrex::Real deltabar,
  amrex::Real& alpha00,
  amrex::Real& alpha01,
  amrex::Real& alpha02,
  amrex::Real& alpha11,
  amrex::Real& alpha12,
  amrex::Real& alpha22,
  amrex::Real& alpha,
  amrex::Real& flux_T0,
  amrex::Real& flux_T1,
  amrex::Real& flux_T2)
{
  // Calculate derivatives at cell centers, second order central difference
  amrex::Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QU) - q(i - 1, j, k, QU));
  dUdx[0][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QU) - q(i, j - 1, k, QU));
  dUdx[0][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QU) - q(i, j, k - 1, QU));
  dUdx[1][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QV) - q(i - 1, j, k, QV));
  dUdx[1][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QV) - q(i, j - 1, k, QV));
  dUdx[1][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QV) - q(i, j, k - 1, QV));
  dUdx[2][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QW) - q(i - 1, j, k, QW));
  dUdx[2][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QW) - q(i, j - 1, k, QW));
  dUdx[2][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QW) - q(i, j, k - 1, QW));

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (dUdx[m][n] + dUdx[n][m]);
      Sijmag += S[m][n] * S[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];
  const amrex::Real mut = q(i, j, k, QRHO) * deltabar * deltabar * Sijmag;

  alpha00 = 2.0 * mut * (S[0][0] - Skk / 3.0);
  alpha01 = 2.0 * mut * (S[0][1]);
  alpha02 = 2.0 * mut * (S[0][2]);
  alpha11 = 2.0 * mut * (S[1][1] - Skk / 3.0);
  alpha12 = 2.0 * mut * (S[1][2]);
  alpha22 = 2.0 * mut * (S[2][2] - Skk / 3.0);
  alpha = 2.0 * mut * Sijmag;

  flux_T0 =
    mut * 0.5 * dxinv[0] * (q(i + 1, j, k, QTEMP) - q(i - 1, j, k, QTEMP));
  flux_T1 =
    mut * 0.5 * dxinv[1] * (q(i, j + 1, k, QTEMP) - q(i, j - 1, k, QTEMP));
  flux_T2 =
    mut * 0.5 * dxinv[2] * (q(i, j, k + 1, QTEMP) - q(i, j, k - 1, QTEMP));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_dynamic_smagorinsky_quantities(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const int fgr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::Array4<amrex::Real>& Kij,
  const amrex::Array4<amrex::Real>& RUT,
  const amrex::Array4<amrex::Real>& alphaij,
  const amrex::Array4<amrex::Real>& alpha,
  const amrex::Array4<amrex::Real>& flux_T)
{
  // indexes
  const int i00 = 0 * AMREX_SPACEDIM + 0; // 0
  const int i01 = 1 * AMREX_SPACEDIM + 0; // 3
  const int i02 = 2 * AMREX_SPACEDIM + 0; // 6
  const int i10 = 0 * AMREX_SPACEDIM + 1; // 1
  const int i11 = 1 * AMREX_SPACEDIM + 1; // 4
  const int i12 = 2 * AMREX_SPACEDIM + 1; // 7
  const int i20 = 0 * AMREX_SPACEDIM + 2; // 2
  const int i21 = 1 * AMREX_SPACEDIM + 2; // 5
  const int i22 = 2 * AMREX_SPACEDIM + 2; // 8

  const amrex::Real deltabar =
    fgr * std::pow(dx[0] * dx[1] * dx[2], static_cast<amrex::Real>(1.0 / 3.0));
  const amrex::Real dxinv[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(1.0 / dx[0]),
    static_cast<amrex::Real>(1.0 / dx[1]),
    static_cast<amrex::Real>(1.0 / dx[2])};

  get_sfs_stresses_cc(
    i, j, k, q, dxinv, deltabar, alphaij(i, j, k, i00), alphaij(i, j, k, i01),
    alphaij(i, j, k, i02), alphaij(i, j, k, i11), alphaij(i, j, k, i12),
    alphaij(i, j, k, i22), alpha(i, j, k), flux_T(i, j, k, 0),
    flux_T(i, j, k, 1), flux_T(i, j, k, 2));
  alphaij(i, j, k, i10) = alphaij(i, j, k, i01);
  alphaij(i, j, k, i20) = alphaij(i, j, k, i02);
  alphaij(i, j, k, i21) = alphaij(i, j, k, i12);

  Kij(i, j, k, 0) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QU);
  Kij(i, j, k, 1) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QV);
  Kij(i, j, k, 2) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QW);
  Kij(i, j, k, 3) = q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QV);
  Kij(i, j, k, 4) = q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QW);
  Kij(i, j, k, 5) = q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QW);

  RUT(i, j, k, 0) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QTEMP);
  RUT(i, j, k, 1) = q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QTEMP);
  RUT(i, j, k, 2) = q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QTEMP);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_dynamic_smagorinsky_coeffs(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const int fgr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::Array4<const amrex::Real>& Kij,
  const amrex::Array4<const amrex::Real>& RUT,
  const amrex::Array4<const amrex::Real>& alphaij,
  const amrex::Array4<const amrex::Real>& alpha,
  const amrex::Array4<const amrex::Real>& flux_T,
  const amrex::Array4<amrex::Real>& coeff)
{
  // indexes
  const int i00 = 0 * AMREX_SPACEDIM + 0;
  const int i01 = 1 * AMREX_SPACEDIM + 0;
  const int i02 = 2 * AMREX_SPACEDIM + 0;
  const int i10 = 0 * AMREX_SPACEDIM + 1;
  const int i11 = 1 * AMREX_SPACEDIM + 1;
  const int i12 = 2 * AMREX_SPACEDIM + 1;
  const int i20 = 0 * AMREX_SPACEDIM + 2;
  const int i21 = 1 * AMREX_SPACEDIM + 2;
  const int i22 = 2 * AMREX_SPACEDIM + 2;

  const amrex::Real deltahat =
    fgr * std::pow(dx[0] * dx[1] * dx[2], static_cast<amrex::Real>(1.0 / 3.0));
  const amrex::Real dxinv[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(1.0 / dx[0]),
    static_cast<amrex::Real>(1.0 / dx[1]),
    static_cast<amrex::Real>(1.0 / dx[2])};

  amrex::Real L[AMREX_SPACEDIM * AMREX_SPACEDIM],
    M[AMREX_SPACEDIM * AMREX_SPACEDIM], betaij[AMREX_SPACEDIM * AMREX_SPACEDIM],
    beta, T[AMREX_SPACEDIM], KE[AMREX_SPACEDIM];

  get_sfs_stresses_cc(
    i, j, k, q, dxinv, deltahat, betaij[i00], betaij[i01], betaij[i02],
    betaij[i11], betaij[i12], betaij[i22], beta, T[0], T[1], T[2]);
  betaij[i10] = betaij[i01];
  betaij[i20] = betaij[i02];
  betaij[i21] = betaij[i12];

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    T[dir] = flux_T(i, j, k, dir) - T[dir];
  }

  // "resolved turbulent stresses" and others
  for (int dir = 0; dir < AMREX_SPACEDIM * AMREX_SPACEDIM; dir++) {
    M[dir] = betaij[dir] - alphaij(i, j, k, dir);
  }

  L[i00] = Kij(i, j, k, 0) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QU);
  L[i01] = Kij(i, j, k, 1) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QV);
  L[i02] = Kij(i, j, k, 2) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QW);
  L[i10] = L[i01];
  L[i11] = Kij(i, j, k, 3) - q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QV);
  L[i12] = Kij(i, j, k, 4) - q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QW);
  L[i20] = L[i02];
  L[i21] = L[i12];
  L[i22] = Kij(i, j, k, 5) - q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QW);

  KE[0] =
    RUT(i, j, k, 0) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QTEMP);
  KE[1] =
    RUT(i, j, k, 1) - q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QTEMP);
  KE[2] =
    RUT(i, j, k, 2) - q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QTEMP);

  // Contractions
  amrex::Real TT = 0.0, KT = 0.0, LM = 0.0, MM = 0.0;
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    TT += T[dir] * T[dir];
    KT += KE[dir] * T[dir];
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      LM += L[dir * AMREX_SPACEDIM + n] * M[dir * AMREX_SPACEDIM + n];
      MM += M[dir * AMREX_SPACEDIM + n] * M[dir * AMREX_SPACEDIM + n];
    }
  }
  const amrex::Real Lkk = L[i00] + L[i11] + L[i22];
  const amrex::Real bma = beta - alpha(i, j, k);

  // Coefficients
  coeff(i, j, k, comp_Cs2) = amrex::max<amrex::Real>(
    LM / (MM + constants::small_num()), constants::small_num());
  coeff(i, j, k, comp_CI) = amrex::max<amrex::Real>(
    Lkk / (bma + constants::small_num()), constants::small_num());
  coeff(i, j, k, comp_Cs2ovPrT) = amrex::max<amrex::Real>(
    KT / (TT + constants::small_num()), constants::small_num());

  // Limit CI to reasonable bounds
  coeff(i, j, k, comp_CI) =
    amrex::min<amrex::Real>(coeff(i, j, k, comp_CI), 1.0);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_dynamic_smagorinsky_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& alphaij,
  const amrex::Array4<const amrex::Real>& alpha,
  const amrex::Array4<const amrex::Real>& flux_T,
  const amrex::Array4<const amrex::Real>& coeff,
  const amrex::Array4<const amrex::Real>& a,
  const int dir,
  const amrex::Array4<amrex::Real>& flx)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real sigmadx =
    coeff(iv, comp_Cs2) * alphaij(iv, 0) -
    coeff(iv, comp_CI) * alpha(iv) / 3.0 * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    coeff(iv, comp_Cs2) * alphaij(iv, 1) -
    coeff(iv, comp_CI) * alpha(iv) / 3.0 * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    coeff(iv, comp_Cs2) * alphaij(iv, 2) -
    coeff(iv, comp_CI) * alpha(iv) / 3.0 * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * coeff(iv, comp_Cs2ovPrT) * flux_T(iv);

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_wale_sfs_stresses(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  const int dir,
  amrex::Real& alpha0,
  amrex::Real& alpha1,
  amrex::Real& alpha2,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  amrex::Real g[3][3] = {{0.0}};

  // Normal velocity gradients
  g[0][dir] = dxinv * (q(iv, QU) - q(ivm, QU));
  g[1][dir] = dxinv * (q(iv, QV) - q(ivm, QV));
  g[2][dir] = dxinv * (q(iv, QW) - q(ivm, QW));
  // Tangential components using precomputed derivatives
  const int c1 = dir == 0 ? 1 : 0;
  const int c2 = dir == 2 ? 1 : 2;

  g[0][c1] = td(iv, 0);
  g[0][c2] = td(iv, 3);
  g[1][c1] = td(iv, 1);
  g[1][c2] = td(iv, 4);
  g[2][c1] = td(iv, 2);
  g[2][c2] = td(iv, 5);

  amrex::Real g2[3][3] = {{0.0}};
  g2[0][0] = g[0][0] * g[0][0] + g[0][1] * g[1][0] + g[0][2] * g[2][0];
  g2[0][1] = g[0][0] * g[0][1] + g[0][1] * g[1][1] + g[0][2] * g[2][1];
  g2[0][2] = g[0][0] * g[0][2] + g[0][1] * g[1][2] + g[0][2] * g[2][2];
  g2[1][0] = g[1][0] * g[0][0] + g[1][1] * g[1][0] + g[1][2] * g[2][0];
  g2[1][1] = g[1][0] * g[0][1] + g[1][1] * g[1][1] + g[1][2] * g[2][1];
  g2[1][2] = g[1][0] * g[0][2] + g[1][1] * g[1][2] + g[1][2] * g[2][2];
  g2[2][0] = g[2][0] * g[0][0] + g[2][1] * g[1][0] + g[2][2] * g[2][0];
  g2[2][1] = g[2][0] * g[0][1] + g[2][1] * g[1][1] + g[2][2] * g[2][1];
  g2[2][2] = g[2][0] * g[0][2] + g[2][1] * g[1][2] + g[2][2] * g[2][2];

  amrex::Real g2kk = g2[0][0] + g2[1][1] + g2[2][2];

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0.0;
  amrex::Real SijSij = 0.0;
  amrex::Real D[3][3] = {{0.0}};
  amrex::Real DijDij = 0.0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (g[m][n] + g[n][m]);
      SijSij += S[m][n] * S[m][n];
      D[m][n] = 0.5 * (g2[m][n] + g2[n][m]);
      if (m == n) {
        D[m][n] = D[m][n] - g2kk / 3.0;
      }
      DijDij += D[m][n] * D[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  amrex::Real DijDij3 = DijDij * DijDij * DijDij;
  amrex::Real DijDij5 = DijDij3 * DijDij * DijDij;
  amrex::Real SijSij5 = SijSij * SijSij * SijSij * SijSij * SijSij;

  const amrex::Real mut =
    0.5 * (q(iv, QRHO) + q(ivm, QRHO)) * deltabar * deltabar *
    (std::sqrt(DijDij3) / (std::sqrt(SijSij5) + std::sqrt(std::sqrt(DijDij5))));

  alpha0 = 2.0 * mut * (S[dir][0] - (dir == 0 ? Skk / 3.0 : 0.0));
  alpha1 = 2.0 * mut * (S[dir][1] - (dir == 1 ? Skk / 3.0 : 0.0));
  alpha2 = 2.0 * mut * (S[dir][2] - (dir == 2 ? Skk / 3.0 : 0.0));
  alpha = 2.0 * 0.5 * (q(iv, QRHO) + q(ivm, QRHO)) * deltabar * deltabar *
          Sijmag * Sijmag;

  const amrex::Real dTdx = dxinv * (q(iv, QTEMP) - q(ivm, QTEMP));
  flux_T = mut * dTdx;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_wale_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const int dir,
  const amrex::Real Cw,
  const amrex::Real CI,
  const amrex::Real PrT,
  const amrex::Array4<amrex::Real>& flx)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real dxinv = 1.0 / dx[dir];
  const amrex::Real deltabar =
    std::pow(dx[0] * dx[1] * dx[2], static_cast<amrex::Real>(1.0 / 3.0));
  const amrex::Real Cw2 = Cw * Cw;
  const amrex::Real third = 1.0 / 3.0;
  amrex::Real alphaij[AMREX_SPACEDIM] = {0.0}, alpha, flux_T;
  get_wale_sfs_stresses(
    i, j, k, q, td, dxinv, deltabar, dir, alphaij[0], alphaij[1], alphaij[2],
    alpha, flux_T);
  const amrex::Real sigmadx =
    Cw2 * alphaij[0] - third * CI * alpha * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    Cw2 * alphaij[1] - third * CI * alpha * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    Cw2 * alphaij[2] - third * CI * alpha * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * Cw2 / PrT * flux_T;

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_vreman_sfs_stresses(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  const int dir,
  amrex::Real& alpha0,
  amrex::Real& alpha1,
  amrex::Real& alpha2,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  amrex::Real alfa[3][3] = {{0.0}};

  // Normal velocity gradients
  alfa[0][dir] = dxinv * (q(iv, QU) - q(ivm, QU));
  alfa[1][dir] = dxinv * (q(iv, QV) - q(ivm, QV));
  alfa[2][dir] = dxinv * (q(iv, QW) - q(ivm, QW));
  // Tangential components using precomputed derivatives
  const int c1 = dir == 0 ? 1 : 0;
  const int c2 = dir == 2 ? 1 : 2;

  alfa[0][c1] = td(iv, 0);
  alfa[0][c2] = td(iv, 3);
  alfa[1][c1] = td(iv, 1);
  alfa[1][c2] = td(iv, 4);
  alfa[2][c1] = td(iv, 2);
  alfa[2][c2] = td(iv, 5);

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real beta[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real aijaij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (alfa[m][n] + alfa[n][m]);
      SijSij += S[m][n] * S[m][n];
      aijaij += alfa[m][n] * alfa[m][n];
      beta[m][n] = deltabar * deltabar *
                   (alfa[0][m] * alfa[0][n] + alfa[1][m] * alfa[1][n] +
                    alfa[2][m] * alfa[2][n]);
    }
  }
  amrex::Real bbeta = 0;
  bbeta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
          beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
          beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];

  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut =
    0.5 * (q(iv, QRHO) + q(ivm, QRHO)) * (std::sqrt(bbeta / aijaij));

  alpha0 = 2.0 * mut * (S[dir][0] - (dir == 0 ? Skk / 3.0 : 0.0));
  alpha1 = 2.0 * mut * (S[dir][1] - (dir == 1 ? Skk / 3.0 : 0.0));
  alpha2 = 2.0 * mut * (S[dir][2] - (dir == 2 ? Skk / 3.0 : 0.0));
  alpha = 2.0 * 0.5 * (q(iv, QRHO) + q(ivm, QRHO)) * deltabar * deltabar *
          Sijmag * Sijmag;

  const amrex::Real dTdx = dxinv * (q(iv, QTEMP) - q(ivm, QTEMP));
  flux_T = mut * dTdx;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_vreman_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const int dir,
  const amrex::Real Cs,
  const amrex::Real CI,
  const amrex::Real PrT,
  const amrex::Array4<amrex::Real>& flx)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real dxinv = 1.0 / dx[dir];
  const amrex::Real deltabar =
    std::pow(dx[0] * dx[1] * dx[2], static_cast<amrex::Real>(1.0 / 3.0));
  const amrex::Real c = 2.5 * Cs * Cs;
  const amrex::Real third = 1.0 / 3.0;
  amrex::Real alphaij[AMREX_SPACEDIM] = {0.0}, alpha, flux_T;
  get_vreman_sfs_stresses(
    i, j, k, q, td, dxinv, deltabar, dir, alphaij[0], alphaij[1], alphaij[2],
    alpha, flux_T);
  const amrex::Real sigmadx =
    c * alphaij[0] - third * CI * alpha * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    c * alphaij[1] - third * CI * alpha * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    c * alphaij[2] - third * CI * alpha * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * c / PrT * flux_T;

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

void computeTangentialVelDerivs(
  const amrex::Array<const amrex::Box, AMREX_SPACEDIM>& eboxes,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& tander_ec,
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& tanders,
  const amrex::Array4<amrex::Real>& q_ar,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx);

void resizeAndSetFlux(
  const amrex::Array<const amrex::Box, AMREX_SPACEDIM>& eboxes,
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  amrex::Array<amrex::FArrayBox, AMREX_SPACEDIM>& flux_ec);

void computeFluxDiv(
  amrex::MultiFab& LESTerm,
  const amrex::MFIter& mfi,
  const amrex::Box& vbox,
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::MultiFab& volume);

#endif // AMREX_SPACEDIM == 3
#endif
