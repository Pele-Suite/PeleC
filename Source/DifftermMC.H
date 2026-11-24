#ifndef DIFFTERMMC_H
#define DIFFTERMMC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>

#include "PelePhysics.H"
#include "IndexDefines.H"

#ifdef PELEC_MULTICOMPONENT_DIFFUSION

// Multicomponent species and energy diffusion fluxes (ideal gas only).
//
// Adds to the face fluxes flx the species mass fluxes
//   F_k = - sum_j Dmat_kj ( d_j + chi_j dlnT )
// with the generalized driving forces d_j = dX_j + (X_j - Y_j) dlnp, the
// flux diffusion matrix Dmat_kj = rho Y_k D_kj and the thermal diffusion
// ratios chi_j = X_j chit_j (EGlib convention, used when do_soret), the
// corresponding species enthalpy flux sum_k h_k F_k and, when do_dufour, the
// Dufour heat flux p sum_k chi_k V_k written as R T sum_k chit_k F_k / W_k.
// Dmat and chit are read from the cell-centered coefficient array and
// averaged arithmetically to the face, which preserves sum_k Dmat_kj = 0 and
// hence sum_k F_k = 0 exactly: no correction velocity is needed.
//
// Must be called after pc_diffusion_flux / pc_diffusion_flux_eb, which have
// set the viscous and Fourier contributions (and zero species fluxes, since
// the mixture-averaged coefficients are zero in this mode) and scaled them by
// the face area; the contributions added here are scaled by the area too.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_mc_species_energy_flux(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& coef_cc,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const amrex::Array4<const amrex::Real>& area,
  const amrex::Array4<amrex::Real>& flx,
  const amrex::Real dxinv,
  const int dir,
  const bool do_soret,
  const bool do_dufour)
{
  amrex::ignore_unused(k); // not part of the IntVect in 2D
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
  if (!(flags(iv).isRegular() || flags(iv).isSingleValued())) {
    return;
  }
  const amrex::Real a = area(iv);
  if (a == 0.0) {
    return;
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
  amrex::Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  amrex::Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    mass1[n] = q(iv, QFS + n);
    mass2[n] = q(ivm, QFS + n);
  }
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);
  eos.T2Hi(q(iv, QTEMP), hi1);
  eos.T2Hi(q(ivm, QTEMP), hi2);

  const amrex::Real pface = 0.5 * (q(iv, QPRES) + q(ivm, QPRES));
  const amrex::Real Tface = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  const amrex::Real dlnp = dxinv * (q(iv, QPRES) - q(ivm, QPRES)) / pface;
  const amrex::Real dlnT = dxinv * (q(iv, QTEMP) - q(ivm, QTEMP)) / Tface;

  // Generalized driving forces (including the Soret part)
  amrex::Real d[NUM_SPECIES];
  amrex::Real chitf[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real Xf = 0.5 * (mole1[n] + mole2[n]);
    const amrex::Real Yf = 0.5 * (mass1[n] + mass2[n]);
    chitf[n] =
      do_soret
        ? 0.5 * (coef_cc(iv, dComp_chit + n) + coef_cc(ivm, dComp_chit + n))
        : 0.0;
    d[n] =
      dxinv * (mole1[n] - mole2[n]) + (Xf - Yf) * dlnp + Xf * chitf[n] * dlnT;
  }

  // F_k = - sum_j Dmat_kj d_j, Dmat averaged arithmetically to the face
  amrex::Real wt[NUM_SPECIES];
  if (do_dufour) {
    eos.molecular_weight(wt);
  }
  amrex::Real dufour = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    amrex::Real Fk = 0.0;
    for (int m = 0; m < NUM_SPECIES; ++m) {
      const int c = dComp_Dmat + n * NUM_SPECIES + m;
      Fk -= 0.5 * (coef_cc(iv, c) + coef_cc(ivm, c)) * d[m];
    }
    flx(iv, UFS + n) += Fk * a;
    flx(iv, UEDEN) += Fk * 0.5 * (hi1[n] + hi2[n]) * a;
    if (do_dufour) {
      dufour += chitf[n] * Fk / wt[n];
    }
  }
  if (do_dufour) {
    flx(iv, UEDEN) += pele::physics::Constants::RU * Tface * dufour * a;
  }
}

#endif // PELEC_MULTICOMPONENT_DIFFUSION
#endif
