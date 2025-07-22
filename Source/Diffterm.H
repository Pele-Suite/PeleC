#ifndef DIFFTERM_H
#define DIFFTERM_H

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>

#include "IndexDefines.H"
#include "EBStencilTypes.H"
#include "PelePhysics.H"
#include "PeleCUtilities.H"
#include "GradUtil.H"
#include "TransCoeff.H"

// This header file contains functions and declarations for diffterm.

namespace {
AMREX_GPU_CONSTANT const amrex::Real weights[3] = {0.0, 1.0, 0.5};
} // namespace

template <typename EOSType>
struct SpeciesEnergyFlux
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect iv,
    const amrex::IntVect ivm,
    const amrex::Real dxinv,
    const amrex::GpuArray<amrex::Real, dComp_lambda + 1>& coef,
    const amrex::Array4<const amrex::Real>& q,
    const amrex::Array4<amrex::Real>& flx)
  {
    auto eos = pele::physics::PhysicsType::eos();

    // Get massfrac, molefrac, enthalpy
    amrex::Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
    amrex::Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
    amrex::Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass1[ns] = q(iv, ns + QFS);
    }
    eos.Y2X(mass1, mole1);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass2[ns] = q(ivm, ns + QFS);
    }
    eos.Y2X(mass2, mole2);

    // Compute species and enthalpy fluxes for ideal EOS
    // Get species/enthalpy diffusion, compute correction vel
    amrex::Real T = q(iv, QTEMP);
    eos.T2Hi(T, hi1);
    T = q(ivm, QTEMP);
    eos.T2Hi(T, hi2);
    amrex::Real Vc = 0.0;
    const amrex::Real dpdx = dxinv * (q(iv, QPRES) - q(ivm, QPRES));
    const amrex::Real dlnp = dpdx / (0.5 * (q(iv, QPRES) + q(ivm, QPRES)));
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      const amrex::Real Xface = 0.5 * (mole1[ns] + mole2[ns]);
      const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
      const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
      const amrex::Real dXdx = dxinv * (mole1[ns] - mole2[ns]);
      const amrex::Real Vd =
        -coef[dComp_rhoD + ns] * (dXdx + (Xface - Yface) * dlnp);
      flx(iv, UFS + ns) = Vd;
      Vc += Vd;
      flx(iv, UEDEN) += Vd * hface;
    }
    // Add correction velocity to fluxes
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
      const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
      flx(iv, UFS + ns) -= Yface * Vc;
      flx(iv, UEDEN) -= Yface * hface * Vc;
    }
  }
};

template <>
struct SpeciesEnergyFlux<pele::physics::eos::SRK>
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect iv,
    const amrex::IntVect ivm,
    const amrex::Real dxinv,
    const amrex::GpuArray<amrex::Real, dComp_lambda + 1>& coef,
    const amrex::Array4<const amrex::Real>& q,
    const amrex::Array4<amrex::Real>& flx)
  {
    pele::physics::eos::SRK eos;

    // Get massfrac, molefrac, enthalpy
    amrex::Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
    amrex::Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
    amrex::Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass1[ns] = q(iv, ns + QFS);
    }
    eos.Y2X(mass1, mole1);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass2[ns] = q(ivm, ns + QFS);
    }
    eos.Y2X(mass2, mole2);

    // Compute species and enthalpy fluxes accounting for nonideal EOS
    // Implementation note: nonideal EOS coeffs are evaluated at cell centers,
    // then inteprolated to faces. Could be slightly more efficient to
    // interpolate state to face, then do the EOS calls.

    const amrex::Real rho1 = q(iv, QRHO);
    const amrex::Real rho2 = q(ivm, QRHO);
    const amrex::Real T1 = q(iv, QTEMP);
    const amrex::Real T2 = q(ivm, QTEMP);
    const amrex::Real dpdx = dxinv * (q(iv, QPRES) - q(ivm, QPRES));

    amrex::Real Vc = 0.0;
    amrex::Real diP1[NUM_SPECIES], dijY1[NUM_SPECIES][NUM_SPECIES];
    eos.RTY2transport(rho1, T1, mass1, diP1, dijY1);
    eos.RTY2Hi(rho1, T1, mass1, hi1);
    amrex::Real diP2[NUM_SPECIES], dijY2[NUM_SPECIES][NUM_SPECIES];
    eos.RTY2transport(rho2, T2, mass2, diP2, dijY2);
    eos.RTY2Hi(rho2, T2, mass2, hi2);
    amrex::Real dYdx[NUM_SPECIES], ddrive[NUM_SPECIES];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      dYdx[ns] = dxinv * (mass1[ns] - mass2[ns]);
      ddrive[ns] = 0.5 * (diP1[ns] + diP2[ns]) * dpdx;
    }
    for (int ns1 = 0; ns1 < NUM_SPECIES; ++ns1) {
      for (int ns2 = 0; ns2 < NUM_SPECIES; ++ns2) {
        ddrive[ns1] += 0.5 * (dijY1[ns1][ns2] + dijY2[ns1][ns2]) * dYdx[ns2];
      }
    }
    amrex::Real dsum = 0.0;
    for (amrex::Real ns : ddrive) {
      dsum += ns;
    }
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
      const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
      ddrive[ns] -= Yface * dsum;
      const amrex::Real Vd = -coef[dComp_rhoD + ns] * ddrive[ns];
      flx(iv, UFS + ns) = Vd;
      Vc += Vd;
      flx(iv, UEDEN) += Vd * hface;
    }

    // Add correction velocity to fluxes
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
      const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
      flx(iv, UFS + ns) -= Yface * Vc;
      flx(iv, UEDEN) -= Yface * hface * Vc;
    }
  }
};

struct FluxTypes
{
  using SpeciesEnergyFluxType = SpeciesEnergyFlux<pele::physics::EosType>;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_diffusion_flux_eb(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::GpuArray<amrex::Real, dComp_lambda + 1>& coef,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const amrex::Array4<const amrex::Real>& area,
  const amrex::Array4<amrex::Real>& flx,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
  const int dir)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
  amrex::Real taux = 0.0, tauy = 0.0, tauz = 0.0;
  const bool update = flags(iv).isRegular() || flags(iv).isSingleValued();

  if (update) {
    if (dir == 0) {
      AMREX_D_TERM(
        const amrex::Real dudx = dxinv[dir] * (q(iv, QU) - q(ivm, QU));
        , const amrex::Real dvdx = dxinv[dir] * (q(iv, QV) - q(ivm, QV));
        , const amrex::Real dwdx = dxinv[dir] * (q(iv, QW) - q(ivm, QW)););

      const int jhip = j + static_cast<int>(flags(iv).isConnected(0, 1, 0));
      const int jhim = j - static_cast<int>(flags(iv).isConnected(0, -1, 0));
      const int jlop = j + static_cast<int>(flags(ivm).isConnected(0, 1, 0));
      const int jlom = j - static_cast<int>(flags(ivm).isConnected(0, -1, 0));
      const amrex::Real wjhi = weights[jhip - jhim];
      const amrex::Real wjlo = weights[jlop - jlom];

      const amrex::Real dudy =
        (0.5 * dxinv[1]) *
        ((q(i, jhip, k, QU) - q(i, jhim, k, QU)) * wjhi +
         (q(i - 1, jlop, k, QU) - q(i - 1, jlom, k, QU)) * wjlo);

      const amrex::Real dvdy =
        (0.5 * dxinv[1]) *
        ((q(i, jhip, k, QV) - q(i, jhim, k, QV)) * wjhi +
         (q(i - 1, jlop, k, QV) - q(i - 1, jlom, k, QV)) * wjlo);

#if AMREX_SPACEDIM == 3
      const int khip = k + static_cast<int>(flags(iv).isConnected(0, 0, 1));
      const int khim = k - static_cast<int>(flags(iv).isConnected(0, 0, -1));
      const int klop = k + static_cast<int>(flags(ivm).isConnected(0, 0, 1));
      const int klom = k - static_cast<int>(flags(ivm).isConnected(0, 0, -1));
      const amrex::Real wkhi = weights[khip - khim];
      const amrex::Real wklo = weights[klop - klom];

      const amrex::Real dudz =
        (0.5 * dxinv[2]) *
        ((q(i, j, khip, QU) - q(i, j, khim, QU)) * wkhi +
         (q(i - 1, j, klop, QU) - q(i - 1, j, klom, QU)) * wklo);
      const amrex::Real dwdz =
        (0.5 * dxinv[2]) *
        ((q(i, j, khip, QW) - q(i, j, khim, QW)) * wkhi +
         (q(i - 1, j, klop, QW) - q(i - 1, j, klom, QW)) * wklo);
#endif
      const amrex::Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

      AMREX_D_TERM(
        taux = coef[dComp_mu] * (2.0 * dudx - 2.0 / 3.0 * divu) +
               coef[dComp_xi] * divu;
        , tauy = coef[dComp_mu] * (dudy + dvdx);
        , tauz = coef[dComp_mu] * (dudz + dwdx););
    } else if (dir == 1) {
      AMREX_D_TERM(
        const amrex::Real dudy = dxinv[dir] * (q(iv, QU) - q(ivm, QU));
        , const amrex::Real dvdy = dxinv[dir] * (q(iv, QV) - q(ivm, QV));
        , const amrex::Real dwdy = dxinv[dir] * (q(iv, QW) - q(ivm, QW));)

      const int ihip = i + static_cast<int>(flags(iv).isConnected(1, 0, 0));
      const int ihim = i - static_cast<int>(flags(iv).isConnected(-1, 0, 0));
      const int ilop = i + static_cast<int>(flags(ivm).isConnected(1, 0, 0));
      const int ilom = i - static_cast<int>(flags(ivm).isConnected(-1, 0, 0));
      const amrex::Real wihi = weights[ihip - ihim];
      const amrex::Real wilo = weights[ilop - ilom];

      const amrex::Real dudx =
        (0.5 * dxinv[0]) *
        ((q(ihip, j, k, QU) - q(ihim, j, k, QU)) * wihi +
         (q(ilop, j - 1, k, QU) - q(ilom, j - 1, k, QU)) * wilo);
      const amrex::Real dvdx =
        (0.5 * dxinv[0]) *
        ((q(ihip, j, k, QV) - q(ihim, j, k, QV)) * wihi +
         (q(ilop, j - 1, k, QV) - q(ilom, j - 1, k, QV)) * wilo);

#if AMREX_SPACEDIM == 3
      const int khip = k + static_cast<int>(flags(iv).isConnected(0, 0, 1));
      const int khim = k - static_cast<int>(flags(iv).isConnected(0, 0, -1));
      const int klop = k + static_cast<int>(flags(ivm).isConnected(0, 0, 1));
      const int klom = k - static_cast<int>(flags(ivm).isConnected(0, 0, -1));
      const amrex::Real wkhi = weights[khip - khim];
      const amrex::Real wklo = weights[klop - klom];

      const amrex::Real dvdz =
        (0.5 * dxinv[2]) *
        ((q(i, j, khip, QV) - q(i, j, khim, QV)) * wkhi +
         (q(i, j - 1, klop, QV) - q(i, j - 1, klom, QV)) * wklo);
      const amrex::Real dwdz =
        (0.5 * dxinv[2]) *
        ((q(i, j, khip, QW) - q(i, j, khim, QW)) * wkhi +
         (q(i, j - 1, klop, QW) - q(i, j - 1, klom, QW)) * wklo);
#endif

      const amrex::Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);
      AMREX_D_TERM(
        taux = coef[dComp_mu] * (dudy + dvdx);
        , tauy = coef[dComp_mu] * (2.0 * dvdy - 2.0 / 3.0 * divu) +
                 coef[dComp_xi] * divu;
        , tauz = coef[dComp_mu] * (dwdy + dvdz);)
    } else if (dir == 2) {
#if AMREX_SPACEDIM == 3
      const amrex::Real dudz = dxinv[dir] * (q(iv, QU) - q(ivm, QU));
      const amrex::Real dvdz = dxinv[dir] * (q(iv, QV) - q(ivm, QV));
      const amrex::Real dwdz = dxinv[dir] * (q(iv, QW) - q(ivm, QW));

      const int ihip = i + static_cast<int>(flags(iv).isConnected(1, 0, 0));
      const int ihim = i - static_cast<int>(flags(iv).isConnected(-1, 0, 0));
      const int ilop = i + static_cast<int>(flags(ivm).isConnected(1, 0, 0));
      const int ilom = i - static_cast<int>(flags(ivm).isConnected(-1, 0, 0));
      const amrex::Real wihi = weights[ihip - ihim];
      const amrex::Real wilo = weights[ilop - ilom];

      const amrex::Real dudx =
        (0.5 * dxinv[0]) *
        ((q(ihip, j, k, QU) - q(ihim, j, k, QU)) * wihi +
         (q(ilop, j, k - 1, QU) - q(ilom, j, k - 1, QU)) * wilo);
      const amrex::Real dwdx =
        (0.5 * dxinv[0]) *
        ((q(ihip, j, k, QW) - q(ihim, j, k, QW)) * wihi +
         (q(ilop, j, k - 1, QW) - q(ilom, j, k - 1, QW)) * wilo);

      const int jhip = j + static_cast<int>(flags(iv).isConnected(0, 1, 0));
      const int jhim = j - static_cast<int>(flags(iv).isConnected(0, -1, 0));
      const int jlop = j + static_cast<int>(flags(ivm).isConnected(0, 1, 0));
      const int jlom = j - static_cast<int>(flags(ivm).isConnected(0, -1, 0));
      const amrex::Real wjhi = weights[jhip - jhim];
      const amrex::Real wjlo = weights[jlop - jlom];

      const amrex::Real dvdy =
        (0.5 * dxinv[1]) *
        ((q(i, jhip, k, QV) - q(i, jhim, k, QV)) * wjhi +
         (q(i, jlop, k - 1, QV) - q(i, jlom, k - 1, QV)) * wjlo);
      const amrex::Real dwdy =
        (0.5 * dxinv[1]) *
        ((q(i, jhip, k, QW) - q(i, jhim, k, QW)) * wjhi +
         (q(i, jlop, k - 1, QW) - q(i, jlom, k - 1, QW)) * wjlo);

      const amrex::Real divu = dudx + dvdy + dwdz;
      taux = coef[dComp_mu] * (dudz + dwdx);
      tauy = coef[dComp_mu] * (dvdz + dwdy);
      tauz = coef[dComp_mu] * (2.0 * dwdz - 2.0 / 3.0 * divu) +
             coef[dComp_xi] * divu;
#endif
    }
  }
  const amrex::Real dTdd =
    update ? dxinv[dir] * (q(iv, QTEMP) - q(ivm, QTEMP)) : 0.0;

  flx(iv, UMX) = -taux;
  flx(iv, UMY) = -tauy;
  flx(iv, UMZ) = -tauz;
  flx(iv, UEDEN) =
    0.5 * (AMREX_D_TERM(
            -taux * (q(iv, QU) + q(ivm, QU)), -tauy * (q(iv, QV) + q(ivm, QV)),
            -tauz * (q(iv, QW) + q(ivm, QW)))) -
    coef[dComp_lambda] * dTdd;

  if (update) {
    FluxTypes::SpeciesEnergyFluxType()(iv, ivm, dxinv[dir], coef, q, flx);
  }

  // Scale by area
  AMREX_D_TERM(flx(iv, UMX) *= area(iv);, flx(iv, UMY) *= area(iv);
               , flx(iv, UMZ) *= area(iv););
  flx(iv, UEDEN) *= area(iv);
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    flx(iv, UFS + ns) *= area(iv);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_diffusion_flux(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::GpuArray<amrex::Real, dComp_lambda + 1>& coef,
  const amrex::Array4<const amrex::Real>& area,
  const amrex::Array4<amrex::Real>& flx,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
  const int dir)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
  amrex::Real taux = 0.0, tauy = 0.0, tauz = 0.0;

  if (dir == 0) {
    AMREX_D_TERM(
      const amrex::Real dudx = dxinv[dir] * (q(iv, QU) - q(ivm, QU));
      , const amrex::Real dvdx = dxinv[dir] * (q(iv, QV) - q(ivm, QV));
      , const amrex::Real dwdx = dxinv[dir] * (q(iv, QW) - q(ivm, QW)););
    const amrex::Real dudy = (q(i, j + 1, k, QU) + q(i - 1, j + 1, k, QU) -
                              q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU)) *
                             (0.25 * dxinv[1]);
    const amrex::Real dvdy = (q(i, j + 1, k, QV) + q(i - 1, j + 1, k, QV) -
                              q(i, j - 1, k, QV) - q(i - 1, j - 1, k, QV)) *
                             (0.25 * dxinv[1]);

#if (AMREX_SPACEDIM == 3)
    const amrex::Real dudz = (q(i, j, k + 1, QU) + q(i - 1, j, k + 1, QU) -
                              q(i, j, k - 1, QU) - q(i - 1, j, k - 1, QU)) *
                             (0.25 * dxinv[2]);
    const amrex::Real dwdz = (q(i, j, k + 1, QW) + q(i - 1, j, k + 1, QW) -
                              q(i, j, k - 1, QW) - q(i - 1, j, k - 1, QW)) *
                             (0.250 * dxinv[2]);
#endif
    const amrex::Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

    taux =
      coef[dComp_mu] * (2.0 * dudx - 2.0 / 3.0 * divu) + coef[dComp_xi] * divu;
    AMREX_D_TERM(
      , tauy = coef[dComp_mu] * (dudy + dvdx);
      , tauz = coef[dComp_mu] * (dudz + dwdx););
  } else if (dir == 1) {

    AMREX_D_TERM(
      const amrex::Real dudy = dxinv[dir] * (q(iv, QU) - q(ivm, QU));
      , const amrex::Real dvdy = dxinv[dir] * (q(iv, QV) - q(ivm, QV));
      , const amrex::Real dwdy = dxinv[dir] * (q(iv, QW) - q(ivm, QW)););

    const amrex::Real dudx = (q(i + 1, j, k, QU) + q(i + 1, j - 1, k, QU) -
                              q(i - 1, j, k, QU) - q(i - 1, j - 1, k, QU)) *
                             (0.25 * dxinv[0]);
    const amrex::Real dvdx = (q(i + 1, j, k, QV) + q(i + 1, j - 1, k, QV) -
                              q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV)) *
                             (0.25 * dxinv[0]);
#if (AMREX_SPACEDIM == 3)
    const amrex::Real dvdz = (q(i, j, k + 1, QV) + q(i, j - 1, k + 1, QV) -
                              q(i, j, k - 1, QV) - q(i, j - 1, k - 1, QV)) *
                             (0.25 * dxinv[2]);
    const amrex::Real dwdz = (q(i, j, k + 1, QW) + q(i, j - 1, k + 1, QW) -
                              q(i, j, k - 1, QW) - q(i, j - 1, k - 1, QW)) *
                             (0.25 * dxinv[2]);
#endif

    const amrex::Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

    taux = coef[dComp_mu] * (dudy + dvdx);
    tauy =
      coef[dComp_mu] * (2.0 * dvdy - 2.0 / 3.0 * divu) + coef[dComp_xi] * divu;
    AMREX_D_TERM(, , tauz = coef[dComp_mu] * (dwdy + dvdz););
  } else if (dir == 2) {
#if AMREX_SPACEDIM == 3
    const amrex::Real dudz = dxinv[dir] * (q(iv, QU) - q(ivm, QU));
    const amrex::Real dvdz = dxinv[dir] * (q(iv, QV) - q(ivm, QV));
    const amrex::Real dwdz = dxinv[dir] * (q(iv, QW) - q(ivm, QW));
    const amrex::Real dudx = (q(i + 1, j, k, QU) + q(i + 1, j, k - 1, QU) -
                              q(i - 1, j, k, QU) - q(i - 1, j, k - 1, QU)) *
                             (0.25 * dxinv[0]);
    const amrex::Real dwdx = (q(i + 1, j, k, QW) + q(i + 1, j, k - 1, QW) -
                              q(i - 1, j, k, QW) - q(i - 1, j, k - 1, QW)) *
                             (0.25 * dxinv[0]);
    const amrex::Real dvdy = (q(i, j + 1, k, QV) + q(i, j + 1, k - 1, QV) -
                              q(i, j - 1, k, QV) - q(i, j - 1, k - 1, QV)) *
                             (0.25 * dxinv[1]);
    const amrex::Real dwdy = (q(i, j + 1, k, QW) + q(i, j + 1, k - 1, QW) -
                              q(i, j - 1, k, QW) - q(i, j - 1, k - 1, QW)) *
                             (0.25 * dxinv[1]);
    const amrex::Real divu = dudx + dvdy + dwdz;

    taux = coef[dComp_mu] * (dudz + dwdx);
    tauy = coef[dComp_mu] * (dvdz + dwdy);
    tauz =
      coef[dComp_mu] * (2.0 * dwdz - 2.0 / 3.0 * divu) + coef[dComp_xi] * divu;
#endif
  }
  const amrex::Real dTdd = dxinv[dir] * (q(iv, QTEMP) - q(ivm, QTEMP));

  flx(iv, UMX) = -taux;
  flx(iv, UMY) = -tauy;
  flx(iv, UMZ) = -tauz;
  flx(iv, UEDEN) =
    0.5 * (AMREX_D_TERM(
            -taux * (q(iv, QU) + q(ivm, QU)), -tauy * (q(iv, QV) + q(ivm, QV)),
            -tauz * (q(iv, QW) + q(ivm, QW)))) -
    coef[dComp_lambda] * dTdd;

  FluxTypes::SpeciesEnergyFluxType()(iv, ivm, dxinv[dir], coef, q, flx);

  // Scale by area
  AMREX_D_TERM(flx(iv, UMX) *= area(iv);, flx(iv, UMY) *= area(iv);
               , flx(iv, UMZ) *= area(iv););
  flx(iv, UEDEN) *= area(iv);
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    flx(iv, UFS + ns) *= area(iv);
  }
}

// This function computes the flux divergence.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_flux_div(
  const int i,
  const int j,
  const int k,
  const int n,
  AMREX_D_DECL(
    const amrex::Array4<const amrex::Real>& fx,
    const amrex::Array4<const amrex::Real>& fy,
    const amrex::Array4<const amrex::Real>& fz),
  const amrex::Array4<const amrex::Real>& V,
  const amrex::Array4<amrex::Real>& D)
{
  D(i, j, k, n) =
    -(AMREX_D_TERM(
      fx(i + 1, j, k, n) - fx(i, j, k, n), +fy(i, j + 1, k, n) - fy(i, j, k, n),
      +fz(i, j, k + 1, n) - fz(i, j, k, n))) /
    V(i, j, k);
}

// This function computes diffusive energy flux through isothermal walls
// and *adds* it to the flx array
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_isothermal_wall_fluxes(
  const int i,
  const int j,
  const int k,
  const int dir,
  const int normal,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& Twall_arr,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const amrex::Array4<const amrex::Real>& area,
  const amrex::Array4<amrex::Real>& flx,
  const amrex::GeometryData& geomdata,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::PhysicsType::transport_type> const* tparm,
  ProbParmDevice const& prob_parm)
{

  // iv: face centered index where flux is computed
  amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  // ivm: cell centered index of interior cell
  amrex::IntVect ivm{AMREX_D_DECL(i, j, k)};
  ivm[dir] += ((normal > 0) ? -1 : 0);
  // ivr: boundary cell relative to interior cell
  amrex::IntVect ivr{0};
  ivr[dir] = normal;

  // Only do anything if the face isn't covered by an EB
  if (flags(ivm).isConnected(ivr)) {

    // First get density at wall temp (only comes into play for SRK)
    amrex::Real rho_wall = 0.0;
    const amrex::Real Twall = Twall_arr(iv);
    const amrex::Real pwall = q(ivm, QPRES);
    amrex::Real Ywall[NUM_SPECIES] = {0.0};
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      Ywall[ns] = q(ivm, ns + QFS);
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(pwall, Ywall, Twall, rho_wall);

    // Now compute lambda at the wall temperature
    const bool wtr_get_xi = false;
    const bool wtr_get_mu = false;
    const bool wtr_get_lam = true;
    const bool wtr_get_Ddiag = false;
    const bool wtr_get_chi = false;
    amrex::Real* dummy_Ddiag = nullptr;
    amrex::Real* dummy_chi_mix = nullptr;
    amrex::Real dummy_mu, dummy_xi;
    amrex::Real lambda;

    const amrex::RealVect x = pc_cmp_loc(ivm, geomdata);
    pc_transcoeff(
      wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, Twall,
      rho_wall, Ywall, dummy_Ddiag, dummy_chi_mix, dummy_mu, dummy_xi, lambda,
      tparm, prob_parm, x);

    // Compute Fourier flux and scale by area
    amrex::Real dTdx =
      2.0 * (q(ivm, QTEMP) - Twall) / geomdata.CellSize(dir) * normal;
    flx(i, j, k, UEDEN) += lambda * dTdx * area(iv);
  }
}

#endif
