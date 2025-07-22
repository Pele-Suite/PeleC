#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

#include "PelePhysics.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PeleCUtilities.H"
#include "TransCoeff.H"

// EstDt routines

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pc_trans4dt(
  const int which_trans,
  amrex::Real T,
  amrex::Real rho,
  amrex::Real massfrac[],
  amrex::Real& D,
  const amrex::RealVect& x,
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type> const* trans_parm,
  ProbParmDevice const& prob_parm)
{
  bool get_xi = false, get_mu = false, get_lam = false, get_Ddiag = false,
       get_chi = false;
  amrex::Real dum1 = 0., dum2 = 0.;

  if (which_trans == 0) {
    get_mu = true;
    pc_transcoeff(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac, nullptr,
      nullptr, D, dum1, dum2, trans_parm, prob_parm, x);
  } else if (which_trans == 1) {
    get_lam = true;
    pc_transcoeff(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac, nullptr,
      nullptr, dum1, dum2, D, trans_parm, prob_parm, x);
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_hydro(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::Array4<const amrex::EBCellFlag>& flags,
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept
{
  amrex::Real dt = std::numeric_limits<amrex::Real>::max();

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
    if (!flags(i, j, k).isCovered()) {
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real massfrac[NUM_SPECIES];
      amrex::Real c;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = u(i, j, k, UFS + n) * rhoInv;
      }
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2Cs(rho, T, massfrac, c);
      AMREX_D_TERM(
        const amrex::Real ux = u(i, j, k, UMX) * rhoInv;
        const amrex::Real dt1 = dx / (c + std::abs(ux));
        dt = amrex::min<amrex::Real>(dt, dt1);
        , const amrex::Real uy = u(i, j, k, UMY) * rhoInv;
        const amrex::Real dt2 = dy / (c + std::abs(uy));
        dt = amrex::min<amrex::Real>(dt, dt2);
        , const amrex::Real uz = u(i, j, k, UMZ) * rhoInv;
        const amrex::Real dt3 = dz / (c + std::abs(uz));
        dt = amrex::min<amrex::Real>(dt, dt3););
    }
  });
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0, "ERROR: dt needs to be positive.");
  return dt;
}

// Diffusion Velocity
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_veldif(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::Array4<const amrex::EBCellFlag>& flags,
  const amrex::GeometryData& geomdata,
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type> const* trans_parm,
  ProbParmDevice const& prob_parm) noexcept
{
  amrex::Real dt = std::numeric_limits<amrex::Real>::max();

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
    if (!flags(i, j, k).isCovered()) {
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = u(i, j, k, n + UFS) * rhoInv;
      }
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real D = 0.0;
      const int which_trans = 0;
      const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, geomdata);
      pc_trans4dt(which_trans, T, rho, massfrac, D, x, trans_parm, prob_parm);
      D *= rhoInv;
      if (D == 0.0) {
        D = constants::small_num();
      }
      AMREX_D_TERM(
        const amrex::Real dt1 = 0.5 * geomdata.CellSize(0) *
                                geomdata.CellSize(0) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt1);
        , const amrex::Real dt2 = 0.5 * geomdata.CellSize(1) *
                                  geomdata.CellSize(1) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt2);
        , const amrex::Real dt3 = 0.5 * geomdata.CellSize(2) *
                                  geomdata.CellSize(2) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt3););
    }
  });
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0, "ERROR: dt needs to be positive.");
  return dt;
}

// Diffusion Temperature
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_tempdif(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::Array4<const amrex::EBCellFlag>& flags,
  const amrex::GeometryData& geomdata,
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type> const* trans_parm,
  ProbParmDevice const& prob_parm) noexcept
{
  amrex::Real dt = std::numeric_limits<amrex::Real>::max();

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
    if (!flags(i, j, k).isCovered()) {
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = u(i, j, k, n + UFS) * rhoInv;
      }
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real D = 0.0;
      const int which_trans = 1;
      const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, geomdata);
      pc_trans4dt(which_trans, T, rho, massfrac, D, x, trans_parm, prob_parm);
      amrex::Real cv;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2Cv(rho, T, massfrac, cv);
      D *= rhoInv / cv;
      if (D == 0.0) {
        D = constants::small_num();
      }
      AMREX_D_TERM(
        const amrex::Real dt1 = 0.5 * geomdata.CellSize(0) *
                                geomdata.CellSize(0) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt1);
        , const amrex::Real dt2 = 0.5 * geomdata.CellSize(1) *
                                  geomdata.CellSize(1) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt2);
        , const amrex::Real dt3 = 0.5 * geomdata.CellSize(2) *
                                  geomdata.CellSize(2) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt3););
    }
  });
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0, "ERROR: dt needs to be positive.");
  return dt;
}

// Diffusion Enthalpy
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
pc_estdt_enthdif(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::Array4<const amrex::EBCellFlag>& flags,
  const amrex::GeometryData& geomdata,
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type> const* trans_parm,
  ProbParmDevice const& prob_parm) noexcept
{
  amrex::Real dt = std::numeric_limits<amrex::Real>::max();

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
    if (!flags(i, j, k).isCovered()) {
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = u(i, j, k, n + UFS) * rhoInv;
      }
      amrex::Real T = u(i, j, k, UTEMP);
      amrex::Real cp;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2Cp(rho, T, massfrac, cp);
      amrex::Real D;
      const int which_trans = 1;
      const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, geomdata);
      pc_trans4dt(which_trans, T, rho, massfrac, D, x, trans_parm, prob_parm);
      D *= rhoInv / cp;
      AMREX_D_TERM(
        const amrex::Real dt1 = 0.5 * geomdata.CellSize(0) *
                                geomdata.CellSize(0) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt1);
        , const amrex::Real dt2 = 0.5 * geomdata.CellSize(1) *
                                  geomdata.CellSize(1) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt2);
        , const amrex::Real dt3 = 0.5 * geomdata.CellSize(2) *
                                  geomdata.CellSize(2) / (AMREX_SPACEDIM * D);
        dt = amrex::min<amrex::Real>(dt, dt3););
    }
  });
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0, "ERROR: dt needs to be positive.");
  return dt;
}

#endif
