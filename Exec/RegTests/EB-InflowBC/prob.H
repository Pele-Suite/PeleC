#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "PeleCUtilities.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = 0.0;
  }

  // Set the states
  state(i, j, k, URHO) = prob_parm.rho_l;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEDEN) = prob_parm.rhoe_l;
  state(i, j, k, UEINT) = prob_parm.rhoe_l;
  state(i, j, k, UTEMP) = prob_parm.T_l;
  state(i, j, k, UFS + prob_parm.left_gas_id) = state(i, j, k, URHO);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
}

void pc_prob_close();

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static bool problem_eb_state(
    amrex::GeometryData const& geomdata,
    const amrex::Real /*vf*/,
    const amrex::IntVect iv,
    AMREX_D_DECL(
      const amrex::Real /*nx*/,
      const amrex::Real /*ny*/,
      const amrex::Real /*nz*/),
    const amrex::Real ql[5 + NUM_SPECIES],
    const amrex::Real* /*spl[NUM_SPECIES]*/,
    const amrex::Real /*rhoe_l*/,
    const amrex::Real /*gamc_l*/,
    amrex::Real qr[5 + NUM_SPECIES],
    amrex::Real spr[NUM_SPECIES],
    amrex::Real& rhoe_r,
    amrex::Real& gamc_r,
    ProbParmDevice const* prob_parm)
  {
    bool do_ebfill = false;
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const amrex::Real max_r = 0.1;
    const amrex::Real radius =
      std::sqrt(AMREX_D_TERM(, x[1] * x[1], +x[2] * x[2]));
    if (radius < max_r) {
      qr[R_RHO] = prob_parm->rho_r;
      qr[R_UN] = -prob_parm->U_r;
      qr[R_UT1] = ql[R_UT1];
      qr[R_UT2] = ql[R_UT2];
      for (int n = 0; n < NUM_SPECIES; n++) {
        spr[n] = 0.0;
      }
      spr[prob_parm->right_gas_id] = 1.0;

      amrex::Real eos_state_rho = qr[R_RHO];
      amrex::Real eos_state_p = qr[R_P];
      amrex::Real eos_state_e;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RYP2E(eos_state_rho, spr, eos_state_p, eos_state_e);
      rhoe_r = eos_state_rho * eos_state_e;
      amrex::Real eos_state_T;
      eos.RYP2T(eos_state_rho, spr, eos_state_p, eos_state_T);
      eos.TY2G(eos_state_T, spr, gamc_r);
      do_ebfill = true;
    }

    return do_ebfill;
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;

#endif
