#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "Utilities.H"
#include "PeleCUtilities.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real x[3] = {
    prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
    prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
    prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

  // Uniform density, temperature, pressure (internal energy)
  state(i, j, k, URHO) = prob_parm.rho0;
  state(i, j, k, UEINT) = prob_parm.rho0 * prob_parm.eint0;
  state(i, j, k, UTEMP) = prob_parm.T0;
  const amrex::Real massfrac[NUM_SPECIES] = {1.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = prob_parm.rho0 * massfrac[n];
  }

  // Fill in the velocities and energy.
  amrex::Real u[3] = {0.0};
  amrex::Real uinterp[3] = {0.0};

  // Interpolation factors
  amrex::Real mod[3] = {0.0};
  int idx[3] = {0};
  int idxp1[3] = {0};
  amrex::Real slp[3] = {0.0};
  for (int cnt = 0; cnt < 3; cnt++) {
    mod[cnt] = std::fmod(x[cnt], prob_parm.Linput);
    pele::physics::utilities::locate(
      prob_parm.d_xarray, prob_parm.inres, mod[cnt], idx[cnt]);
    idxp1[cnt] = (idx[cnt] + 1) % prob_parm.inres;
    slp[cnt] =
      (mod[cnt] - prob_parm.d_xarray[idx[cnt]]) / prob_parm.d_xdiff[idx[cnt]];
  }

  const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
  const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
  const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
  const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
  const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
  const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
  const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
  const amrex::Real f7 = slp[0] * slp[1] * slp[2];

  uinterp[0] =
    prob_parm.d_uinput
        [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
      f0 +
    prob_parm.d_uinput
        [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
      f1 +
    prob_parm.d_uinput
        [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
      f2 +
    prob_parm.d_uinput
        [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
      f3 +
    prob_parm.d_uinput
        [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
      f4 +
    prob_parm.d_uinput
        [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
      f5 +
    prob_parm.d_uinput
        [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
      f6 +
    prob_parm.d_uinput
        [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
      f7;
  uinterp[1] =
    prob_parm.d_vinput
        [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
      f0 +
    prob_parm.d_vinput
        [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
      f1 +
    prob_parm.d_vinput
        [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
      f2 +
    prob_parm.d_vinput
        [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
      f3 +
    prob_parm.d_vinput
        [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
      f4 +
    prob_parm.d_vinput
        [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
      f5 +
    prob_parm.d_vinput
        [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
      f6 +
    prob_parm.d_vinput
        [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
      f7;
  uinterp[2] =
    prob_parm.d_winput
        [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
      f0 +
    prob_parm.d_winput
        [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
      f1 +
    prob_parm.d_winput
        [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
      f2 +
    prob_parm.d_winput
        [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
      f3 +
    prob_parm.d_winput
        [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
      f4 +
    prob_parm.d_winput
        [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
      f5 +
    prob_parm.d_winput
        [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
      f6 +
    prob_parm.d_winput
        [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
      f7;

  u[0] = uinterp[0] + prob_parm.forcing_u0;
  u[1] = uinterp[1] + prob_parm.forcing_v0;
  u[2] = uinterp[2] + prob_parm.forcing_w0;
  state(i, j, k, UMX) = prob_parm.rho0 * u[0];
  state(i, j, k, UMY) = prob_parm.rho0 * u[1];
  state(i, j, k, UMZ) = prob_parm.rho0 * u[2];
  state(i, j, k, UEDEN) =
    prob_parm.rho0 *
    (prob_parm.eint0 + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
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

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
