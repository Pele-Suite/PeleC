#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"

#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "PeleCUtilities.H"
#include "turbinflow.H"

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
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.Y_ox[n];
  }

  amrex::Real x_velocity = 0.0;
  amrex::Real y_velocity = 0.0;
  amrex::Real z_velocity = 0.0;

  amrex::Real rho = 0.0, eint = 0;
  eos.PYT2R(prob_parm.P_mean, massfrac, prob_parm.T_ox, rho);
  eos.RTY2E(rho, prob_parm.T_ox, massfrac, eint);

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * x_velocity;
  state(i, j, k, UMY) = rho * y_velocity;
  state(i, j, k, UMZ) = rho * z_velocity;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (x_velocity * x_velocity + y_velocity * y_velocity +
                         z_velocity * z_velocity));
  state(i, j, k, UTEMP) = prob_parm.T_ox;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_eta(
  amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau)
{
  amrex::Real eta_time;
  eta_time = 0.5 * (std::tanh((time - start) / (duration * tau)) *
                      std::tanh((-time + start + duration) / (duration * tau)) +
                    1.0);
  eta_time =
    amrex::min<amrex::Real>(amrex::max<amrex::Real>(eta_time, 0.0), 1.0);
  return eta_time;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& turb_fluc)
{
  if (sgn == 1 && idir == 2) {

    // Save off the inflow velocity if specified
    amrex::Real velx = 0.0, vely = 0.0, velz = 0.0;
    if (prob_parm.turbulence) {
      velx = turb_fluc[0];
      vely = turb_fluc[1];
      velz = turb_fluc[2] + prob_parm.v_in;
    } else {
      velz += prob_parm.v_in;
    }

    // default noslip wall condition is prepopulated in s_ext
    // only need to edit for inflows

    const amrex::Real radius = std::sqrt(
      (x[0] - prob_parm.center_xy[0]) * (x[0] - prob_parm.center_xy[0]) +
      (x[1] - prob_parm.center_xy[1]) * (x[1] - prob_parm.center_xy[1]));

    if (radius < 1.1 * 0.5 * prob_parm.D) {
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      auto eos = pele::physics::PhysicsType::eos();

      // interior pressure
      amrex::Real p_int;
      amrex::Real massfrac_int[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_int[n] = s_int[UFS + n] / s_int[URHO];
      }
      eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac_int, p_int);

      amrex::Real eta = 0.5 * (1.0 - std::tanh(
                                       2.0 * (radius - prob_parm.D * 0.5) /
                                       (0.1 * prob_parm.D * 0.5)));
      eta *=
        time_eta(time, prob_parm.inj_start, prob_parm.inj_dur, prob_parm.tau);
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] =
          eta * prob_parm.Y_fuel[n] + (1.0 - eta) * prob_parm.Y_ox[n];
      }
      amrex::Real int_energy =
        eta * prob_parm.U_fuel + (1.0 - eta) * prob_parm.U_ox;

      s_ext[UTEMP] = 500.0; // initial guess to help EY2T converge
      eos.EY2T(int_energy, massfrac, s_ext[UTEMP]);
      //    GET_T_GIVEN_HY(&enthalpy,massfrac,&s_ext[DEF_Temp],&ierr);

      amrex::Real rho, eint;
      eos.PYT2RE(p_int, massfrac, s_ext[UTEMP], rho, eint);

      s_ext[URHO] = rho;
      s_ext[UEINT] = rho * eint;

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[UFS + n] = rho * massfrac[n];
      }
      // Velocity
      s_ext[UMX] = rho * velx * eta;
      s_ext[UMY] = rho * vely * eta;
      s_ext[UMZ] = rho * velz * eta;
      s_ext[UEDEN] = s_ext[UEINT] + rho * eta * eta * 0.5 *
                                      (velx * velx + vely * vely + velz * velz);
    }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
