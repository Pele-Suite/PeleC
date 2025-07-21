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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real Cs;
  massfrac[O2_ID] = prob_parm.Yox_amb;
  massfrac[N2_ID] = 1.0 - prob_parm.Yox_amb;

  amrex::Real rho = 0.0, eint = 0;
  eos.PYT2R(prob_parm.Pres_amb, massfrac, prob_parm.Temp_amb, rho);
  eos.RTY2E(rho, prob_parm.Temp_amb, massfrac, eint);
  eos.RTY2Cs(rho, prob_parm.Temp_amb, massfrac, Cs);

  // Init options
  // 0 : uniform x_vel everywhere
  // 1 : zero velocity in cavity region
  // 2 : linear velocity profile in cavity region
  amrex::Real x_velocity = 0.0;
  amrex::Real y_velocity = 0.0;
  amrex::Real z_velocity = 0.0;
  const amrex::Real cavity_height = prob_parm.cavity_depth;
  if (prob_parm.init_type == 0) {
    x_velocity = Cs * prob_parm.Mach_channel;
    // Do nothing - constant velocity
  } else if (prob_parm.init_type == 1) {
    // Zero velocity in cavity
    amrex::Real yval =
      geomdata.ProbLo(1) + (amrex::Real(j) + 0.5) * geomdata.CellSize(1);
    x_velocity = yval <= cavity_height ? 0.0 : Cs * prob_parm.Mach_channel;
  } else if (prob_parm.init_type == 2) {
    // Velocity gradient in cavity
    amrex::Real yval =
      geomdata.ProbLo(1) + (amrex::Real(j) + 0.5) * geomdata.CellSize(1);
    amrex::Real interp = amrex::min(
      amrex::max((yval - 0.5 * cavity_height) / (0.5 * cavity_height), -1.0),
      1.0);
    x_velocity = interp * Cs * prob_parm.Mach_channel;
  } else {
    amrex::Abort(
      "prob.init_type must be 0 (uniform velocity), 1 (quiescent "
      "cavity), or 2 (linear velocity in cavity)");
  }

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * x_velocity;
  state(i, j, k, UMY) = rho * y_velocity;
  state(i, j, k, UMZ) = rho * z_velocity;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (x_velocity * x_velocity + y_velocity * y_velocity +
                         z_velocity * z_velocity));
  state(i, j, k, UTEMP) = prob_parm.Temp_amb;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
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
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{

  if (idir == 1) {

    // default noslip wall condition is prepopulated in s_ext
    // only need to edit for inflows

    if (prob_parm.inject_fuel) {
      amrex::Real hole_cx = prob_parm.centx;
      amrex::Real hole_cz = prob_parm.centz;

      const amrex::Real radius = std::sqrt(
        (x[0] - hole_cx) * (x[0] - hole_cx) +
        (x[2] - hole_cz) * (x[2] - hole_cz));

      if (radius < prob_parm.r_hole) {
        auto eos = pele::physics::PhysicsType::eos();

        // interior pressure
        amrex::Real p_int;
        amrex::Real massfrac_int[NUM_SPECIES] = {0.0};
        amrex::Real massfrac_ext[NUM_SPECIES] = {0.0};
        for (int n = 0; n < NUM_SPECIES; n++) {
          massfrac_int[n] = s_int[UFS + n] / s_int[URHO];
        }
        eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac_int, p_int);

        massfrac_ext[H2_ID] = prob_parm.Yfuel_jet;
        massfrac_ext[O2_ID] =
          (1.0 - massfrac_ext[H2_ID]) * 0.23; // making the rest air
        massfrac_ext[N2_ID] = 1.0 - massfrac_ext[H2_ID] - massfrac_ext[O2_ID];
        amrex::Real vjet_pipeflow =
          prob_parm.vel_jet * (1.0 - std::pow(radius / prob_parm.r_hole, 2.0));
        amrex::Real vx_in = 0.0;
        amrex::Real vy_in = vjet_pipeflow;
        amrex::Real vz_in = 0.0;

        s_ext[UTEMP] = prob_parm.temp_jet;
        amrex::Real rho, eint;
        eos.PYT2RE(p_int, massfrac_ext, s_ext[UTEMP], rho, eint);

        s_ext[URHO] = rho;
        s_ext[UEINT] = rho * eint;

        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[UFS + n] = rho * massfrac_ext[n];
        }
        // Velocity
        s_ext[UMX] = rho * vx_in;
        s_ext[UMY] = rho * vy_in;
        s_ext[UMZ] = rho * vz_in;
        s_ext[UEDEN] =
          s_ext[UEINT] +
          0.5 * rho * (vx_in * vx_in + vy_in * vy_in + vz_in * vz_in);
      }
    }
  }
  if (idir == 0) {

    // supersonic inflow
    amrex::Real rho, eint, Cs;
    s_ext[UTEMP] = prob_parm.Temp_amb;
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac_ext[NUM_SPECIES] = {0.0};
    massfrac_ext[O2_ID] = prob_parm.Yox_channel;
    massfrac_ext[H2_ID] = 0.0;
    massfrac_ext[N2_ID] = 1 - prob_parm.Yox_channel - massfrac_ext[H2_ID];
    eos.PYT2RE(prob_parm.Pres_amb, massfrac_ext, s_ext[UTEMP], rho, eint);
    eos.RTY2Cs(rho, prob_parm.Temp_amb, massfrac_ext, Cs);

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * Cs * prob_parm.Mach_channel;
    s_ext[UMY] = 0.0;
    s_ext[UMZ] = 0.0;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = s_ext[UEINT] + 0.5 * s_ext[UMX] * s_ext[UMX] / rho;

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * massfrac_ext[n];
    }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
