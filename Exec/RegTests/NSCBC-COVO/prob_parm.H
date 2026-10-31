#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  // 0 = convected isentropic vortex (COVO)
  // 1 = circular Gaussian pressure pulse
  int probtype = 0;

  // COVO only: apply the characteristic treatment at the inflow as well as at
  // the outflow.  Setting this to 0 leaves the inflow to the ordinary
  // bcnormal(), which isolates the outflow -- the two runs then differ in
  // exactly one boundary.
  int nscbc_inflow = 1;

  // Ambient / reference state (CGS)
  amrex::Real p_amb = 1.01325e6; // dyn/cm^2
  amrex::Real T_amb = 300.0;     // K
  amrex::Real rho_amb = 0.0;     // filled by amrex_probinit
  amrex::Real cs_amb = 0.0;      // filled by amrex_probinit
  amrex::Real gamma = 1.4;       // filled by amrex_probinit

  // Mean convection, as a Mach number (COVO: along x; pulse: quiescent)
  amrex::Real mach = 0.2;
  amrex::Real u_mean = 0.0; // filled by amrex_probinit

  // COVO: vortex strength (as a fraction of the mean speed at the core
  // radius) and core radius as a fraction of the domain height.
  amrex::Real vortex_strength = 0.5;
  amrex::Real vortex_radius = 0.075;

  // Pulse: amplitude relative to p_amb, and width as a fraction of the
  // domain half-width.
  amrex::Real pulse_amp = 1.0e-3;
  amrex::Real pulse_width = 0.06;

  // Centre of the vortex / pulse, as a fraction of the domain in each
  // direction.  The COVO default sits upstream so the vortex has room to
  // convect; the pulse default is the domain centre so the wave reaches all
  // four faces at once, including the corners.
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> centre = {0.25, 0.5};

  // Filled by amrex_probinit from geometry.prob_lo/hi.
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xc = {0.0};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> Lbox = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
