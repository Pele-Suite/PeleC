#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0;

  // Ignition kernel: a half-disc of burnt products centred on the closed-end
  // wall at (problo[0] + kernel_x, kernel_y * Ly), radius kernel_r [cm].  The
  // 1-D flame profile is wrapped around its edge by radial signed distance.
  // kernel_x = 0 is the domain wall (the duct variants); the EB box variants
  // put the closed end at the chamber's EB back wall.
  amrex::Real kernel_r = 0.1;
  amrex::Real kernel_x = 0.0;
  amrex::Real kernel_y = 0.5;

  // Coordinate of the flame in the 1-D PMF profile (see NSCBC-FlameOutflow).
  amrex::Real pmf_flame_loc = 3.0;

  // Laminar flame speed, from the profile's inlet velocity; banner only.
  amrex::Real s_L = 0.0; // filled by amrex_probinit

  // Box-variant ignition clip: when the chamber is interior EB
  // (eb2.geom_type = nscbc-chamber-box), the kernel's radial profile must
  // not exist outside the chamber -- the closed-end slab is thinner than
  // the kernel radius, so an unclipped disc pokes through it and ignites
  // the plenum charge behind the chamber, and the PMF preheat tail leaks
  // past the side walls.  Cells outside [clip_lo, clip_hi] initialise as
  // pure fresh charge.  Off (kernel_clip = 0) for the duct variants:
  // bit-identical to the historical IC.
  int kernel_clip = 0;
  amrex::GpuArray<amrex::Real, 2> clip_lo = {{0.0, 0.0}};
  amrex::GpuArray<amrex::Real, 2> clip_hi = {{0.0, 0.0}};

  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> L = {{1.0}};
  int pmf_N = 0;
  int pmf_M = 0;
  bool pmf_do_average = false;
  amrex::Real* d_pmf_X = nullptr;
  amrex::Real* d_pmf_Y = nullptr;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_pmf_X;
  amrex::Vector<amrex::Real> h_pmf_Y;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_X;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_Y;
  ProbParmHost() : pmf_X(0), pmf_Y(0) {}
};

#endif
