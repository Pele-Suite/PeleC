#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0 * 100.0;
  amrex::Real phi_in = -0.2;
  amrex::Real T_in = 298.0;
  amrex::Real vn_in = 0.2;
  amrex::Real pertmag = 0.0;
  amrex::Real standoff = 0.0;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> L = {{1.0}};
  int pmf_N = 0;
  int pmf_M = 0;
  bool pmf_do_average = false;

  amrex::GpuArray<amrex::Real, NVAR> fuel_state = {{0.0}};
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
