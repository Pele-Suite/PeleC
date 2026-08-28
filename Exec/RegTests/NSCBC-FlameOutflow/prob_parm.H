#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0;
  amrex::Real phi_in = -0.5;
  amrex::Real T_in = 298.0;
  amrex::Real vn_in = 0.2;
  amrex::Real pertmag = 0.0;

  // Mean position of the flame sheet as a fraction of the domain
  // length, and the wrinkle amplitude and wavenumber.  The default parks the
  // mean sheet position exactly ON the outflow plane, so that half the boundary
  // is burnt, half is fresh, and the reaction zone sits inside the boundary
  // cells at the two crossings from the first step.  At U = 2 S_L the sheet
  // then drifts out at S_L -- which is slow compared with the boundary
  // relaxation time, so it is effectively a stationary straddling flame over
  // the length of a measurement (see README).
  amrex::Real sheet_x = 1.0;
  amrex::Real wrinkle_amp = 0.08;
  amrex::Real wrinkle_k = 1.0;
  amrex::Real sheet_x0 = 0.0; // absolute, filled by amrex_probinit
  amrex::Real ylo = 0.0;      // problo[1], so the wrinkle phase is anchored

  // Treat the INLET characteristically as well.  Setting this to 0 leaves the
  // inlet a hard Dirichlet and puts the characteristic treatment on the
  // outflow alone, which is what isolates the outflow when comparing against a
  // longer-domain reference: the two domains then differ in nothing but where
  // the outflow is, and the inflow relaxation rate K = relax_u c / L_ref -- a
  // rate that depends on the domain length -- cannot contaminate the result.
  int nscbc_inflow = 1;

  // Mean inlet speed as a multiple of the laminar flame speed.  The only
  // requirement is U > S_L so the sheet drifts downstream and needs no
  // anchor; 2 is the natural choice.
  amrex::Real u_ratio = 2.0;

  // Laminar flame speed, taken from the inlet velocity of the PMF profile.
  amrex::Real s_L = 0.0;  // filled by amrex_probinit
  amrex::Real u_in = 0.0; // = u_ratio * s_L

  // Coordinate of the flame in the 1-D PMF profile: the initial condition maps
  // that profile onto the sheet through the signed normal distance s,
  // evaluating it at (pmf_flame_loc + s).
  amrex::Real pmf_flame_loc = 3.0;

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
