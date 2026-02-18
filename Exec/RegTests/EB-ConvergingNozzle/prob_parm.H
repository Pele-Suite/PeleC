#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

AMREX_ENUM(TurbInflowType, None, TurbulentFluctuations, TotalVelocity);

struct ProbParmDevice
{
  amrex::Real d_inlet = 8;
  amrex::Real l_inlet = 24;
  amrex::Real M_inlet = 0.1;
  amrex::Real T_inlet = 300.0;
  amrex::Real p_inlet = 1.013e6;
  amrex::Real l_nozzle = 5;
  amrex::Real d_exit = 5;
  amrex::Real p_exit = 1.013e6;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  int inlet_type = 0;
  int outlet_type = 1;
  bool do_sponge_zones = true;
  TurbInflowType turb_inflow_type = TurbInflowType::None;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
