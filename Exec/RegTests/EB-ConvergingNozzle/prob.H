#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "EB.H"
#include <AMReX_EB2_IF_Rotation.H>
#include "PeleCUtilities.H"
#include "Geometry.H"

void parse_params(std::unique_ptr<ProbParmDevice>& prob_parm_device);

class EBConvergingNozzle
  : public pele::pelec::Geometry::Register<EBConvergingNozzle>
{
public:
  static std::string identifier() { return "converging-nozzle"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

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

  const amrex::Real u = 0.0;
  const amrex::Real v = 0.0;
  const amrex::Real w = 0.0;
  const amrex::Real p = prob_parm.p0;
  const amrex::Real T = prob_parm.T0;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real rho, eint;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, T, rho, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = rho * v;
  state(i, j, k, UMZ) = rho * w;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * (u * u + v * v + w * w));
  state(i, j, k, UTEMP) = T;
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
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& turb_fluc)
{

  if ((idir == 0) && (sgn == 1)) {
    // Inlet types:
    //    0 - FOExtrap pressure (recommended)
    //    1 - Enforce constant pressure
    //    2 - Characteristic extrapolation of pressure (Blazek Eqn 8.32 - 8.36)

    // Inlet conditions - use interior mass fraction
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = s_int[UFS + n] / s_int[URHO];
    }

    // Apply different inlet types for pressure
    amrex::Real p_in = 0.0;
    if (prob_parm.inlet_type == 0) {
      eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac, p_in);
    } else if (prob_parm.inlet_type == 1) {
      p_in = prob_parm.p_inlet;
    } else {
      amrex::Abort("Invalid inlet boundary type specified");
    }

    // Inlet conditions
    amrex::Real T_in = prob_parm.T_inlet;
    amrex::Real rho_in, eint_in;
    eos.PYT2RE(p_in, massfrac, T_in, rho_in, eint_in);
    amrex::Real Cs_in;
    eos.RTY2Cs(rho_in, T_in, massfrac, Cs_in);
    amrex::Real u_in = prob_parm.M_inlet * Cs_in;
    amrex::Real v_in = 0.0;
    amrex::Real w_in = 0.0;
    // Use TurbInflow stuff if needed
    if (prob_parm.turb_inflow_type == TurbInflowType::TurbulentFluctuations) {
      AMREX_D_TERM(
        u_in += turb_fluc[0];, v_in += turb_fluc[1];, w_in += turb_fluc[2];);
    } else if (prob_parm.turb_inflow_type == TurbInflowType::TotalVelocity) {
      AMREX_D_TERM(
        u_in = turb_fluc[0];, v_in = turb_fluc[1];, w_in = turb_fluc[2];);
    }

    // Set boundary cells
    s_ext[URHO] = rho_in;
    s_ext[UMX] = rho_in * u_in;
    s_ext[UMY] = rho_in * v_in;
    s_ext[UMZ] = rho_in * w_in;
    s_ext[UEINT] = rho_in * eint_in;
    s_ext[UEDEN] =
      rho_in * (eint_in + 0.5 * (u_in * u_in + v_in * v_in + w_in * w_in));
    s_ext[UTEMP] = T_in;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho_in * massfrac[n];
    }
  } else if ((idir == 0) && (sgn == -1)) {
    // outlet types:
    //    0 - FOExtrap everything
    //    1 - Characteristic Pressure Outlet (recommended)
    //    2 - Fixed Pressure

    // Initialize all quantities with FOExtrap values
    amrex::Real rho = s_int[URHO];
    amrex::Real u = s_int[UMX] / rho;
    amrex::Real v = s_int[UMY] / rho;
    amrex::Real w = s_int[UMZ] / rho;
    amrex::Real eint = s_int[UEINT] / rho;
    amrex::Real T = s_int[UTEMP];
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = s_int[UFS + n] / rho;
    }

    // If using a pressure-based outlet, compute that
    if (prob_parm.outlet_type == 0) {
      // Nothing further to be done
    } else if (prob_parm.outlet_type == 1) {
      // Following Blazek p 279, eq. 8.23 - reimplemented from EB-C10
      // NOTE: it is assumed that mass fractions are constant across the
      // boundary Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_inter = rho;
      const amrex::Real u_inter = u;
      const amrex::Real v_inter = v;
      const amrex::Real w_inter = w;
      const amrex::Real T_inter = T;
      amrex::Real p_inter = 0.0, cs_inter = 0.0;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2P(rho_inter, T_inter, massfrac, p_inter);
      eos.RTY2Cs(rho_inter, T_inter, massfrac, cs_inter);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p0;
      const amrex::Real rhob =
        s_int[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
      const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);
      const amrex::Real vb = v_inter;
      const amrex::Real wb = w_inter;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
      const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;

      eos.RYP2E(rho, massfrac, p, eint);
      eos.REY2T(rho, eint, massfrac, T);

      u = (ub - u_inter) / (xb - xd) * (x[0] - xd) + u_inter;
      v = (vb - v_inter) / (xb - xd) * (x[0] - xd) + v_inter;
      w = (wb - w_inter) / (xb - xd) * (x[0] - xd) + w_inter;
    } else if (prob_parm.outlet_type == 2) {
      amrex::Real p = prob_parm.p0;
      auto eos = pele::physics::PhysicsType::eos();
      eos.PYT2RE(p, massfrac, T, rho, eint);
    } else {
      amrex::Abort("Invalid outlet boundary type specified");
    }

    // Write to s_ext
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * massfrac[n];
    }
  }
}

void pc_prob_close();

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void problem_modify_transport_coeffs(
    const bool /*get_xi*/,
    const bool get_mu,
    const bool get_lam,
    const bool /*get_Ddiag*/,
    const bool /*get_chi*/,
    const amrex::Real /*Tloc*/,
    const amrex::Real /*rholoc*/,
    amrex::Real* /*Yloc*/,
    amrex::Real* /*Ddiag*/,
    amrex::Real* /*chi_mix*/,
    amrex::Real& mu,
    amrex::Real& /*xi*/,
    amrex::Real& lam,
    pele::physics::transport::TransParm<
      pele::physics::PhysicsType::eos_type,
      pele::physics::PhysicsType::transport_type> const* /*tparm*/,
    ProbParmDevice const& prob_parm,
    const amrex::RealVect& x)
  {

    // Create sponge zones at inlet and outlet
    if (prob_parm.do_sponge_zones) {
      amrex::Real stretch = 2.0;
      amrex::Real buffer = 3.0;
      amrex::Real end = 18.0;
      amrex::Real beg = 0.0;
      amrex::Real factor =
        0.5 * (1 + std::tanh(stretch * (beg + buffer - x[0]))) +
        0.5 * (1 + std::tanh(stretch * (x[0] - end + buffer)));

      amrex::Real visc = 1.0;
      if (get_mu) {
        mu = factor * visc;
      }
      amrex::Real cond = 1.0;
      if (get_lam) {
        lam = factor * cond;
      }
    }
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;

#endif
