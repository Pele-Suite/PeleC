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
#include "Constants.H"
#include "PelePhysics.H"
#include "EOS.H"

#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "EB.H"
#include <AMReX_EB2_IF_Rotation.H>
#include "PeleCUtilities.H"
#include "Geometry.H"

class EBConcCylinders : public pele::pelec::Geometry::Register<EBConcCylinders>
{
public:
  static std::string identifier() { return "conc_cylinders"; }

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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  amrex::Real u = 0.0;
  amrex::Real v = 0.0;
  amrex::Real w = 0.0;
  const amrex::Real p = prob_parm.p0;
  const amrex::Real T = prob_parm.T0;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real rho, eint;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, T, rho, eint);

  if (prob_parm.do_rf) {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      prob_parm.rf_axis_x, prob_parm.rf_axis_y, prob_parm.rf_axis_z};
    amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    amrex::RealVect r = get_rotaxis_vec(iv, axis_loc, geomdata);
    amrex::RealVect omgavec(0.0, 0.0, 0.0);
    omgavec[prob_parm.rf_axis] = prob_parm.rf_omega;
    amrex::RealVect w_cross_r = omgavec.crossProduct(r);
    u = -w_cross_r[0];
    v = -w_cross_r[1];
    w = -w_cross_r[2];
  }

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

  if (prob_parm.do_rf) {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      prob_parm.rf_axis_x, prob_parm.rf_axis_y, prob_parm.rf_axis_z};
    amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    amrex::Real rotenrg = get_rot_energy(
      iv, prob_parm.rf_omega, prob_parm.rf_axis, axis_loc, geomdata);
    state(i, j, k, UEDEN) -= rho * rotenrg;
  }
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
