#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

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

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
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
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  AMREX_D_TERM(
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

  AMREX_D_TERM(
    const amrex::Real Lx = prob_hi[0] - prob_lo[0];
    , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
    , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

  AMREX_D_TERM(
    const amrex::Real xc = prob_lo[0] + 0.5 * Lx;
    , const amrex::Real yc = prob_lo[1] + 0.5 * Ly;
    , const amrex::Real zc = prob_lo[2] + 0.5 * Lz;);

  const amrex::Real deltax = x - prob_parm.xvort;
  const amrex::Real deltay = y - prob_parm.yvort;
  const amrex::Real d_sq = deltax * deltax + deltay * deltay;
  const amrex::Real r_sq = prob_parm.rvort * prob_parm.rvort;
  const amrex::Real u_vort =
    -prob_parm.forcevort * deltay / r_sq * exp(-d_sq / r_sq / 2.);
  const amrex::Real v_vort =
    prob_parm.forcevort * deltax / r_sq * exp(-d_sq / r_sq / 2.);
  const amrex::Real w_vort = 0.;

  const amrex::Real u = u_vort;
  const amrex::Real v = v_vort;
  const amrex::Real w = w_vort;
  const amrex::Real p = prob_parm.P_mean;
  const amrex::Real T = prob_parm.T_mean;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = 0.233;
  massfrac[N2_ID] = 0.767;
  if (std::abs(x - xc) < 0.1 * Lx) {
    massfrac[O2_ID] = 0.283;
    massfrac[N2_ID] = 0.717;
  }

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
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
