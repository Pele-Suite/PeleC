#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

// ---------------------------------------------------------------------------
//  NSCBC-COVO -- two multi-dimensional tests of the characteristic boundary
//  treatment, sharing one executable.
//
//  probtype = 0  Convected isentropic vortex (COVO).  A Yee/Shu vortex is
//                carried by a uniform stream from a subsonic inflow at x-lo to
//                a subsonic outflow at x-hi, with y periodic.  The vortex is a
//                pure vorticity/entropy structure: an ideal boundary lets it
//                leave without generating any acoustic signal at all, so the
//                pressure left behind is a direct measure of what the boundary
//                got wrong.  It is also the only test here that exercises the
//                INFLOW path end to end.
//
//  probtype = 1  Circular Gaussian pressure pulse, after Motheau et al. (2017).
//                A quiescent square with all four faces non-reflecting.  The
//                pulse splits into an outgoing circular wave that meets every
//                face, and every corner, at once.  The diagnostic is the SHAPE
//                of the wavefront: with a good boundary the pressure contours
//                stay circular right through the boundary; with a poor one they
//                flatten and buckle where the wave meets the face, and worst of
//                all at the corners, where the wave arrives at 45 degrees and
//                the neglected transverse terms matter most.  This is the test
//                that shows honestly what a 1-D-normal LODI boundary cannot do.
// ---------------------------------------------------------------------------
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
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x =
    prob_lo[0] + (static_cast<amrex::Real>(i) + 0.5) * dx[0];
  const amrex::Real y =
    prob_lo[1] + (static_cast<amrex::Real>(j) + 0.5) * dx[1];

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  amrex::Real rho = prob_parm.rho_amb;
  amrex::Real p = prob_parm.p_amb;
  amrex::Real T = prob_parm.T_amb;
  amrex::Real u = prob_parm.u_mean;
  amrex::Real v = 0.0;

  const amrex::Real dxv = x - prob_parm.xc[0];
  const amrex::Real dyv = y - prob_parm.xc[1];

  if (prob_parm.probtype == 0) {
    // Isentropic vortex.  beta is the circulation parameter; the velocity
    // perturbation peaks at r = R with magnitude beta*u_mean/sqrt(e).
    const amrex::Real R = prob_parm.vortex_radius * prob_parm.Lbox[1];
    const amrex::Real r2 = (dxv * dxv + dyv * dyv) / (R * R);
    const amrex::Real beta = prob_parm.vortex_strength * prob_parm.u_mean;
    const amrex::Real fac = beta * std::exp(0.5 * (1.0 - r2));

    u += -fac * dyv / R;
    v += fac * dxv / R;

    // Isentropic closure: the temperature deficit balances the swirl.
    const amrex::Real g = prob_parm.gamma;
    const amrex::Real dT = -(g - 1.0) * beta * beta * std::exp(1.0 - r2) /
                           (2.0 * g * prob_parm.p_amb / prob_parm.rho_amb) *
                           prob_parm.T_amb;
    T = prob_parm.T_amb + dT;
    p = prob_parm.p_amb * std::pow(T / prob_parm.T_amb, g / (g - 1.0));
    eos.PYT2R(p, massfrac, T, rho);
  } else {
    // Circular pressure pulse, isentropic and at rest, so that it splits
    // symmetrically into an ingoing and an outgoing circular acoustic wave
    // and carries no entropy content of its own.
    const amrex::Real w = prob_parm.pulse_width * prob_parm.Lbox[0];
    const amrex::Real r2 = (dxv * dxv + dyv * dyv) / (w * w);
    const amrex::Real dp =
      prob_parm.pulse_amp * prob_parm.p_amb * std::exp(-r2);
    p = prob_parm.p_amb + dp;
    rho = prob_parm.rho_amb + dp / (prob_parm.cs_amb * prob_parm.cs_amb);
    eos.RYP2T(rho, massfrac, p, T);
    u = prob_parm.u_mean;
    v = 0.0;
  }

  amrex::Real eint = 0.0;
  eos.RTY2E(rho, T, massfrac, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = rho * v;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * (u * u + v * v));
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

// The reference "bad" boundary, used when pelec.bc_nscbc = 0.  Ambient
// pressure and temperature imposed hard, velocity extrapolated -- the
// straightforward thing to write, and the thing the characteristic treatment
// is meant to replace.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* s_int,
  amrex::Real* s_ext,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  const amrex::Real u = s_int[UMX] / s_int[URHO];
  const amrex::Real v = s_int[UMY] / s_int[URHO];

  amrex::Real rho = 0.0, eint = 0.0;
  eos.PYT2RE(prob_parm.p_amb, massfrac, prob_parm.T_amb, rho, eint);

  s_ext[URHO] = rho;
  s_ext[UMX] = rho * u;
  s_ext[UMY] = rho * v;
  s_ext[UMZ] = 0.0;
  s_ext[UEINT] = rho * eint;
  s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v));
  s_ext[UTEMP] = prob_parm.T_amb;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * massfrac[n];
  }
}

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static pc_nscbc::Target bcnormal_nscbc(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::Real* /*s_int[NVAR]*/,
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
  {
    pc_nscbc::Target t;
    if (prob_parm.probtype == 0) {
      if ((idir == 0) && (sgn == +1) && (prob_parm.nscbc_inflow != 0)) {
        // Subsonic inflow: velocity, temperature and composition are the
        // incoming characteristics and so are what gets specified.
        t.type = pc_nscbc::Type::inflow;
        t.u[0] = prob_parm.u_mean;
        t.u[1] = 0.0;
        t.T = prob_parm.T_amb;
        t.Y[0] = 1.0;
      } else if ((idir == 0) && (sgn == -1)) {
        // Subsonic outflow: pressure is the only incoming characteristic.
        t.type = pc_nscbc::Type::outflow;
        t.p = prob_parm.p_amb;
      }
    } else {
      // All four faces non-reflecting.
      t.type = pc_nscbc::Type::outflow;
      t.p = prob_parm.p_amb;
    }
    return t;
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;

void pc_prob_close();

#endif
