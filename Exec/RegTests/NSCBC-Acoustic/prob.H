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
//  NSCBC-Acoustic
//
//  Two pulse geometries in one dimension-agnostic setup.  The same sources
//  build in 1-D, 2-D and 3-D; only the inputs file changes.
//
//  pulse_type = 0  PLANAR.  A right-running isentropic pulse in a quiescent
//                  gas, launched along +x toward a single non-reflecting
//                  outflow at x-hi, everything else periodic or wall.  This is
//                  the AMReX-side counterpart of check C4 in
//                  Verification/NSCBC1D: the pulse should leave and the
//                  residual left behind should be a small fraction of the
//                  incident amplitude.  Historical configuration; unchanged.
//
//  pulse_type = 1  RADIAL.  A Gaussian pressure pulse at the box centre,
//                  isentropic and AT REST, so it splits into an outgoing
//                  spherical (2-D: circular) wave and an imploding one.  Every
//                  face is a characteristic outflow.
//
//  Geometry 1 is the reason this case is worth building in 3-D at all.  A
//  planar pulse loads one face at normal incidence, which is what the 2-D run
//  already tests; a radial pulse reaches the six faces at normal incidence,
//  the twelve edges at 45 degrees, and the eight corners along the body
//  diagonal, in a single run.  Those are the places where a ghost cell lies
//  outside the domain in more than one direction at once, so they test the
//  corner-ownership rule in BCfill.cpp rather than the 1-D algebra in NSCBC.H.
//  Because the initial condition is exactly isotropic about the box centre,
//  ANY departure from isotropy in the departing front is boundary-generated,
//  and no reference solution is needed to say so.
//
//  Run either with pelec.bc_nscbc = 0 to see what a hard boundary does.
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

  const amrex::Real xv[AMREX_SPACEDIM] = {AMREX_D_DECL(
    prob_lo[0] + (static_cast<amrex::Real>(i) + 0.5) * dx[0],
    prob_lo[1] + (static_cast<amrex::Real>(j) + 0.5) * dx[1],
    prob_lo[2] + (static_cast<amrex::Real>(k) + 0.5) * dx[2])};

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  // Distance from the pulse, and the direction the pulse runs in.  For the
  // planar pulse that direction is +x everywhere; for the radial pulse the
  // wave is at rest initially and no direction is needed at all.  The duct
  // (pulse_type = 2) has no pulse: a uniform stream at u0, forced only
  // through the boundary target.
  amrex::Real s = 0.0;
  if (prob_parm.pulse_type == 0) {
    s = xv[0] - prob_parm.xc[0];
  } else {
    const amrex::Real d[AMREX_SPACEDIM] = {AMREX_D_DECL(
      xv[0] - prob_parm.xc[0], xv[1] - prob_parm.xc[1],
      xv[2] - prob_parm.xc[2])};
    s = std::sqrt(AMREX_D_TERM(d[0] * d[0], +d[1] * d[1], +d[2] * d[2]));
  }

  const amrex::Real dp =
    (prob_parm.pulse_type == 2)
      ? 0.0
      : prob_parm.amp * prob_parm.p_amb *
          std::exp(-(s / prob_parm.w) * (s / prob_parm.w));

  // Isentropic in both cases: drho = dp / c^2.  The planar pulse is given the
  // right-running velocity dp = rho c u so that it is a single wave; the
  // radial pulse is left at rest so that it splits symmetrically, which is
  // what makes the outgoing front isotropic by construction.
  const amrex::Real p = prob_parm.p_amb + dp;
  const amrex::Real rho =
    prob_parm.rho_amb + dp / (prob_parm.cs_amb * prob_parm.cs_amb);

  amrex::Real vel[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};
  if (prob_parm.pulse_type == 0) {
    vel[0] = dp / (prob_parm.rho_amb * prob_parm.cs_amb);
  } else if (prob_parm.pulse_type == 2) {
    vel[0] = prob_parm.u0;
  }
  const amrex::Real ke =
    0.5 * AMREX_D_TERM(vel[0] * vel[0], +vel[1] * vel[1], +vel[2] * vel[2]);

  amrex::Real T = 0.0, eint = 0.0;
  eos.RYP2T(rho, massfrac, p, T);
  eos.RTY2E(rho, T, massfrac, eint);

  state(i, j, k, URHO) = rho;
  // PeleC always carries three momentum slots; the ones this build does not
  // use stay at zero.
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  AMREX_D_TERM(state(i, j, k, UMX) = rho * vel[0];
               , state(i, j, k, UMY) = rho * vel[1];
               , state(i, j, k, UMZ) = rho * vel[2];)
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + ke);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

// The ordinary ghost fill.  Only reached when pelec.bc_nscbc = 0, or for
// points where bcnormal_nscbc returns `off`.  Deliberately crude -- a hard
// imposition of the ambient density and pressure, keeping the interior
// velocity -- so that the comparison against the characteristic treatment is
// stark.
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

  const amrex::Real rho = prob_parm.rho_amb;
  const amrex::Real rho_inv_int = 1.0 / s_int[URHO];
  const amrex::Real vel[AMREX_SPACEDIM] = {AMREX_D_DECL(
    s_int[UMX] * rho_inv_int, s_int[UMY] * rho_inv_int,
    s_int[UMZ] * rho_inv_int)};
  const amrex::Real ke =
    0.5 * AMREX_D_TERM(vel[0] * vel[0], +vel[1] * vel[1], +vel[2] * vel[2]);

  amrex::Real T = 0.0, eint = 0.0;
  eos.RYP2T(rho, massfrac, prob_parm.p_amb, T);
  eos.RTY2E(rho, T, massfrac, eint);

  s_ext[URHO] = rho;
  s_ext[UMX] = 0.0;
  s_ext[UMY] = 0.0;
  s_ext[UMZ] = 0.0;
  AMREX_D_EXPR(
    s_ext[UMX] = rho * vel[0], s_ext[UMY] = rho * vel[1],
    s_ext[UMZ] = rho * vel[2]);
  s_ext[UEINT] = rho * eint;
  s_ext[UEDEN] = rho * (eint + ke);
  s_ext[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * massfrac[n];
  }
}

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  // The characteristic treatment.  For a pure non-reflecting outflow the
  // target pressure is the only thing that needs to be specified.
  //
  // Which faces are characteristic is a property of the pulse geometry, not of
  // the inputs file: the planar pulse only ever reaches x-hi, so making the
  // other faces characteristic would test nothing and would remove the wall
  // that keeps the historical configuration what it was.  The radial pulse
  // reaches all of them.
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static pc_nscbc::Target bcnormal_nscbc(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::Real* /*s_int[NVAR]*/,
    const int idir,
    const int sgn,
    const amrex::Real time,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
  {
    pc_nscbc::Target t;
    if (prob_parm.pulse_type == 2) {
      // The duct: a characteristic INFLOW at x-lo whose velocity target
      // carries the harmonic forcing -- the time argument this hook always
      // received, finally consumed.  x-hi stays `off`, so the plumbing falls
      // back to bcnormal's hard-pressure ghost: the fully reflecting far end.
      if ((idir == 0) && (sgn == +1)) {
        t.type = pc_nscbc::Type::inflow;
        const amrex::Real om =
          2.0 * amrex::Math::pi<amrex::Real>() * prob_parm.force_freq;
        t.u[0] = prob_parm.u0 + prob_parm.force_amp * std::sin(om * time);
        // The stateless injection slot (Target::dudt): the analytic rate of
        // the forcing, so the signal does not ride the relaxation.  Opt-in
        // (prob.force_feedforward) so the classical-inlet tables stay
        // reproducible.
        if (prob_parm.force_feedforward != 0) {
          t.dudt = prob_parm.force_amp * om * std::cos(om * time);
        }
        t.T = prob_parm.T_amb;
        t.Y[0] = 1.0;
        for (int n = 1; n < NUM_SPECIES; n++) {
          t.Y[n] = 0.0;
        }
      }
      return t;
    }
    const bool live =
      (prob_parm.pulse_type == 0) ? ((idir == 0) && (sgn == -1)) : true;
    if (live) {
      t.type = pc_nscbc::Type::outflow;
      t.p = prob_parm.p_amb;
    }
    return t;
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;

void pc_prob_close();

#endif
