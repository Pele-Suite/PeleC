#ifndef PROB_H
#define PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "Constants.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  const amrex::Real xlo,
  const amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector,
  const ProbParmDevice& prob_parm)
{
  if (prob_parm.pmf_do_average) {
    int lo_loside = 0;
    int lo_hiside = 0;
    int hi_loside = 0;
    int hi_hiside = 0;
    if (xlo < prob_parm.d_pmf_X[0]) {
      lo_loside = 0;
      lo_hiside = 0;
    }
    if (xhi < prob_parm.d_pmf_X[0]) {
      hi_loside = 0;
      hi_hiside = 0;
    }
    if (xlo > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      lo_loside = prob_parm.pmf_N - 1;
      lo_hiside = prob_parm.pmf_N - 1;
    }
    if (xhi > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      hi_loside = prob_parm.pmf_N - 1;
      hi_hiside = prob_parm.pmf_N - 1;
    }
    if (lo_loside == 0) {
      for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if ((xlo > prob_parm.d_pmf_X[i]) && (xlo < prob_parm.d_pmf_X[i + 1])) {
          lo_loside = i;
          lo_hiside = i + 1;
        }
      }
    }
    if (hi_loside == 0) {
      for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if ((xhi > prob_parm.d_pmf_X[i]) && (xhi < prob_parm.d_pmf_X[i + 1])) {
          hi_loside = i;
          hi_hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < prob_parm.pmf_M; j++) {
      amrex::Real x1 = prob_parm.d_pmf_X[lo_loside];
      amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + lo_loside];
      amrex::Real x2 = prob_parm.d_pmf_X[lo_hiside];
      amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + lo_hiside];
      amrex::Real dydx = 0.0;
      if (lo_loside == lo_hiside) {
        dydx = 0.0;
      } else {
        dydx = (y2 - y1) / (x2 - x1);
      }
      amrex::Real ylo = y1 + dydx * (xlo - x1);
      amrex::Real yhi = 0.0;
      if (lo_loside == hi_loside) {
        yhi = y1 + dydx * (xhi - x1);
        y_vector[j] = 0.5 * (ylo + yhi);
      } else {
        amrex::Real sum = (x2 - xlo) * 0.5 * (ylo + y2);
        x1 = prob_parm.d_pmf_X[hi_loside];
        y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hi_loside];
        x2 = prob_parm.d_pmf_X[hi_hiside];
        y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hi_hiside];
        if (hi_loside == hi_hiside) {
          dydx = 0.0;
        } else {
          dydx = (y2 - y1) / (x2 - x1);
        }
        yhi = y1 + dydx * (xhi - x1);
        sum += (xhi - x1) * 0.5 * (yhi + y1);
        for (int k = lo_hiside; k < hi_loside - 1; k++) {
          sum += (prob_parm.d_pmf_X[k + 1] - prob_parm.d_pmf_X[k]) * 0.5 *
                 (prob_parm.d_pmf_Y[prob_parm.pmf_N * j + k] +
                  prob_parm.d_pmf_Y[prob_parm.pmf_N * j + k + 1]);
        }
        y_vector[j] = sum / (xhi - xlo);
      }
    }
  } else {
    amrex::Real xmid = 0.5 * (xlo + xhi);
    int loside = -1;
    int hiside = -1;
    if (xmid < prob_parm.d_pmf_X[0]) {
      loside = 0;
      hiside = 0;
    }
    if (xmid > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      loside = prob_parm.pmf_N - 1;
      hiside = prob_parm.pmf_N - 1;
    }
    if (loside == -1) {
      for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if (
          (xmid >= prob_parm.d_pmf_X[i]) &&
          (xmid <= prob_parm.d_pmf_X[i + 1])) {
          loside = i;
          hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < prob_parm.pmf_M; j++) {
      const amrex::Real x1 = prob_parm.d_pmf_X[loside];
      const amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + loside];
      const amrex::Real x2 = prob_parm.d_pmf_X[hiside];
      const amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hiside];
      amrex::Real dydx = 0.0;
      if (loside == hiside) {
        dydx = 0.0;
      } else {
        dydx = (y2 - y1) / (x2 - x1);
      }
      y_vector[j] = y1 + dydx * (xmid - x1);
    }
  }
}

// ---------------------------------------------------------------------------
//  NSCBC-FlameOutflow -- a premixed flame whose reaction zone lies ACROSS a
//  subsonic characteristic outflow: a wrinkled planar sheet, unanchored, in a
//  uniform stream at U = 2 S_L, straddling the outflow plane.
//
//  It needs no anchor: the stream carries the sheet downstream at
//  U - S_L = S_L, so the flame simply leaves, and there is no ignition kernel,
//  no apex, and no steady-state to wait for.  The
//  wrinkle x_f(y) = x0 + A cos(2 pi k y / Ly) with x0 on the outflow plane puts
//  half the boundary in burnt gas and half in fresh, with the reaction zone
//  passing through the boundary cells at the two crossings and the flame normal
//  tilted up to atan(A k) ~ 52 deg away from the boundary normal there.  So a
//  single configuration loads all three of the things the other cases cannot:
//  heat release inside the boundary cell (bc_nscbc_beta_s), strong tangential
//  gradients along the boundary (bc_nscbc_beta), and a multi-species diffusive
//  structure crossing it.
//
//  What this case does NOT do is advect out within a measurement.  The drift
//  Mach number is S_L/c ~ 6e-4, so a run long enough for the flame to clear the
//  boundary is some 400 acoustic transits.  The sheet is therefore initialised
//  already straddling, and is quasi-frozen over the few relaxation times the
//  boundary needs to equilibrate.  That is not a compromise of the test -- the
//  quantity under test is the boundary's response to a live reaction zone in
//  its own cells, and that is present from step 0.
// ---------------------------------------------------------------------------

// Signed normal distance to a WRINKLED PLANAR SHEET at x = xf(y), positive on
// the burnt (downstream) side, with the unit normal pointing fresh -> burnt.
//
//     xf(y) = x0 + A cos(2 pi k y / Ly)
//     phi   = x - xf(y)          (level set; zero on the sheet)
//     n_hat = grad(phi)/|grad(phi)| = (1, -xf'(y)) / sqrt(1 + xf'^2)
//
// The gradient magnitude divides out so that phi/|grad phi| is the true
// distance to first order, which is what the 1-D profile should be sampled at.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
sheet_signed_distance(
  const amrex::Real x,
  const amrex::Real y,
  ProbParmDevice const& prob_parm,
  amrex::Real& nx,
  amrex::Real& ny)
{
  const amrex::Real Ly = prob_parm.L[1];
  const amrex::Real kk = 2.0 * constants::PI() * prob_parm.wrinkle_k / Ly;
  const amrex::Real yr = y - prob_parm.ylo;
  const amrex::Real xf =
    prob_parm.sheet_x0 + prob_parm.wrinkle_amp * std::cos(kk * yr);
  const amrex::Real dxf = -prob_parm.wrinkle_amp * kk * std::sin(kk * yr);
  const amrex::Real g = std::sqrt(1.0 + dxf * dxf);
  nx = 1.0 / g;
  ny = -dxf / g;
  return (x - xf) / g;
}

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
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  // Wrinkled sheet: sb is the distance on the BURNT side.
  amrex::Real nx = 0.0, ny = 0.0;
  const amrex::Real sb = sheet_signed_distance(x, y, prob_parm, nx, ny);

  // Sample the 1-D flame.  The profile runs fresh (low) to burnt (high), so
  // the burnt-side distance adds to the flame location.
  const amrex::Real xp = prob_parm.pmf_flame_loc + sb;
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  pmf(xp, xp, pmf_vals, prob_parm);

  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    molefrac[n] = pmf_vals[3 + n];
  }
  const amrex::Real T = pmf_vals[0];
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  amrex::Real rho = 0.0, energy = 0.0;
  eos.PYT2RE(prob_parm.pamb, massfrac, T, rho, energy);

  // Uniform stream plus the flame's own dilatation along its normal.  On the
  // fresh side un = S_L and this is exactly (U, 0); across the flame the gas
  // is accelerated normal to the sheet by the density ratio, as it must be.
  const amrex::Real un = pmf_vals[1]; // flame-normal speed from the profile
  const amrex::Real u = prob_parm.u_in + (un - prob_parm.s_L) * nx;
  const amrex::Real v = (un - prob_parm.s_L) * ny;

  const amrex::Real ke = 0.5 * (u * u + v * v);
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = rho * v;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * energy;
  state(i, j, k, UEDEN) = rho * (energy + ke);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

// Fresh-stream inlet state.  Shared by the reference bcnormal and by the
// characteristic inflow target so the two cannot drift apart.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
inlet_state(
  ProbParmDevice const& prob_parm,
  amrex::Real& T,
  amrex::Real massfrac[NUM_SPECIES])
{
  const amrex::Real xp = prob_parm.d_pmf_X[0]; // fresh end of the profile
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  pmf(xp, xp, pmf_vals, prob_parm);
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    molefrac[n] = pmf_vals[3 + n];
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  T = pmf_vals[0];
}

// The reference boundary condition, used when pelec.bc_nscbc = 0: the inlet
// state imposed hard, and ambient pressure imposed hard at the outflow.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real* s_int,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if (idir != 0) {
    return; // transverse faces are periodic
  }
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real T = 0.0, rho = 0.0, energy = 0.0, u = 0.0, v = 0.0;

  if (sgn == 1) {
    inlet_state(prob_parm, T, massfrac);
    u = prob_parm.u_in;
    v = 0.0;
    eos.PYT2RE(prob_parm.pamb, massfrac, T, rho, energy);
  } else {
    // Outflow: ambient pressure imposed, everything else extrapolated.
    rho = s_int[URHO];
    u = s_int[UMX] / rho;
    v = s_int[UMY] / rho;
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = s_int[UFS + n] / rho;
    }
    eos.RYP2T(rho, massfrac, prob_parm.pamb, T);
    eos.RTY2E(rho, T, massfrac, energy);
  }

  const amrex::Real ke = 0.5 * (u * u + v * v);
  s_ext[URHO] = rho;
  s_ext[UMX] = rho * u;
  s_ext[UMY] = rho * v;
  s_ext[UMZ] = 0.0;
  s_ext[UEINT] = rho * energy;
  s_ext[UEDEN] = rho * (energy + ke);
  s_ext[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * massfrac[n];
  }
}

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static pc_nscbc::Target bcnormal_nscbc(
    const amrex::Real x[AMREX_SPACEDIM],
    const amrex::Real* /*s_int[NVAR]*/,
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
  {
    pc_nscbc::Target t;
    if (idir != 0) {
      return t; // periodic transverse
    }
    if (sgn == -1) {
      t.type = pc_nscbc::Type::outflow;
      t.p = prob_parm.pamb;
      return t;
    }
    if (prob_parm.nscbc_inflow == 0) {
      return t; // Type::off -> the inlet falls through to the hard bcnormal
    }
    amrex::Real T = 0.0, massfrac[NUM_SPECIES] = {0.0};
    inlet_state(prob_parm, T, massfrac);
    t.type = pc_nscbc::Type::inflow;
    t.T = T;
    t.u[0] = prob_parm.u_in;
    t.u[1] = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      t.Y[n] = massfrac[n];
    }
    return t;
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;

void pc_prob_close();

#endif
