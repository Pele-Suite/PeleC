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
#include "Geometry.H"

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
//  NSCBC-Chamber -- the T7 mini-SydGex capstone (Lesson 9): a closed-end
//  premixed chamber venting through its open end, run two ways --
//
//    chamber-vent.inp     the characteristic boundary IS the vent plane
//    chamber-plenum.inp   the duct continues past the vent plane and the
//                         boundary sits two chamber lengths downstream
//
//  so that the difference of the two overpressure traces is the boundary's
//  contribution to the application-level answer.  A half-disc of burnt
//  products on the closed-end wall ignites a quiescent H2/air phi = 0.4
//  charge; the flame fingers down the duct, drives flow out of the open end,
//  crosses the vent plane, and the emptied chamber rings down at its
//  longitudinal modes -- every earlier lesson in one trace.
//
//  Scaling, stated honestly: true SydGex is a 25 cm chamber and a multi-ms
//  turbulent burn.  This chamber is 1.2 x 0.3 cm (the Sydney L/H = 4) with
//  the burn laminar and the "plenum" a straight duct extension -- the
//  boundary-placement question T7 poses, without EB geometry.  A laterally
//  expanding plenum (chamber inside a box) is the EB follow-up.
// ---------------------------------------------------------------------------

// Signed distance to the kernel edge, positive on the BURNT (inside) side,
// measured from the half-disc centred at (problo[0], kernel_y * Ly).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
kernel_burnt_distance(
  const amrex::Real x,
  const amrex::Real y,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real xc = geomdata.ProbLo()[0] + prob_parm.kernel_x;
  const amrex::Real yc =
    geomdata.ProbLo()[1] + prob_parm.kernel_y * prob_parm.L[1];
  const amrex::Real r =
    std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
  return prob_parm.kernel_r - r;
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

  // Sample the 1-D flame profile at the burnt-side distance from the kernel
  // edge: products inside, fresh charge outside, the laminar structure across
  // the interface.  The charge is QUIESCENT -- the profile's velocity column
  // is deliberately not used; the kernel starts the flow itself.  In the box
  // variants the ignition structure is clipped to the chamber interior (see
  // ProbParmDevice::kernel_clip): outside it, pure fresh charge.
  amrex::Real sb = kernel_burnt_distance(x, y, geomdata, prob_parm);
  if (
    (prob_parm.kernel_clip != 0) &&
    ((x < prob_parm.clip_lo[0]) || (x > prob_parm.clip_hi[0]) ||
     (y < prob_parm.clip_lo[1]) || (y > prob_parm.clip_hi[1]))) {
    sb = -1.0e10; // clamps to the profile's fresh end
  }
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

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * energy;
  state(i, j, k, UEDEN) = rho * energy;
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

// The reference boundary at the open end, used when pelec.bc_nscbc = 0:
// ambient pressure imposed, everything else extrapolated.  The other faces
// are walls and never reach here.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x*/,
  const amrex::Real* s_int,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if ((idir != 0) || (sgn != -1)) {
    return;
  }
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  const amrex::Real rho = s_int[URHO];
  const amrex::Real u = s_int[UMX] / rho;
  const amrex::Real v = s_int[UMY] / rho;
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = s_int[UFS + n] / rho;
  }
  amrex::Real T = 0.0, energy = 0.0;
  eos.RYP2T(rho, massfrac, prob_parm.pamb, T);
  eos.RTY2E(rho, T, massfrac, energy);
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

// ---------------------------------------------------------------------------
//  The EB follow-ups (Lesson 9's deliberate omissions): the chamber as
//  embedded walls INSIDE a laterally larger box, so the vent opens into a
//  genuine plenum rather than a straight duct continuation; and, on top of
//  that, the Sydney baffle -- two plates across the chamber leaving a
//  central gap.  The whole chamber (closed end included) is interior EB:
//  nothing touches a domain face, which is the documented
//  EB-at-domain-boundary limitation.  Three slabs (bottom, top, back) plus
//  optionally two baffle plates, as a union of boxes.  Requires
//  pelec.eb_zero_body_state = 1 like every EB+NSCBC configuration.
// ---------------------------------------------------------------------------
class NSCBCChamberBox : public pele::pelec::Geometry::Register<NSCBCChamberBox>
{
public:
  static std::string identifier() { return "nscbc-chamber-box"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static pc_nscbc::Target bcnormal_nscbc(
    const amrex::Real* /*x*/,
    const amrex::Real* /*s_int*/,
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
  {
    pc_nscbc::Target t;
    if ((idir == 0) && (sgn == -1)) {
      t.type = pc_nscbc::Type::outflow;
      t.p = prob_parm.pamb;
    }
    return t;
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;

void pc_prob_close();

#endif
