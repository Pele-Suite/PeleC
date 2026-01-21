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

#include "SootModel.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  const amrex::Real xlo,
  const amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector,
  const ProbParmDevice& prob_parm)
{
  if (prob_parm.pmf_do_average == 1) {
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
      y_vector[j] = y1 + dydx * (xlo - x1);
    }
  }
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
  const int tempCol = prob_parm.tempCol;
  const int velCol = prob_parm.velCol;
  const int specCol = prob_parm.specCol;
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real x1 = (x - prob_parm.standoff - 0.5 * dx[0]);
  amrex::Real x2 = (x - prob_parm.standoff + 0.5 * dx[0]);
  pmf(x1, x2, pmf_vals, prob_parm);
  const amrex::Real T = pmf_vals[tempCol];
  const amrex::Real pres = prob_parm.pamb;
  amrex::Real u[3] = {0.0};
  u[0] = pmf_vals[velCol];
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real norm = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = amrex::Clamp(pmf_vals[specCol + n], 0.0, 1.0);
    norm += massfrac[n];
  }
  for (auto& mf : massfrac) {
    mf /= norm;
  }

  amrex::Real rho = 0.0;
  amrex::Real energy = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(pres, massfrac, T, rho, energy);
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEINT) = rho * energy;
  state(i, j, k, UEDEN) =
    rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; n++) {
    state(i, j, k, UFSOOT + n) = prob_parm.soot_vals[n];
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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  const int tempCol = prob_parm.tempCol;
  const int velCol = prob_parm.velCol;
  const int specCol = prob_parm.specCol;
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real u[3] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real rho, energy, T, pres;
  auto eos = pele::physics::PhysicsType::eos();

  // Outflow, extrapolate density and Y, set T, p, u
  if (sgn == -1) {
    rho = s_int[URHO];
    pmf(prob_hi[idir], prob_hi[idir], pmf_vals, prob_parm);
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = s_int[UFS + n] / rho;
    }
    T = pmf_vals[tempCol];
    pres = prob_parm.pamb;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[0] = pmf_vals[velCol];
    eos.PYT2RE(pres, massfrac, T, rho, energy);
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u[0];
    s_ext[UMY] = rho * u[1];
    s_ext[UMZ] = rho * u[2];
    s_ext[UEINT] = rho * energy;
    s_ext[UEDEN] =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * massfrac[n];
    }
  } else {
    pmf(prob_lo[idir], prob_lo[idir], pmf_vals, prob_parm);
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = pmf_vals[specCol + n];
    }
    T = pmf_vals[tempCol];
    pres = prob_parm.pamb;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[0] = pmf_vals[velCol];
    eos.PYT2RE(pres, massfrac, T, rho, energy);
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u[0];
    s_ext[UMY] = rho * u[1];
    s_ext[UMZ] = rho * u[2];
    s_ext[UEINT] = rho * energy;
    s_ext[UEDEN] =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * massfrac[n];
    }
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; n++) {
    s_ext[UFSOOT + n] = s_int[UFSOOT + n];
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
