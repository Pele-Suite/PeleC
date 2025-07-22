#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "Utilities.H"
#include "PeleCUtilities.H"
#include "Constants.H"

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

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  const amrex::Real p_si2cgs = 10.0;

  // TG functions
  amrex::Real u[3] = {0.0};
  u[0] = prob_parm.u0 * sin(prob_parm.omega_x * x / prob_parm.L_x) *
         cos(prob_parm.omega_y * y / prob_parm.L_y) *
         cos(prob_parm.omega_z * z / prob_parm.L_z);
  u[1] = -prob_parm.u0 * cos(prob_parm.omega_x * x / prob_parm.L_x) *
         sin(prob_parm.omega_y * y / prob_parm.L_y) *
         cos(prob_parm.omega_z * z / prob_parm.L_z);
  amrex::Real rho = 0.0, eint = 0.0, p, T;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();

  if ((prob_parm.case_type_int == 1) || (prob_parm.case_type_int == 2)) {
    amrex::Real xmod = std::fmod(x, prob_parm.L_x);
    int m = 0;
    pele::physics::utilities::locate(prob_parm.d_xarray, prob_parm.nx, xmod, m);
    const int mp1 = (m + 1) % prob_parm.nx;
    const amrex::Real fact =
      (x - prob_parm.d_input[prob_parm.input_x + m * prob_parm.nvars]) /
      prob_parm.d_dxinput[m];
    const amrex::Real f0 = (1 - fact);
    const amrex::Real f1 = fact;

    p = (prob_parm.d_input[prob_parm.input_p + m * prob_parm.nvars] * f0 +
         prob_parm.d_input[prob_parm.input_p + mp1 * prob_parm.nvars] * f1) *
        p_si2cgs;
    T = prob_parm.d_input[prob_parm.input_Tad + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_Tad + mp1 * prob_parm.nvars] * f1;
    if (prob_parm.case_type_int == 1) {
#ifdef USE_FUEGO_EOS
      massfrac[H2_ID] =
        prob_parm.d_input[prob_parm.input_H2_init + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_H2_init + mp1 * prob_parm.nvars] * f1;
      massfrac[O2_ID] =
        prob_parm.d_input[prob_parm.input_O2_init + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_O2_init + mp1 * prob_parm.nvars] * f1;
      massfrac[N2_ID] =
        prob_parm.d_input[prob_parm.input_N2_init + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_N2_init + mp1 * prob_parm.nvars] * f1;
      eos.PYT2RE(p, massfrac, T, rho, eint);
#endif
    } else if (prob_parm.case_type_int == 2) {
#ifdef USE_FUEGO_EOS
      massfrac[H2_ID] =
        prob_parm.d_input[prob_parm.input_H2 + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_H2 + mp1 * prob_parm.nvars] * f1;
      massfrac[O2_ID] =
        prob_parm.d_input[prob_parm.input_O2 + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_O2 + mp1 * prob_parm.nvars] * f1;
      massfrac[N2_ID] =
        prob_parm.d_input[prob_parm.input_N2 + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_N2 + mp1 * prob_parm.nvars] * f1;
      massfrac[H2O_ID] =
        prob_parm.d_input[prob_parm.input_H2O + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_H2O + mp1 * prob_parm.nvars] * f1;
      massfrac[H_ID] =
        prob_parm.d_input[prob_parm.input_H + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_H + mp1 * prob_parm.nvars] * f1;
      massfrac[O_ID] =
        prob_parm.d_input[prob_parm.input_O + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_O + mp1 * prob_parm.nvars] * f1;
      massfrac[OH_ID] =
        prob_parm.d_input[prob_parm.input_OH + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_OH + mp1 * prob_parm.nvars] * f1;
      massfrac[HO2_ID] =
        prob_parm.d_input[prob_parm.input_HO2 + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_HO2 + mp1 * prob_parm.nvars] * f1;
      massfrac[H2O2_ID] =
        prob_parm.d_input[prob_parm.input_H2O2 + m * prob_parm.nvars] * f0 +
        prob_parm.d_input[prob_parm.input_H2O2 + mp1 * prob_parm.nvars] * f1;
      eos.PYT2RE(p, massfrac, T, rho, eint);
#endif
    }
  } else {
    massfrac[0] = 1.0;
    rho = prob_parm.rho0;
    T = prob_parm.T0;
    p = prob_parm.p0;
    eos.RYP2E(rho, massfrac, p, eint);
  }

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
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

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  static void add_problem_drive(
    amrex::DeriveList& derive_lst, amrex::DescriptorList& desc_lst)
  {
    derive_lst.add(
      "z_vorticity", amrex::IndexType::TheCellType(), 1, pc_dervortz,
      amrex::DeriveRec::GrowBoxByOne);
    derive_lst.addComponent("z_vorticity", desc_lst, State_Type, 0, NVAR);

#ifdef PELE_USE_REACTIONS
    derive_lst.add(
      "heat_release", amrex::IndexType::TheCellType(), 1, pc_derheatrelease,
      the_same_box);
    derive_lst.addComponent("heat_release", desc_lst, State_Type, 0, NVAR);
#endif
  }

  static void pc_dervortz(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int /*level*/)
  {
    auto const dat = datfab.array();
    auto vort = derfab.array();

    const amrex::Box& gbx = amrex::grow(bx, 1);

    amrex::FArrayBox local(gbx, 2, amrex::The_Async_Arena());
    auto larr = local.array();

    // Convert momentum to velocity.
    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
      larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
      larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
    });

    const amrex::Real dx = geomdata.CellSize(0);
    const amrex::Real dy = geomdata.CellSize(1);

    // Calculate vorticity.
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real vx =
        0.5 * (larr(i + 1, j, k, 1) - larr(i - 1, j, k, 1)) / dx;
      const amrex::Real uy =
        0.5 * (larr(i, j + 1, k, 0) - larr(i, j - 1, k, 0)) / dy;
      vort(i, j, k) = vx - uy;
    });
  }

#ifdef PELE_USE_REACTIONS
  static void pc_derheatrelease(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int /*level*/)
  {
    auto const dat = datfab.array();
    auto hr = derfab.array();

    // Calculate heat release.
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real rho = dat(i, j, k, URHO);
      amrex::Real T = dat(i, j, k, UTEMP);
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      const amrex::Real rhoInv = 1.0 / rho;
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
      amrex::Real wdot[NUM_SPECIES] = {0.0};
      amrex::Real hi[NUM_SPECIES] = {0.0};
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2WDOT(rho, T, massfrac, wdot);
      eos.T2Hi(T, hi);

      hr(i, j, k) = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++)
        hr(i, j, k) += wdot[n] * hi[n];
    });
  }
#endif
};

void pc_prob_close();

using ProblemSpecificFunctions = MyProblemSpecificFunctions;

#endif
