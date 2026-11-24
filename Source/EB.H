#ifndef EB_H
#define EB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Extrusion.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_IF_Difference.H>
#include <AMReX_EBFluxRegister_C.H>
#include <AMReX_EB_utils.H>

#include "PeleC.H"
#include "IndexDefines.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int
argmax(const amrex::Real array[AMREX_SPACEDIM])
{
  const auto max_num =
    amrex::max<amrex::Real>(AMREX_D_DECL(array[0], array[1], array[2]));
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    if (std::abs(max_num - array[i]) < constants::very_small_num()) {
      return i;
    }
  }
  return 0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
idxsort(const amrex::Real array[AMREX_SPACEDIM], int idx[AMREX_SPACEDIM])
{
  amrex::Real absarray[AMREX_SPACEDIM] = {0.0};
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    absarray[dir] = std::abs(array[dir]);
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    const int i = argmax(absarray);
    idx[dir] = i;
    absarray[i] = -1.0;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_set_body_state(
  const int i,
  const int j,
  const int k,
  const int n,
  const amrex::Array4<const amrex::EBCellFlag>& flags,
  const amrex::GpuArray<amrex::Real, NVAR>& body_state,
  const amrex::Array4<amrex::Real>& S)
{
  if (flags(i, j, k).isCovered()) {
    S(i, j, k, n) = body_state[n];
  }
}

// least-squares specific functions
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_qmat(
  AMREX_D_DECL(
    amrex::Real delta_x_i[NLSPTS],
    amrex::Real delta_y_i[NLSPTS],
    amrex::Real delta_z_i[NLSPTS]),
  int N,
  amrex::Real qmat[NEL_TRIMAT])
{

  // note: all distances are non-dimensional ~ 1, does not depend on dx
  // therefore this choice of eps is good
  const amrex::Real eps = 1e-10;
  amrex::Real r11 = 0.0;
  amrex::Real r12 = 0.0;
  amrex::Real r22 = 0.0;

#if AMREX_SPACEDIM > 2
  amrex::Real r13 = 0.0;
  amrex::Real r23 = 0.0;
  amrex::Real r33 = 0.0;
#endif

  for (int i = 0; i < N; i++) {
    r11 += pow(delta_x_i[i], 2.0);
    r12 += delta_x_i[i] * delta_y_i[i];
    r22 += pow(delta_y_i[i], 2.0);

#if AMREX_SPACEDIM > 2
    r13 += delta_x_i[i] * delta_z_i[i];
    r23 += delta_y_i[i] * delta_z_i[i];
    r33 += pow(delta_z_i[i], 2.0);
#endif
  }

  r11 = sqrt(r11) + eps; // eps added to prevent division issues
  r12 = r12 / r11;
  r22 = sqrt(r22 - r12 * r12) + eps;

#if AMREX_SPACEDIM > 2
  r13 = r13 / r11;
  r23 = r23 / r22 - r12 * r13 / r22;
  r33 = sqrt(r33 - r13 * r13 - r23 * r23) + eps;
#endif

  qmat[0] = r11;
  qmat[1] = r12;
  qmat[2] = r22;

#if AMREX_SPACEDIM > 2
  qmat[3] = r13;
  qmat[4] = r23;
  qmat[5] = r33;
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_weightvec(
  AMREX_D_DECL(
    const amrex::Real delta_x_i[NLSPTS],
    const amrex::Real delta_y_i[NLSPTS],
    const amrex::Real delta_z_i[NLSPTS]),
  int N,
  const amrex::Real qmat[NEL_TRIMAT],
  amrex::Real wvec[NLSPTS][AMREX_SPACEDIM])
{
  amrex::Real r11 = qmat[0];
  amrex::Real r12 = qmat[1];
  amrex::Real r22 = qmat[2];

#if AMREX_SPACEDIM > 2
  amrex::Real r13 = qmat[3];
  amrex::Real r23 = qmat[4];
  amrex::Real r33 = qmat[5];
  amrex::Real beta = (r12 * r23 - r13 * r22) / (r11 * r22);
#endif

  for (int i = 0; i < N; i++) {
    amrex::Real alpha_ij_1 = delta_x_i[i] / (r11 * r11);
    amrex::Real alpha_ij_2 =
      (delta_y_i[i] - r12 / r11 * delta_x_i[i]) / (r22 * r22);

#if AMREX_SPACEDIM > 2
    amrex::Real alpha_ij_3 =
      (delta_z_i[i] - r23 / r22 * delta_y_i[i] + beta * delta_x_i[i]) /
      (r33 * r33);
#endif

#if AMREX_SPACEDIM == 2
    wvec[i][0] = alpha_ij_1 - r12 / r11 * alpha_ij_2;
    wvec[i][1] = alpha_ij_2;
#else
    wvec[i][0] = alpha_ij_1 - r12 / r11 * alpha_ij_2 + beta * alpha_ij_3;
    wvec[i][1] = alpha_ij_2 - r23 / r22 * alpha_ij_3;
    wvec[i][2] = alpha_ij_3;
#endif
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
is_cut_neighborhood(
  const amrex::IntVect& iv,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  const amrex::IntVect& ng = amrex::IntVect(1))
{
  // Check if there's a cut cell in an ng neighborhood around iv
  bool has_cut_cell = flags(iv).isSingleValued();
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    const auto dimvec = amrex::IntVect::TheDimensionVector(idir);
    for (int n = 1; n <= ng[idir]; n++) {
      has_cut_cell |= flags(iv - n * dimvec).isSingleValued() ||
                      flags(iv + n * dimvec).isSingleValued();
    }
  }
  return has_cut_cell;
}

void pc_fill_sv_ebg(
  const amrex::Box& /*bx*/,
  const int /*Nebg*/,
  const amrex::Array4<const amrex::Real>& /*vfrac*/,
  const amrex::Array4<const amrex::Real>& /*bcent*/,
  AMREX_D_DECL(
    const amrex::Array4<const amrex::Real>&,
    const amrex::Array4<const amrex::Real>&,
    const amrex::Array4<const amrex::Real>&),
  EBBndryGeom* /*ebg*/);

void pc_fill_bndry_grad_stencil_ls(
  const amrex::Box& /*bx*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*dx*/,
  const int /*unused*/,
  const EBBndryGeom* /*ebg*/,
  const int /*Nsten*/,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  EBBndrySten* /*grad_stencil*/);

void pc_check_bndry_grad_stencil(
  const amrex::Box& /*bx*/,
  const int /*Nsten*/,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const EBBndrySten* /*grad_stencil*/);

void pc_fill_bndry_grad_stencil_quadratic(
  const amrex::Box& /*bx*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*dx*/,
  const int /*unused*/,
  const EBBndryGeom* /*ebg*/,
  const int /*Nsten*/,
  EBBndrySten* /*grad_stencil*/);

void pc_fill_flux_interp_stencil(
  const amrex::Box& /*bx*/,
  const int /*Nsten*/,
  const amrex::Array4<const amrex::Real>& /*fc*/,
  const amrex::Array4<const amrex::Real>& /*fa*/,
  FaceSten* /*sten*/);

void pc_apply_face_stencil(
  const amrex::Box& /*bx*/,
  const amrex::Box /*unused*/,
  const FaceSten* /*sten*/,
  const int /*Nsten*/,
  const int /*dir*/,
  const int /*nc*/,
  const amrex::Array4<amrex::Real>& /*vout*/);

void pc_eb_div(
  const amrex::Box& /*bx*/,
  const amrex::Real /*vol*/,
  const int /*nc*/,
  const EBBndryGeom* /*sv_ebg*/,
  const int /*Ncut*/,
  AMREX_D_DECL(
    const amrex::Array4<const amrex::Real>&,
    const amrex::Array4<const amrex::Real>&,
    const amrex::Array4<const amrex::Real>&),
  const amrex::Real* /*ebflux*/,
  const amrex::Array4<const amrex::Real>& /*vf*/,
  const amrex::Array4<amrex::Real>& /*DC*/);

void pc_apply_eb_boundry_visc_flux_stencil(
  const amrex::Box& /*bx*/,
  const EBBndrySten* /*sten*/,
  const int /*Nsten*/,
  const EBBndryGeom* /*ebg*/,
  const int /*unused*/,
  amrex::Array4<const amrex::Real> const& /*q*/,
  amrex::Array4<const amrex::Real> const& /*coeff*/,
  const amrex::Real* /*bcval*/,
  const int /*unused*/,
  amrex::Real* /*bcflux*/,
  const int /*Nflux*/);

void pc_apply_eb_boundry_flux_stencil(
  const amrex::Box& /*bx*/,
  const EBBndrySten* /*sten*/,
  const int /*Nsten*/,
  amrex::Array4<const amrex::Real> const& /*s*/,
  const int /*scomp*/,
  amrex::Array4<const amrex::Real> const& /*D*/,
  const int /*Dcomp*/,
  const amrex::Real* /*bcval*/,
  const int /*unused*/,
  amrex::Real* /*bcflux*/,
  const int /*Nflux*/,
  const int /*nc*/);

void pc_eb_clean_massfrac(
  const amrex::Box& /*bx*/,
  const amrex::Real /*dt*/,
  const amrex::Real /*threshold*/,
  amrex::Array4<const amrex::Real> const& /*state*/,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<amrex::Real> const& /*scratch*/,
  amrex::Array4<amrex::Real> const& /*div*/);

void pc_post_eb_redistribution(
  const amrex::Box& /*bx*/,
  const amrex::Real /*dt*/,
  const bool /*eb_clean_massfrac*/,
  const amrex::Real /*threshold*/,
  amrex::Array4<const amrex::Real> const& /*state*/,
  const amrex::FabType& /*typ*/,
  const amrex::Array4<amrex::EBCellFlag const>& /*flags*/,
  amrex::Array4<amrex::Real> const& /*scratch*/,
  amrex::Array4<amrex::Real> const& /*div*/);

#endif
