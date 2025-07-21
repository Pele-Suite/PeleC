#ifndef DERIVE_H
#define DERIVE_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#ifdef PELE_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "EB.H"
#include "PeleCUtilities.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_idx(const int i, int& im, int& ip)
{
  im = i - 1;
  ip = i + 1;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_idx(
  const int i,
  const int dir,
  const bool all_regular,
  const amrex::EBCellFlag& flag,
  int& im,
  int& ip)
{
  if (all_regular) {
    get_idx(i, im, ip);
  } else {
    if (flag.isCovered()) {
      im = i;
      ip = i;
    } else {
      const amrex::IntVect ivm = -amrex::IntVect::TheDimensionVector(dir);
      const amrex::IntVect ivp = amrex::IntVect::TheDimensionVector(dir);
      im = i - static_cast<int>(flag.isConnected(ivm));
      ip = i + static_cast<int>(flag.isConnected(ivp));
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
get_weight(const int im, const int ip)
{
  const int diff = ip - im;
  return diff == 0 ? 0.0 : (diff == 1 ? 1.0 : 0.5);
}

template <int state_idx, int num_comp>
void
pc_derdividebyrho(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto derived = derfab.array();

  amrex::ParallelFor(
    bx, num_comp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      derived(i, j, k, n) = dat(i, j, k, state_idx + n) / dat(i, j, k, URHO);
    });
}

template <int dir>
void
PeleC::pc_dervel_if(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
#if AMREX_SPACEDIM > 2
  const amrex::Geometry& geomdata,
#else
  const amrex::Geometry& /*geomdata*/,
#endif
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.const_array();
  auto vel = derfab.array();
  int rotframeflag = do_rf;
  if (rotframeflag) {
#if AMREX_SPACEDIM > 2
    int axis = rf_axis;
    amrex::Real omega = rf_omega;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};
    auto const gdata = geomdata.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::IntVect iv(AMREX_D_DECL(i, j, k));
      amrex::RealVect r = get_rotaxis_vec(iv, axis_loc, gdata);
      amrex::RealVect w(AMREX_D_DECL(0.0, 0.0, 0.0));
      w[axis] = omega;
      amrex::RealVect w_cross_r = w.crossProduct(r);
      vel(i, j, k) =
        dat(i, j, k, UMX + dir) / dat(i, j, k, URHO) + w_cross_r[dir];
    });
#endif
  } else {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      vel(i, j, k) = dat(i, j, k, UMX + dir) / dat(i, j, k, URHO);
    });
  }
}

void pc_dermagvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermagmom(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derkineng(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derlogden(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermagvort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derdivu(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derenstrophy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dernull(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermolefrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derentropy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dermachnumber(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dertemp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derspectrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level,
  const int idx);

void pc_derradialvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dercp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dercv(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dercoord(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

#ifdef PELE_USE_MASA
void pc_derrhommserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derummserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_dervmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derwmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void pc_derpmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);
#endif
#endif
