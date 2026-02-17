#ifndef PELEC_UTILITIES_H
#define PELEC_UTILITIES_H

#include <AMReX_IArrayBox.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include "Constants.H"
#include "IndexDefines.H"
#include "PelePhysics.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::RealVect
pc_cmp_loc(const amrex::IntVect& idx, const amrex::GeometryData& geomdata)
{
  const amrex::RealVect x(AMREX_D_DECL(
    geomdata.ProbLo(0) + (idx[0] + 0.5) * geomdata.CellSize(0),
    geomdata.ProbLo(1) + (idx[1] + 0.5) * geomdata.CellSize(1),
    geomdata.ProbLo(2) + (idx[2] + 0.5) * geomdata.CellSize(2)));
  return x;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::RealVect
get_rotaxis_vec(
  const amrex::IntVect& iv,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc,
  amrex::GeometryData const& geomdata)
{
  amrex::RealVect axloc(AMREX_D_DECL(axis_loc[0], axis_loc[1], axis_loc[2]));
  amrex::RealVect r(AMREX_D_DECL(0.0, 0.0, 0.0));
  amrex::RealVect x = pc_cmp_loc(iv, geomdata);
  r = x - axloc;
  return (r);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
get_rotaxis_dist(
  const amrex::IntVect& iv,
  int axis,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc,
  amrex::GeometryData const& geomdata)
{
  amrex::RealVect r = get_rotaxis_vec(iv, axis_loc, geomdata);
  amrex::Real rad2 = r.radSquared();
  rad2 -= r[axis] * r[axis]; // only in plane
  return (std::sqrt(rad2));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
get_rot_energy(
  const amrex::IntVect& iv,
  amrex::Real omega,
  int axis,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc,
  amrex::GeometryData const& geomdata)
{
  amrex::Real rad = get_rotaxis_dist(iv, axis, axis_loc, geomdata);
  return (0.5 * omega * omega * rad * rad);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_cmpTemp(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& S)
{
  amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
  amrex::Real T = S(i, j, k, UTEMP);
  amrex::Real e = S(i, j, k, UEINT) * rhoInv;
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = S(i, j, k, UFS + n) * rhoInv;
  }
  amrex::Real rho = S(i, j, k, URHO);
  auto eos = pele::physics::PhysicsType::eos();
  eos.REY2T(rho, e, massfrac, T);
  S(i, j, k, UTEMP) = T;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_rst_int_e(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& S,
  const bool allow_small_energy,
  const bool allow_negative_energy,
  const bool dual_energy_update_E_from_e,
  const amrex::Real dual_energy_eta2,
  amrex::Real omega = 0.0,
  amrex::Real rad = 0.0)
{
  if (std::is_same_v<
        pele::physics::PhysicsType::eos_type, pele::physics::eos::SRK>) {
    if ((!allow_small_energy) || (!allow_negative_energy)) {
      amrex::Error("Small/negative energy must be allowed for nonideal EOS");
    }
  }
  auto eos = pele::physics::PhysicsType::eos();
  if (!allow_small_energy) {
    const amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke =
      0.5 * (Up * Up + Vp * Vp + Wp * Wp) - 0.5 * omega * omega * rad * rad;
    const amrex::Real eden = S(i, j, k, UEDEN) * rhoInv;
    const amrex::Real eos_state_T = constants::very_small_num();
    amrex::Real eos_state_massfrac[NUM_SPECIES] = {0.0};
    amrex::Real eos_state_ei[NUM_SPECIES] = {0.0};
    eos.T2Ei(eos_state_T, eos_state_ei);
    amrex::Real eos_state_e = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
      eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
    }
    const amrex::Real small_e = eos_state_e;
    if (eden < small_e) {
      if (S(i, j, k, UEINT) * rhoInv < small_e) {
        const amrex::Real eos_state_T_loc = amrex::max<amrex::Real>(
          S(i, j, k, UTEMP), std::numeric_limits<amrex::Real>::min());
        eos.T2Ei(eos_state_T_loc, eos_state_ei);
        eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
      S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
    } else {
      const amrex::Real rho_eint = S(i, j, k, UEDEN) - S(i, j, k, URHO) * ke;
      if (
        (rho_eint >
         (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) &&
        ((rho_eint -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min()) /
         (S(i, j, k, UEDEN) -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) >
          dual_energy_eta2) {
        S(i, j, k, UEINT) = rho_eint;
      }
      if (S(i, j, k, UEINT) * rhoInv < small_e) {
        const amrex::Real eos_state_T_loc = amrex::max<amrex::Real>(
          S(i, j, k, UTEMP), std::numeric_limits<amrex::Real>::min());
        eos.T2Ei(eos_state_T_loc, eos_state_ei);
        eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        if (dual_energy_update_E_from_e) {
          S(i, j, k, UEDEN) =
            S(i, j, k, UEDEN) +
            (S(i, j, k, URHO) * eos_state_e - S(i, j, k, UEINT));
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
    }
  } else if (!allow_negative_energy) {
    const amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke =
      0.5 * (Up * Up + Vp * Vp + Wp * Wp) - 0.5 * omega * omega * rad * rad;
    if (
      S(i, j, k, UEDEN) <
      (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) {
      if (
        S(i, j, k, UEINT) <
        (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) {
        const amrex::Real eos_state_T = constants::very_small_num();
        amrex::Real eos_state_massfrac[NUM_SPECIES] = {0.0};
        amrex::Real eos_state_ei[NUM_SPECIES] = {0.0};
        eos.T2Ei(eos_state_T, eos_state_ei);
        amrex::Real eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
      S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
    } else {
      const amrex::Real rho_eint = S(i, j, k, UEDEN) - S(i, j, k, URHO) * ke;
      if (
        (rho_eint >
         (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) &&
        ((rho_eint -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min()) /
         (S(i, j, k, UEDEN) -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) >
          dual_energy_eta2) {
        S(i, j, k, UEINT) = rho_eint;
      } else if (
        S(i, j, k, UEINT) >
          (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min()) &&
        (dual_energy_update_E_from_e)) {
        S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
      } else if (
        S(i, j, k, UEINT) <=
        (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) {
        const amrex::Real eos_state_T = constants::very_small_num();
        amrex::Real eos_state_massfrac[NUM_SPECIES] = {0.0};
        amrex::Real eos_state_ei[NUM_SPECIES] = {0.0};
        eos.T2Ei(eos_state_T, eos_state_ei);
        amrex::Real eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        const amrex::Real eint_new = eos_state_e;
        if (dual_energy_update_E_from_e) {
          S(i, j, k, UEDEN) = S(i, j, k, UEDEN) +
                              (S(i, j, k, URHO) * eint_new - S(i, j, k, UEINT));
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eint_new;
      }
    }
  } else {
    const amrex::Real rho = S(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke =
      0.5 * (Up * Up + Vp * Vp + Wp * Wp) - 0.5 * omega * omega * rad * rad;
    S(i, j, k, UEINT) = S(i, j, k, UEDEN) - rho * ke;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& qa,
  amrex::Real omega = 0.0,
  amrex::Real rad = 0.0)
{
  auto eos = pele::physics::PhysicsType::eos();
  const amrex::Real rho = u(i, j, k, URHO);
  const amrex::Real rhoinv = 1.0 / rho;
  const amrex::Real vx = u(i, j, k, UMX) * rhoinv;
  const amrex::Real vy = u(i, j, k, UMY) * rhoinv;
  const amrex::Real vz = u(i, j, k, UMZ) * rhoinv;
  const amrex::Real kineng = 0.5 * rho * (vx * vx + vy * vy + vz * vz) -
                             0.5 * rho * omega * omega * rad * rad;
  q(i, j, k, QRHO) = rho;
  q(i, j, k, QU) = vx;
  q(i, j, k, QV) = vy;
  q(i, j, k, QW) = vz;
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    q(i, j, k, QFA + n) = u(i, j, k, UFA + n) / rho;
  }
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    q(i, j, k, QFS + n) = u(i, j, k, UFS + n) / rho;
  }
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    q(i, j, k, QFX + n) = u(i, j, k, UFX + n);
  }
#endif
  // Handle passive variables that are not per unit mass
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    q(i, j, k, QLIN + n) = u(i, j, k, ULIN + n);
  }
#endif

  const amrex::Real e = (u(i, j, k, UEDEN) - kineng) * rhoinv;
  amrex::Real T = u(i, j, k, UTEMP);
  amrex::Real massfrac[NUM_SPECIES];
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    if (
      (-1e-4 * std::numeric_limits<amrex::Real>::epsilon() <
       q(i, j, k, sp + QFS)) &&
      (q(i, j, k, sp + QFS) <
       1e-4 * std::numeric_limits<amrex::Real>::epsilon())) {
      q(i, j, k, sp + QFS) = 0.0;
    }
    massfrac[sp] = q(i, j, k, sp + QFS);
  }

  amrex::Real dpdr_e, dpde, gam1, cs, wbar, p;

  eos.Y2WBAR(massfrac, wbar);
  eos.REY2T(rho, e, massfrac, T);
  eos.RTY2P(rho, T, massfrac, p);
  eos.RTY2Cs(rho, T, massfrac, cs);
  eos.RTY2G(rho, T, massfrac, gam1);
  eos.RTY2dpde_dpdre(rho, T, massfrac, dpde, dpdr_e);

  q(i, j, k, QTEMP) = T;
  q(i, j, k, QREINT) = e * rho;
  q(i, j, k, QPRES) = p;
  q(i, j, k, QGAME) = p / (e * rho) + 1.0;

  // Auxiliary Fab
  qa(i, j, k, QDPDR) = dpdr_e;
  qa(i, j, k, QDPDE) = dpde;
  qa(i, j, k, QGAMC) = gam1;
  qa(i, j, k, QC) = cs;
  qa(i, j, k, QCSML) = amrex::max<amrex::Real>(
    constants::small_num(), constants::small_num() * cs);
  qa(i, j, k, QRSPEC) = pele::physics::Constants::RU / wbar;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_move_transcoefs_to_ec(
  AMREX_D_DECL(const int i, const int j, const int k),
  const int n,
  const amrex::Array4<const amrex::Real>& carr,
  amrex::Real* earr,
  const int dir,
  const bool harmonic_mean)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
  if (harmonic_mean) {
    const amrex::Real a = carr(iv, n);
    const amrex::Real b = carr(ivm, n);
    earr[n] = (a * b > 0.0) ? 2.0 * (a * b) / (a + b) : 0.0;
  } else {
    earr[n] = 0.5 * (carr(iv, n) + carr(ivm, n));
  }
}

AMREX_FORCE_INLINE
void
setV(
  const amrex::Box& box,
  const int Ncomp,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      array(i, j, k, n) = val;
    });
}

AMREX_FORCE_INLINE
void
setC(
  const amrex::Box& box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++) {
      array(i, j, k, n) = val;
    }
  });
}

AMREX_FORCE_INLINE
void
copy_array4(
  const amrex::Box& box,
  const int Ncomp,
  const amrex::Array4<const amrex::Real>& in,
  const amrex::Array4<amrex::Real>& out)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      out(i, j, k, n) = in(i, j, k, n);
    });
}

AMREX_FORCE_INLINE
void
lincomb_array4(
  const amrex::Box& box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
  const amrex::Real alpha,
  const amrex::Real beta,
  const amrex::Array4<amrex::Real>& d)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++) {
      d(i, j, k, n) = alpha * f1(i, j, k, n) + beta * f2(i, j, k, n);
    }
  });
}

AMREX_FORCE_INLINE
std::string
read_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

void read_binary(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<double>& data); // needs to be double

void read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<double>& data); // needs to be double

// Find position of element in vector
template <typename T>
int
find_position(const std::vector<T>& vec, const T& element)
{
  auto it = std::find(vec.begin(), vec.end(), element);
  if (it != vec.end()) {
    return distance(vec.begin(), it);
  }
  return -1;
}

std::string convertIntGG(int number);

// Clean the mass fractions on state, given a mask
void clean_massfrac(
  const amrex::Box& /*bx*/,
  const amrex::Real /*threshold*/,
  amrex::Array4<const int> const& /*mask*/,
  amrex::Array4<amrex::Real> const& /*state*/);

#endif
