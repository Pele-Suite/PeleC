// ============================================================================
//  nscbc1d -- standalone 1-D verification driver for Source/NSCBC.H
//
//  This driver compiles the production NSCBC kernel unmodified against a
//  minimal 1-D finite-volume Euler solver.  Its purpose is to separate two
//  failure modes that are otherwise indistinguishable in a full PeleC run:
//
//      "the boundary condition is wrong"   vs   "the AMReX plumbing is wrong"
//
//  It is cheap enough to sweep parameters, so it also produces the reference
//  curves (reflection coefficient vs sigma, relaxation rate vs resolution)
//  that the AMReX-side regression tests are checked against.
//
//  Build:  see CMakeLists.txt in this directory.
//  Run:    ./nscbc1d              (runs every check, prints a summary)
//          ./nscbc1d sweep        (also dumps the sigma sweep as CSV)
//
//  Solver: MUSCL-Hancock, minmod-limited primitive reconstruction, HLLC flux,
//  SSP-RK2.  Deliberately simple -- it is the boundary condition under test,
//  not the interior scheme.
// ============================================================================

#include <AMReX.H>
#include <AMReX_Print.H>

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "NSCBC.H"

using amrex::Real;

namespace {

constexpr int NG = 4; // ghost layers, matching PeleC::numGrow() without EB

// ---------------------------------------------------------------------------
//  A 1-D state vector, laid out exactly like PeleC's conserved state so that
//  the kernel sees the indices it expects.
// ---------------------------------------------------------------------------
struct Field
{
  int n;
  std::vector<Real> s; // (n + 2*NG) * NVAR

  explicit Field(int n_) : n(n_), s(static_cast<size_t>(n + 2 * NG) * NVAR, 0.0)
  {
  }
  Real* at(int i) { return &s[static_cast<size_t>(i + NG) * NVAR]; }
  const Real* at(int i) const { return &s[static_cast<size_t>(i + NG) * NVAR]; }
};

struct Prim
{
  Real rho, u, p, T;
  Real Y[NUM_SPECIES];
};

void
set_state(Real* s, Real rho, Real u, Real T, const Real Y[NUM_SPECIES])
{
  auto eos = pele::physics::PhysicsType::eos();
  Real e = 0.0;
  eos.RTY2E(rho, T, Y, e);
  s[URHO] = rho;
  s[UMX] = rho * u;
  s[UMY] = 0.0;
  s[UMZ] = 0.0;
  s[UEINT] = rho * e;
  s[UEDEN] = rho * (e + 0.5 * u * u);
  s[UTEMP] = T;
  for (int k = 0; k < NUM_SPECIES; k++) {
    s[UFS + k] = rho * Y[k];
  }
#if NUM_ADV > 0
  // Nonzero, distinct, per-unit-mass values, so that a build with passive
  // scalars actually exercises pack_ghost's UFA handling: C1's uniform-state
  // loop runs over ALL NVAR components and catches any slot the fill drops
  // or mis-scales.
  for (int k = 0; k < NUM_ADV; k++) {
    s[UFA + k] = rho * 0.1 * static_cast<Real>(k + 1);
  }
#endif
}

Prim
get_prim(const Real* s)
{
  auto eos = pele::physics::PhysicsType::eos();
  Prim q{};
  q.rho = s[URHO];
  q.u = s[UMX] / q.rho;
  q.T = s[UTEMP];
  Real ysum = 0.0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    q.Y[k] = s[UFS + k] / q.rho;
    ysum += q.Y[k];
  }
  for (int k = 0; k < NUM_SPECIES; k++) {
    q.Y[k] /= ysum;
  }
  eos.RTY2P(q.rho, q.T, q.Y, q.p);
  return q;
}

Real
sound_speed(const Prim& q)
{
  auto eos = pele::physics::PhysicsType::eos();
  Real c = 0.0;
  eos.RTY2Cs(q.rho, q.T, q.Y, c);
  return c;
}

// ---------------------------------------------------------------------------
//  Boundary fill: exactly the dispatch the AMReX-side BCfill.cpp will do.
//  layer runs 1..NG outward; the stencil walks inward from the boundary cell.
// ---------------------------------------------------------------------------
void
fill_bcs(
  Field& f,
  const pc_nscbc::Target& lo_tgt,
  const pc_nscbc::Target& hi_tgt,
  const pc_nscbc::Params& prm,
  Real dx,
  amrex::Long* diag = nullptr)
{
  const int n = f.n;
  for (int layer = 1; layer <= NG; layer++) {
    if (lo_tgt.type != pc_nscbc::Type::off) {
      pc_nscbc::apply(
        f.at(0), f.at(1), f.at(2), 3, dx, /*idir=*/0, /*sgn=*/+1, layer, lo_tgt,
        prm, f.at(-layer), diag);
    } else {
      for (int v = 0; v < NVAR; v++) {
        f.at(-layer)[v] = f.at(0)[v];
      }
    }
    if (hi_tgt.type != pc_nscbc::Type::off) {
      pc_nscbc::apply(
        f.at(n - 1), f.at(n - 2), f.at(n - 3), 3, dx, /*idir=*/0, /*sgn=*/-1,
        layer, hi_tgt, prm, f.at(n - 1 + layer), diag);
    } else {
      for (int v = 0; v < NVAR; v++) {
        f.at(n - 1 + layer)[v] = f.at(n - 1)[v];
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  HLLC flux for the multi-species Euler equations.
// ---------------------------------------------------------------------------
void
hllc(const Prim& L, const Prim& R, Real* flx)
{
  auto eos = pele::physics::PhysicsType::eos();
  const Real cL = sound_speed(L), cR = sound_speed(R);
  const Real sL = std::min(L.u - cL, R.u - cR);
  const Real sR = std::max(L.u + cL, R.u + cR);
  const Real sM =
    (R.p - L.p + L.rho * L.u * (sL - L.u) - R.rho * R.u * (sR - R.u)) /
    (L.rho * (sL - L.u) - R.rho * (sR - R.u));

  auto cons_flux = [&](const Prim& q, Real* F) {
    Real e = 0.0;
    eos.RTY2E(q.rho, q.T, q.Y, e);
    const Real E = q.rho * (e + 0.5 * q.u * q.u);
    F[URHO] = q.rho * q.u;
    F[UMX] = q.rho * q.u * q.u + q.p;
    F[UMY] = 0.0;
    F[UMZ] = 0.0;
    F[UEDEN] = (E + q.p) * q.u;
    F[UEINT] = q.rho * e * q.u;
    F[UTEMP] = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      F[UFS + k] = q.rho * q.Y[k] * q.u;
    }
  };
  auto star_state = [&](const Prim& q, Real s, Real* U) {
    Real e = 0.0;
    eos.RTY2E(q.rho, q.T, q.Y, e);
    const Real E = q.rho * (e + 0.5 * q.u * q.u);
    const Real fac = (s - q.u) / (s - sM);
    const Real rs = q.rho * fac;
    U[URHO] = rs;
    U[UMX] = rs * sM;
    U[UMY] = 0.0;
    U[UMZ] = 0.0;
    U[UEDEN] = rs * (E / q.rho + (sM - q.u) * (sM + q.p / (q.rho * (s - q.u))));
    U[UEINT] = rs * e;
    U[UTEMP] = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      U[UFS + k] = rs * q.Y[k];
    }
  };

  Real FL[NVAR], FR[NVAR], UL[NVAR], UR[NVAR], US[NVAR];
  cons_flux(L, FL);
  cons_flux(R, FR);
  auto cons_state = [&](const Prim& q, Real* U) {
    Real e = 0.0;
    eos.RTY2E(q.rho, q.T, q.Y, e);
    U[URHO] = q.rho;
    U[UMX] = q.rho * q.u;
    U[UMY] = 0.0;
    U[UMZ] = 0.0;
    U[UEDEN] = q.rho * (e + 0.5 * q.u * q.u);
    U[UEINT] = q.rho * e;
    U[UTEMP] = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      U[UFS + k] = q.rho * q.Y[k];
    }
  };
  cons_state(L, UL);
  cons_state(R, UR);

  if (sL >= 0.0) {
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FL[v];
    }
  } else if (sR <= 0.0) {
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FR[v];
    }
  } else if (sM >= 0.0) {
    star_state(L, sL, US);
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FL[v] + sL * (US[v] - UL[v]);
    }
  } else {
    star_state(R, sR, US);
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FR[v] + sR * (US[v] - UR[v]);
    }
  }
}

Real
mm(Real a, Real b)
{
  return (a * b <= 0.0) ? 0.0 : ((std::abs(a) < std::abs(b)) ? a : b);
}

// One SSP-RK2 stage: conserved-variable MUSCL with minmod slopes.
// Optional constant-conductivity conduction in the mini solver, for C12: the
// energy flux gains -lambda dT/dx at every face.  Zero (off) everywhere else.
Real g_lambda = 0.0;

void
stage(const Field& in, Field& out, Real dx, Real dt)
{
  const int n = in.n;
  std::vector<Real> flux(static_cast<size_t>(n + 1) * NVAR, 0.0);
  for (int i = 0; i <= n; i++) {
    Real sl[NVAR], sr[NVAR];
    for (int v = 0; v < NVAR; v++) {
      const Real dL =
        mm(in.at(i - 1)[v] - in.at(i - 2)[v], in.at(i)[v] - in.at(i - 1)[v]);
      const Real dR =
        mm(in.at(i)[v] - in.at(i - 1)[v], in.at(i + 1)[v] - in.at(i)[v]);
      sl[v] = in.at(i - 1)[v] + 0.5 * dL;
      sr[v] = in.at(i)[v] - 0.5 * dR;
    }
    // Recover T for each face state so the EOS is self-consistent.
    auto eos = pele::physics::PhysicsType::eos();
    auto to_prim = [&](Real* s) {
      Prim q{};
      q.rho = std::max(s[URHO], 1e-12);
      q.u = s[UMX] / q.rho;
      Real ys = 0.0;
      for (int k = 0; k < NUM_SPECIES; k++) {
        q.Y[k] = std::max(s[UFS + k] / q.rho, 0.0);
        ys += q.Y[k];
      }
      for (int k = 0; k < NUM_SPECIES; k++) {
        q.Y[k] /= ys;
      }
      const Real e = s[UEDEN] / q.rho - 0.5 * q.u * q.u;
      q.T = s[UTEMP] > 0.0 ? s[UTEMP] : 300.0;
      eos.REY2T(q.rho, e, q.Y, q.T);
      eos.RTY2P(q.rho, q.T, q.Y, q.p);
      return q;
    };
    hllc(to_prim(sl), to_prim(sr), &flux[static_cast<size_t>(i) * NVAR]);
    if (g_lambda > 0.0) {
      // Conduction between the adjacent CELL CENTRES, like PeleC's diffusion
      // operator: at i = 0 and i = n this reads a ghost temperature, so the
      // outflow's ghost closure IS the boundary heat flux here too.
      // UEDEN only: stage() recomputes UEINT from UEDEN afterwards.
      flux[static_cast<size_t>(i) * NVAR + UEDEN] -=
        g_lambda * (in.at(i)[UTEMP] - in.at(i - 1)[UTEMP]) / dx;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int v = 0; v < NVAR; v++) {
      out.at(i)[v] =
        in.at(i)[v] - dt / dx *
                        (flux[static_cast<size_t>(i + 1) * NVAR + v] -
                         flux[static_cast<size_t>(i) * NVAR + v]);
    }
    // keep UTEMP consistent
    auto eos = pele::physics::PhysicsType::eos();
    Real* s = out.at(i);
    const Real rho = s[URHO];
    const Real u = s[UMX] / rho;
    Real Y[NUM_SPECIES], ys = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      Y[k] = std::max(s[UFS + k] / rho, 0.0);
      ys += Y[k];
    }
    for (int k = 0; k < NUM_SPECIES; k++) {
      Y[k] /= ys;
    }
    const Real e = s[UEDEN] / rho - 0.5 * u * u;
    Real T = s[UTEMP] > 0.0 ? s[UTEMP] : 300.0;
    eos.REY2T(rho, e, Y, T);
    s[UTEMP] = T;
    s[UEINT] = rho * e;
  }
}

// PeleC and PelePhysics work in CGS: lengths in cm, pressures in dyn/cm^2,
// velocities in cm/s, densities in g/cm^3.  Constants::PATM = 1.01325e6.
struct Case
{
  int n = 400;
  Real L = 10.0;       // cm
  Real p0 = 1.01325e6; // dyn/cm^2 (1 atm)
  Real T0 = 300.0;     // K
  Real u0 = 0.0;       // cm/s
  Real cfl = 0.4;
};

// ---------------------------------------------------------------------------
//  Reporting
// ---------------------------------------------------------------------------
int n_pass = 0, n_fail = 0;

void
check(bool ok, const std::string& name, const std::string& detail)
{
  if (ok) {
    n_pass++;
    std::printf("  PASS  %-46s %s\n", name.c_str(), detail.c_str());
  } else {
    n_fail++;
    std::printf("  FAIL  %-46s %s\n", name.c_str(), detail.c_str());
  }
}

// A measurement worth recording that is not a pass/fail criterion: either the
// prediction it would test cannot be isolated in the configuration available,
// or the number is diagnostic rather than normative.
void
report(const std::string& name, const std::string& detail)
{
  std::printf("  ....  %-46s %s\n", name.c_str(), detail.c_str());
}

// Air, for whatever mechanism this was built against.  Must not assume a
// two-species mechanism: with LiDryer's nine species the old form
// ("0.233 if O2 else 0.767") sums to well over one and every check downstream
// fails for reasons that have nothing to do with the boundary condition.
Real
air_Y(int k)
{
  if (k == O2_ID) {
    return 0.233;
  }
  if (k == N2_ID) {
    return 0.767;
  }
  return 0.0;
}

} // namespace

// ===========================================================================
//  Checks
// ===========================================================================

// C1: a uniform state must be reproduced exactly in every ghost layer, for
//     both inflow and outflow, at any sigma.  If this fails, the kernel is
//     manufacturing a gradient out of nothing and everything downstream is
//     meaningless.
void
check_uniform()
{
  Case cs;
  Field f(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  for (int i = -NG; i < cs.n + NG; i++) {
    set_state(f.at(i), rho0, cs.u0 + 2.0e3, cs.T0, Y); // 20 m/s
  }
  const Real dx = cs.L / cs.n;

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  pc_nscbc::Target off, out, in;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;
  in.type = pc_nscbc::Type::inflow;
  in.u[0] = 2.0e3;
  in.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    in.Y[k] = Y[k];
  }

  for (bool mat : {false, true}) {
    prm.extrap_material = mat;
    for (Real sig : {0.0, 0.25, 1.0}) {
      prm.sigma = sig;
      prm.relax_u = sig;
      prm.relax_t = sig;
      Field g = f;
      fill_bcs(g, in, out, prm, dx);
      Real worst = 0.0;
      for (int layer = 1; layer <= NG; layer++) {
        for (int v = 0; v < NVAR; v++) {
          const Real ref = std::max(std::abs(f.at(0)[v]), 1.0);
          worst = std::max(worst, std::abs(g.at(-layer)[v] - f.at(0)[v]) / ref);
          worst = std::max(
            worst,
            std::abs(g.at(cs.n - 1 + layer)[v] - f.at(cs.n - 1)[v]) / ref);
        }
      }
      char buf[160];
      std::snprintf(
        buf, sizeof(buf), "sigma=%.2f mat=%d  max rel ghost error = %.3e", sig,
        static_cast<int>(mat), worst);
#ifdef USE_SRK_EOS
      // SRK's (rho,e,Y)->T and (rho,Y,p)->T are Newton solves; round-trip
      // convergence is ~1e-11, not machine epsilon.  The check's meaning --
      // no manufactured gradients -- survives at that floor.
      check(worst < 1e-9, "uniform state is reproduced exactly", buf);
#else
      check(worst < 1e-12, "uniform state is reproduced exactly", buf);
#endif
    }
  }
}

// C2: the relaxation must move the boundary TOWARD the target, for every
//     target quantity and on both faces.  This is the assertion that catches
//     the sign errors the legacy Fortran leaked to its users as "relax_T must
//     be negative".
void
check_relaxation_signs()
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 0.25;
  prm.relax_u = 2.0;
  prm.relax_t = 0.2;
  prm.order = 1; // isolate the relaxation from the extrapolation

  // --- outflow: interior pressure above target -> ghost pressure below the
  //     non-relaxed value (i.e. pulled toward the target)
  for (int sgn : {+1, -1}) {
    Real rho = 0.0, e = 0.0;
    const Real p_int = 1.05 * cs.p0;
    eos.PYT2RE(p_int, Y, cs.T0, rho, e);
    Real sN[NVAR], sg_relaxed[NVAR], sg_free[NVAR];
    set_state(sN, rho, (sgn == +1) ? -3.0e3 : 3.0e3, cs.T0, Y); // outflow
    pc_nscbc::Target t;
    t.type = pc_nscbc::Type::outflow;
    t.p = cs.p0;
    pc_nscbc::apply(sN, sN, sN, 1, dx, 0, sgn, 1, t, prm, sg_relaxed);
    pc_nscbc::Params p0prm = prm;
    p0prm.sigma = 0.0;
    pc_nscbc::apply(sN, sN, sN, 1, dx, 0, sgn, 1, t, p0prm, sg_free);
    const Prim qr = get_prim(sg_relaxed), qf = get_prim(sg_free);
    char buf[160];
    std::snprintf(
      buf, sizeof(buf), "%s face: p_ghost %.3f vs %.3f dyn/cm2 unrelaxed",
      (sgn == +1 ? "lo" : "hi"), qr.p, qf.p);
    check(qr.p < qf.p, "outflow sigma pulls pressure toward target", buf);
  }

  // --- inflow: interior normal velocity above target -> ghost moves down
  //     toward it; interior T above target -> ghost T moves down.
  for (int sgn : {+1, -1}) {
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(cs.p0, Y, 400.0, rho, e);
    const Real u_int = (sgn == +1) ? 4.0e3 : -4.0e3; // inflow through the face
    Real sN[NVAR], sg[NVAR];
    set_state(sN, rho, u_int, 400.0, Y);
    pc_nscbc::Target t;
    t.type = pc_nscbc::Type::inflow;
    t.u[0] = (sgn == +1) ? 2.0e3 : -2.0e3; // want less inflow
    t.T = 300.0;                           // want colder
    for (int k = 0; k < NUM_SPECIES; k++) {
      t.Y[k] = Y[k];
    }
    pc_nscbc::apply(sN, sN, sN, 1, dx, 0, sgn, 1, t, prm, sg);
    const Prim q = get_prim(sg);
    // "toward the target" in the outward frame
    const Real n_sgn = -static_cast<Real>(sgn);
    const Real uo_int = n_sgn * u_int, uo_g = n_sgn * q.u,
               uo_t = n_sgn * t.u[0];
    char buf[200];
    std::snprintf(
      buf, sizeof(buf), "%s face: u_out %.2f -> %.2f cm/s (target %.2f)",
      (sgn == +1 ? "lo" : "hi"), uo_int, uo_g, uo_t);
    check(
      (uo_g - uo_t) * (uo_int - uo_t) >= 0.0 &&
        std::abs(uo_g - uo_t) < std::abs(uo_int - uo_t),
      "inflow relax_u moves velocity toward target", buf);
    std::snprintf(
      buf, sizeof(buf), "%s face: T %.3f -> %.3f (target %.1f)",
      (sgn == +1 ? "lo" : "hi"), 400.0, q.T, t.T);
    check(
      q.T < 400.0 && q.T > t.T,
      "inflow relax_t moves temperature toward target", buf);
  }
}

// C3: species handling.  At an outflow a composition gradient must be
//     extrapolated, NOT clamped to any target -- the legacy code imposed the
//     target composition at outflows, which over-specifies the problem.  At an
//     inflow the target composition must be imposed exactly.
void
check_species()
{
  Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 0.25;

  // A linear composition ramp running out through a hi outflow face.
  Real sN[NVAR], sNm1[NVAR], sNm2[NVAR], sg[NVAR];
  auto mk = [&](Real* s, Real yO2) {
    // Zero-initialised: with a mechanism larger than air's two species, the
    // remaining entries would otherwise be stack garbage -- which is exactly
    // the "nothing here may assume a particular mechanism" trap the README
    // documents, and it made this check fail under LiDryer while passing
    // under air.
    Real Y[NUM_SPECIES] = {0.0};
    Y[O2_ID] = yO2;
    Y[N2_ID] = 1.0 - yO2;
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(cs.p0, Y, cs.T0, rho, e);
    set_state(s, rho, 3.0e3, cs.T0, Y);
  };
  mk(sN, 0.30);
  mk(sNm1, 0.28);
  mk(sNm2, 0.26);
  pc_nscbc::Target t;
  t.type = pc_nscbc::Type::outflow;
  t.p = cs.p0;
  pc_nscbc::apply(sN, sNm1, sNm2, 3, dx, 0, -1, 1, t, prm, sg);
  const Prim q = get_prim(sg);
  char buf[200];
  std::snprintf(
    buf, sizeof(buf), "Y(O2) 0.26,0.28,0.30 -> ghost %.6f (expect 0.32)",
    q.Y[O2_ID]);
  check(
    std::abs(q.Y[O2_ID] - 0.32) < 1e-6,
    "outflow extrapolates composition (not imposed)", buf);

  // Sum of mass fractions must be exactly one after extrapolation+clipping.
  Real ysum = 0.0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    ysum += sg[UFS + k];
  }
  std::snprintf(
    buf, sizeof(buf), "|sum(rhoY)/rho - 1| = %.3e",
    std::abs(ysum / sg[URHO] - 1.0));
  check(
    std::abs(ysum / sg[URHO] - 1.0) < 1e-14, "sum(Y) == 1 after outflow fill",
    buf);

  // Inflow imposes the target exactly.
  pc_nscbc::Target ti;
  ti.type = pc_nscbc::Type::inflow;
  ti.T = cs.T0;
  ti.u[0] = -3.0e3;
  ti.Y[O2_ID] = 1.0;
  ti.Y[N2_ID] = 0.0;
  Real sIn[NVAR];
  mk(sIn, 0.30);
  set_state(sIn, get_prim(sIn).rho, -3.0e3, cs.T0, [&] {
    static Real Y[NUM_SPECIES];
    Y[O2_ID] = 0.30;
    Y[N2_ID] = 0.70;
    return Y;
  }());
  pc_nscbc::apply(sIn, sIn, sIn, 1, dx, 0, -1, 1, ti, prm, sg);
  const Prim qi = get_prim(sg);
  std::snprintf(
    buf, sizeof(buf), "interior Y(O2)=0.30, target 1.0 -> ghost %.6f",
    qi.Y[O2_ID]);
  check(
    std::abs(qi.Y[O2_ID] - 1.0) < 1e-12, "inflow imposes composition exactly",
    buf);

  // Energy identity must hold exactly.
  const Real ke = 0.5 * sg[UMX] * sg[UMX] / sg[URHO];
  const Real resid =
    std::abs(sg[UEDEN] - sg[UEINT] - ke) / std::max(std::abs(sg[UEDEN]), 1.0);
  std::snprintf(buf, sizeof(buf), "|UEDEN - UEINT - KE|/UEDEN = %.3e", resid);
  check(resid < 1e-14, "UEDEN == UEINT + KE exactly", buf);
}

// C4: acoustic reflection.  Launch a Gaussian pressure pulse at an outflow and
//     measure the amplitude that comes back.  Reports the sigma sweep, which
//     is the curve the AMReX regression test is checked against.
Real
reflection_coefficient(
  Real sigma, int n, int order, bool pin, Real* p_drift, bool mat = false,
  bool ndnr = false)
{
  Case cs;
  cs.n = n;
  Field f(cs.n), g(cs.n), h(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;

  const Real amp = 1.0e-3; // 0.1% pressure pulse -- linear regime
  const Real x0 = 0.35 * cs.L, w = 0.04 * cs.L;
  Real rho_ref = 0.0, e_ref = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho_ref, e_ref);
  const Prim qref = [&] {
    Prim q{};
    q.rho = rho_ref;
    q.u = 0.0;
    q.p = cs.p0;
    q.T = cs.T0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      q.Y[k] = Y[k];
    }
    return q;
  }();
  const Real c0 = sound_speed(qref);

  for (int i = -NG; i < cs.n + NG; i++) {
    const Real x = (i + 0.5) * dx;
    const Real dp = amp * cs.p0 * std::exp(-std::pow((x - x0) / w, 2));
    // Right-running isentropic acoustic pulse.
    const Real p = cs.p0 + dp;
    const Real rho = rho_ref + dp / (c0 * c0);
    const Real u = dp / (rho_ref * c0);
    Real T = 0.0;
    eos.RYP2T(rho, Y, p, T);
    set_state(f.at(i), rho, u, T, Y);
  }

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = sigma;
  prm.order = order;
  prm.pin_farfield = pin;
  prm.extrap_material = mat;
  pc_nscbc::Target off, out;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;

  const Real t_end = 1.6 * cs.L / c0; // pulse out, any reflection back in
  Real t = 0.0;
  Real peak_in = amp * cs.p0, peak_out = 0.0;
  // measurement window: the left half, after the pulse has left
  const Real t_measure = 0.9 * cs.L / c0;

  // NDNR register (queue item 4, phase B): a driver-held EMA of the
  // boundary-cell pressure, tau = 3 t_a, updated once per step outside the
  // stages.  The kernel sees a composed target p_eff = p_b - EMA + p0, so
  // its relaxation K (p_b - p_eff) = K (EMA - p0) acts on the slow mean
  // only -- the acoustic fluctuation is stripped before sigma touches it.
  Real ema_p = cs.p0;
  const Real tau_ema = 3.0 * (2.0 * cs.L / c0);

  while (t < t_end) {
    Real cmax = 0.0;
    for (int i = 0; i < cs.n; i++) {
      const Prim q = get_prim(f.at(i));
      cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
    }
    Real dt = cs.cfl * dx / cmax;
    dt = std::min(dt, t_end - t);
    if (dt <= 0.0) {
      break;
    }

    pc_nscbc::Target out_eff = out;
    if (ndnr) {
      const Prim qb = get_prim(f.at(cs.n - 1));
      const Real w = dt / (tau_ema + dt);
      ema_p += w * (qb.p - ema_p);
      out_eff.p = qb.p - ema_p + cs.p0;
    }

    fill_bcs(f, off, out_eff, prm, dx);
    stage(f, g, dx, dt);
    fill_bcs(g, off, out_eff, prm, dx);
    stage(g, h, dx, dt);
    for (int i = 0; i < cs.n; i++) {
      for (int v = 0; v < NVAR; v++) {
        f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
      }
    }
    t += dt;

    if (t > t_measure) {
      // Measure the wave content only: the deviation from the instantaneous
      // domain mean.  Measuring against p0 instead would count the sigma
      // anchoring transient (a uniform, non-propagating pressure adjustment)
      // as if it were a reflected wave, which it is not.
      Real pbar = 0.0;
      for (int i = 0; i < cs.n; i++) {
        pbar += get_prim(f.at(i)).p;
      }
      pbar /= cs.n;
      for (int i = 0; i < cs.n / 2; i++) {
        peak_out = std::max(peak_out, std::abs(get_prim(f.at(i)).p - pbar));
      }
    }
  }
  if (p_drift != nullptr) {
    Real pm = 0.0;
    for (int i = 0; i < cs.n; i++) {
      pm += get_prim(f.at(i)).p;
    }
    *p_drift = pm / cs.n - cs.p0;
  }
  return peak_out / peak_in;
}

// The inflow-side analogue of the outflow measurement above: a LEFT-running
// pulse into a characteristic inflow holding a quiescent target, reflection
// measured in the right half after the bounce.  A soft inflow (small relax_u)
// lets the pulse push the boundary and swallows most of it; a stiff one is a
// Dirichlet condition in disguise and reflects like a wall.  This curve is
// how relax_u should be chosen, and it did not exist before: C2 checks only
// the relaxation DIRECTIONS at an inflow.
Real
inflow_reflection(Real relax_u, int n)
{
  Case cs;
  cs.n = n;
  Field f(cs.n), g(cs.n), h(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;

  const Real amp = 1.0e-3;
  const Real x0 = 0.65 * cs.L, w = 0.04 * cs.L;
  Real rho_ref = 0.0, e_ref = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho_ref, e_ref);
  const Prim qref = [&] {
    Prim q{};
    q.rho = rho_ref;
    q.u = 0.0;
    q.p = cs.p0;
    q.T = cs.T0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      q.Y[k] = Y[k];
    }
    return q;
  }();
  const Real c0 = sound_speed(qref);

  // A base INFLOW is required: with a quiescent target the incident pulse
  // pushes gas out through the inflow face, the kernel's reversal guard
  // (correctly) drops to a zero-gradient copy, and the "inflow" being
  // measured is not the inflow model at all -- every relax_u then returns
  // R = 0.  The pulse's velocity perturbation (~24 cm/s) rides on u0 = 2000
  // cm/s, so the face stays an inflow throughout.
  const Real u0 = 2.0e3;
  for (int i = -NG; i < cs.n + NG; i++) {
    const Real x = (i + 0.5) * dx;
    const Real dp = amp * cs.p0 * std::exp(-std::pow((x - x0) / w, 2));
    // LEFT-running isentropic pulse: u' = -dp/(rho c).
    const Real p = cs.p0 + dp;
    const Real rho = rho_ref + dp / (c0 * c0);
    const Real u = u0 - dp / (rho_ref * c0);
    Real T = 0.0;
    eos.RYP2T(rho, Y, p, T);
    set_state(f.at(i), rho, u, T, Y);
  }

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 0.0; // hi outflow: perfectly non-reflecting, out of the way
  prm.relax_u = relax_u;
  prm.relax_t = 0.2;
  pc_nscbc::Target in, out;
  in.type = pc_nscbc::Type::inflow;
  in.u[0] = u0;
  in.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    in.Y[k] = Y[k];
  }
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;

  const Real t_end = 1.6 * cs.L / c0;
  const Real t_measure = 0.9 * cs.L / c0; // pulse hits the inflow at ~0.65
  Real t = 0.0;
  Real peak_in = amp * cs.p0, peak_out = 0.0;

  while (t < t_end) {
    Real cmax = 0.0;
    for (int i = 0; i < cs.n; i++) {
      const Prim q = get_prim(f.at(i));
      cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
    }
    Real dt = cs.cfl * dx / cmax;
    dt = std::min(dt, t_end - t);
    if (dt <= 0.0) {
      break;
    }
    fill_bcs(f, in, out, prm, dx);
    stage(f, g, dx, dt);
    fill_bcs(g, in, out, prm, dx);
    stage(g, h, dx, dt);
    for (int i = 0; i < cs.n; i++) {
      for (int v = 0; v < NVAR; v++) {
        f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
      }
    }
    t += dt;

    if (t > t_measure) {
      Real pbar = 0.0;
      for (int i = 0; i < cs.n; i++) {
        pbar += get_prim(f.at(i)).p;
      }
      pbar /= cs.n;
      for (int i = cs.n / 2; i < cs.n; i++) {
        peak_out = std::max(peak_out, std::abs(get_prim(f.at(i)).p - pbar));
      }
    }
  }
  return peak_out / peak_in;
}

void
check_reflection(bool sweep)
{
  char buf[220];
  Real drift = 0.0;
  const Real R2 = reflection_coefficient(0.25, 400, 2, false, &drift);
  std::snprintf(
    buf, sizeof(buf), "R = %.4f %% (sigma=0.25, n=400, order=2)", 100.0 * R2);
  check(R2 < 0.01, "outflow reflection below 1%", buf);

  // Judge the extrapolation order at sigma = 0.  At sigma = 0.25 the residual
  // in the domain is dominated by the anchoring transient, not by reflection,
  // so an order comparison there measures the wrong thing.
  const Real R0o1 = reflection_coefficient(0.0, 400, 1, false, nullptr);
  const Real R0o2 = reflection_coefficient(0.0, 400, 2, false, nullptr);
  std::snprintf(
    buf, sizeof(buf), "sigma=0: order1 R = %.5f %%, order2 R = %.5f %%",
    100.0 * R0o1, 100.0 * R0o2);
  check(
    R0o2 <= R0o1 * 1.05, "2nd-order extrapolation is not worse than 1st", buf);

  std::snprintf(
    buf, sizeof(buf), "sigma=0 R = %.5f %%, sigma=0.25 R = %.5f %%",
    100.0 * R0o2, 100.0 * R2);
  check(R0o2 <= R2 * 1.05, "sigma=0 is the least reflecting", buf);

  // The material-slope continuation must not disturb the ACOUSTIC behaviour:
  // a right-running pulse has dR_- = 0 to linear order, so both the
  // reflection and the sigma anchoring must come out essentially unchanged.
  // This is the acoustic half of the extrap_material contract; C9(a) gates
  // the material half.
  Real drift_mat = 0.0;
  const Real R2m = reflection_coefficient(0.25, 400, 2, false, &drift_mat, true);
  std::snprintf(
    buf, sizeof(buf),
    "R = %.4f %% (entropy %.4f %%), drift %.2f (entropy %.2f)", 100.0 * R2m,
    100.0 * R2, drift_mat, drift);
  check(R2m < 0.01, "extrap_material keeps reflection below 1%", buf);
  check(
    std::abs(drift_mat) < 2.0 * std::abs(drift) + 0.5,
    "extrap_material keeps the sigma anchoring", buf);

  // The inflow curve.  R must rise monotonically with relax_u -- softer
  // swallows more, stiffer walls more -- and the two ends must actually
  // differ, or relax_u is a dial connected to nothing.
  const Real Ri[4] = {
    inflow_reflection(0.5, 400), inflow_reflection(2.0, 400),
    inflow_reflection(10.0, 400), inflow_reflection(50.0, 400)};
  std::snprintf(
    buf, sizeof(buf),
    "R = %.3f / %.3f / %.3f / %.3f at relax_u = 0.5 / 2 / 10 / 50",
    Ri[0], Ri[1], Ri[2], Ri[3]);
  check(
    (Ri[0] <= Ri[1] * 1.05) && (Ri[1] <= Ri[2] * 1.05) &&
      (Ri[2] <= Ri[3] * 1.05),
    "inflow reflection rises monotonically with relax_u", buf);
  check(
    Ri[3] > 2.0 * Ri[0],
    "relax_u spans soft to stiff (the ends differ)", buf);

  if (sweep) {
    std::printf("\n  sigma sweep (n=400, order=2)\n");
    std::printf(
      "  %8s  %12s  %16s  %14s\n", "sigma", "R [%]", "p drift [dyn/cm2]",
      "tau_relax [s]");
    for (Real s :
         {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0,
          16.0}) {
      Real d = 0.0;
      const Real R = reflection_coefficient(s, 400, 2, false, &d);
      const Real tau = (s > 0.0) ? Case().L / (s * 34783.7) : 0.0;
      std::printf("  %8.3f  %12.5f  %16.6e  %14.4e\n", s, 100.0 * R, d, tau);
    }
    Real d = 0.0;
    const Real Rp = reflection_coefficient(0.0, 400, 2, true, &d);
    std::printf(
      "  %8s  %12.5f  %16.6e  %14s   (pin_farfield)\n", "--", 100.0 * Rp, d,
      "0 (value pin)");
  }
}

// C5: the relaxation rate must be a RATE -- grid-independent, and equal to
//     K = sigma (1-M^2) c / L.  This is the check that distinguishes the
//     Poinsot-Lele parameterisation adopted here, as against a value-blend,
//     whose effective rate is c/dx and therefore doubles when the mesh does.
void
check_relaxation_rate()
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real sigma = 0.5;

  auto measure = [&](int n, bool mat = false) {
    Case c2;
    c2.n = n;
    Field f(n), g(n), h(n);
    const Real dx = c2.L / n;
    const Real dp0 = 0.02 * c2.p0;
    Real rho = 0.0, e = 0.0, T = 0.0;
    eos.PYT2RE(c2.p0 + dp0, Y, c2.T0, rho, e);
    eos.RYP2T(rho, Y, c2.p0 + dp0, T);
    for (int i = -NG; i < n + NG; i++) {
      set_state(f.at(i), rho, 0.0, T, Y);
    }
    const Prim q0 = get_prim(f.at(0));
    const Real c0 = sound_speed(q0);

    pc_nscbc::Params prm;
    prm.L_ref = c2.L;
    prm.sigma = sigma;
    prm.extrap_material = mat;
    pc_nscbc::Target out;
    out.type = pc_nscbc::Type::outflow;
    out.p = c2.p0;
    // Both faces outflow so the box simply depressurises.
    pc_nscbc::Target out2 = out;

    const Real t_end = 3.0 * c2.L / c0;
    Real t = 0.0;
    while (t < t_end) {
      Real cmax = 0.0;
      for (int i = 0; i < n; i++) {
        const Prim q = get_prim(f.at(i));
        cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
      }
      Real dt = std::min(c2.cfl * dx / cmax, t_end - t);
      if (dt <= 0.0) {
        break;
      }
      fill_bcs(f, out2, out, prm, dx);
      stage(f, g, dx, dt);
      fill_bcs(g, out2, out, prm, dx);
      stage(g, h, dx, dt);
      for (int i = 0; i < n; i++) {
        for (int v = 0; v < NVAR; v++) {
          f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
        }
      }
      t += dt;
    }
    Real pm = 0.0;
    for (int i = 0; i < n; i++) {
      pm += get_prim(f.at(i)).p;
    }
    pm /= n;
    // exp decay over t_end
    const Real ratio = std::max((pm - c2.p0) / dp0, 1e-12);
    return std::make_pair(-std::log(ratio) / t_end, c0);
  };

  const auto r200 = measure(200);
  const auto r800 = measure(800);
  const Real K_expected = sigma * r200.second / Case().L;
  char buf[240];
  std::snprintf(
    buf, sizeof(buf), "K(n=200)=%.2f, K(n=800)=%.2f 1/s  (ratio %.3f)",
    r200.first, r800.first, r800.first / r200.first);
  check(
    std::abs(r800.first / r200.first - 1.0) < 0.15,
    "relaxation rate is grid-independent", buf);
  std::snprintf(
    buf, sizeof(buf), "measured %.2f vs sigma*c/L = %.2f 1/s (ratio %.3f)",
    r800.first, K_expected, r800.first / K_expected);
  check(
    r800.first / K_expected > 0.2 && r800.first / K_expected < 5.0,
    "relaxation rate is within an order of K=sigma*c/L", buf);

  // The material-slope continuation adds a slope, not a rate: the relaxation
  // must decay the same offset at the same K.
  const auto r200m = measure(200, true);
  std::snprintf(
    buf, sizeof(buf), "K = %.2f with extrap_material, %.2f without (ratio %.3f)",
    r200m.first, r200.first, r200m.first / r200.first);
  check(
    std::abs(r200m.first / r200.first - 1.0) < 0.1,
    "extrap_material leaves the relaxation rate unchanged", buf);
}

// C6: robustness.  Every fallback path must return a finite, physical state
//     and must be counted.
void
check_fallbacks()
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  pc_nscbc::Target out;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;

  Real rho = 0.0, e = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho, e);

  amrex::Long diag[pc_nscbc::Diag::count] = {0};
  Real sN[NVAR], sg[NVAR];

  // supersonic outflow
  set_state(sN, rho, 8.0e4, cs.T0, Y);
  pc_nscbc::apply(sN, sN, sN, 1, dx, 0, -1, 1, out, prm, sg, diag);
  bool finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(
    finite && diag[pc_nscbc::Diag::supersonic] == 1 && sg[URHO] == sN[URHO],
    "supersonic outflow -> exact zero-gradient copy", "counted, state finite");

  // reversed flow at an outflow face
  set_state(sN, rho, -3.0e3, cs.T0, Y);
  pc_nscbc::apply(sN, sN, sN, 1, dx, 0, -1, 1, out, prm, sg, diag);
  finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(
    finite && diag[pc_nscbc::Diag::reversed] == 1 && sg[URHO] > 0.0,
    "reversed flow at outflow -> counted, standard closure",
    "counted, state finite (continuity and relaxation gated by C13)");

  // supersonic inflow: every characteristic is incoming, the FULL state is
  // required.  A Target without a pressure is counted (target_incomplete)
  // and the interior pressure substituted -- visible, not silent; with the
  // pressure supplied the imposition is exact.
  {
    pc_nscbc::Target sin_t;
    sin_t.type = pc_nscbc::Type::inflow;
    sin_t.T = 400.0;
    sin_t.u[0] = -8.0e4; // into the domain through the hi face, supersonic
    for (int k = 0; k < NUM_SPECIES; k++) {
      sin_t.Y[k] = Y[k];
    }
    set_state(sN, rho, -8.0e4, cs.T0, Y);
    char buf[220];

    pc_nscbc::apply(sN, sN, sN, 3, dx, 0, -1, 1, sin_t, prm, sg, diag);
    const pc_nscbc::CellPrim qg0 = pc_nscbc::cell_primitives(sg);
    std::snprintf(
      buf, sizeof(buf),
      "counted %lld; ghost T %.1f (target 400), p %.4g (interior %.4g "
      "substituted)",
      static_cast<long long>(diag[pc_nscbc::Diag::target_incomplete]), qg0.T,
      qg0.p, cs.p0);
    check(
      (diag[pc_nscbc::Diag::target_incomplete] == 1) &&
        (std::abs(qg0.T - 400.0) < 1.0e-6 * 400.0) &&
        (std::abs(qg0.p - cs.p0) < 1.0e-6 * cs.p0),
      "supersonic inflow without Target.p -> counted substitution", buf);

    sin_t.p = 1.2 * cs.p0;
    pc_nscbc::apply(sN, sN, sN, 3, dx, 0, -1, 1, sin_t, prm, sg, diag);
    const pc_nscbc::CellPrim qg1 = pc_nscbc::cell_primitives(sg);
    std::snprintf(
      buf, sizeof(buf), "ghost p %.6g vs target %.6g, T %.2f, counter still "
      "%lld",
      qg1.p, 1.2 * cs.p0, qg1.T,
      static_cast<long long>(diag[pc_nscbc::Diag::target_incomplete]));
    check(
      (diag[pc_nscbc::Diag::target_incomplete] == 1) &&
        (std::abs(qg1.p - 1.2 * cs.p0) < 1.0e-6 * cs.p0) &&
        (std::abs(qg1.T - 400.0) < 1.0e-6 * 400.0),
      "supersonic inflow with Target.p -> exact full-state imposition", buf);
  }

  // EB body state in the stencil
  Real body[NVAR];
  for (int v = 0; v < NVAR; v++) {
    body[v] = -1.0;
  }
  set_state(sN, rho, 3.0e3, cs.T0, Y);
  const amrex::Long before = diag[pc_nscbc::Diag::body_state];
  pc_nscbc::apply(sN, body, body, 3, dx, 0, -1, 1, out, prm, sg, diag);
  finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(
    finite && diag[pc_nscbc::Diag::body_state] > before && sg[URHO] > 0.0,
    "covered cells in stencil -> order degraded, no FPE",
    "counted, state finite");

  // fully covered boundary cell
  pc_nscbc::apply(body, body, body, 3, dx, 0, -1, 1, out, prm, sg, diag);
  finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(finite, "covered boundary cell -> finite fallback", "state finite");
}

// C13: the outflow closure must be CONTINUOUS through u_out = 0, and a
//      transient reversal must still be RELAXED -- the ghost normal velocity
//      has to respond outward when the interior is over-pressured, because
//      that response is the only feedback resisting counter-gradient inflow.
//
//      This is the NSCBC-Chamber defect gate
//      (Docs/NSCBC-reversal-branch-defect.md): the dedicated reversal branch
//      produced a ghost state that dropped dR+, S_p and T_in relative to the
//      forward branch -- an O(1) discontinuity at u_out = 0 exactly where a
//      flame finishes its transit -- and froze the ghost velocity, so a
//      breathing vent saw no restoring push.  The production A/B showed the
//      flow dithering across the branch boundary with growing amplitude
//      (-185 -> -962 cm/s in 0.25 ms) and a spurious 0.3 atm chamber spike.
//
//      (a) sweeps the interior normal velocity through zero over a flame-like
//      graded stencil (live dR+, and in reacting builds live S_p via
//      beta_s = 0) and asserts the crossing jump in ghost p and ghost u_out
//      is no outlier against the sweep's own smooth variation.
//      (b) holds a uniform over-pressured reversed state and asserts the
//      ghost pressure moves toward the target AND the ghost normal velocity
//      moves outward from the interior's -- a frozen ghost velocity fails.
void
check_reversal_continuity()
{
  Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 0.25;
  prm.beta_s = 0.0; // the flame-crossing recipe: the reaction source is live
  pc_nscbc::Target out;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;

  // A flame-like face: hot boundary cell, mild inward T and u gradients, 8%
  // overpressure so the relaxation has a defined restoring direction.
  Real Y[NUM_SPECIES] = {0.0};
#if NUM_REACTIONS == 0
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
#else
  Y[H2_ID] = 0.0285;
  Y[O2_ID] = 0.2265;
  Y[N2_ID] = 0.7450;
#endif
  const Real T_face = 1400.0;
  const Real p_in = 1.08 * cs.p0;
  const Real du = 500.0; // cm/s of normal-velocity gradient per cell

  // hi face: sgn = -1, so u_out = +u and positive u is outflow.
  auto fill = [&](Real u_face, bool graded, Real* sg) {
    Real sN[NVAR], sNm1[NVAR], sNm2[NVAR];
    const Real gT = graded ? 150.0 : 0.0;
    const Real gu = graded ? du : 0.0;
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(p_in, Y, T_face, rho, e);
    set_state(sN, rho, u_face, T_face, Y);
    eos.PYT2RE(p_in, Y, T_face - gT, rho, e);
    set_state(sNm1, rho, u_face - gu, T_face - gT, Y);
    eos.PYT2RE(p_in, Y, T_face - 2.0 * gT, rho, e);
    set_state(sNm2, rho, u_face - 2.0 * gu, T_face - 2.0 * gT, Y);
    pc_nscbc::apply(sN, sNm1, sNm2, 3, dx, 0, -1, 1, out, prm, sg, nullptr);
  };

  // (a) the sweep.  10 cm/s steps across [-1200, 1200] cm/s -- fine enough to
  // resolve the material-upwinding band below u_out = 0 -- and the crossing
  // pair's jump in ghost p, u AND T must sit inside the sweep's own
  // variation, not an order above it.
  const int ns = 241;
  const Real u_lo = -1200.0, u_hi = 1200.0;
  std::vector<Real> pg(ns), ug(ns), tg(ns), uu(ns);
  for (int i = 0; i < ns; i++) {
    uu[i] = u_lo + (u_hi - u_lo) * static_cast<Real>(i) / (ns - 1);
    Real sg[NVAR];
    fill(uu[i], true, sg);
    const pc_nscbc::CellPrim qg = pc_nscbc::cell_primitives(sg);
    pg[i] = qg.p;
    ug[i] = qg.u[0]; // n_sgn = +1 on the hi face
    tg[i] = qg.T;
  }
  Real jp_cross = 0.0, ju_cross = 0.0, jt_cross = 0.0;
  Real jp_other = 0.0, ju_other = 0.0, jt_other = 0.0;
  for (int i = 0; i + 1 < ns; i++) {
    const Real jp = std::abs(pg[i + 1] - pg[i]);
    const Real ju = std::abs(ug[i + 1] - ug[i]);
    const Real jt = std::abs(tg[i + 1] - tg[i]);
    const bool crossing = (uu[i] < 0.0) && (uu[i + 1] >= 0.0);
    if (crossing) {
      jp_cross = jp;
      ju_cross = ju;
      jt_cross = jt;
    } else {
      jp_other = std::max(jp_other, jp);
      ju_other = std::max(ju_other, ju);
      jt_other = std::max(jt_other, jt);
    }
  }
  char buf[220];
  std::snprintf(
    buf, sizeof(buf),
    "crossing jump: p %.3e (<= %.3e elsewhere), u %.3e (<= %.3e), "
    "T %.3e (<= %.3e)",
    jp_cross, jp_other, ju_cross, ju_other, jt_cross, jt_other);
  check(
    (jp_cross <= 3.0 * jp_other + 1.0e-9 * cs.p0) &&
      (ju_cross <= 3.0 * ju_other + 1.0e-9 * 1.0e5) &&
      (jt_cross <= 3.0 * jt_other + 1.0e-6 * cs.T0),
    "ghost state is continuous through u_out = 0", buf);

  // (b) reversal is still relaxed.  Uniform over-pressured stencil, firmly
  // reversed: no extrapolation content, so the only ghost response is the
  // incoming-wave relaxation, and it must move p toward the target and u_out
  // outward.  A ghost velocity frozen at the interior's value is exactly the
  // no-feedback closure that let the chamber runaway feed.
  {
    const Real u_rev = -600.0;
    Real sg[NVAR], sN[NVAR];
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(p_in, Y, T_face, rho, e);
    set_state(sN, rho, u_rev, T_face, Y);
    fill(u_rev, false, sg);
    const pc_nscbc::CellPrim qN = pc_nscbc::cell_primitives(sN);
    const pc_nscbc::CellPrim qg = pc_nscbc::cell_primitives(sg);
    std::snprintf(
      buf, sizeof(buf), "ghost du_out = %+.4e cm/s, ghost dp = %+.4e dyn/cm^2",
      qg.u[0] - qN.u[0], qg.p - qN.p);
    check(
      (qg.u[0] > qN.u[0] + 1.0e-10 * qN.c) && (qg.p < qN.p),
      "over-pressured reversal relaxes p AND pushes u_out outward", buf);
  }

  // (c) reversal must NOT extrapolate material content.  The stencil here is
  // the chamber's fatal configuration: temperature falling TOWARD the face
  // (cold gas at the boundary).  Under firm backflow an extrapolated ghost T
  // continues that ramp downward, the inflow advects the colder ghost gas
  // back in, and the loop refrigerates the boundary cell to cryogenic
  // temperatures while 1/(rho c) amplifies the relaxation -- the production
  // vent run went 385 -> 241 -> 89 K in 42 us and NaN'd.  The ghost
  // temperature under firm reversal must be the interior cell's own.
  {
    const Real u_rev = -600.0;
    Real sN[NVAR], sNm1[NVAR], sNm2[NVAR], sg[NVAR];
    Real rho = 0.0, e = 0.0;
    const Real T_cold = 400.0; // face cell, coldest; T RISES inward
    eos.PYT2RE(p_in, Y, T_cold, rho, e);
    set_state(sN, rho, u_rev, T_cold, Y);
    eos.PYT2RE(p_in, Y, T_cold + 300.0, rho, e);
    set_state(sNm1, rho, u_rev, T_cold + 300.0, Y);
    eos.PYT2RE(p_in, Y, T_cold + 600.0, rho, e);
    set_state(sNm2, rho, u_rev, T_cold + 600.0, Y);
    pc_nscbc::apply(sN, sNm1, sNm2, 3, dx, 0, -1, 1, out, prm, sg, nullptr);
    const pc_nscbc::CellPrim qN = pc_nscbc::cell_primitives(sN);
    const pc_nscbc::CellPrim qg = pc_nscbc::cell_primitives(sg);
    std::snprintf(
      buf, sizeof(buf), "ghost T = %.2f K vs interior %.2f K", qg.T, qN.T);
    check(
      std::abs(qg.T - qN.T) < 1.0e-3 * qN.T,
      "firm reversal freezes material: ghost T is the interior's", buf);
  }
}

// C14: sustained recirculation.  A duct held in steady INFLOW through a face
//      configured as an outflow: the lo face relaxes hard toward a low
//      pressure (the sink), the hi face is an outflow at ambient whose
//      Target also carries the reservoir state (T, Y).  Gas enters through
//      the hi face for the whole run -- not a breath but a regime, the
//      configuration queue item 2 was written for.  The material question is
//      what the entering gas carries: the frozen-material reversal closure
//      recycles the interior state, so a hot domain drawing from a cold
//      reservoir stays hot on its own exhaust forever (the measured negative
//      control); with backflow_material the lambda_0 ghost content ramps to
//      the reservoir state under firm reversal and the domain flushes cold
//      on the advective clock.  A static continuity probe confirms the ramp
//      is inert at breathing amplitudes (|M| < 1e-3).
Real
c14_run_flush(bool backflow, Real T_hot, Real T_res, Real dp_sink, Real t_end)
{
  Case cs;
  cs.n = 100;
  cs.L = 2.5;
  Field f(cs.n), g(cs.n), h(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  Real rho_h = 0.0, e_h = 0.0;
  eos.PYT2RE(cs.p0, Y, T_hot, rho_h, e_h);
  for (int i = -NG; i < cs.n + NG; i++) {
    set_state(f.at(i), rho_h, 0.0, T_hot, Y);
  }

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 2.0;
  prm.backflow_material = backflow;
  pc_nscbc::Target sink, res;
  sink.type = pc_nscbc::Type::outflow;
  sink.p = cs.p0 - dp_sink;
  res.type = pc_nscbc::Type::outflow;
  res.p = cs.p0;
  res.T = T_res; // the reservoir state the backflow may draw on
  for (int k = 0; k < NUM_SPECIES; k++) {
    res.Y[k] = Y[k];
  }

  Real t = 0.0;
  while (t < t_end) {
    Real cmax = 0.0;
    for (int i = 0; i < cs.n; i++) {
      const Prim q = get_prim(f.at(i));
      cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
    }
    Real dt = cs.cfl * dx / cmax;
    dt = std::min(dt, t_end - t);
    if (dt <= 0.0) {
      break;
    }
    fill_bcs(f, sink, res, prm, dx);
    stage(f, g, dx, dt);
    fill_bcs(g, sink, res, prm, dx);
    stage(g, h, dx, dt);
    for (int i = 0; i < cs.n; i++) {
      for (int v = 0; v < NVAR; v++) {
        f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
      }
    }
    t += dt;
  }
  // Mean temperature of the half nearest the reversed (hi) face.
  Real Tm = 0.0;
  for (int i = cs.n / 2; i < cs.n; i++) {
    Tm += get_prim(f.at(i)).T;
  }
  return Tm / (cs.n - cs.n / 2);
}

void
check_sustained_recirculation()
{
  const Real T_hot = 600.0, T_res = 300.0, dp_sink = 5.0e4;
  const Real t_end = 4.0e-3;
  char buf[220];

  const Real T_frozen = c14_run_flush(false, T_hot, T_res, dp_sink, t_end);
  const Real T_flush = c14_run_flush(true, T_hot, T_res, dp_sink, t_end);
  std::snprintf(
    buf, sizeof(buf),
    "hi-half mean T: frozen closure %.1f K, backflow_material %.1f K "
    "(interior 600, reservoir 300)",
    T_frozen, T_flush);
  check(
    T_frozen > 0.8 * T_hot,
    "frozen material recycles: the domain feeds on its own exhaust", buf);
  check(
    std::abs(T_flush - T_res) < 0.15 * T_res,
    "backflow_material flushes to the reservoir state", buf);

  // Continuity at breathing amplitude: at |M| ~ 1e-4 the ramp must be inert
  // -- the ghost with the flag on is the ghost with it off.
  {
    Case cs;
    Real Y[NUM_SPECIES];
    for (int k = 0; k < NUM_SPECIES; k++) {
      Y[k] = air_Y(k);
    }
    auto eos = pele::physics::PhysicsType::eos();
    const Real dx = cs.L / cs.n;
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(1.05 * cs.p0, Y, 500.0, rho, e);
    Real sN[NVAR], sg_off[NVAR], sg_on[NVAR];
    set_state(sN, rho, -4.0, 500.0, Y); // |M| ~ 1e-4: a breath
    pc_nscbc::Params prm;
    prm.L_ref = cs.L;
    pc_nscbc::Target res;
    res.type = pc_nscbc::Type::outflow;
    res.p = cs.p0;
    res.T = 300.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      res.Y[k] = Y[k];
    }
    prm.backflow_material = false;
    pc_nscbc::apply(sN, sN, sN, 3, dx, 0, -1, 1, res, prm, sg_off);
    prm.backflow_material = true;
    pc_nscbc::apply(sN, sN, sN, 3, dx, 0, -1, 1, res, prm, sg_on);
    Real dmax = 0.0;
    for (int v = 0; v < NVAR; v++) {
      dmax = std::max(
        dmax, std::abs(sg_on[v] - sg_off[v]) /
                std::max(std::abs(sg_off[v]), 1.0e-30));
    }
    std::snprintf(buf, sizeof(buf), "max relative ghost difference %.2e", dmax);
    check(
      dmax < 1.0e-12,
      "the ramp is inert at breathing amplitudes (|M| ~ 1e-4)", buf);
  }
}

// C8: does the ghost carry a sensible TEMPERATURE gradient?
//
//     This is not a hyperbolic question.  PeleC's diffusion operator forms the
//     conductive and species fluxes at a physical boundary face from these same
//     ghost cells, so whatever normal temperature gradient the ghost happens to
//     carry IS the heat flux leaving the domain.  Nothing in the characteristic
//     algebra is chosen with that in mind: at an outflow the ghost density
//     comes from the extrapolated entropy invariant and the ghost pressure from
//     the acoustic pair, and T is then whatever the EOS returns.
//
//     The test puts a flame-like temperature ramp on the boundary stencil at
//     uniform pressure -- 300 K rising at 2e4 K/cm, which is what a hydrocarbon
//     flame does -- and asks how far the ghost temperature is from the linear
//     continuation of the interior profile, as a fraction of one cell's dT.
//     A boundary with a controlled diffusive flux should sit near zero.
void
check_diffusive_gradient()
{
  const Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }

  const Real dx = 1.0e-2;  // 0.01 cm, ~10 cells through a flame
  const Real dTdx = 2.0e4; // K/cm
  const Real p0 = cs.p0;
  const Real u0 = 3.0e3; // outflow, subsonic

  // Interior stencil N, N-1, N-2 at uniform pressure with a linear T ramp,
  // T increasing toward the boundary as it does on the burnt side of a flame.
  Real sN[NVAR], sM1[NVAR], sM2[NVAR], sg[NVAR];
  const Real T_N = 1400.0;
  const Real T_M1 = T_N - dTdx * dx;
  const Real T_M2 = T_N - 2.0 * dTdx * dx;
  auto rho_of = [&](const Real T) {
    Real r = 0.0, e = 0.0;
    eos.PYT2RE(p0, Y, T, r, e);
    return r;
  };
  set_state(sN, rho_of(T_N), u0, T_N, Y);
  set_state(sM1, rho_of(T_M1), u0, T_M1, Y);
  set_state(sM2, rho_of(T_M2), u0, T_M2, Y);

  pc_nscbc::Params prm;
  prm.sigma = 0.25;
  prm.L_ref = 10.0;
  prm.order = 2;
  pc_nscbc::Target tgt;
  tgt.type = pc_nscbc::Type::outflow;
  tgt.p = p0;

  char buf[256];
  // Layer 1 is the one that sets the face flux.
  pc_nscbc::apply(sN, sM1, sM2, 3, dx, 0, -1, 1, tgt, prm, sg);
  const Real T_g = get_prim(sg).T;
  const Real T_lin = T_N + dTdx * dx; // linear continuation
  const Real dT_cell = dTdx * dx;     // one cell's worth
  const Real err = (T_g - T_lin) / dT_cell;

  std::snprintf(
    buf, sizeof(buf),
    "ghost T = %.2f K, linear continuation %.2f K, error %.3f of a cell dT",
    T_g, T_lin, err);
  // A tolerance of one full cell dT is deliberately loose: this check exists to
  // MEASURE the discrepancy and put a number in the log, not to gate on a
  // tight bound the current formulation was never designed to meet.
  check(
    std::abs(err) < 1.0, "ghost T within one cell dT of the interior ramp",
    buf);

  // The implied conductive flux error, in the only units that matter.  Reported
  // rather than gated: lambda is problem-dependent, so the fraction is the
  // transferable number.
  std::snprintf(
    buf, sizeof(buf), "implied face dT/dx is %.1f%% of the interior value",
    100.0 * (T_g - get_prim(sN).T) / dT_cell);
  check(true, "  (reported) face temperature gradient", buf);

  // Species: Y is minmod-extrapolated, so it should be a clean continuation.
  // Uniform composition here, so the ghost must reproduce it exactly.
  const Prim qg = get_prim(sg);
  Real dYmax = 0.0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    dYmax = std::max(dYmax, std::abs(qg.Y[k] - Y[k]));
  }
  std::snprintf(buf, sizeof(buf), "max |dY| = %.3e", dYmax);
  check(dYmax < 1.0e-12, "uniform composition passes through unchanged", buf);

  // Now the temperature closure, which exists precisely to fix the above.
  prm.extrap_temperature = true;
  pc_nscbc::apply(sN, sM1, sM2, 3, dx, 0, -1, 1, tgt, prm, sg);
  const Real T_gT = get_prim(sg).T;
  const Real errT = (T_gT - T_lin) / dT_cell;
  std::snprintf(
    buf, sizeof(buf), "ghost T = %.2f K vs %.2f K linear, error %.2e of a cell",
    T_gT, T_lin, errT);
  check(
    std::abs(errT) < 1.0e-9, "extrap_temperature reproduces the ramp exactly",
    buf);

  std::snprintf(
    buf, sizeof(buf),
    "face dT/dx: entropy closure %.1f%%, temperature closure "
    "%.1f%% of the interior value",
    100.0 * (T_g - get_prim(sN).T) / dT_cell,
    100.0 * (T_gT - get_prim(sN).T) / dT_cell);
  check(true, "  (reported) the two closures side by side", buf);

  // The point of the closure is the DIFFUSIVE flux, so it must not have cost
  // anything on the hyperbolic side: a uniform state must still come back
  // exactly, at every layer.
  Real su[NVAR], sgu[NVAR];
  set_state(su, rho_of(cs.T0), u0, cs.T0, Y);
  Real worst = 0.0;
  for (int layer = 1; layer <= 4; layer++) {
    pc_nscbc::apply(su, su, su, 3, dx, 0, -1, layer, tgt, prm, sgu);
    for (int v = 0; v < NVAR; v++) {
      const Real ref = std::abs(su[v]) > 1.0e-30 ? std::abs(su[v]) : 1.0;
      worst = std::max(worst, std::abs(sgu[v] - su[v]) / ref);
    }
  }
  std::snprintf(buf, sizeof(buf), "max rel ghost error = %.3e", worst);
  check(
    worst < 1.0e-12, "extrap_temperature still reproduces a uniform state",
    buf);
}

// C9: the ghost-pressure bias -- the mechanism blamed in
//     Exec/RegTests/NSCBC-FlameOutflow for the mean-pressure error at a
//     front-crossing outflow, reproduced here where nothing else can be
//     responsible.
//
//     The claim is that extrapolating the OUTGOING invariant R_+ = u_out +
//     p/(rho c) across a region with a normal velocity gradient manufactures a
//     ghost pressure that has nothing to do with any acoustic wave, because
//     dR_+/dn there is dominated by dilatation.  Two parts:
//
//     (a) STATIC.  With p_N already at the target the relaxation contributes
//         nothing, so the algebra predicts exactly
//
//             p_ghost - p_N = 1/2 rho c * layer * du_out
//
//         and predicts zero at order = 1, where the slope is discarded.  No
//         fitted quantity: if the mechanism is what it is claimed to be, this
//         is an identity.
//
//     (b) DYNAMIC, and NOT an isolation of (a) -- see the note at the order
//         control below.  A prescribed heat band straddling the outflow,
//         against a matched-heat band in the interior as a control.  The
//         difference is what the boundary added, it falls with sigma, and it
//         has the same shape as the sigma sweep in
//         Exec/RegTests/NSCBC-FlameOutflow.  But the order = 1 control shows it
//         is dominated by the unmodelled energy source in the boundary cells
//         rather than by the extrapolation bias of part (a).  Reported, not
//         gated.
void
check_ghost_pressure_bias()
{
  const Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  char buf[256];

  // ---- (a) static -------------------------------------------------------
  {
    const Real dx = 2.5e-2;
    const Real u_N = 4.0e3; // 40 m/s at the boundary cell
    const Real du = 5.0e2;  // 5 m/s per cell of ramp -> du/dn = 2e4 1/s
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(cs.p0, Y, cs.T0, rho, e);
    const Prim qref = [&] {
      Prim q{};
      q.rho = rho;
      q.u = u_N;
      q.p = cs.p0;
      q.T = cs.T0;
      for (int k = 0; k < NUM_SPECIES; k++) {
        q.Y[k] = Y[k];
      }
      return q;
    }();
    const Real c0 = sound_speed(qref);

    Real sN[NVAR], sM1[NVAR], sM2[NVAR], sg[NVAR];
    set_state(sN, rho, u_N, cs.T0, Y);
    set_state(sM1, rho, u_N - du, cs.T0, Y);
    set_state(sM2, rho, u_N - 2.0 * du, cs.T0, Y);

    pc_nscbc::Params prm;
    prm.L_ref = cs.L;
    prm.sigma = 0.25;
    prm.order = 2;
    pc_nscbc::Target out;
    out.type = pc_nscbc::Type::outflow;
    out.p = cs.p0; // p_N is already on target

    for (int layer = 1; layer <= 2; layer++) {
      pc_nscbc::apply(sN, sM1, sM2, 3, dx, 0, -1, layer, out, prm, sg);
      const Real p_g = get_prim(sg).p;
      const Real predicted = 0.5 * rho * c0 * layer * du;
      const Real rel = (p_g - cs.p0 - predicted) / predicted;
      std::snprintf(
        buf, sizeof(buf),
        "layer %d: p_ghost - p_N = %.2f, predicted 1/2 rho c l du = %.2f "
        "(rel %.2e)",
        layer, p_g - cs.p0, predicted, rel);
      check(
        std::abs(rel) < 1.0e-6, "ghost pressure bias matches the algebra", buf);
    }

    prm.order = 1;
    pc_nscbc::apply(sN, sM1, sM2, 3, dx, 0, -1, 1, out, prm, sg);
    const Real p_g1 = get_prim(sg).p;
    std::snprintf(
      buf, sizeof(buf), "order 1: p_ghost - p_N = %.3e (bias term discarded)",
      p_g1 - cs.p0);
    check(
      std::abs(p_g1 - cs.p0) < 1.0e-6 * cs.p0,
      "the bias is carried entirely by the extrapolation", buf);

    // The transit guard must stay quiet here: a velocity ramp at uniform
    // density has dS = 0, and quiet-on-acoustics is the guard's whole value.
    {
      amrex::Long diag[pc_nscbc::Diag::count] = {0};
      prm.order = 2;
      pc_nscbc::apply(sN, sM1, sM2, 3, dx, 0, -1, 1, out, prm, sg, diag);
      std::snprintf(
        buf, sizeof(buf), "structure count = %lld on a dS = 0 ramp",
        static_cast<long long>(diag[pc_nscbc::Diag::structure]));
      check(
        diag[pc_nscbc::Diag::structure] == 0,
        "transit guard is quiet without entropy structure", buf);
    }

    // The fix, on the structure it exists for.  The uniform-density ramp
    // above is a synthetic state -- steady continuity does not admit it -- so
    // the entropy-family bound in extrap_material correctly sees nothing
    // there.  Rebuild the ramp as the mass-conserving structure of C10
    // (rho u uniform, pressure carrying the momentum flux): the measured
    // slope of R_- and the entropy bound then agree, the ghost continues the
    // interior's u and p slopes, and the 1/2 rho c du bias is gone while
    // order = 2 keeps the full structure that order = 1 throws away.
    {
      const Real mdot = rho * u_N;
      auto ramp_state = [&](Real* s, const Real u) {
        const Real rr = mdot / u;
        const Real pp = cs.p0 + mdot * (u_N - u); // p_N lands on the target
        Real TT = 0.0;
        eos.RYP2T(rr, Y, pp, TT);
        set_state(s, rr, u, TT, Y);
      };
      ramp_state(sN, u_N);
      ramp_state(sM1, u_N - du);
      ramp_state(sM2, u_N - 2.0 * du);
      const Real dp_cell = -mdot * du; // exact per-cell momentum-flux slope
      const Real bias0 = 0.5 * rho * c0 * du; // what the entropy closure adds

      prm.order = 2;
      prm.extrap_material = true;
      amrex::Long diag[pc_nscbc::Diag::count] = {0};
      for (int layer = 1; layer <= 2; layer++) {
        pc_nscbc::apply(sN, sM1, sM2, 3, dx, 0, -1, layer, out, prm, sg, diag);
        const Prim qg = get_prim(sg);
        const Real p_resid = qg.p - (cs.p0 + layer * dp_cell);
        std::snprintf(
          buf, sizeof(buf),
          "layer %d: p_ghost - p_expected = %.2f (entropy closure bias %.1f), "
          "u_ghost = %.1f (slope continued: %.1f)",
          layer, p_resid, layer * bias0, qg.u, u_N + layer * du);
        check(
          std::abs(p_resid) < 0.05 * layer * bias0,
          "extrap_material removes the ghost-pressure bias", buf);
        check(
          std::abs(qg.u - (u_N + layer * du)) < 0.05 * layer * du,
          "extrap_material keeps the full du/dn in the ghost", buf);
      }
      prm.extrap_material = false;

      // ... and the transit guard must FIRE here: this ramp's per-cell
      // density change is far past the 5% threshold, and it is exactly the
      // structure whose crossing the sigma = 0.25 default does not survive.
      std::snprintf(
        buf, sizeof(buf), "structure count = %lld on the mass-conserving ramp",
        static_cast<long long>(diag[pc_nscbc::Diag::structure]));
      check(
        diag[pc_nscbc::Diag::structure] > 0,
        "transit guard fires on material structure", buf);
    }
  }

  // ---- (b) dynamic ------------------------------------------------------
  // A prescribed heat band sustains a velocity gradient.  Two placements at
  // matched total heat release:
  //
  //   INTERIOR  band at L/2  -- the flow has finished accelerating long before
  //                            the outflow, so du/dn at the boundary is ~0
  //   BOUNDARY  band at L    -- the gradient is IN the boundary cells, exactly
  //                            as the flame straddles the outflow in
  //                            Exec/RegTests/NSCBC-FlameOutflow
  //
  // Heating a duct at fixed inflow raises the mean pressure whatever the
  // boundary does -- that is real physics, set by mass and energy balance, and
  // it is the same in both placements once the total heat is matched.  The
  // DIFFERENCE is what the boundary added, and it is the only quantity here
  // that the bias can be responsible for.
  auto run_heated = [&](
                      const Real sigma, const Real Qmax, const int n,
                      const Real xq_frac, const int order, Real& dudn_out) {
    const Real dx = cs.L / n;
    Field f(n), g(n), h(n);
    Real rho0 = 0.0, e0 = 0.0;
    eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
    const Real u_in = 3.0e3;
    for (int i = -NG; i < n + NG; i++) {
      set_state(f.at(i), rho0, u_in, cs.T0, Y);
    }

    pc_nscbc::Params prm;
    prm.L_ref = cs.L;
    prm.sigma = sigma;
    prm.order = order;
    pc_nscbc::Target off, out;
    out.type = pc_nscbc::Type::outflow;
    out.p = cs.p0;

    const Real xq = xq_frac * cs.L, wq = 0.4;
    std::vector<Real> qd(n, 0.0);
    Real qtot = 0.0;
    for (int i = 0; i < n; i++) {
      qd[i] = std::exp(-std::pow(((i + 0.5) * dx - xq) / wq, 2));
      qtot += qd[i] * dx;
    }
    // Normalise so both placements deliver the same integrated heat, despite
    // the boundary band being half outside the domain.
    for (int i = 0; i < n; i++) {
      qd[i] *= Qmax * wq * std::sqrt(M_PI) / qtot;
    }

    const Real t_end = 4.0 * cs.L / u_in;
    Real t = 0.0;
    while (t < t_end) {
      Real cmax = 0.0;
      for (int i = 0; i < n; i++) {
        const Prim q = get_prim(f.at(i));
        cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
      }
      Real dt = std::min(cs.cfl * dx / cmax, t_end - t);
      if (dt <= 0.0) {
        break;
      }

      for (int layer = 1; layer <= NG; layer++) {
        set_state(f.at(-layer), rho0, u_in, cs.T0, Y);
      }
      fill_bcs(f, off, out, prm, dx);
      stage(f, g, dx, dt);
      for (int layer = 1; layer <= NG; layer++) {
        set_state(g.at(-layer), rho0, u_in, cs.T0, Y);
      }
      fill_bcs(g, off, out, prm, dx);
      stage(g, h, dx, dt);
      for (int i = 0; i < n; i++) {
        for (int v = 0; v < NVAR; v++) {
          f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
        }
        f.at(i)[UEDEN] += dt * qd[i];
        f.at(i)[UEINT] += dt * qd[i];
        Real* sp = f.at(i);
        const Real rr = sp[URHO];
        const Real uu = sp[UMX] / rr;
        Real Yl[NUM_SPECIES], ys = 0.0;
        for (int k = 0; k < NUM_SPECIES; k++) {
          Yl[k] = sp[UFS + k] / rr;
          ys += Yl[k];
        }
        for (int k = 0; k < NUM_SPECIES; k++) {
          Yl[k] /= ys;
        }
        Real Tl = sp[UTEMP];
        eos.REY2T(rr, sp[UEDEN] / rr - 0.5 * uu * uu, Yl, Tl);
        sp[UTEMP] = Tl;
      }
      t += dt;
    }

    dudn_out = (get_prim(f.at(n - 1)).u - get_prim(f.at(n - 2)).u) / dx;
    Real psum = 0.0;
    for (int i = 0; i < n; i++) {
      psum += get_prim(f.at(i)).p;
    }
    return psum / n - cs.p0;
  };

  const Real Q = 2.0e9; // erg/cm^3/s
  const int NN = 200;
  std::printf(
    "\n    %7s %9s %11s %11s %11s %11s\n", "sigma", "order", "du/dn int",
    "du/dn bnd", "<p> int", "<p> bnd");
  Real bias[4];
  const Real sigs[4] = {0.25, 1.0, 4.0, 16.0};
  for (int k = 0; k < 4; k++) {
    Real di = 0.0, db = 0.0;
    const Real oi = run_heated(sigs[k], Q, NN, 0.5, 2, di);
    const Real ob = run_heated(sigs[k], Q, NN, 1.0, 2, db);
    bias[k] = ob - oi;
    std::printf(
      "    %7.2f %9d %11.0f %11.0f %11.1f %11.1f\n", sigs[k], 2, di, db, oi,
      ob);
  }
  std::printf("      bias (bnd - int): ");
  for (int k = 0; k < 4; k++) {
    std::printf("%10.1f", bias[k]);
  }
  std::printf("\n");

  std::snprintf(
    buf, sizeof(buf), "bias %.1f -> %.1f as sigma 0.25 -> 16 (x%.2f)", bias[0],
    bias[3], bias[0] / (std::abs(bias[3]) > 1e-30 ? bias[3] : 1.0));
  check(
    std::abs(bias[0]) > 2.0 * std::abs(bias[3]),
    "a source AT the boundary adds an offset that sigma suppresses", buf);

  // The order control, and the reason the dynamic part of this check reports
  // rather than gates.
  //
  // Statically (part a) order = 1 gives EXACTLY zero ghost-pressure bias.  If
  // the dynamic offset above were the extrapolation bias, dropping to first
  // order would remove most of it.  It does not -- the two agree to ~2%.  So
  // what the heat band actually measures is an unmodelled ENERGY SOURCE in the
  // boundary cells, which is the beta_s situation, not the beta one, and at
  // ~12% of ambient it swamps the extrapolation term entirely.
  //
  // A heat band cannot do better: there is no way to hold a velocity gradient
  // at the boundary with a local source without also putting that source in the
  // boundary cells.  Isolating the extrapolation dynamically needs a
  // source-free sustainer -- an established velocity/density ramp advected out
  // through the face -- which is not built yet.  Until it is, part (a) is the
  // isolation and this part is context.
  Real d1i = 0.0, d1b = 0.0;
  const Real o1i = run_heated(1.0, Q, NN, 0.5, 1, d1i);
  const Real o1b = run_heated(1.0, Q, NN, 1.0, 1, d1b);
  std::snprintf(
    buf, sizeof(buf),
    "order 2 %.1f vs order 1 %.1f -- agree, so this is NOT the extrapolation",
    bias[1], o1b - o1i);
  check(true, "  (reported) order control on the dynamic offset", buf);
}

// ---------------------------------------------------------------------------
//  fit_profile_ghosts -- the C11x profile-fit ghost closure, shared by the
//  C10 (release) and C11 (sustained front) experiments.
//
//  The closure knows a tanh profile FAMILY -- end states (u0, ratio_a*u0),
//  rho = mdot/u, p = p0 + mdot*(u0 - u) -- but not its position or thickness;
//  both are fitted per call by inverting T at the last two interior cells
//  through the family (a closed-form value-and-slope match, stateless).  It
//  overwrites only the MATERIAL content of the kernel-filled ghosts: T from
//  the fitted profile, rho from the EOS at the kernel's ghost pressure, and
//  (with_u) the profile's u -- the dilatation structure.  p always stays the
//  kernel's.
//
//  q_src >= 0 applies the SOURCE-CONSISTENCY bound: a steady front obeys
//  du/dn = dp/dt|_src/(rho c^2) with dp/dt|_src = (gamma-1) q, so the
//  continuation blends toward the plain kernel by
//  w = min(1, du/dn_sustainable / du/dn_measured) -- inert on a front the
//  source sustains, a full release on one it cannot.
// ---------------------------------------------------------------------------
void
fit_profile_ghosts(
  Field& w,
  const Real dx,
  const Real u0,
  const Real ratio_a,
  const Real mdot,
  const Real p0,
  const Real Y[NUM_SPECIES],
  const bool with_u,
  const Real q_src)
{
  auto eos = pele::physics::PhysicsType::eos();
  auto T_of_u = [&](const Real u) {
    const Real rho = mdot / u;
    const Real p = p0 + mdot * (u0 - u);
    Real T = 0.0;
    eos.RYP2T(rho, Y, p, T);
    return T;
  };
  const Real ua = u0 * (1.0 + 1.0e-9);
  const Real ub = u0 * ratio_a * (1.0 - 1.0e-9);
  auto u_of_T = [&](const Real T) {
    if (T <= T_of_u(ua)) {
      return ua;
    }
    if (T >= T_of_u(ub)) {
      return ub;
    }
    Real a = ua, b = ub;
    for (int it = 0; it < 60; it++) {
      const Real m = 0.5 * (a + b);
      (T_of_u(m) < T ? a : b) = m;
    }
    return 0.5 * (a + b);
  };
  const int N = w.n - 1;
  const Prim qN = get_prim(w.at(N));
  const Prim qM = get_prim(w.at(N - 1));
  auto arg = [&](const Real T) {
    Real g = (u_of_T(T) / u0 - 1.0) / (ratio_a - 1.0);
    g = std::min(std::max(g, 1.0e-9), 1.0 - 1.0e-9);
    return std::atanh(2.0 * g - 1.0);
  };
  const Real aN = arg(qN.T);
  const Real aM = arg(qM.T);
  if (!(aN > aM) || !std::isfinite(aN) || !std::isfinite(aM)) {
    return; // no usable structure in T; the kernel's ghosts stand
  }
  Real wgt = 1.0;
  if (q_src >= 0.0) {
    const Real c_N = sound_speed(qN);
    const Real gam = qN.rho * c_N * c_N / qN.p;
    const Real g_sus = (gam - 1.0) * q_src / (qN.rho * c_N * c_N);
    const Real g_meas = (qN.u - qM.u) / dx;
    wgt =
      (g_meas > 1.0e-12 * u0 / dx) ? std::min(1.0, g_sus / g_meas) : 0.0;
    if (wgt <= 0.0) {
      return; // nothing sustainable; released to the plain kernel
    }
  }
  const Real lam = dx / (aN - aM);
  const Real sh = (N + 0.5) * dx - lam * aN;
  for (int layer = 1; layer <= NG; layer++) {
    const int i = N + layer;
    const Real x = (i + 0.5) * dx;
    const Real gg = 0.5 * (1.0 + std::tanh((x - sh) / lam));
    const Real uf = u0 * (1.0 + (ratio_a - 1.0) * gg);
    const Real rf = mdot / uf;
    const Real pf = p0 + mdot * (u0 - uf);
    Real Tf = 0.0;
    eos.RYP2T(rf, Y, pf, Tf);
    const Prim qg = get_prim(w.at(i)); // the kernel's fill: p always kept
    const Real Tb = wgt * Tf + (1.0 - wgt) * qg.T;
    const Real ub2 = with_u ? (wgt * uf + (1.0 - wgt) * qg.u) : qg.u;
    Real rho_g = 0.0, e_g = 0.0;
    eos.PYT2RE(qg.p, Y, Tb, rho_g, e_g);
    set_state(w.at(i), rho_g, ub2, Tb, Y);
  }
}

// C10: does the ghost-pressure bias of C9(a) actually DRIVE the solution?
//
//      C9(a) establishes the bias exactly, but statically.  C9(b) tried to make
//      it dynamic with a heat source and failed: a source strong enough to hold
//      a velocity gradient at the boundary also deposits energy there, and the
//      order control showed the source, not the extrapolation, setting the
//      offset.  C10 removes every source instead.
//
//      The sustainer is a steady-flame structure with the chemistry taken out:
//      uniform mass flux rho*u, velocity rising through a tanh ramp by an
//      expansion ratio, density falling to match, and the pressure carrying the
//      momentum flux so that rho*u^2 + p is uniform as well.  It straddles the
//      outflow, exactly as the sheet does in Exec/RegTests/NSCBC-FlameOutflow.
//      Truth is the identical initial condition on a domain five times longer,
//      whose own outflow is 40 cm from the common region -- far enough that
//      going from three times to five times longer moved the answer by 0.1 in
//      9754, which is the shielding argument checked rather than asserted.
//
//      THE RESULT IS NEGATIVE, AND IT IS THE POINT OF THE CHECK.  A source-free
//      expansion cannot be sustained in a constant-area duct.  Mass, momentum
//      and energy together admit no smooth steady state with du/dx != 0: with
//      rho*u and rho*u^2 + p held uniform the energy flux still varies as
//      (cp/R) p du/dx ~ 3.5e6 * du/dx, which is ~20% of rho*E per 3e-4 s.  The
//      run confirms it -- du/dn at the measurement point falls from 449 to -2
//      in the REFERENCE, whose boundary is 40 cm away and cannot be blamed.  A
//      flame's expansion is held up by its heat release; take the chemistry out
//      and the expansion does not survive long enough to be advected out.
//
//      So the two quantitative predictions from C9(a) -- error proportional to
//      du_out, error proportional to 1/sigma -- cannot be tested here, and are
//      reported rather than gated.  C9(b) and C10 are the two horns: hold the
//      gradient with a source and the source dominates; remove the source and
//      the gradient does not persist.
//
//      What survives is the order control, which changes ONLY the outgoing
//      extrapolation and leaves everything else identical.  Source-free it
//      moves the accumulated pressure error by 5x (1477 -> 7591), so the
//      extrapolation does drive the solution.  But the sign is the opposite of
//      the one assumed: order 1, which has no slope and therefore none of the
//      C9(a) bias, is five times WORSE.  Removing the bias term is not a fix.
//      Likewise sigma: relaxing harder toward p_inf makes the error grow, not
//      shrink, because p_inf is not the correct pressure while a structure is
//      crossing the boundary.
void
check_extrapolation_drives_solution()
{
  const Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  char buf[256];

  const int n = 200;      // cells over cs.L = 10 cm, so dx = 0.05 cm
  const Real dx = cs.L / n;
  const Real u0 = 3.0e2;  // cm/s, M ~ 0.01: quasi-frozen
  const Real wr = 1.0;    // ramp half-width, cm -- 20 cells, and wide enough
                          // that advecting it 0.2 cm leaves the gradient at the
                          // boundary essentially unchanged over the run
  const Real t_end = 3.0e-4; // ~1 relaxation time at sigma = 1
  const int NS = 6;          // history samples

  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  const Real mdot = rho0 * u0;

  // The ramp, centred on the outflow face at x = cs.L.
  auto uofx = [&](const Real x, const Real ratio) {
    const Real g = 0.5 * (1.0 + std::tanh((x - cs.L) / wr));
    return u0 * (1.0 + (ratio - 1.0) * g);
  };

  // Fill `nc` cells.  rho*u is uniform, and the pressure carries the momentum
  // flux so that rho*u^2 + p is uniform too -- otherwise the initial state is
  // out of momentum balance by mdot*u0*(ratio-1) and rings acoustically from
  // the first step for reasons that have nothing to do with the boundary.
  auto init = [&](Field& f, const int nc, const Real ratio) {
    for (int i = -NG; i < nc + NG; i++) {
      const Real x = (i + 0.5) * dx;
      const Real u = uofx(x, ratio);
      const Real rho = mdot / u;
      const Real p = cs.p0 + mdot * (u0 - u);
      Real T = 0.0;
      eos.RYP2T(rho, Y, p, T);
      set_state(f.at(i), rho, u, T, Y);
    }
  };

  // Advance `nc` cells to t_end with a characteristic outflow, sampling the
  // mean pressure over the FIRST n cells (the region the short and the long
  // domain share) and du/dn at the outflow face at NS times.
  // bmode: 0 = plain kernel fill, 1 = extrap_material, 2 = profile-fitU --
  // the C11x closure (fit the tanh family's position and thickness to T at
  // the last two interior cells each fill; overwrite ghost T and u from the
  // fitted profile, rho from the EOS at the kernel's ghost pressure; p stays
  // the relaxation's).  Here the family is the run's own initial ramp, but
  // the STRUCTURE IS DECAYING: the release question is whether the stateless
  // re-fit lets it go (tracking the weakening interior slope, and bailing to
  // the plain kernel when the structure is gone) or holds it alive at the
  // face the way extrap_material measurably does.
  auto run = [&](const int nc, const Real sigma, const int order,
                 const Real ratio, const int bmode, Real* pm, Real* gr) {
    Field f(nc), g(nc), h(nc);
    init(f, nc, ratio);

    pc_nscbc::Params prm;
    prm.L_ref = nc * dx;
    prm.sigma = sigma;
    prm.order = order;
    prm.extrap_material = (bmode == 1);
    pc_nscbc::Target off, out;
    out.type = pc_nscbc::Type::outflow;
    out.p = cs.p0;

    auto sample = [&](const int k) {
      Real psum = 0.0;
      for (int i = 0; i < n; i++) {
        psum += get_prim(f.at(i)).p;
      }
      pm[k] = psum / n;
      gr[k] = (get_prim(f.at(n - 1)).u - get_prim(f.at(n - 2)).u) / dx;
    };

    Real t = 0.0;
    int k = 0;
    sample(k++);
    while (k < NS) {
      const Real t_next = t_end * static_cast<Real>(k) / (NS - 1);
      while (t < t_next) {
        Real cmax = 0.0;
        for (int i = 0; i < nc; i++) {
          const Prim q = get_prim(f.at(i));
          cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
        }
        Real dt = std::min(cs.cfl * dx / cmax, t_next - t);
        if (dt <= 0.0) {
          break;
        }
        // Fixed inflow: this test is about the outflow only.
        for (int layer = 1; layer <= NG; layer++) {
          set_state(f.at(-layer), rho0, u0, cs.T0, Y);
        }
        fill_bcs(f, off, out, prm, dx);
        if (bmode >= 2) {
          fit_profile_ghosts(
            f, dx, u0, ratio, mdot, cs.p0, Y, true,
            (bmode == 3) ? 0.0 : -1.0);
        }
        stage(f, g, dx, dt);
        for (int layer = 1; layer <= NG; layer++) {
          set_state(g.at(-layer), rho0, u0, cs.T0, Y);
        }
        fill_bcs(g, off, out, prm, dx);
        if (bmode >= 2) {
          fit_profile_ghosts(
            g, dx, u0, ratio, mdot, cs.p0, Y, true,
            (bmode == 3) ? 0.0 : -1.0);
        }
        stage(g, h, dx, dt);
        for (int i = 0; i < nc; i++) {
          for (int v = 0; v < NVAR; v++) {
            f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
          }
        }
        t += dt;
      }
      sample(k++);
    }
  };

  // One shielded reference per (sigma, order, ratio).  Its own outflow is 40 cm
  // from the common region: the expanded gas reaches ~1200 K, where c ~ 7e4
  // cm/s, so 20 cm would only just be shielded over t_end and 40 cm is not.
  auto err = [&](const Real sigma, const int order, const Real ratio,
                 const int bmode, Real* eh, Real* gh, Real* gr) {
    Real pr[NS], pt[NS];
    run(5 * n, sigma, order, ratio, bmode, pr, gr);
    run(n, sigma, order, ratio, bmode, pt, gh);
    for (int k = 0; k < NS; k++) {
      eh[k] = pt[k] - pr[k];
    }
  };

  Real eh[NS], gh[NS], gr[NS];
  auto row = [&](const Real sigma, const int order, const Real ratio,
                 const int bmode = 0) {
    const char* bl[4] = {"  ", " m", " f", " b"};
    err(sigma, order, ratio, bmode, eh, gh, gr);
    std::printf("    %6.2f %6d %6.1f%s", sigma, order, ratio, bl[bmode]);
    for (int k = 1; k < NS; k++) {
      std::printf(" %9.1f", eh[k]);
    }
    std::printf(
      "   | du/dn %5.0f ->%5.0f (ref %5.0f ->%5.0f)\n", gh[0], gh[NS - 1],
      gr[0], gr[NS - 1]);
    return eh[NS - 1];
  };

  std::printf("\n    %6s %6s %6s   %s\n", "sigma", "order", "ratio",
              "<p>-<p>_ref at t_end/5 .. t_end");

  // --- prediction 1: error proportional to du_out ------------------------
  const Real e_r2 = row(1.0, 2, 2.0);
  const Real g_r2 = gh[0];
  const Real e_r4 = row(1.0, 2, 4.0);
  const Real g_r4 = gh[0];
  const Real g_ratio = g_r4 / (std::abs(g_r2) > 1e-30 ? g_r2 : 1.0);
  const Real e_ratio = e_r4 / (std::abs(e_r2) > 1e-30 ? e_r2 : 1.0);
  std::snprintf(
    buf, sizeof(buf),
    "du/dn(0) x%.2f -> error x%.2f -- not a law, the ramp is gone by t_end",
    g_ratio, e_ratio);
  report("error vs the velocity gradient", buf);

  // --- prediction 2: error proportional to 1/sigma -----------------------
  const Real e_lo = row(0.5, 2, 4.0);
  row(2.0, 2, 4.0);
  const Real e_hi = row(8.0, 2, 4.0);
  std::snprintf(
    buf, sizeof(buf),
    "sigma 0.5 -> 8 (x16): error %.1f -> %.1f -- stronger relaxation is WORSE",
    e_lo, e_hi);
  report("error vs the relaxation strength", buf);

  // --- the order control, which is the whole point -----------------------
  const Real e_o1 = row(1.0, 1, 4.0);
  std::snprintf(
    buf, sizeof(buf), "order 2 error %.1f, order 1 error %.1f (order 1 worse)",
    e_r4, e_o1);
  check(
    std::abs(e_o1 - e_r4) > 0.25 * std::abs(e_r4),
    "the outgoing extrapolation drives the solution", buf);

  // --- extrap_material on the same configuration -------------------------
  // Reported, not gated, and the sign of the result matters more than its
  // size: the continuation HOLDS the ramp at the boundary (du/dn stays high
  // where the reference decays to zero), so the short domain keeps venting a
  // structure whose exact solution is busy dying.  That is not the fix
  // misbehaving -- it is C10's own negative result seen from the other side:
  // a source-free ramp is not sustainable, so a boundary condition that
  // faithfully continues the structure it sees disagrees with a reference in
  // which that structure decays.  The gate for the fix is C11, where the
  // ramp is SUSTAINED and the exact answer is to hold it.  What this row
  // establishes is the honest cost: do not leave extrap_material on at a
  // boundary whose structure is transient and should be allowed to die out.
  const Real e_mat = row(1.0, 2, 4.0, 1);
  std::snprintf(
    buf, sizeof(buf),
    "error at t_end: %.1f without, %.1f with -- the continuation holds a "
    "DECAYING ramp alive at the face",
    e_r4, e_mat);
  report("extrap_material on a decaying (unsustainable) ramp", buf);

  // --- C11x release test: the profile-fit on the same dying ramp ---------
  // The C11 result (fitU holds a SUSTAINED front at oracle level) is only
  // half a qualification: extrap_material also helps there and then fails
  // HERE, holding this unsustainable ramp alive at the face while the
  // reference lets it die.  A usable structure closure must pass both.  The
  // hoped-for release mechanism is in the fit itself: value-and-slope
  // matching tracks the weakening interior gradient (lambda grows as the
  // ramp dies), and when no monotone structure remains the fit declines and
  // the plain kernel fill stands.  Reported, not gated.
  const Real e_fit = row(1.0, 2, 4.0, 2);
  const Real g_fit_end = gh[NS - 1];
  std::snprintf(
    buf, sizeof(buf),
    "error at t_end: plain %.1f, extrap_material %.1f, profile-fitU %.1f; "
    "du/dn at face -> %.0f (ref -> %.0f)",
    e_r4, e_mat, e_fit, g_fit_end, gr[NS - 1]);
  report("C11x release: the fit on a ramp that must be let go", buf);

  // --- and the same fit behind the source-consistency bound ---------------
  // The closure measures the local source (here exactly zero), computes the
  // sustainable dilatation, and refuses any continuation beyond it.  On a
  // source-free ramp that is a full refusal, so this row should reproduce
  // the plain-kernel row: the release is complete, by construction rather
  // than by decay-tracking.
  const Real e_fitB = row(1.0, 2, 4.0, 3);
  std::snprintf(
    buf, sizeof(buf),
    "error at t_end: plain %.1f, unbounded fitU %.1f, source-bounded %.1f "
    "-- the bound releases what the source cannot sustain",
    e_r4, e_fit, e_fitB);
  report("C11x release under the source-consistency bound", buf);
}

// C11: the sustained ramp -- a front crossing the outflow, with an exact
//      steady solution.  This is the test C9(b) and C10 both said was
//      missing: C9(b)'s heat band could not hold a velocity gradient at the
//      boundary without dumping unmodelled energy into the boundary cells,
//      and C10's source-free ramp could not hold it at all.  The resolution
//      is manufactured: take the mass- and momentum-consistent ramp of C10
//      (rho u uniform, p carrying the momentum flux) and add the energy
//      source S_E(x) = mdot dH/dx that makes it an EXACT steady solution of
//      the sourced Euler equations -- which is precisely a flame's mechanical
//      structure sustained by its heat release, minus the chemistry.
//
//      The exact solution is the initial condition, indefinitely.  A perfect
//      boundary holds it; every departure is boundary error.  The entropy
//      closure converts the ramp's du/dn into ghost pressure (C9(a)) and the
//      domain drifts to the equilibrium offset that NSCBC-FlameOutflow
//      measures; extrap_material continues the structure and must hold both
//      the mean pressure AND the du/dn at the face.
void
check_sustained_ramp()
{
  const Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  char buf[256];

  const int n = 200;
  const Real dx = cs.L / n;
  const Real u0 = 3.0e2;
  const Real ratio = 4.0;
  const Real wr = 1.0;
  const Real t_end = 6.0e-4; // ~2 relaxation times at sigma = 1
  const int NS = 4;

  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  const Real mdot = rho0 * u0;

  // The TRUE front shape.  0 = tanh (the family the fit closure assumes);
  // 1 = a Richards curve with k = 3, centred so g = 1/2 at the boundary:
  // the same end states and width scale but genuinely outside the tanh
  // family -- e^{3 xi} approach on the fresh side, a slow e^{-xi} tail on
  // the burnt side, the way a real flame's preheat structure is one-sided.
  int true_shape = 0;
  auto uofx = [&](const Real x) {
    const Real xi = (x - cs.L) / wr;
    Real g;
    if (true_shape == 0) {
      g = 0.5 * (1.0 + std::tanh(xi));
    } else {
      const Real x0 = 1.3475805944; // -ln(2^{1/3} - 1): g(0) = 1/2
      g = std::pow(1.0 + std::exp(-(xi + x0)), -3.0);
    }
    return u0 * (1.0 + (ratio - 1.0) * g);
  };
  // Total enthalpy of the exact profile at x.
  auto Hofx = [&](const Real x) {
    const Real u = uofx(x);
    const Real rho = mdot / u;
    const Real p = cs.p0 + mdot * (u0 - u);
    Real T = 0.0, e = 0.0;
    eos.RYP2T(rho, Y, p, T);
    eos.RTY2E(rho, T, Y, e);
    return e + p / rho + 0.5 * u * u;
  };

  auto init = [&](Field& f, const int nc) {
    for (int i = -NG; i < nc + NG; i++) {
      const Real x = (i + 0.5) * dx;
      const Real u = uofx(x);
      const Real rho = mdot / u;
      const Real p = cs.p0 + mdot * (u0 - u);
      Real T = 0.0;
      eos.RYP2T(rho, Y, p, T);
      set_state(f.at(i), rho, u, T, Y);
    }
  };

  // The manufactured energy source, cell-centred, frozen in time.
  auto make_source = [&](const int nc) {
    std::vector<Real> q(nc, 0.0);
    for (int i = 0; i < nc; i++) {
      const Real x = (i + 0.5) * dx;
      q[i] = mdot * (Hofx(x + 0.5 * dx) - Hofx(x - 0.5 * dx)) / dx;
    }
    return q;
  };

  auto add_source = [&](Field& f, const std::vector<Real>& q, const Real dt) {
    for (int i = 0; i < f.n; i++) {
      Real* s = f.at(i);
      s[UEDEN] += dt * q[i];
      s[UEINT] += dt * q[i];
      const Real rho = s[URHO];
      Real Yl[NUM_SPECIES], ys = 0.0;
      for (int k = 0; k < NUM_SPECIES; k++) {
        Yl[k] = std::max(s[UFS + k] / rho, 0.0);
        ys += Yl[k];
      }
      for (int k = 0; k < NUM_SPECIES; k++) {
        Yl[k] /= ys;
      }
      const Real e = s[UEINT] / rho;
      Real T = s[UTEMP] > 0.0 ? s[UTEMP] : 300.0;
      eos.REY2T(rho, e, Yl, T);
      s[UTEMP] = T;
    }
  };

  // Advance nc cells to t_end; sample mean p over the first n cells and
  // du/dn at x = cs.L (the short domain's outflow face) at NS times.
  // mode: 0 = entropy closure, 1 = extrap_material, 2 = ORACLE -- the hi
  // ghosts are overwritten with the exact profile every fill, i.e. the best
  // any ghost-cell closure can possibly do.  The oracle is the yardstick the
  // closures are gated against: the truncated discrete problem has its own
  // attractor (see below), and no ghost fill can beat the oracle's.
  auto run = [&](const int nc, const Real sigma, const int mode, Real* pm,
                 Real* gr) {
    Field f(nc), g(nc), h(nc);
    init(f, nc);
    const auto q = make_source(nc);

    auto oracle_ghosts = [&](Field& w) {
      auto eos2 = pele::physics::PhysicsType::eos();
      for (int layer = 1; layer <= NG; layer++) {
        const int i = w.n - 1 + layer;
        const Real x = (i + 0.5) * dx;
        const Real u = uofx(x);
        const Real rho = mdot / u;
        const Real p = cs.p0 + mdot * (u0 - u);
        Real T = 0.0;
        eos2.RYP2T(rho, Y, p, T);
        set_state(w.at(i), rho, u, T, Y);
      }
    };

    pc_nscbc::Params prm;
    prm.L_ref = nc * dx;
    prm.sigma = sigma;
    prm.extrap_material = (mode == 1);
    prm.extrap_temperature = (mode == 3);
    pc_nscbc::Target off, out;
    out.type = pc_nscbc::Type::outflow;
    out.p = cs.p0 + mdot * (u0 - uofx(cs.L)); // the exact face pressure

    auto sample = [&](const int k) {
      Real psum = 0.0;
      for (int i = 0; i < n; i++) {
        psum += get_prim(f.at(i)).p;
      }
      pm[k] = psum / n;
      gr[k] = (get_prim(f.at(n - 1)).u - get_prim(f.at(n - 2)).u) / dx;
    };

    Real t = 0.0;
    int k = 0;
    sample(k++);
    while (k < NS) {
      const Real t_next = t_end * static_cast<Real>(k) / (NS - 1);
      while (t < t_next) {
        Real cmax = 0.0;
        for (int i = 0; i < nc; i++) {
          const Prim qq = get_prim(f.at(i));
          cmax = std::max(cmax, std::abs(qq.u) + sound_speed(qq));
        }
        Real dt = std::min(cs.cfl * dx / cmax, t_next - t);
        if (dt <= 0.0) {
          break;
        }
        for (int layer = 1; layer <= NG; layer++) {
          set_state(f.at(-layer), rho0, u0, cs.T0, Y);
        }
        fill_bcs(f, off, out, prm, dx);
        if (mode == 2) {
          oracle_ghosts(f);
        } else if (mode >= 4) {
          fit_profile_ghosts(
            f, dx, u0, (mode == 5 || mode == 7) ? 0.85 * ratio : ratio,
            mdot, cs.p0, Y, mode >= 6, (mode == 8) ? q[nc - 1] : -1.0);
        }
        stage(f, g, dx, dt);
        add_source(g, q, dt);
        for (int layer = 1; layer <= NG; layer++) {
          set_state(g.at(-layer), rho0, u0, cs.T0, Y);
        }
        fill_bcs(g, off, out, prm, dx);
        if (mode == 2) {
          oracle_ghosts(g);
        } else if (mode >= 4) {
          fit_profile_ghosts(
            g, dx, u0, (mode == 5 || mode == 7) ? 0.85 * ratio : ratio,
            mdot, cs.p0, Y, mode >= 6, (mode == 8) ? q[nc - 1] : -1.0);
        }
        stage(g, h, dx, dt);
        add_source(h, q, dt);
        for (int i = 0; i < nc; i++) {
          for (int v = 0; v < NVAR; v++) {
            f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
          }
        }
        t += dt;
      }
      sample(k++);
    }
  };

  const Real g0 =
    (uofx(cs.L - 0.5 * dx) - uofx(cs.L - 1.5 * dx)) / dx; // initial du/dn

  // ---- The face flux, statically -----------------------------------------
  // Before any dynamics: fill the ghosts from the exact interior state with
  // each closure, reconstruct the boundary face exactly as stage() does, and
  // compare the HLLC flux against the exact steady flux the face must carry
  // (mdot, mdot u + p, mdot H).  This is the flux-level form of C9(a), and it
  // has no translational mode to hide behind (see below).
  {
    auto face_flux = [&](const bool mat, Real* flx) {
      Field f(n);
      init(f, n);
      pc_nscbc::Params prm;
      prm.L_ref = cs.L;
      prm.sigma = 1.0;
      prm.extrap_material = mat;
      pc_nscbc::Target off, out;
      out.type = pc_nscbc::Type::outflow;
      out.p = cs.p0 + mdot * (u0 - uofx(cs.L));
      fill_bcs(f, off, out, prm, dx);
      Real sl[NVAR], sr[NVAR];
      for (int v = 0; v < NVAR; v++) {
        const Real dL = mm(
          f.at(n - 1)[v] - f.at(n - 2)[v], f.at(n)[v] - f.at(n - 1)[v]);
        const Real dR =
          mm(f.at(n)[v] - f.at(n - 1)[v], f.at(n + 1)[v] - f.at(n)[v]);
        sl[v] = f.at(n - 1)[v] + 0.5 * dL;
        sr[v] = f.at(n)[v] - 0.5 * dR;
      }
      auto to_prim_face = [&](Real* s) {
        Prim q{};
        q.rho = s[URHO];
        q.u = s[UMX] / q.rho;
        Real ys = 0.0;
        for (int k = 0; k < NUM_SPECIES; k++) {
          q.Y[k] = std::max(s[UFS + k] / q.rho, 0.0);
          ys += q.Y[k];
        }
        for (int k = 0; k < NUM_SPECIES; k++) {
          q.Y[k] /= ys;
        }
        const Real e = s[UEDEN] / q.rho - 0.5 * q.u * q.u;
        q.T = 300.0;
        eos.REY2T(q.rho, e, q.Y, q.T);
        eos.RTY2P(q.rho, q.T, q.Y, q.p);
        return q;
      };
      hllc(to_prim_face(sl), to_prim_face(sr), flx);
    };
    const Real F_mass = mdot;
    const Real F_ener = mdot * Hofx(cs.L);
    Real fe[NVAR], fm[NVAR];
    face_flux(false, fe);
    face_flux(true, fm);
    const Real ee_ent = fe[UEDEN] - F_ener;
    const Real ee_mat = fm[UEDEN] - F_ener;
    const Real em_ent = fe[URHO] - F_mass;
    const Real em_mat = fm[URHO] - F_mass;
    std::snprintf(
      buf, sizeof(buf),
      "energy flux error: entropy %.3e, extrap_material %.3e (exact %.3e); "
      "mass: %.2e vs %.2e (exact %.2e)",
      ee_ent, ee_mat, F_ener, em_ent, em_mat, F_mass);
    check(
      std::abs(ee_mat) < 0.35 * std::abs(ee_ent) &&
        std::abs(em_mat) < 0.5 * std::abs(em_ent),
      "extrap_material corrects the boundary-face flux", buf);
  }

  // The reference: same sourced problem, outflow 4 L downstream.  It must
  // HOLD the steady state; its residual drift is the discretisation floor.
  Real pr[NS], gref[NS];
  run(5 * n, 1.0, 0, pr, gref);
  std::snprintf(
    buf, sizeof(buf),
    "reference <p> drift %.1f dyn/cm2 over %.0e s, du/dn %.0f -> %.0f (of "
    "%.0f)",
    pr[NS - 1] - pr[0], t_end, gref[0], gref[NS - 1], g0);
  check(
    std::abs(pr[NS - 1] - pr[0]) < 150.0,
    "the manufactured steady state holds in the long domain", buf);

  // The truncated domain has its OWN discrete attractor: cutting the ramp
  // mid-structure and representing its continuation with 4 ghost cells
  // shifts the balance point, for ANY ghost fill.  The oracle row measures
  // that attractor -- it is the floor no ghost-cell closure can beat -- so
  // the closures are judged by their distance from the oracle, not from the
  // long reference.
  std::printf(
    "\n     sigma  ghost   <p>-<p>_ref at t_end/3 .. t_end   | du/dn at the "
    "face\n");
  Real pt[NS], gt[NS];
  const char* label[9] = {"ent",  "mat",  "orc",   "e_T", "fit",
                          "fitX", "fitU", "fitUX", "fitB"};
  auto row = [&](const Real sigma, const int mode) {
    run(n, sigma, mode, pt, gt);
    std::printf("    %6.2f  %5s  ", sigma, label[mode]);
    for (int k = 1; k < NS; k++) {
      std::printf(" %9.1f", pt[k] - pr[k]);
    }
    std::printf("   | %5.0f -> %5.0f (exact %.0f)\n", gt[0], gt[NS - 1], g0);
    return pt[NS - 1] - pr[NS - 1];
  };

  const Real e_orc = row(1.0, 2);
  const Real g_orc = gt[NS - 1];
  const Real g1_orc = gt[1];
  const Real e1_orc = pt[1] - pr[1];
  const Real e_ent1 = row(1.0, 0);
  const Real g_ent = gt[NS - 1];
  const Real g1_ent = gt[1];
  const Real e1_ent = pt[1] - pr[1];
  const Real e_ent4 = row(4.0, 0);
  const Real e_mat1 = row(1.0, 1);
  const Real g_mat = gt[NS - 1];
  const Real g1_mat = gt[1];
  const Real e1_mat = pt[1] - pr[1];
  row(1.0, 3);
  const Real g_eT = gt[NS - 1];
  const Real e1_eT = pt[1] - pr[1];
  row(1.0, 4);
  const Real g_fit = gt[NS - 1];
  const Real e1_fit = pt[1] - pr[1];
  row(1.0, 5);
  const Real g_fitX = gt[NS - 1];
  const Real e1_fitX = pt[1] - pr[1];
  const Real eE_fitU = row(1.0, 6);
  const Real g_fitU = gt[NS - 1];
  const Real e1_fitU = pt[1] - pr[1];
  row(1.0, 7);
  const Real g_fitUX = gt[NS - 1];
  const Real e1_fitUX = pt[1] - pr[1];
  const Real eE_fitB = row(1.0, 8);
  const Real g_fitB = gt[NS - 1];
  const Real e1_fitB = pt[1] - pr[1];

  // What the oracle row establishes: a ghost fill that carries the exact
  // continuation HOLDS the sustained front, in this same truncated domain,
  // with this same solver -- so nothing below can be blamed on the
  // architecture.  What the late-time columns then measure is an artefact of
  // the MANUFACTURED source: q(x) is frozen in space, so once a closure's
  // early flux error has nudged the structure off its source the mismatch
  // feeds itself and every non-oracle run walks to the same shifted
  // equilibrium (~+27000 here, sigma-independent -- note sigma 1 vs 4).  A
  // real flame carries its source with its front and has no such mode, which
  // is why the gates below sit inside the boundary-equilibration window
  // (t_end/3 ~ 0.7 relaxation times), where the columns still measure the
  // boundary.
  std::snprintf(
    buf, sizeof(buf),
    "late-time error: ent %.1f, mat %.1f, oracle %.1f -- the frozen source, "
    "not the boundary",
    e_ent1, e_mat1, e_orc);
  report("the truncated MMS walks off its source at late time", buf);
  std::snprintf(
    buf, sizeof(buf), "sigma 1 -> 4 at late time: %.1f -> %.1f (the C9(a) "
    "bias would fall x4)",
    e_ent1, e_ent4);
  report("the late-time offset is not the C9(a) bias", buf);

  std::snprintf(
    buf, sizeof(buf),
    "error at 0.7 tau: %.1f with extrap_material, %.1f without, oracle %.1f",
    e1_mat, e1_ent, e1_orc);
  check(
    std::abs(e1_mat - e1_orc) < 0.35 * std::abs(e1_ent - e1_orc),
    "extrap_material holds the front while the boundary equilibrates", buf);

  // Gated at 0.7 tau, the same window as the pressure gate and for the same
  // reason (see the frozen-source caveat above); the t_end values are
  // reported for context but sit in the artifact-dominated region.
  std::snprintf(
    buf, sizeof(buf),
    "du/dn at 0.7 tau: %.0f with, %.0f without, oracle %.0f (exact %.0f; "
    "t_end: %.0f / %.0f / %.0f)",
    g1_mat, g1_ent, g1_orc, g0, g_mat, g_ent, g_orc);
  check(
    std::abs(g1_mat - g1_orc) < 0.5 * std::abs(g1_ent - g1_orc),
    "extrap_material preserves the structure the oracle preserves", buf);

  // ---- C11x: the profile-fit experiment (reported, not gated) ------------
  // The question: how much of the oracle's advantage does a FITTED library
  // profile recover, when it supplies material structure only?  The fraction
  // below is (e_T - fit) / (e_T - oracle) in the early window -- 1.0 means
  // the fit is as good as knowing the exact answer, 0.0 means it bought
  // nothing over linear-in-T ghosts.  fitX is the same fit through a family
  // whose end state is 15% wrong: the stretch/curvature analogue.
  {
    const Real denom = e1_eT - e1_orc;
    const Real frac =
      (std::abs(denom) > 1.0e-30) ? (e1_eT - e1_fit) / denom : 0.0;
    const Real fracX =
      (std::abs(denom) > 1.0e-30) ? (e1_eT - e1_fitX) / denom : 0.0;
    std::snprintf(
      buf, sizeof(buf),
      "at 0.7 tau: ent %.1f, e_T %.1f, fit %.1f (fitX %.1f), mat %.1f, "
      "fitU %.1f (fitUX %.1f), oracle %.1f",
      e1_ent, e1_eT, e1_fit, e1_fitX, e1_mat, e1_fitU, e1_fitUX, e1_orc);
    report("C11x the ladder from line to oracle", buf);
    const Real denomU = e1_mat - e1_orc;
    const Real fracU =
      (std::abs(denomU) > 1.0e-30) ? (e1_mat - e1_fitU) / denomU : 0.0;
    std::snprintf(
      buf, sizeof(buf),
      "material-only fit recovers %.0f%% of the e_T -> oracle gap (%.0f%% "
      "with the 15%%-wrong family); with the profile's u, %.0f%% of the "
      "mat -> oracle gap",
      100.0 * frac, 100.0 * fracX, 100.0 * fracU);
    report("C11x profile-fit recovery fractions", buf);
    std::snprintf(
      buf, sizeof(buf),
      "du/dn at t_end: e_T %.0f, fit %.0f, mat %.0f, fitU %.0f (fitUX "
      "%.0f), oracle %.0f (exact %.0f)",
      g_eT, g_fit, g_mat, g_fitU, g_fitUX, g_orc, g0);
    report("C11x structure at the face", buf);
    std::snprintf(
      buf, sizeof(buf),
      "source-bounded fitU on the SUSTAINED front: %.1f at 0.7 tau, %.1f at "
      "t_end, du/dn -> %.0f (unbounded fitU %.1f / %.1f / %.0f)",
      e1_fitB, eE_fitB, g_fitB, e1_fitU, eE_fitU, g_fitU);
    report("C11x the bound must not tax the sustained front", buf);
  }

  // ---- C11x shape distortion: the truth is not in the family -------------
  // The fitX row above distorted the family's END STATE and the fit absorbed
  // it; this block distorts the SHAPE.  The truth becomes the Richards front
  // (asymmetric, outside any tanh), the manufactured source and the oracle
  // follow it automatically, and the fit still assumes tanh.  What survives
  // of the 97% is the measure of how much the closure leans on the library
  // profile actually matching the flame.
  {
    true_shape = 1;
    Real prd[NS], grd[NS], ptd[NS], gtd[NS];
    run(5 * n, 1.0, 0, prd, grd); // fresh shielded reference, distorted truth
    const Real g0d =
      (uofx(cs.L - 0.5 * dx) - uofx(cs.L - 1.5 * dx)) / dx;
    std::printf(
      "\n     distorted truth (Richards k=3), tanh-family fit; exact du/dn "
      "%.0f\n",
      g0d);
    auto rowd = [&](const int mode) {
      run(n, 1.0, mode, ptd, gtd);
      std::printf("    %6.2f  %5s  ", 1.0, label[mode]);
      for (int k = 1; k < NS; k++) {
        std::printf(" %9.1f", ptd[k] - prd[k]);
      }
      std::printf(
        "   | %5.0f -> %5.0f (exact %.0f)\n", gtd[0], gtd[NS - 1], g0d);
      return ptd[1] - prd[1];
    };
    const Real d1_orc = rowd(2);
    const Real d1_ent = rowd(0);
    const Real d1_fitU = rowd(6);
    const Real dE_fitU = ptd[NS - 1] - prd[NS - 1];
    const Real d1_fitB = rowd(8);
    const Real denom = d1_ent - d1_orc;
    const Real frac =
      (std::abs(denom) > 1.0e-30) ? (d1_ent - d1_fitU) / denom : 0.0;
    std::snprintf(
      buf, sizeof(buf),
      "at 0.7 tau: ent %.1f, fitU %.1f (t_end %.1f), fitB %.1f, oracle %.1f "
      "-- the tanh fit recovers %.0f%% of the gap on a non-tanh truth",
      d1_ent, d1_fitU, dE_fitU, d1_fitB, d1_orc, 100.0 * frac);
    report("C11x shape-distorted truth", buf);
    true_shape = 0;
  }
}

// C7: the reaction source term.  Chemistry changes the pressure only through
//     the composition, so the closed-form expression in reaction_dpdt() is
//     checked against a DIRECTIONAL finite difference along the reaction path
//     at fixed rho and e -- refining tau, which should show first-order
//     convergence toward the analytic value.  That FD is also exactly what the
//     real-gas branch of reaction_dpdt() computes, so this validates both
//     paths at once.
void
check_reaction_source()
{
#if NUM_REACTIONS == 0
  std::printf(
    "  SKIP  reaction source: %s has no reactions\n",
    pele::physics::PhysicsType::identifier().c_str());
#else
  auto eos = pele::physics::PhysicsType::eos();
  // Stoichiometric H2/air, hot enough to react briskly.
  Real Y[NUM_SPECIES] = {0.0};
  Y[H2_ID] = 0.0285;
  Y[O2_ID] = 0.2265;
  Y[N2_ID] = 0.7450;
  const Real p0 = 1.01325e6, T0 = 1400.0;
  Real rho = 0.0, e0 = 0.0;
  eos.PYT2RE(p0, Y, T0, rho, e0);

  Real s[NVAR];
  set_state(s, rho, 0.0, T0, Y);
  const pc_nscbc::CellPrim q = pc_nscbc::cell_primitives(s);
  bool ok = true;
  const Real analytic = pc_nscbc::reaction_dpdt(q, ok);

  // Directional FD: Y'(tau) = Y + tau * wdot / rho, which preserves sum Y = 1
  // because sum wdot = 0; then T' from (rho, e) and p' from (rho, T', Y').
  Real wdot[NUM_SPECIES];
  eos.RTY2WDOT(q.rho, q.T, q.Y, wdot);
  Real wsum = 0.0, wmax = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    wsum += wdot[n];
    wmax = std::max(wmax, std::abs(wdot[n]));
  }
  auto fd = [&](Real tau) {
    Real Yp[NUM_SPECIES], sum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Yp[n] = std::max(q.Y[n] + tau * wdot[n] / q.rho, 0.0);
      sum += Yp[n];
    }
    for (int n = 0; n < NUM_SPECIES; n++) {
      Yp[n] /= sum;
    }
    Real Tp = q.T, pp = 0.0;
    eos.REY2T(q.rho, q.e, Yp, Tp);
    eos.RTY2P(q.rho, Tp, Yp, pp);
    return (pp - q.p) / tau;
  };

  char buf[220];
  std::snprintf(
    buf, sizeof(buf), "|sum wdot| / max|wdot| = %.3e",
    std::abs(wsum) / std::max(wmax, 1e-300));
  check(
    std::abs(wsum) / std::max(wmax, 1e-300) < 1e-10,
    "chemistry conserves mass (sum wdot = 0)", buf);

  const Real tau0 = 1.0e-6 / (wmax / q.rho);
#ifdef USE_SRK_EOS
  // Under SRK the kernel's reaction_dpdt IS a directional FD, so refining a
  // second FD toward it measures nothing but the difference of two step
  // sizes.  The meaningful gate is agreement at a matching step.
  const Real v0 = fd(tau0);
  std::snprintf(
    buf, sizeof(buf), "kernel %.6e vs test FD %.6e (rel %.2e)", analytic, v0,
    std::abs(analytic - v0) / std::abs(v0));
  check(
    ok && std::abs(analytic - v0) / std::abs(v0) < 1e-2,
    "real-gas FD path agrees with an independent FD", buf);
#else
  // A refining FD trades truncation error for cancellation error, and where
  // the crossover lands depends on the toolchain (llvm/arm64 hits the floor
  // one refinement earlier than gcc/x86).  So the gate is floor-aware:
  // errors must decrease monotonically UNTIL the minimum, and the minimum
  // must sit at the round-off floor -- not monotone-to-the-end, which gates
  // the machine rather than the mathematics.
  Real err_k[4], min_err = 1e300;
  int argmin = 0;
  std::printf("      analytic dp/dt|_react = %.6e dyn/(cm^2 s)\n", analytic);
  for (int k = 0; k < 4; k++) {
    const Real tau = tau0 / std::pow(4.0, k);
    const Real v = fd(tau);
    err_k[k] = std::abs(v - analytic) / std::abs(analytic);
    std::printf(
      "      tau = %.3e  FD = %.6e  rel err = %.3e\n", tau, v, err_k[k]);
    if (err_k[k] < min_err) {
      min_err = err_k[k];
      argmin = k;
    }
  }
  bool converging = true;
  for (int k = 1; k <= argmin; k++) {
    converging = converging && (err_k[k] < err_k[k - 1] * 1.05);
  }
  std::snprintf(buf, sizeof(buf), "final relative error %.3e", err_k[3]);
  check(
    ok && err_k[3] < 1e-3, "closed-form dp/dt matches the directional FD",
    buf);
  std::snprintf(
    buf, sizeof(buf),
    "error falls monotonically to %.2e at tau/%.0f, then sits at the "
    "round-off floor",
    min_err, std::pow(4.0, argmin));
  check(
    converging && (min_err < 1e-6),
    "FD converges to the closed form's round-off floor", buf);
#endif

  // A frozen (cold) state must give exactly zero, so beta_s is a no-op there.
  Real s_cold[NVAR];
  Real rho_c = 0.0, e_c = 0.0;
  eos.PYT2RE(p0, Y, 300.0, rho_c, e_c);
  set_state(s_cold, rho_c, 0.0, 300.0, Y);
  const pc_nscbc::CellPrim qc = pc_nscbc::cell_primitives(s_cold);
  bool ok_c = true;
  const Real cold = pc_nscbc::reaction_dpdt(qc, ok_c);
  std::snprintf(
    buf, sizeof(buf), "dp/dt = %.3e at 300 K vs %.3e at 1400 K", cold,
    analytic);
  check(
    ok_c && std::abs(cold) < 1e-6 * std::abs(analytic),
    "cold, unreacting state gives a negligible source", buf);
#endif
}

// C12: where the diffusive capability gap actually lives.
//
//      Real conduction is switched on in the mini solver (g_lambda), so the
//      boundary heat flux is formed from the ghost temperatures exactly as
//      PeleC's diffusion operator forms it, and a temperature ramp is parked
//      with its high-curvature flank in the short domain's outflow cells.
//      Against a 5x shielded reference (the C10/C11 protocol), this gates
//      DYNAMICALLY what C8 gates statically: the conductive boundary error
//      belongs to the ghost TEMPERATURE CLOSURE, and extrap_temperature
//      removes most of it.
//
//      It is also where a would-be "viscous condition" for the incoming wave
//      went to die, and the numbers are worth keeping: a modelled
//      dp/dt|_diffusion of the boundary cell (normal conduction + species
//      diffusion of the resolved fields, EOS-exact on quadratic profiles to
//      1e-9) moved the error measured here from +104 to -911 dyn/cm2, and
//      PeleC's flame-outflow at sigma = 1 from +1200 to +1771.  In the
//      ghost-cell form the diffusion operator READS the ghosts, so a correct
//      ghost closure already carries the diffusive physics and an
//      amplitude-side term double-counts it.  The term was removed; this
//      check keeps the closure honest instead.
void
check_diffusive_dynamics(const Real lam_cond)
{
  const Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  char buf[256];

  const int n = 200;
  const Real dxs = cs.L / n;
  const Real u0 = 3.0e2;
  const Real wr = 1.0;
  const Real Tratio = 4.0;
  const Real xc = cs.L - wr; // high-curvature flank in the outflow cells
  const Real t_end = 3.0e-4;
  const int NS = 4;

  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);

  auto Tofx = [&](const Real x) {
    const Real g = 0.5 * (1.0 + std::tanh((x - xc) / wr));
    return cs.T0 * (1.0 + (Tratio - 1.0) * g);
  };
  auto init = [&](Field& f, const int nc) {
    for (int i = -NG; i < nc + NG; i++) {
      const Real x = (i + 0.5) * dxs;
      const Real T = Tofx(x);
      Real rho = 0.0, e = 0.0;
      eos.PYT2RE(cs.p0, Y, T, rho, e);
      set_state(f.at(i), rho, u0, T, Y);
    }
  };

  auto run = [&](const int nc, const bool extrap_T, Real* pm) {
    Field f(nc), g(nc), h(nc);
    init(f, nc);
    pc_nscbc::Params prm;
    prm.L_ref = nc * dxs;
    prm.sigma = 1.0;
    prm.extrap_temperature = extrap_T;
    pc_nscbc::Target off, out;
    out.type = pc_nscbc::Type::outflow;
    out.p = cs.p0;
    Real t = 0.0;
    int k = 0;
    auto sample = [&](const int kk) {
      Real psum = 0.0;
      for (int i = 0; i < n; i++) {
        psum += get_prim(f.at(i)).p;
      }
      pm[kk] = psum / n;
    };
    sample(k++);
    while (k < NS) {
      const Real t_next = t_end * static_cast<Real>(k) / (NS - 1);
      while (t < t_next) {
        Real cmax = 0.0;
        for (int i = 0; i < nc; i++) {
          const Prim q = get_prim(f.at(i));
          cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
        }
        Real dt = std::min(cs.cfl * dxs / cmax, t_next - t);
        if (dt <= 0.0) {
          break;
        }
        for (int layer = 1; layer <= NG; layer++) {
          set_state(f.at(-layer), rho0, u0, cs.T0, Y);
        }
        fill_bcs(f, off, out, prm, dxs);
        stage(f, g, dxs, dt);
        for (int layer = 1; layer <= NG; layer++) {
          set_state(g.at(-layer), rho0, u0, cs.T0, Y);
        }
        fill_bcs(g, off, out, prm, dxs);
        stage(g, h, dxs, dt);
        for (int i = 0; i < nc; i++) {
          for (int v = 0; v < NVAR; v++) {
            f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
          }
        }
        t += dt;
      }
      sample(k++);
    }
  };

  g_lambda = lam_cond; // conduction ON, in both domains
  Real pr[NS], p_ent[NS], p_tmp[NS];
  run(5 * n, false, pr);
  run(n, false, p_ent);
  run(n, true, p_tmp);
  g_lambda = 0.0;

  const Real e_ent = p_ent[NS - 1] - pr[NS - 1];
  const Real e_tmp = p_tmp[NS - 1] - pr[NS - 1];
  std::snprintf(
    buf, sizeof(buf),
    "error vs shielded reference: entropy closure %.1f, extrap_temperature "
    "%.1f (x%.2f)",
    e_ent, e_tmp, e_tmp / (std::abs(e_ent) > 1e-30 ? e_ent : 1.0));
  check(
    std::abs(e_tmp) < 0.5 * std::abs(e_ent),
    "resolved conduction is handled by the ghost T closure", buf);
}

// ===========================================================================
//  Duct modes (CLI: t1 / t3 / t5) -- the injection-fidelity bed.
//
//  One duct, two experiments, run only when asked (not part of the default
//  suite): an NSCBC value-relaxation INFLOW at the lo face and a hard
//  pressure ghost at the hi face (FOExtrap with p pinned to p0 -- the fully
//  reflecting open end, R ~ -1).
//
//  T1 (Dupuy 2026 Figs 5a/8a): step the inlet velocity target at t = 0 and
//  measure the convergence time versus relax_u.  The CLR theory predicts a
//  narrow viable band with an interior optimum near 0.3 and growth on both
//  sides (drift below, reflection-delayed convergence above); that
//  relationship is the gate.
//
//  T3 (Daviller 2019 sec. 4-5): force the inlet target harmonically,
//  u^t = u0 + A sin(2 pi f t), at frequencies straddling the duct's
//  quarter-wave resonance, and measure what actually enters: the achieved
//  velocity amplitude at the boundary cell over A (I_u), and the incoming
//  invariant's amplitude over the ideal injector's 2A (I_in).  Our inflow
//  has no amplitude slot (design doc I.3e) -- a time-varying target rides
//  the same relaxation that holds the mean -- so injection fidelity is
//  relax_u- and frequency-dependent BY CONSTRUCTION, and relax_u = 0
//  injects nothing at all.  Those two properties are the gates; the
//  deterioration near resonance is the reported measurement, with the
//  stiff-limit analytic response 1/|1 + e^{i theta}| (theta the Doppler-
//  corrected round-trip phase; divergent at the quarter-wave) as the
//  reference column.
//
//  T5 (Daviller 2019 sec. 7, laminar core): same runs, but the deliverable
//  is P_RMS(x) against the analytic standing-wave envelope |sin| shape --
//  node/antinode geometry gated by shape correlation, antinode amplitude
//  reported.
// ===========================================================================

// FOExtrap with pressure pinned: the fully reflecting hi end of the duct.
void
hard_p_fill_hi(Field& f, Real p0)
{
  auto eos = pele::physics::PhysicsType::eos();
  const int n = f.n;
  const Prim q = get_prim(f.at(n - 1));
  Real rho = 0.0, e = 0.0;
  eos.PYT2RE(p0, q.Y, q.T, rho, e);
  for (int layer = 1; layer <= NG; layer++) {
    set_state(f.at(n - 1 + layer), rho, q.u, q.T, q.Y);
  }
}

struct DuctForcedResult
{
  Real I_u = 0.0;  // achieved u' amplitude at the inlet cell / target A
  Real I_in = 0.0; // incoming-invariant amplitude / ideal injector's 2A
  std::vector<Real> prms; // P_RMS(x) per cell over the measurement window
};

// Harmonic forcing of the inflow target against the reflecting far end.
// Settles for n_settle acoustic round trips (forcing on throughout), then
// projects n_per integer periods.  The Fourier projection at f ignores DC,
// so the mean-flow offset never needs subtracting.
//
// nri = true prototypes the BOUNDARY-REGISTER architecture (design doc II.2,
// queue item 4) without touching the kernel: the driver -- standing in for
// the once-per-advance register update PeleC would own -- keeps an EMA of
// the outgoing invariant R- = u - p/(rho c) at the boundary cell
// (tau = 3 t_a, the NDNR prescription), splits the instantaneous outgoing
// wave off it, u_minus = (R- - EMA)/2, and hands the kernel a target
// u^t + u_minus.  The relaxation then fights only the incoming content and
// the mean drift; the wave the duct sends back is granted passage instead
// of being mistaken for error.  The kernel remains a pure function of
// (state, Target): all the state lives in the caller, which is the entire
// architectural claim being measured.
// mode: 0 = classical, 1 = driver-held NRI register, 2 = stateless
// feed-forward (Target::dudt), 3 = both.
DuctForcedResult
duct_forced(
  Real relax_u, Real freq, int n, int n_settle, int n_per, int mode = 0)
{
  Case cs;
  cs.n = n;
  Field f(cs.n), g(cs.n), h(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  Prim q0{};
  q0.rho = rho0;
  q0.u = 0.0;
  q0.p = cs.p0;
  q0.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    q0.Y[k] = Y[k];
  }
  const Real c0 = sound_speed(q0);
  const Real u0 = 2.0e3; // mean inflow; forcing amplitude rides well below it
  const Real A = 1.0e-3 * c0; // ~35 cm/s: linear regime

  for (int i = -NG; i < cs.n + NG; i++) {
    set_state(f.at(i), rho0, u0, cs.T0, Y);
  }

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.relax_u = relax_u;
  prm.relax_t = 0.2;
  pc_nscbc::Target in, off;
  in.type = pc_nscbc::Type::inflow;
  in.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    in.Y[k] = Y[k];
  }

  const Real t_a = 2.0 * cs.L / c0;
  const Real t_meas0 = n_settle * t_a;
  const Real T_w = n_per / freq;
  const Real t_end = t_meas0 + T_w;
  const Real om = 2.0 * M_PI * freq;

  DuctForcedResult r;
  r.prms.assign(cs.n, 0.0);
  std::vector<Real> psum(cs.n, 0.0), p2sum(cs.n, 0.0);
  Real us = 0.0, uc = 0.0, rs = 0.0, rc = 0.0, wsum = 0.0;

  // The register: an EMA of the outgoing invariant at the boundary cell.
  // Seeded with the quiescent value; tau = 3 t_a per the NDNR prescription.
  const Real rhoc0 = rho0 * c0;
  Real ema_Rm = u0 - cs.p0 / rhoc0;
  const Real tau_ema = 3.0 * (2.0 * cs.L / c0);

  Real t = 0.0;
  while (t < t_end) {
    Real cmax = 0.0;
    for (int i = 0; i < cs.n; i++) {
      const Prim q = get_prim(f.at(i));
      cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
    }
    Real dt = cs.cfl * dx / cmax;
    dt = std::min(dt, t_end - t);
    if (dt <= 0.0) {
      break;
    }
    // The register update: ONCE per step, outside the stages, from the
    // completed state -- exactly the cadence PeleC's once-per-advance
    // update would have.  The stages below read it frozen.
    Real u_minus = 0.0;
    if ((mode == 1) || (mode == 3)) {
      const Prim qb = get_prim(f.at(0));
      // mode >= 4 flag reserved; local-impedance probe via env is clunky --
      // hardwire the experiment: LOCAL rho c instead of the frozen ambient.
      const Real rhoc_loc = qb.rho * sound_speed(qb);
      const bool local_rc = (std::getenv("NSCBC1D_LOCAL_RC") != nullptr);
      const Real rc = local_rc ? rhoc_loc : rhoc0;
      const Real Rm = qb.u - qb.p / rc;
      const Real w = dt / (tau_ema + dt);
      ema_Rm += w * (Rm - ema_Rm);
      u_minus = 0.5 * (Rm - ema_Rm);
    }
    // The RK2 stages see the target at their own stage times.  The hard-p
    // fill runs AFTER fill_bcs: the hi target is `off` there, whose
    // zero-gradient copy would otherwise overwrite the reflecting ghost.
    const bool ff = (mode == 2) || (mode == 3);
    in.u[0] = u0 + A * std::sin(om * t) + u_minus;
    in.dudt = ff ? A * om * std::cos(om * t) : 0.0;
    fill_bcs(f, in, off, prm, dx);
    hard_p_fill_hi(f, cs.p0);
    stage(f, g, dx, dt);
    in.u[0] = u0 + A * std::sin(om * (t + dt)) + u_minus;
    in.dudt = ff ? A * om * std::cos(om * (t + dt)) : 0.0;
    fill_bcs(g, in, off, prm, dx);
    hard_p_fill_hi(g, cs.p0);
    stage(g, h, dx, dt);
    for (int i = 0; i < cs.n; i++) {
      for (int v = 0; v < NVAR; v++) {
        f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
      }
    }
    t += dt;

    if (t > t_meas0) {
      const Prim qb = get_prim(f.at(0));
      const Real sn = std::sin(om * t), cn = std::cos(om * t);
      us += qb.u * sn * dt;
      uc += qb.u * cn * dt;
      const Real Rin = qb.u + qb.p / (rho0 * c0); // enters through the lo face
      rs += Rin * sn * dt;
      rc += Rin * cn * dt;
      wsum += dt;
      for (int i = 0; i < cs.n; i++) {
        const Real p = get_prim(f.at(i)).p;
        psum[i] += p * dt;
        p2sum[i] += p * p * dt;
      }
    }
  }

  const Real amp_u = 2.0 / wsum * std::sqrt(us * us + uc * uc);
  const Real amp_R = 2.0 / wsum * std::sqrt(rs * rs + rc * rc);
  r.I_u = amp_u / A;
  r.I_in = amp_R / (2.0 * A);
  for (int i = 0; i < cs.n; i++) {
    const Real pm = psum[i] / wsum;
    const Real var = std::max(p2sum[i] / wsum - pm * pm, 0.0);
    r.prms[i] = std::sqrt(var);
  }
  return r;
}

// Doppler-corrected round-trip phase and the duct's quarter-wave frequency.
Real
duct_roundtrip_phase(Real freq, Real L, Real c, Real u0)
{
  return 2.0 * (2.0 * M_PI * freq) * L * c / (c * c - u0 * u0);
}

// T1: step the inlet target, measure convergence time vs relax_u.
Real
duct_step_tconv(Real relax_u, int n, Real du, Real t_end_ta)
{
  Case cs;
  cs.n = n;
  Field f(cs.n), g(cs.n), h(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  Prim q0{};
  q0.rho = rho0;
  q0.p = cs.p0;
  q0.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    q0.Y[k] = Y[k];
  }
  const Real c0 = sound_speed(q0);
  const Real t_a = 2.0 * cs.L / c0;

  // Quiescent duct; at t = 0 the inlet target is already at du (the step).
  for (int i = -NG; i < cs.n + NG; i++) {
    set_state(f.at(i), rho0, 0.0, cs.T0, Y);
  }
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.relax_u = relax_u;
  prm.relax_t = 0.2;
  pc_nscbc::Target in, off;
  in.type = pc_nscbc::Type::inflow;
  in.u[0] = du;
  in.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    in.Y[k] = Y[k];
  }

  const Real t_end = t_end_ta * t_a;
  Real t = 0.0, t_last = 0.0;
  while (t < t_end) {
    Real cmax = 0.0;
    for (int i = 0; i < cs.n; i++) {
      const Prim q = get_prim(f.at(i));
      cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
    }
    Real dt = cs.cfl * dx / cmax;
    dt = std::min(dt, t_end - t);
    if (dt <= 0.0) {
      break;
    }
    fill_bcs(f, in, off, prm, dx);
    hard_p_fill_hi(f, cs.p0);
    stage(f, g, dx, dt);
    fill_bcs(g, in, off, prm, dx);
    hard_p_fill_hi(g, cs.p0);
    stage(g, h, dx, dt);
    for (int i = 0; i < cs.n; i++) {
      for (int v = 0; v < NVAR; v++) {
        f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
      }
    }
    t += dt;

    Real ubar = 0.0;
    for (int i = 0; i < cs.n; i++) {
      ubar += get_prim(f.at(i)).u;
    }
    ubar /= cs.n;
    if (std::abs(ubar - du) > 1.0e-3 * du) {
      t_last = t; // still outside the band: convergence not yet held
    }
  }
  return t_last / t_a;
}

void
run_t1(int n)
{
  const Real ks[] = {0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0};
  const int nk = static_cast<int>(sizeof(ks) / sizeof(ks[0]));
  std::printf(
    "\nT1  step-change convergence vs relax_u (duct, hard-p far end, "
    "n = %d)\n", n);
  std::printf(
    "    t_conv = last time |<u> - u^t| > 0.1%% of the step, in units of "
    "t_a = 2L/c\n\n");
  std::printf("    %8s %12s\n", "relax_u", "t_conv/t_a");
  Real tc[nk];
  int imin = 0;
  for (int i = 0; i < nk; i++) {
    tc[i] = duct_step_tconv(ks[i], n, 2.0e3, 40.0);
    std::printf("    %8.2f %12.2f\n", ks[i], tc[i]);
    if (tc[i] < tc[imin]) {
      imin = i;
    }
  }
  char buf[220];
  std::snprintf(
    buf, sizeof(buf), "argmin relax_u = %.2f, t_conv %.2f t_a (ends: %.2f, %.2f)",
    ks[imin], tc[imin], tc[0], tc[nk - 1]);
  check(
    (imin > 0) && (imin < nk - 1) && (tc[0] > 1.2 * tc[imin]) &&
      (tc[nk - 1] > 1.2 * tc[imin]),
    "t_conv has an interior minimum (CLR band)", buf);
  check(
    ks[imin] >= 0.1 && ks[imin] <= 1.0,
    "the optimum sits near the CLR prediction (~0.3)", buf);
}

void
run_t3_t5(bool profiles, int n, Real relax_cli)
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  Prim q0{};
  q0.rho = rho0;
  q0.p = cs.p0;
  q0.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    q0.Y[k] = Y[k];
  }
  const Real c0 = sound_speed(q0);
  const Real u0 = 2.0e3;
  // Quarter-wave frequency, Doppler-corrected: round-trip phase = pi.
  const Real f0 = (c0 * c0 - u0 * u0) / (4.0 * cs.L * c0);
  const Real fr[] = {0.8 * f0, f0, 1.2 * f0};

  if (!profiles) {
    const Real kv[] = {0.0, 0.5, 2.0, 5.0};
    std::printf(
      "\nT3  forced-inlet injection fidelity (duct, hard-p far end, "
      "f0 = %.1f Hz)\n", f0);
    std::printf(
      "    I_u = achieved u' at the boundary cell / target amplitude "
      "(1 = faithful);\n"
      "    I_in = incoming-invariant amplitude / ideal injector's; its\n"
      "    velocity-Dirichlet limit is the duct gain 1/|1+e^{i theta}| "
      "(last column,\n    divergent at f0).\n\n");
    std::printf(
      "    %8s %10s %8s %8s %10s\n", "relax_u", "f/f0", "I_u", "I_in",
      "I_in stiff");
    Real Iu_row[4] = {0.0};
    for (int m = 0; m < 4; m++) {
      for (int j = 0; j < 3; j++) {
        const DuctForcedResult r = duct_forced(kv[m], fr[j], n, 6, 4);
        const Real th = duct_roundtrip_phase(fr[j], cs.L, c0, u0);
        const Real stiff =
          1.0 /
          std::max(std::hypot(1.0 + std::cos(th), std::sin(th)), 1.0e-12);
        if (stiff < 999.0) {
          std::printf(
            "    %8.2f %10.2f %8.3f %8.3f %10.2f\n", kv[m], fr[j] / f0,
            r.I_u, r.I_in, stiff);
        } else {
          std::printf(
            "    %8.2f %10.2f %8.3f %8.3f %10s\n", kv[m], fr[j] / f0, r.I_u,
            r.I_in, "inf");
        }
        if (j == 0) { // off-resonance column used for the monotonicity gate
          Iu_row[m] = r.I_u;
        }
      }
    }
    char buf[220];
    std::snprintf(
      buf, sizeof(buf), "I_u(relax_u=0) = %.4f at f = 0.8 f0", Iu_row[0]);
    check(
      Iu_row[0] < 0.05,
      "relax_u = 0 injects nothing (no amplitude slot: I.3e)", buf);
    std::snprintf(
      buf, sizeof(buf), "I_u = %.3f -> %.3f -> %.3f for relax_u 0.5 -> 2 -> 5",
      Iu_row[1], Iu_row[2], Iu_row[3]);
    check(
      (Iu_row[1] < Iu_row[2]) && (Iu_row[2] < Iu_row[3]),
      "off-resonance injection stiffens with relax_u", buf);

    // The queue-item-4 prototypes (reported, not gated -- gates belong to
    // real implementations).  Mode 1: driver-held NRI register (Daviller's
    // I = 1 claim).  Mode 2: the stateless feed-forward amplitude slot
    // (Target::dudt).  Mode 3: both together.
    for (int md = 1; md <= 3; md++) {
      const char* names[] = {"", "NRI register", "dudt feed-forward",
                             "register + feed-forward"};
      std::printf("\n    same matrix, %s:\n", names[md]);
      Real imin = 1.0e30, imax = 0.0;
      for (int m2 = 1; m2 < 4; m2++) {
        for (int j = 0; j < 3; j++) {
          const DuctForcedResult rn = duct_forced(kv[m2], fr[j], n, 6, 4, md);
          imin = std::min(imin, rn.I_in);
          imax = std::max(imax, rn.I_in);
          std::snprintf(
            buf, sizeof(buf),
            "relax_u %.2f, f/f0 %.2f: I_in = %.3f  (I_u = %.3f)", kv[m2],
            fr[j] / f0, rn.I_in, rn.I_u);
          report(names[md], buf);
        }
      }
      if (md == 3) {
        // Phase-A driver gate: the full NRI property, I_in ~ 1 at every
        // stiffness and frequency measured, resonance included.
        std::snprintf(
          buf, sizeof(buf), "I_in in [%.3f, %.3f] over the matrix", imin,
          imax);
        check(
          (imin > 0.85) && (imax < 1.15),
          "register + feed-forward: I_in ~ 1 at any K and f", buf);
      }
    }
  } else {
    const Real k = relax_cli > 0.0 ? relax_cli : 2.0;
    std::printf(
      "\nT5  standing-wave pattern under forcing (relax_u = %.2f, "
      "f0 = %.1f Hz)\n", k, f0);
    std::printf(
      "    P_RMS(x) in dyn/cm^2 at 17 stations, against the analytic\n"
      "    envelope shape |sin(k_eff (L-x))|, k_eff = round-trip phase / "
      "2L.\n\n");
    for (int j = 0; j < 3; j++) {
      const DuctForcedResult r = duct_forced(k, fr[j], n, 6, 4);
      // Analytic envelope shape for the reflecting end: |sin(k_eff (L-x))|
      // with k_eff from the Doppler-mean wavenumber.
      const Real keff = duct_roundtrip_phase(fr[j], cs.L, c0, u0) / (2.0 * cs.L);
      Real num = 0.0, mag_m = 0.0, mag_a = 0.0, pk = 0.0;
      std::printf("    f/f0 = %.2f\n      x/L :", fr[j] / f0);
      for (int s = 0; s < 17; s++) {
        std::printf(" %6.2f", s / 16.0);
      }
      std::printf("\n      Prms:");
      for (int s = 0; s < 17; s++) {
        const int i = std::min(
          static_cast<int>((s / 16.0) * (n - 1) + 0.5), n - 1);
        std::printf(" %6.1f", r.prms[i]);
      }
      std::printf("\n      anly:");
      for (int s = 0; s < 17; s++) {
        const int i = std::min(
          static_cast<int>((s / 16.0) * (n - 1) + 0.5), n - 1);
        const Real x = (i + 0.5) * (cs.L / n);
        const Real a = std::abs(std::sin(keff * (cs.L - x)));
        const Real m = r.prms[i];
        num += a * m;
        mag_m += m * m;
        mag_a += a * a;
        pk = std::max(pk, m);
        std::printf(" %6.2f", a);
      }
      const Real corr = num / std::max(
        std::sqrt(mag_m * mag_a), 1.0e-300);
      char buf[220];
      std::snprintf(
        buf, sizeof(buf), "shape correlation %.3f, antinode P_RMS %.1f",
        corr, pk);
      if (std::abs(fr[j] / f0 - 1.0) > 0.05) {
        check(corr > 0.95, "standing-wave geometry matches the envelope", buf);
      } else {
        report("on-resonance profile (amplitude is the finding)", buf);
      }
    }
  }
}

// Phase-B driver gate: the NDNR register must collapse the reflection at
// strong anchoring while keeping the anchoring itself.  sigma = 16 is the
// measured 28%-reflection price; NDNR's claim is that the price vanishes
// because the acoustics never reach the relaxation.
void
run_ndnr()
{
  std::printf(
    "\nNDNR  EMA-mean relaxation vs classical (pulse reflection, "
    "n = 400)\n\n");
  std::printf(
    "    %8s %12s %12s %14s %14s\n", "sigma", "R_cl [%]", "R_ndnr [%]",
    "drift_cl", "drift_ndnr");
  Real Rc16 = 0.0, Rn16 = 0.0, dc16 = 0.0, dn16 = 0.0;
  for (const Real sg : {0.25, 4.0, 16.0}) {
    Real dc = 0.0, dn = 0.0;
    const Real Rc = reflection_coefficient(sg, 400, 2, false, &dc);
    const Real Rn =
      reflection_coefficient(sg, 400, 2, false, &dn, false, true);
    std::printf(
      "    %8.2f %12.3f %12.3f %14.2f %14.2f\n", sg, 100.0 * Rc, 100.0 * Rn,
      dc, dn);
    if (sg == 16.0) {
      Rc16 = Rc; Rn16 = Rn; dc16 = dc; dn16 = dn;
    }
  }
  char buf[220];
  std::snprintf(
    buf, sizeof(buf), "R 16: %.2f%% -> %.2f%%; drift %.1f -> %.1f",
    100.0 * Rc16, 100.0 * Rn16, dc16, dn16);
  check(
    Rn16 < 0.2 * Rc16,
    "NDNR collapses the sigma = 16 reflection (>5x)", buf);
  check(
    std::abs(dn16) < 3.0 * std::abs(dc16) + 5.0,
    "NDNR keeps the anchoring (drift comparable)", buf);
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, false);
  {
    // CLI: [mode] [key=value ...].  Modes: (none) = the check suite,
    // "sweep" = suite + sigma table, "t1"/"t3"/"t5" = the duct
    // injection-fidelity modes (Part III of the design doc), which run
    // alone.  Keys: n= (duct resolution), relax_u= (T5 inlet stiffness).
    const std::string mode = (argc > 1) ? argv[1] : "";
    int duct_n = 200;
    Real relax_cli = -1.0;
    for (int a = 2; a < argc; a++) {
      const std::string s(argv[a]);
      const auto eq = s.find('=');
      if (eq == std::string::npos) {
        continue;
      }
      const std::string key = s.substr(0, eq);
      const Real val = std::atof(s.c_str() + eq + 1);
      if (key == "n") {
        duct_n = static_cast<int>(val);
      } else if (key == "relax_u") {
        relax_cli = val;
      }
    }
    if (mode == "t1" || mode == "t3" || mode == "t5" || mode == "ndnr") {
      std::printf("\nnscbc1d duct mode: %s\n", mode.c_str());
      if (mode == "t1") {
        run_t1(std::min(duct_n, 200) / 2); // step response resolves fine coarse
      } else if (mode == "ndnr") {
        run_ndnr();
      } else {
        run_t3_t5(mode == "t5", duct_n, relax_cli);
      }
      std::printf("\n%d passed, %d failed\n\n", n_pass, n_fail);
      amrex::Finalize();
      return (n_fail == 0) ? 0 : 1;
    }

    const bool sweep = (mode == "sweep");
    std::printf("\nnscbc1d -- standalone verification of Source/NSCBC.H\n");
    std::printf(
      "EOS: %s,  NUM_SPECIES = %d,  NVAR = %d\n\n",
      pele::physics::PhysicsType::identifier().c_str(), NUM_SPECIES, NVAR);

    std::printf("C1  uniform-state consistency\n");
    check_uniform();
    std::printf("\nC2  relaxation directions\n");
    check_relaxation_signs();
    std::printf("\nC3  species and state identities\n");
    check_species();
#ifndef USE_SRK_EOS
    std::printf("\nC4  acoustic reflection\n");
    check_reflection(sweep);
    std::printf("\nC5  relaxation rate is a rate\n");
    check_relaxation_rate();
#else
    // The dynamic checks (C4, C5, C9b, C10, C11, C12) integrate the mini
    // solver for thousands of steps; under SRK every step pays several
    // Newton solves per cell and the full suite runs for the better part of
    // an hour to re-verify algebra that is EOS-independent and already
    // gated under Fuego.  The SRK build's purpose is the kernel's
    // EOS-portability -- every fill is an algebraic function of EOS calls --
    // which the static checks exercise completely.
    amrex::ignore_unused(sweep);
    std::printf(
      "\nC4/C5 dynamic acoustic checks: SKIPPED under SRK (EOS-independent; "
      "gated under Fuego)\n");
#endif
    std::printf("\nC6  fallbacks\n");
    check_fallbacks();
    std::printf("\nC13 outflow reversal: continuity and relaxation\n");
    check_reversal_continuity();
#ifndef USE_SRK_EOS
    std::printf("\nC14 sustained recirculation: the backflow material\n");
    check_sustained_recirculation();
#else
    std::printf(
      "\nC14 sustained recirculation: dynamic halves SKIPPED under SRK "
      "(EOS-independent; gated under Fuego)\n");
#endif
    std::printf("\nC7  reaction source term\n");
    check_reaction_source();
    std::printf("\nC8  diffusive behaviour of the ghost\n");
    check_diffusive_gradient();
    std::printf("\nC9  ghost-pressure bias from the outgoing extrapolation\n");
    check_ghost_pressure_bias();
#ifndef USE_SRK_EOS
    std::printf("\nC10 does that bias drive the solution? (source-free)\n");
    check_extrapolation_drives_solution();
    std::printf(
      "\nC11 the sustained ramp: a front on the outflow with an exact steady "
      "solution\n");
    check_sustained_ramp();
#else
    std::printf(
      "\nC10/C11/C12 dynamic checks: SKIPPED under SRK (EOS-independent; "
      "gated under Fuego)\n");
#endif

#ifndef USE_SRK_EOS
    // Conductivity for C12, boosted ~100x above air so the conductive
    // boundary error is well above every other error in the test.
    std::printf("\nC12 the diffusive gap lives in the ghost T closure\n");
    check_diffusive_dynamics(2.6e5);
#endif

    std::printf("\n%d passed, %d failed\n\n", n_pass, n_fail);
  }
  amrex::Finalize();
  return (n_fail == 0) ? 0 : 1;
}
