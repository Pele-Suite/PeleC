#ifndef NSCBC_H
#define NSCBC_H
// ============================================================================
//  NSCBC.H -- characteristic (Navier-Stokes Characteristic Boundary Condition)
//             ghost-cell fill for subsonic inflow and outflow faces.
//
//  Deliberately free of any dependency on PeleC.H so the standalone driver in
//  Verification/NSCBC1D compiles it unmodified; depends only on
//  IndexDefines.H, the PelePhysics EOS, and AMReX.
//
//  FORMULATION.  Everything is written in an outward-normal frame,
//  u_out = n_sgn * u_n (n_sgn = +1 on a hi face, -1 on a lo face), so one
//  code path serves all six faces and an outflow always has u_out > 0.  The
//  outward wave speeds are u_out - c (always incoming subsonically), u_out
//  (the lambda_0 family: entropy, tangential velocity, species, passives --
//  outgoing at an outflow, incoming at an inflow) and u_out + c (always
//  outgoing).  Hence a subsonic outflow specifies one quantity (p); a
//  subsonic inflow specifies u, T and Y.
//
//  We work with linearised invariants, the impedance rho*c FROZEN at the
//  boundary cell N (a per-cell impedance would make R differences of two
//  different variables):
//
//      R_+/- = u_out +/- p/(rho c),     S = rho - p/c^2 .
//
//  Outgoing quantities extrapolate outward on a minmod-limited slope; the
//  incoming one takes the Poinsot-Lele relaxation as a gradient increment,
//
//      L_in = K (phi - phi_target),  K = coeff (1 - M^2) c / L_ref  [1/s]
//      R_-,ghost = R_-,N + layer dx L_in / ((c - u_out) rho c) ,
//
//  so K is a mesh-independent RATE and literature coefficients transfer
//  directly.  The gradient form gives L_in an additive slot for the
//  transverse terms (beta) and the reaction source (beta_s); a hard pin of
//  R_- has none.  The fill is a pure algebraic function of the current
//  interior state: idempotent under repeated FillPatch calls,
//  bit-reproducible across restarts and decompositions.
//
//  SIGNS.  Every user-facing coefficient (sigma, relax_u, relax_t) is
//  positive; all sign handling is internal and Params::validate() rejects
//  negatives.
//
//  GPU.  Finiteness guards use amrex::Math::isfinite, never std::isfinite
//  (SYCL fast-math makes the latter constant-false, which would silently
//  degrade every fill to a zero-gradient copy).  No allocation, no virtual
//  dispatch, no printf/Abort on the device path; fallbacks are counted via
//  atomics (see Diag).  Stack footprint is order 10*NVAR doubles per thread
//  and will spill past the register budget for large mechanisms; tolerable
//  because the launch covers only a thin ghost region.
//
//  EOS.  Only density-carrying entry points are used (RTY2P, RTY2Cs, RYP2T,
//  RTY2E, REY2T, PYT2RE), which keeps the kernel real-gas capable: the
//  shortcut forms error out under SRK.
// ============================================================================

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>
#include <cmath>

#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"

namespace pc_nscbc {

// ---------------------------------------------------------------------------
//  Boundary type, decided per boundary POINT by the problem's bcnormal_nscbc().
//  A single face may carry a mixture (a jet inlet surrounded by outflow).
// ---------------------------------------------------------------------------
enum struct Type : int { off = 0, outflow = 1, inflow = 2 };

// ---------------------------------------------------------------------------
//  Numerical parameters.  POD, captured host-side into the boundary functor;
//  never read ParmParse or a class static from inside a device kernel.
// ---------------------------------------------------------------------------
struct Params
{
  // Outflow pressure relaxation, classical Poinsot-Lele sigma:
  //     K = sigma * (1 - M^2) * c / L_ref      [1/s]
  // 0 is perfectly non-reflecting and leaves the mean pressure unanchored;
  // the literature band is 0.15-0.3 (Rudy & Strikwerda 1980 give ~0.27 as
  // optimal for a 1-D pulse).
  amrex::Real sigma = 0.25;

  // Inflow normal-velocity relaxation.  Larger imposes the target velocity
  // more strongly and reflects more; >~10 is a hard Dirichlet in disguise.
  amrex::Real relax_u = 2.0;

  // Inflow temperature and tangential-velocity relaxation.
  amrex::Real relax_t = 0.2;

  // Boundary-normal reference length, probhi[idir] - problo[idir].  Filled per
  // face by the caller.  Only the ratio sigma/L_ref is physical, so this is
  // not a user knob; it exists so that sigma keeps its literature meaning.
  amrex::Real L_ref = 1.0;

  // Order of the outgoing-invariant extrapolation:
  //   1 = zeroth-order copy (constant across ghost layers)
  //   2 = minmod-limited linear in the ghost layer index (default)
  // Exposed so that a bit-comparison against the standalone 1-D driver can be
  // done without recompiling.  Verification knob, not a physics knob.
  int order = 2;

  // Replace the relaxation by a hard pin of R_- to the far-field value
  // (u_out = 0, p = p_target): non-reflecting AND anchored, but it anchors to
  // p_target + rho c u_out, and as a value constraint its effective rate is
  // c/dx -- it does not converge under refinement.  For a boundary onto a
  // quiescent reservoir; wrong for a through-flow duct.  sigma is ignored.
  bool pin_farfield = false;

  AMREX_GPU_HOST bool validate(std::string& why) const
  {
    if (sigma < 0.0) {
      why = "bc_nscbc_sigma must be >= 0 (all NSCBC relaxation coefficients "
            "are positive in PeleC; internal signs are handled by the kernel)";
      return false;
    }
    if (relax_u < 0.0) {
      why = "bc_nscbc_relax_u must be >= 0";
      return false;
    }
    if (relax_t < 0.0) {
      why = "bc_nscbc_relax_t must be >= 0";
      return false;
    }
    if ((order != 1) && (order != 2)) {
      why = "bc_nscbc_order must be 1 or 2";
      return false;
    }
    if (L_ref <= 0.0) {
      why = "NSCBC reference length must be > 0";
      return false;
    }
    return true;
  }
};

// ---------------------------------------------------------------------------
//  Target state returned per boundary point by the problem hook.  Only the
//  entries relevant to `type` are read: outflow -> p; inflow -> u[], T, Y[].
//  Deliberately no relaxation for composition: at an inflow every species
//  characteristic is incoming (hard imposition is well posed), at an outflow
//  all are outgoing (imposing would over-specify).
// ---------------------------------------------------------------------------
struct Target
{
  Type type = Type::off;
  amrex::Real p = 0.0;                // outflow: target pressure [dyn/cm^2]
  amrex::Real u[3] = {0.0, 0.0, 0.0}; // inflow: target velocity [cm/s]
  amrex::Real T = 0.0;                // inflow: target temperature [K]
  amrex::Real Y[NUM_SPECIES] = {0.0}; // inflow: target mass fractions

  // Inflow acoustic-injection feed-forward: the time derivative of the
  // NORMAL velocity target, supplied analytically by the problem hook when
  // it is forcing the inlet.  The value-relaxation alone has no amplitude
  // slot (design doc I.3e) -- an injected signal rides the same relaxation
  // that holds the mean, with the measured deterioration of driver t3 --
  // and a register-based NRI transplant does not create one (measured:
  // subtracting the outgoing wave from the deviation only exposes the bare
  // injection bandwidth limit).  This term is that slot, stateless: the
  // incoming amplitude gains 2 rho c du^t/dt, the GC transplant of the
  // flux form's acoustic forcing, and the relaxation keeps its one job of
  // holding the mean.  Zero (default) recovers the classical inlet
  // identically.
  amrex::Real dudt = 0.0; // d(u_normal target)/dt, LAB frame [cm/s^2]
};

// ---------------------------------------------------------------------------
//  Diagnostic counters.  Every fallback is counted; a silent fallback is a
//  bug that will not be found.  Storage is amrex::Long, not int: with the
//  default sum_interval = -1 the counters are never reset, and the supersonic
//  counter on a large 3-D run walks past 2^31 in ~1e4 steps.
// ---------------------------------------------------------------------------
namespace Diag {
enum : int {
  supersonic = 0,    // |M| >= 1: zero-gradient copy is exact, not a failure
  reversed,          // flow ran the wrong way through the face
  body_state,        // an EB-covered cell was found in the stencil
  eos_failure,       // an EOS call returned a non-finite or non-positive value
  floored,           // a positivity floor was hit
  transverse_drop,   // a transverse derivative was unavailable or unusable
  source_drop,       // a modelled source term was unavailable or unusable
  target_incomplete, // a supersonic inflow whose Target carries no pressure:
                     // every characteristic is incoming, so the FULL state is
                     // required; the interior pressure was substituted, which
                     // makes the ghost depend on a domain it should be
                     // causally upstream of.  Set Target.p in bcnormal_nscbc.
  structure,         // ADVISORY, not a fallback: material structure sat in an
                     // outflow boundary cell (|dS| > 5% of rho).  A front at
                     // the boundary wants the flame closures --
                     // extrap_temperature and beta_s = 0; see
                     // NSCBC-FlameOutflow/README.md.
  count
};
}

// ===========================================================================
//  Small helpers
// ===========================================================================

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
finite_or(const amrex::Real x, const amrex::Real fallback) noexcept
{
  return amrex::Math::isfinite(x) ? x : fallback;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
minmod(const amrex::Real a, const amrex::Real b) noexcept
{
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (std::abs(a) < std::abs(b)) ? a : b;
}

// True if a state looks like EB body state rather than fluid, detected by
// non-positive density.  This is why read_params makes eb_zero_body_state = 1
// mandatory under EB: PeleC's default body state is a sampled fluid state the
// stencil cannot distinguish from interior data.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
is_body_state(const amrex::Real s[NVAR]) noexcept
{
  return !(s[URHO] > constants::very_small_num()) ||
         !amrex::Math::isfinite(s[URHO]) || !amrex::Math::isfinite(s[UEDEN]);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
copy_state(const amrex::Real src[NVAR], amrex::Real dst[NVAR]) noexcept
{
  for (int n = 0; n < NVAR; n++) {
    dst[n] = src[n];
  }
}

// ---------------------------------------------------------------------------
//  Primitive decomposition of one stencil cell.  T is taken from the state
//  (a maintained component; computeTemp runs after every update), keeping
//  every EOS call on this path non-iterative -- which matters under SRK --
//  with the iterative recovery as fallback.
// ---------------------------------------------------------------------------
struct CellPrim
{
  amrex::Real rho, u[3], T, p, c, e;
  amrex::Real Y[NUM_SPECIES];
  bool ok;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE CellPrim
cell_primitives(const amrex::Real s[NVAR]) noexcept
{
  auto eos = pele::physics::PhysicsType::eos();
  CellPrim q;
  q.ok = true;

  q.rho = s[URHO];
  const amrex::Real rhoinv = 1.0 / q.rho;
  q.u[0] = s[UMX] * rhoinv;
  q.u[1] = s[UMY] * rhoinv;
  q.u[2] = s[UMZ] * rhoinv;

  amrex::Real Ysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    q.Y[n] = amrex::max<amrex::Real>(s[UFS + n] * rhoinv, 0.0);
    Ysum += q.Y[n];
  }
  if (Ysum > constants::very_small_num()) {
    const amrex::Real inv = 1.0 / Ysum;
    for (amrex::Real& yv : q.Y) {
      yv *= inv;
    }
  } else {
    q.ok = false;
  }

  const amrex::Real ke =
    0.5 * (q.u[0] * q.u[0] + q.u[1] * q.u[1] + q.u[2] * q.u[2]);
  q.e = s[UEDEN] * rhoinv - ke;

  q.T = s[UTEMP];
  if (!(q.T > 0.0) || !amrex::Math::isfinite(q.T)) {
    q.T = 300.0; // initial guess for the iterative fallback
    eos.REY2T(q.rho, q.e, q.Y, q.T);
  }

  eos.RTY2P(q.rho, q.T, q.Y, q.p);
  eos.RTY2Cs(q.rho, q.T, q.Y, q.c);

  if (
    !amrex::Math::isfinite(q.p) || !(q.p > 0.0) ||
    !amrex::Math::isfinite(q.c) || !(q.c > 0.0) ||
    !amrex::Math::isfinite(q.T) || !(q.T > 0.0)) {
    q.ok = false;
  }
  return q;
}

// ---------------------------------------------------------------------------
//  Pack a set of ghost primitives into a conserved state.
//
//  Every NVAR slot is written, and the algebraic identities that the rest of
//  the code assumes are enforced exactly rather than approximately:
//      sum_k (rho Y)_k == rho          (species renormalised before packing)
//      UEDEN == UEINT + 0.5 rho |u|^2
//  The passive families (advected, auxiliary, linear/soot) ride the lambda_0
//  characteristic and are handled by the caller through `pass`.
// ---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
pack_ghost(
  const amrex::Real rho,
  const amrex::Real u[3],
  const amrex::Real e,
  const amrex::Real T,
  const amrex::Real Y[NUM_SPECIES],
  const amrex::Real s_N[NVAR],
  amrex::Real s_ghost[NVAR]) noexcept
{
  // s_N feeds only the passive families below; in a build with no passive
  // fields every #if block is empty and the parameter would warn as unused.
  amrex::ignore_unused(s_N);
  s_ghost[URHO] = rho;
  s_ghost[UMX] = rho * u[0];
  s_ghost[UMY] = rho * u[1];
  s_ghost[UMZ] = rho * u[2];
  s_ghost[UEINT] = rho * e;
  s_ghost[UEDEN] =
    rho * e + 0.5 * rho * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  s_ghost[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ghost[UFS + n] = rho * Y[n];
  }
  // Passive families default to a mass-scaled zero-gradient; the caller
  // overwrites them where a target applies.
#if NUM_ADV > 0
  const amrex::Real ratio =
    (s_N[URHO] > constants::very_small_num()) ? rho / s_N[URHO] : 1.0;
  for (int n = 0; n < NUM_ADV; n++) {
    s_ghost[UFA + n] = s_N[UFA + n] * ratio; // per unit mass in ctoprim
  }
#endif
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    s_ghost[UFX + n] = s_N[UFX + n]; // NOT per unit mass in ctoprim
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    s_ghost[ULIN + n] = s_N[ULIN + n]; // NOT per unit mass in ctoprim
  }
#endif
}

// ===========================================================================
//  apply() -- fill one ghost cell.
//    s_N, s_Nm1, s_Nm2: interior stencil walking inward; tangential indexing
//      is the caller's (BCfill.cpp), keeping the fill a deterministic pure
//      function of data written before the launch.
//    n_stencil: usable cells (3/2/1), reduced further on EB body state.
//    idir, sgn: boundary normal and side (+1 lo, -1 hi, as bcnormal).
//    layer: 1 = nearest the domain.  diag: optional Long[Diag::count].
// ===========================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
apply(
  const amrex::Real s_N[NVAR],
  const amrex::Real s_Nm1[NVAR],
  const amrex::Real s_Nm2[NVAR],
  int n_stencil,
  const amrex::Real dx_normal,
  const int idir,
  const int sgn,
  const int layer,
  Target const& tgt,
  Params const& prm,
  amrex::Real s_ghost[NVAR],
  amrex::Long* diag = nullptr) noexcept
{
  auto eos = pele::physics::PhysicsType::eos();

  auto bump = [=](const int which) {
    if (diag != nullptr) {
      amrex::HostDevice::Atomic::Add(&diag[which], amrex::Long(1));
    }
  };

  // ---- 0. Stencil hygiene ------------------------------------------------
  if (is_body_state(s_N)) {
    // Nothing usable to build a characteristic state from.  Leave whatever the
    // caller pre-populated (the reflected NoSlipWall state) in place.
    bump(Diag::body_state);
    copy_state(s_N, s_ghost);
    return;
  }
  if ((n_stencil > 2) && is_body_state(s_Nm2)) {
    n_stencil = 2;
    bump(Diag::body_state);
  }
  if ((n_stencil > 1) && is_body_state(s_Nm1)) {
    n_stencil = 1;
    bump(Diag::body_state);
  }
  // Extrapolation order can never exceed what the stencil supports.  A
  // minmod-limited slope needs three cells; with two it would be unlimited.
  const int order =
    amrex::min<int>(prm.order, amrex::max<int>(n_stencil - 1, 1));
  const bool linear = (order >= 2) && (n_stencil >= 3);

  const CellPrim qN = cell_primitives(s_N);
  if (!qN.ok) {
    bump(Diag::eos_failure);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 1. Outward-normal frame ------------------------------------------
  const auto n_sgn = static_cast<amrex::Real>(-sgn);
  const amrex::Real u_out_N = n_sgn * qN.u[idir];
  const amrex::Real c_N = qN.c;
  const amrex::Real mach = u_out_N / c_N;

  // Frozen acoustic impedance (see the header note: this must be a single
  // value shared by every stencil cell).
  const amrex::Real rho_c = qN.rho * c_N;

  // ---- 2. Regimes with an exact answer ----------------------------------
  if (std::abs(mach) >= 1.0) {
    if (mach >= 1.0) {
      // Supersonic outflow: zero-gradient is exact, not an approximation.
      bump(Diag::supersonic);
      copy_state(s_N, s_ghost);
      return;
    }
    // Supersonic inflow: every characteristic is incoming, so the full state
    // is specified.  Only meaningful with an inflow target.  A Target
    // without a pressure under-specifies it: the interior pressure is
    // substituted and COUNTED (target_incomplete), because the ghost then
    // depends on a domain it should be causally upstream of.
    if (tgt.type == Type::inflow) {
      if (!(tgt.p > 0.0)) {
        bump(Diag::target_incomplete);
      }
      amrex::Real rho_g = 0.0, e_g = 0.0;
      eos.PYT2RE(tgt.p > 0.0 ? tgt.p : qN.p, tgt.Y, tgt.T, rho_g, e_g);
      if (
        amrex::Math::isfinite(rho_g) && (rho_g > 0.0) &&
        amrex::Math::isfinite(e_g)) {
        amrex::Real u_g[3] = {tgt.u[0], tgt.u[1], tgt.u[2]};
        pack_ghost(rho_g, u_g, e_g, tgt.T, tgt.Y, s_N, s_ghost);
        bump(Diag::supersonic);
        return;
      }
      bump(Diag::eos_failure);
    }
    bump(Diag::supersonic);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 3. Direction sanity ----------------------------------------------
  const bool outflow_face = (tgt.type == Type::outflow);
  if (outflow_face && (u_out_N < 0.0)) {
    // Transient reversal at an outflow: COUNTED, but no longer a separate
    // closure.  The acoustic side of the standard outflow path is well
    // defined and restoring for u_out < 0 (lambda_in = c - u_out grows,
    // 1 - M^2 is even in M, K (p - p_tgt) pulls the ghost pressure toward
    // the target AND pushes the ghost normal velocity outward), so one
    // closure runs on both sides of u_out = 0; only the lambda_0 material
    // SLOPES are upwinded off during reversal (see w_mat below), leaving the
    // ghost the interior's material VALUES.  Both dedicated branches tried
    // here were defects the chamber found: a hard pin of tgt.p is a step
    // discontinuity that NaN'd under vent breathing (backflow -144 -> -2628
    // cm/s in 45 steps, a 138 K cell), and the soft pressure-only relaxation
    // that replaced it dropped dR+, S_p and T_in relative to the forward
    // branch and froze the ghost velocity -- an O(1) discontinuity at
    // u_out = 0 that a flame finishing its transit turned into a growing
    // dither and a spurious 0.3 atm chamber spike
    // (Docs/NSCBC-reversal-branch-defect.md; gates: driver C13 a/b/c).
    // Sustained recirculation as a local inflow is Params::backflow_material.
    //
    // The COUNTER carries a roundoff deadband (the closure itself has no
    // branch left to care): a face holding a quiescent charge sits at
    // exactly u_out = 0 and the solver's noise dithers the sign, which read
    // as 314k "reversals" in the first 100 steps of a chamber whose ignition
    // wave was still a centimetre from the boundary.  Counting only
    // |u_out| > 1e-9 c keeps "reversal count zero" meaning NO PHYSICAL
    // BACKFLOW rather than "no roundoff", which is the semantics the
    // README's advice ("sustained nonzero means a misplaced outflow")
    // depends on.
    if (u_out_N < -1.0e-9 * c_N) {
      bump(Diag::reversed);
    }
  }
  if (!outflow_face && (u_out_N > 0.0)) {
    // Outflow through a configured inflow would impose D+NUM_SPECIES+1
    // conditions where one is admitted; zero-gradient and count it.
    bump(Diag::reversed);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 4. Outgoing acoustic invariant, R_+ ------------------------------
  // R_+ leaves the domain at both subsonic inflow and subsonic outflow
  // (lambda_+ = u_out + c > 0 in both cases), so it is always extrapolated.
  const amrex::Real Rp_N = u_out_N + qN.p / rho_c;
  amrex::Real dRp = 0.0;
  CellPrim qNm1, qNm2;
  if (linear) {
    qNm1 = cell_primitives(s_Nm1);
    qNm2 = cell_primitives(s_Nm2);
    if (qNm1.ok && qNm2.ok) {
      const amrex::Real Rp_Nm1 = n_sgn * qNm1.u[idir] + qNm1.p / rho_c;
      const amrex::Real Rp_Nm2 = n_sgn * qNm2.u[idir] + qNm2.p / rho_c;
      dRp = minmod(Rp_N - Rp_Nm1, Rp_Nm1 - Rp_Nm2);
    } else {
      bump(Diag::eos_failure);
    }
  }
  const amrex::Real Rp_g = Rp_N + static_cast<amrex::Real>(layer) * dRp;

  // ---- 5. Incoming acoustic invariant, R_- ------------------------------
  // lambda_- < 0 always (subsonic), so R_- is always modelled.  Sign checks:
  // outflow, p_N > p_target => dR_-/dn > 0 and (dp/dR_- = -rho c/2) the
  // ghost pressure falls toward target; inflow, u above target => smaller
  // R_-,ghost and the ghost normal velocity falls toward target.
  const amrex::Real Rm_N = u_out_N - qN.p / rho_c;
  const amrex::Real lambda_in = c_N - u_out_N; // > 0, magnitude of the speed
  const amrex::Real one_m_M2 =
    amrex::max<amrex::Real>(1.0 - mach * mach, constants::very_small_num());

  // Upwinding of the lambda_0 (material) family.  Its characteristics leave
  // the domain only for u_out > 0; during a reversal an EXTRAPOLATED
  // material ghost is advected back in and feeds on itself.  That loop is
  // the chamber's cold runaway: the outward-cooling T ramp refrigerates the
  // boundary cell, the steepened ramp extrapolates colder still, and
  // 1/(rho c) amplifies the relaxation as T falls -- 385 -> 241 -> 89 K in
  // 42 us, then NaN (Docs/NSCBC-reversal-branch-defect.md).  So the material
  // SLOPES are used only where their upwind justification holds: full for
  // u_out >= 0 (every forward-flow result is untouched to the bit), ramped
  // off across a narrow band below zero (1e-3 c -- a continuity
  // regularisation of the fill, not a physics knob; C13 gates it), and zero
  // under firm reversal, where the ghost keeps the interior VALUES -- which
  // is what the surviving closure always did.
  const amrex::Real w_mat =
    (u_out_N >= 0.0)
      ? 1.0
      : amrex::max<amrex::Real>(0.0, 1.0 + u_out_N / (1.0e-3 * c_N));

  const amrex::Real dRm = 0.0;

  amrex::Real Rm_g;
  if (outflow_face && prm.pin_farfield) {
    // Hard far-field pin: R_- takes its far-field value (u_out = 0,
    // p = p_target).  Value constraint, not a rate, and there is nowhere in it
    // to put a transverse term -- so pin_farfield ignores beta.
    Rm_g = -tgt.p / rho_c;
  } else {
    amrex::Real L_in;
    if (outflow_face) {
      const amrex::Real K = prm.sigma * one_m_M2 * c_N / prm.L_ref;
      L_in = K * (qN.p - tgt.p);
    } else {
      const amrex::Real u_out_tgt = n_sgn * tgt.u[idir];
      const amrex::Real K = prm.relax_u * one_m_M2 * c_N / prm.L_ref;
      L_in = -K * qN.rho * c_N * (u_out_N - u_out_tgt);
      // Acoustic-injection feed-forward (Target::dudt): the incoming wave
      // carries the target's rate directly, so a forced signal does not
      // have to fight its way through the relaxation.  n_sgn converts the
      // lab-frame rate to the outward frame; the sign follows the
      // relaxation term's convention (positive L_in raises u_out through
      // the Rm inversion), and the factor 2 is u = (R+ + R-)/2: imposing a
      // velocity rate through the incoming invariant alone needs twice it.
      L_in += 2.0 * rho_c * (n_sgn * tgt.dudt);
    }
    // There is deliberately NO diffusive source term on L_in.  The
    // ghost-cell form does not have the flux-form's viscous-condition gap:
    // the diffusion operator reads these ghost cells, so a correct ghost
    // closure (extrap_temperature, gated by C8 and C12) already carries the
    // diffusive physics, and an amplitude-side term counts it twice
    // (measured: it moves the C12 driver error 104 -> -911).  Re-verify
    // against C12 before revisiting.
    Rm_g =
      Rm_N + static_cast<amrex::Real>(layer) * dRm +
      static_cast<amrex::Real>(layer) * dx_normal * L_in / (lambda_in * rho_c);
  }

  // ---- 6. Invert the acoustic pair --------------------------------------
  amrex::Real u_out_g = 0.5 * (Rp_g + Rm_g);
  amrex::Real p_g = 0.5 * rho_c * (Rp_g - Rm_g);

  const amrex::Real p_floor = amrex::max<amrex::Real>(
    1.0e-3 * (tgt.p > 0.0 ? tgt.p : qN.p), constants::very_small_num());
  if (!amrex::Math::isfinite(p_g) || (p_g < p_floor)) {
    p_g = p_floor;
    bump(Diag::floored);
  }
  if (!amrex::Math::isfinite(u_out_g)) {
    u_out_g = u_out_N;
    bump(Diag::floored);
  }

  amrex::Real u_g[3] = {qN.u[0], qN.u[1], qN.u[2]};
  u_g[idir] = n_sgn * u_out_g;

  // ---- 7. The lambda_0 family -------------------------------------------
  // Entropy, tangential velocities, species and passive scalars all convect at
  // u_out.  At an OUTFLOW they leave the domain and are extrapolated; at an
  // INFLOW they enter and are set from the target -- composition hard (all
  // NUM_SPECIES characteristics are incoming, so it is well posed), and
  // temperature and tangential velocity through a relaxation.
  amrex::Real rho_g = 0.0, T_g = 0.0, e_g = 0.0;
  amrex::Real Y_g[NUM_SPECIES];

  if (outflow_face) {
    // Entropy invariant extrapolated like R_+; zero slope reduces to the
    // linearised isentrope rho_g = rho_N + (p_g - p_N)/c^2.
    const amrex::Real inv_c2 = 1.0 / (c_N * c_N);
    const amrex::Real S_N = qN.rho - qN.p * inv_c2;
    amrex::Real dS = 0.0;
    amrex::Real dY[NUM_SPECIES] = {0.0};
    amrex::Real dut[3] = {0.0, 0.0, 0.0};
    if (linear && qNm1.ok && qNm2.ok) {
      // All lambda_0 slopes carry the material upwinding factor: extrapolated
      // material content is only valid where it advects OUT (see w_mat).
      const amrex::Real S_Nm1 = qNm1.rho - qNm1.p * inv_c2;
      const amrex::Real S_Nm2 = qNm2.rho - qNm2.p * inv_c2;
      dS = w_mat * minmod(S_N - S_Nm1, S_Nm1 - S_Nm2);
      for (int n = 0; n < NUM_SPECIES; n++) {
        dY[n] = w_mat * minmod(qN.Y[n] - qNm1.Y[n], qNm1.Y[n] - qNm2.Y[n]);
      }
      for (int d = 0; d < 3; d++) {
        if (d != idir) {
          dut[d] = w_mat * minmod(qN.u[d] - qNm1.u[d], qNm1.u[d] - qNm2.u[d]);
        }
      }
    }
    const auto fl = static_cast<amrex::Real>(layer);
    const amrex::Real S_g = S_N + fl * dS;
    rho_g = S_g + p_g * inv_c2;

    amrex::Real Ysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_g[n] = amrex::max<amrex::Real>(qN.Y[n] + fl * dY[n], 0.0);
      Ysum += Y_g[n];
    }
    if (Ysum > constants::very_small_num()) {
      const amrex::Real inv = 1.0 / Ysum;
      for (amrex::Real& yg : Y_g) {
        yg *= inv;
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_g[n] = qN.Y[n];
      }
      bump(Diag::floored);
    }

    for (int d = 0; d < 3; d++) {
      if (d != idir) {
        u_g[d] = qN.u[d] + fl * dut[d];
      }
    }

    const amrex::Real rho_floor =
      amrex::max<amrex::Real>(1.0e-6 * qN.rho, constants::very_small_num());
    if (!amrex::Math::isfinite(rho_g) || (rho_g < rho_floor)) {
      rho_g = rho_floor;
      bump(Diag::floored);
    }
    eos.RYP2T(rho_g, Y_g, p_g, T_g);
    eos.RTY2E(rho_g, T_g, Y_g, e_g);

  } else {
    // Inflow.  One dimensionless nudge factor shared by temperature and the
    // tangential velocities -- both ride lambda_0, magnitude |u_out|:
    //     nu = layer dx relax_t c / (L_ref |u_out|), clamped to [0, 1],
    // applied as an explicit move toward the target so no sign convention
    // can go wrong.
    const amrex::Real u_mag =
      amrex::max<amrex::Real>(std::abs(u_out_N), constants::smallu());
    const amrex::Real nu = amrex::min<amrex::Real>(
      static_cast<amrex::Real>(layer) * dx_normal * prm.relax_t * c_N /
        (prm.L_ref * u_mag),
      1.0);

    T_g = qN.T - nu * (qN.T - tgt.T);
    for (int d = 0; d < 3; d++) {
      if (d != idir) {
        u_g[d] = qN.u[d] - nu * (qN.u[d] - tgt.u[d]);
      }
    }
    // Composition: hard.  Every species characteristic is incoming.
    amrex::Real Ysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_g[n] = amrex::max<amrex::Real>(tgt.Y[n], 0.0);
      Ysum += Y_g[n];
    }
    if (Ysum > constants::very_small_num()) {
      const amrex::Real inv = 1.0 / Ysum;
      for (amrex::Real& yg : Y_g) {
        yg *= inv;
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_g[n] = qN.Y[n];
      }
      bump(Diag::floored);
    }
    if (!amrex::Math::isfinite(T_g) || !(T_g > 0.0)) {
      T_g = qN.T;
      bump(Diag::floored);
    }
    // Close on (p, T, Y): p from the acoustic pair, T and Y from the target.
    eos.PYT2RE(p_g, Y_g, T_g, rho_g, e_g);
  }

  if (
    !amrex::Math::isfinite(rho_g) || !(rho_g > 0.0) ||
    !amrex::Math::isfinite(e_g) || !amrex::Math::isfinite(T_g) ||
    !(T_g > 0.0)) {
    bump(Diag::eos_failure);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 8. Pack and sanitise ---------------------------------------------
  pack_ghost(rho_g, u_g, e_g, T_g, Y_g, s_N, s_ghost);

  for (int n = 0; n < NVAR; n++) {
    if (!amrex::Math::isfinite(s_ghost[n])) {
      s_ghost[n] = s_N[n];
      bump(Diag::floored);
    }
  }
}

} // namespace pc_nscbc

#endif
