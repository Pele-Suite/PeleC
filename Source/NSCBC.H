#ifndef NSCBC_H
#define NSCBC_H
// ============================================================================
//  NSCBC.H -- characteristic (Navier-Stokes Characteristic Boundary Condition)
//             ghost-cell fill for subsonic inflow and outflow faces.
//
//  Deliberately free of any dependency on PeleC.H, so that a standalone 1-D
//  driver (a later PR in this series) can compile it unmodified; depends
//  only on IndexDefines.H, the PelePhysics EOS, and AMReX.
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

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>
#include <cmath>
#include <limits>
#include <string>

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
  // Exposed so that a bit-comparison against a standalone 1-D driver can be
  // done without recompiling.  Verification knob, not a physics knob.
  int order = 2;

  // Replace the relaxation by a hard pin of R_- to the far-field value
  // (u_out = 0, p = p_target): non-reflecting AND anchored, but it anchors to
  // p_target + rho c u_out, and as a value constraint its effective rate is
  // c/dx -- it does not converge under refinement.  For a boundary onto a
  // quiescent reservoir; wrong for a through-flow duct.  sigma is ignored.
  bool pin_farfield = false;

  // Every coefficient must be finite as well as non-negative: a `< 0` test
  // alone lets NaN and inf through (ParmParse accepts both), and an infinite
  // sigma or relax_u drives every ghost onto the pressure floor.
  AMREX_GPU_HOST bool validate(std::string& why) const
  {
    auto bad = [](const amrex::Real x) {
      return !(x >= 0.0) || !amrex::Math::isfinite(x);
    };
    if (bad(sigma)) {
      why = "bc_nscbc_sigma must be finite and >= 0 (all NSCBC relaxation "
            "coefficients are positive in PeleC; internal signs are handled by "
            "the kernel; for a hard pressure outlet use "
            "bc_nscbc_pin_farfield)";
      return false;
    }
    if (bad(relax_u)) {
      why = "bc_nscbc_relax_u must be finite and >= 0";
      return false;
    }
    if (bad(relax_t)) {
      why = "bc_nscbc_relax_t must be finite and >= 0";
      return false;
    }
    if ((order != 1) && (order != 2)) {
      why = "bc_nscbc_order must be 1 or 2";
      return false;
    }
    if (!(L_ref > 0.0) || !amrex::Math::isfinite(L_ref)) {
      why = "NSCBC reference length must be finite and > 0";
      return false;
    }
    return true;
  }
};

// ---------------------------------------------------------------------------
//  Target state returned per boundary point by the problem hook.  Only the
//  entries relevant to `type` are read: outflow -> p; inflow -> u[], T, Y[]
//  (and p where the inflow is supersonic; dudt and the passives are
//  optional).  A missing p, T or Y is counted (Diag::target_incomplete) and
//  the interior value stands in for it.
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
  // slot: an injected signal rides the same relaxation that holds the mean
  // and is attenuated by it (measured on the forced-inlet case of a later
  // PR in this series), and subtracting the outgoing wave from the deviation
  // does not create one either (it only exposes the bare injection
  // bandwidth limit).  This term is that slot, stateless: the
  // incoming amplitude gains 2 rho c du^t/dt, the GC transplant of the
  // flux form's acoustic forcing, and the relaxation keeps its one job of
  // holding the mean.  Zero (default) recovers the classical inlet
  // identically.
  amrex::Real dudt = 0.0; // d(u_normal target)/dt, LAB frame [cm/s^2]

  // Inflow: passive scalars entering with the stream.  Their characteristics
  // are incoming at an inflow, so like Y they must be specified: set
  // `passives = true` and fill adv (per unit mass, UFA/rho), aux and lin (as
  // stored in UFX and ULIN).  Left unset, the ghost keeps the interior values
  // and the fill counts passive_unset.  Ignored at an outflow.
  bool passives = false;
#if NUM_ADV > 0
  amrex::Real adv[NUM_ADV] = {0.0};
#endif
#if NUM_AUX > 0
  amrex::Real aux[NUM_AUX] = {0.0};
#endif
#if NUM_LIN > 0
  amrex::Real lin[NUM_LIN] = {0.0};
#endif
};

// ---------------------------------------------------------------------------
//  Diagnostic counters.  Every fallback is counted; a silent fallback is a
//  bug that will not be found.  Storage is amrex::Long, not int: the
//  counters are reset only when reported (at least every 100 coarse steps,
//  and once more at the end of the run), and the supersonic counter on a
//  large 3-D run can pass 2^31 between two reports.
//
//  The counters count FILL EVENTS, not distinct boundary cells: a ghost cell
//  is counted every time a FillPatch fills it, including the fills done for
//  tagging, derived quantities and plotfiles, and a periodic-seam ghost is
//  counted once on each FAB that holds it.  Totals therefore scale with the
//  number of fills; what carries meaning is whether a counter is zero, and
//  whether it stays nonzero report after report.
// ---------------------------------------------------------------------------
namespace Diag {
enum : int {
  supersonic = 0,    // |M| >= 1: zero-gradient copy is exact, not a failure
  reversed,          // flow ran the wrong way through the face (subsonic or
                     // supersonic; beyond a roundoff deadband)
  body_state,        // an EB-covered cell was found in the stencil
  eos_failure,       // an EOS call returned a non-finite or non-positive value
  floored,           // a positivity floor was hit
  stencil_short,     // the FAB held fewer interior cells along the normal
                     // than the extrapolation order wants (the coarse patch
                     // under a fine grid on the face): the order dropped to
                     // what fits, and the ghost depends on the box layout.
  target_incomplete, // the Target lacks a value the boundary type needs.
                     // Outflow without p: nothing to relax toward, so the
                     // interior pressure was substituted and the ghost
                     // pressure is not anchored at all.  Inflow without T
                     // (T <= 0) or without a composition (Y summing to
                     // zero): the interior value was substituted, so nothing
                     // holds the inlet temperature/composition.  Supersonic
                     // inflow without p: every characteristic is incoming, so
                     // the FULL state is required; the interior pressure was
                     // substituted, which makes the ghost depend on a domain
                     // it should be causally upstream of.  Set the missing
                     // field in bcnormal_nscbc.
  passive_unset,     // an inflow Target left its passive scalars unset
                     // (Target::passives): the interior values were kept, so
                     // nothing about them enters through the inflow.
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

// Roundoff deadband of the flow-reversal COUNTER, as a fraction of c.  1e-9
// in double precision; in single precision one ulp of pressure across a cell
// already drives |u| ~ 1e-8 c to 1e-6 c, so the band is widened with the
// precision (100 epsilon) or a quiescent face would count roundoff as
// reversals.  Only the counter reads this; the closure itself is unaffected.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr amrex::Real
reversal_deadband() noexcept
{
  constexpr amrex::Real eps100 =
    amrex::Real(100.0) * std::numeric_limits<amrex::Real>::epsilon();
  return (eps100 > amrex::Real(1.0e-9)) ? eps100 : amrex::Real(1.0e-9);
}

// True if a state cannot be fluid: non-positive or non-finite density or
// energy (the zero EB body state of eb_zero_body_state = 1, say).  This is a
// second line of defence only: under EB the caller (BCfill.cpp) identifies
// covered stencil cells from the EB geometry itself, since PeleC's DEFAULT
// body state is a sampled fluid state that no test on the data can tell from
// interior data.
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

// Clamp a target composition at zero and renormalise it to unit sum, as
// cell_primitives does for an interior state, so that a packed ghost keeps
// sum_k (rho Y)_k == rho whatever the problem hook returned.  False, with
// Y_out untouched, if the target carries no mass at all.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
normalize_Y(
  const amrex::Real Y_in[NUM_SPECIES], amrex::Real Y_out[NUM_SPECIES]) noexcept
{
  amrex::Real Ysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    Ysum += amrex::max<amrex::Real>(Y_in[n], 0.0);
  }
  if (!(Ysum > constants::very_small_num())) {
    return false;
  }
  const amrex::Real inv = 1.0 / Ysum;
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y_out[n] = amrex::max<amrex::Real>(Y_in[n], 0.0) * inv;
  }
  return true;
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
  // Passive families default to a mass-scaled zero-gradient; at an inflow
  // apply() then blends them toward the Target's (inflow_passives).
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

// ---------------------------------------------------------------------------
//  Passive scalars at an inflow: blend the packed (interior) values toward the
//  Target's with weight w (1 = full inflow; less across the reversal band).
// ---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
inflow_passives(
  Target const& tgt, const amrex::Real w, amrex::Real s_ghost[NVAR]) noexcept
{
  amrex::ignore_unused(tgt, w, s_ghost);
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    s_ghost[UFA + n] =
      (1.0 - w) * s_ghost[UFA + n] + w * s_ghost[URHO] * tgt.adv[n];
  }
#endif
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    s_ghost[UFX + n] = (1.0 - w) * s_ghost[UFX + n] + w * tgt.aux[n];
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    s_ghost[ULIN + n] = (1.0 - w) * s_ghost[ULIN + n] + w * tgt.lin[n];
  }
#endif
}

// ---------------------------------------------------------------------------
//  Last-resort sanitiser shared by every path that packs a ghost: a
//  non-finite slot takes the boundary cell's value, counted as floored.
// ---------------------------------------------------------------------------
template <typename Bump>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sanitize_ghost(
  const amrex::Real s_N[NVAR], amrex::Real s_ghost[NVAR], Bump const& bump)
{
  for (int n = 0; n < NVAR; n++) {
    if (!amrex::Math::isfinite(s_ghost[n])) {
      s_ghost[n] = s_N[n];
      bump(Diag::floored);
    }
  }
}

// ===========================================================================
//  apply() -- fill one ghost cell.
//    s_N, s_Nm1, s_Nm2: interior stencil walking inward; tangential indexing
//      is the caller's (BCfill.cpp), keeping the fill a deterministic pure
//      function of data written before the launch.
//    n_stencil: usable cells (3/2/1), reduced further on EB body state.
//    idir, sgn: boundary normal and side (+1 lo, -1 hi, as bcnormal).
//    layer: 1 = nearest the domain.  diag: optional Long[Diag::count].
//  Returns true when s_ghost was written.  Returns false, with s_ghost
//  untouched, when s_N is EB body state: there is nothing to build a
//  characteristic state from, and the caller must fill the cell some other
//  way (copying the body state out would put a zero-density ghost next to a
//  covered cell).
// ===========================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
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
    // Nothing usable to build a characteristic state from, and the body
    // state itself must not be copied out as a ghost.  Decline; the caller
    // leaves this cell to its ordinary fill.
    bump(Diag::body_state);
    return false;
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
    return true;
  }

  // ---- 1. Outward-normal frame ------------------------------------------
  const auto n_sgn = static_cast<amrex::Real>(-sgn);
  const amrex::Real u_out_N = n_sgn * qN.u[idir];
  const amrex::Real c_N = qN.c;
  const amrex::Real mach = u_out_N / c_N;

  // Frozen acoustic impedance (see the header note: this must be a single
  // value shared by every stencil cell).
  const amrex::Real rho_c = qN.rho * c_N;

  // An inflow Target's temperature, checked once for both inflow paths: an
  // unset T (default 0) would relax the subsonic inlet toward 0 K -- far
  // before any floor fires -- and make PYT2RE divide by zero on the
  // supersonic one.  The interior temperature stands in, counted.
  const bool inflow_tgt = (tgt.type == Type::inflow);
  const bool T_tgt_ok = (tgt.T > 0.0) && amrex::Math::isfinite(tgt.T);
  const amrex::Real T_tgt = T_tgt_ok ? tgt.T : qN.T;

  // ---- 2. Regimes with an exact answer ----------------------------------
  if (std::abs(mach) >= 1.0) {
    if (mach >= 1.0) {
      // Supersonic outflow: zero-gradient is exact, not an approximation.
      // Through a face the problem declared an INFLOW it is also a flow
      // reversal, counted as such (the subsonic path counts it too).
      if (inflow_tgt) {
        bump(Diag::reversed);
      }
      bump(Diag::supersonic);
      copy_state(s_N, s_ghost);
      return true;
    }
    // Supersonic inflow: every characteristic is incoming, so the full state
    // is specified.  Only meaningful with an inflow target.  A Target
    // without a pressure under-specifies it: the interior pressure is
    // substituted and COUNTED (target_incomplete), because the ghost then
    // depends on a domain it should be causally upstream of.
    if (inflow_tgt) {
      if (!(tgt.p > 0.0) || !T_tgt_ok) {
        bump(Diag::target_incomplete);
      }
      // The target composition is cleaned as the subsonic inflow cleans
      // its own, so the packed ghost keeps sum(rho Y) == rho.
      amrex::Real Y_g[NUM_SPECIES];
      if (!normalize_Y(tgt.Y, Y_g)) {
        for (int n = 0; n < NUM_SPECIES; n++) {
          Y_g[n] = qN.Y[n];
        }
        bump(Diag::target_incomplete);
      }
      amrex::Real rho_g = 0.0, e_g = 0.0;
      eos.PYT2RE(tgt.p > 0.0 ? tgt.p : qN.p, Y_g, T_tgt, rho_g, e_g);
      if (
        amrex::Math::isfinite(rho_g) && (rho_g > 0.0) &&
        amrex::Math::isfinite(e_g)) {
        amrex::Real u_g[3] = {tgt.u[0], tgt.u[1], tgt.u[2]};
        pack_ghost(rho_g, u_g, e_g, T_tgt, Y_g, s_N, s_ghost);
        if (tgt.passives) {
          inflow_passives(tgt, 1.0, s_ghost);
        } else if (NUM_ADV + NUM_AUX + NUM_LIN > 0) {
          bump(Diag::passive_unset);
        }
        // The Target's velocity and passives go straight into the ghost on
        // this path; sanitise them as the subsonic path does.
        sanitize_ghost(s_N, s_ghost, bump);
        bump(Diag::supersonic);
        return true;
      }
      bump(Diag::eos_failure);
    } else {
      // Supersonic flow INTO the domain through a face the problem declared
      // an outflow: every characteristic enters and the copy below is
      // ill-posed.  A blast or detonation driven back through an outlet is a
      // flow reversal and must show up as one.
      bump(Diag::reversed);
    }
    bump(Diag::supersonic);
    copy_state(s_N, s_ghost);
    return true;
  }

  // ---- 3. Direction sanity ----------------------------------------------
  const bool outflow_face = (tgt.type == Type::outflow);
  // An outflow Target without a pressure has nothing to relax toward.  The
  // interior pressure stands in (the relaxation then vanishes, and the
  // far-field pin pins to the interior) and the fill is COUNTED: a
  // default-constructed p = 0 would otherwise pull the ghost pressure toward
  // zero on every fill until the floor caught it, visible only as a
  // "floored" count.  An inflow reads p_tgt only through the floor below.
  const amrex::Real p_tgt = (tgt.p > 0.0) ? tgt.p : qN.p;
  if (outflow_face && !(tgt.p > 0.0)) {
    bump(Diag::target_incomplete);
  }
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
    // dither and a spurious 0.3 atm chamber spike (the chamber case and
    // the 1-D driver's reversal tests, both later PRs in this series).
    //
    // The COUNTER carries a roundoff deadband (reversal_deadband()) (the
    // closure itself has no branch left to care): a face holding a quiescent
    // charge sits at exactly u_out = 0 and the solver's noise dithers the sign,
    // which read as 314k "reversals" in the first 100 steps of a chamber whose
    // ignition wave was still a centimetre from the boundary.  Counting only
    // |u_out| > 1e-9 c (wider in single precision) keeps "reversal count
    // zero" meaning NO PHYSICAL
    // BACKFLOW rather than "no roundoff", which is the semantics the
    // advice "a sustained nonzero count means a misplaced outflow"
    // depends on.
    if (u_out_N < -reversal_deadband() * c_N) {
      bump(Diag::reversed);
    }
  }
  // Reversal at an INFLOW (a chamber pushing back through its inlet) runs
  // the same closure, mirrored: the acoustic side keeps relaxing u toward
  // the target -- restoring on both sides of u_out = 0, as at an outflow --
  // and only the lambda_0 content changes hands, from the Target's (entering)
  // to the interior's values (leaving) across a narrow band of outward Mach
  // (1e-3, the width of w_mat's band below).  w_tgt is the Target's share;
  // a hard switch here would be the O(1) discontinuity at u_out = 0 that the
  // outflow side had to remove.  Counted with the same roundoff deadband.
  const amrex::Real w_tgt =
    (u_out_N <= 0.0)
      ? 1.0
      : amrex::max<amrex::Real>(0.0, 1.0 - u_out_N / (1.0e-3 * c_N));
  if (!outflow_face && (u_out_N > reversal_deadband() * c_N)) {
    bump(Diag::reversed);
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
  // 42 us, then NaN (the chamber case of a later PR).  So the material
  // SLOPES are used only where their upwind justification holds: full for
  // u_out >= 0 (every forward-flow result is untouched to the bit), ramped
  // off across a narrow band below zero (1e-3 c -- a continuity
  // regularisation of the fill, not a physics knob; the 1-D driver's
  // reversal tests gate it), and zero
  // under firm reversal, where the ghost keeps the interior VALUES -- which
  // is what the surviving closure always did.
  const amrex::Real w_mat =
    (u_out_N >= 0.0)
      ? 1.0
      : amrex::max<amrex::Real>(0.0, 1.0 + u_out_N / (1.0e-3 * c_N));

  amrex::Real Rm_g;
  if (outflow_face && prm.pin_farfield) {
    // Hard far-field pin: R_- takes its far-field value (u_out = 0,
    // p = p_target).  Value constraint, not a rate, and there is nowhere in it
    // to put a transverse term -- so pin_farfield ignores beta.
    Rm_g = -p_tgt / rho_c;
  } else {
    amrex::Real L_in;
    if (outflow_face) {
      const amrex::Real K = prm.sigma * one_m_M2 * c_N / prm.L_ref;
      L_in = K * (qN.p - p_tgt);
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
    // the diffusion operator reads these ghost cells, so the ghost closure
    // already carries the diffusive physics, and an amplitude-side term
    // counts it twice (measured on the viscous tests of the 1-D driver, a
    // later PR: it moves that error 104 -> -911).  Re-verify there before
    // revisiting.
    Rm_g = Rm_N + static_cast<amrex::Real>(layer) * dx_normal * L_in /
                    (lambda_in * rho_c);
  }

  // ---- 6. Invert the acoustic pair --------------------------------------
  amrex::Real u_out_g = 0.5 * (Rp_g + Rm_g);
  amrex::Real p_g = 0.5 * rho_c * (Rp_g - Rm_g);

  const amrex::Real p_floor =
    amrex::max<amrex::Real>(1.0e-3 * p_tgt, constants::very_small_num());
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

    // Across the reversal band the Target's share w_tgt falls to zero and
    // the interior values take over (w_tgt = 1 for any true inflow).
    const amrex::Real wt = w_tgt * nu;
    // T_tgt is the interior temperature, counted, when the Target left T
    // unset; counted only where the Target still has a share.
    if (!T_tgt_ok && (w_tgt > 0.0)) {
      bump(Diag::target_incomplete);
    }
    T_g = qN.T - wt * (qN.T - T_tgt);
    for (int d = 0; d < 3; d++) {
      if (d != idir) {
        u_g[d] = qN.u[d] - wt * (qN.u[d] - tgt.u[d]);
      }
    }
    // Composition: hard.  Every species characteristic is incoming.  A
    // Target carrying no composition at all is counted, and the interior
    // composition stands in (without this, an unset Y would be swamped by
    // the (1 - w_tgt) interior share inside the reversal band and land on
    // the floored path at w_tgt = 1, under the wrong name).
    amrex::Real Y_t[NUM_SPECIES];
    if (!normalize_Y(tgt.Y, Y_t)) {
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_t[n] = qN.Y[n];
      }
      if (w_tgt > 0.0) {
        bump(Diag::target_incomplete);
      }
    }
    amrex::Real Ysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_g[n] = (1.0 - w_tgt) * qN.Y[n] + w_tgt * Y_t[n];
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
    return true;
  }

  // ---- 8. Pack and sanitise ---------------------------------------------
  pack_ghost(rho_g, u_g, e_g, T_g, Y_g, s_N, s_ghost);
  if (!outflow_face) {
    if (tgt.passives) {
      inflow_passives(tgt, w_tgt, s_ghost);
    } else if ((NUM_ADV + NUM_AUX + NUM_LIN > 0) && (w_tgt > 0.0)) {
      bump(Diag::passive_unset);
    }
  }

  sanitize_ghost(s_N, s_ghost, bump);
  return true;
}

} // namespace pc_nscbc

#endif
