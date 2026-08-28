# NSCBC-Chamber — the T7 mini-SydGex capstone (Lesson 9)

A closed-end premixed chamber venting through its open end, run **two ways**
so that the difference of the two overpressure traces is the boundary
condition's contribution to an application-level answer:

* **`chamber-vent.inp`** — the characteristic boundary IS the vent plane
  (x = 1.2 cm), carrying the flame closures; the front crosses it live.
* **`chamber-plenum.inp`** — the same chamber, but the duct continues to
  x = 3.6 cm and the boundary sits two chamber lengths downstream. Identical
  boundary settings: *placement* is the only difference.

A half-disc of burnt products (r = 1 mm) on the closed-end wall ignites a
quiescent H₂/air φ = 0.4 charge (the shipped PMF profile; S_L = 22.8 cm/s).
The flame fingers down the duct, drives flow out of the open end, crosses the
vent plane, and the emptied chamber rings at its longitudinal modes — the
laminar-phase dilatation (Lesson 5), the crossing (Lesson 6), and the
boundary-set ring-down (Lessons 1 and 3) in one trace.

## Scaling, stated honestly

True SydGex is a 25 cm chamber and a multi-ms turbulent burn. This chamber is
**1.2 × 0.3 cm** (the Sydney L/H = 4) with a laminar burn, ~17 flame
thicknesses long and ~4 across, and the "plenum" is a straight duct
extension — the boundary-placement question T7 poses, without EB geometry.
Two follow-ups are deliberate omissions: a laterally-expanding plenum
(chamber inside a box) needs EB walls, and the baffle of the Sydney rig needs
either EB or an embedded thin wall. Scale up freely on a bigger machine: keep
dx ≤ 0.0125 cm (≥ 5.5 cells per thermal thickness; 0.00625 for the measured
11), grow `prob_hi`/`amr.n_cell` together, and grow `stop_time` with the
front-travel estimate below.

## Running

```sh
make -j TPL && make -j             # or the CMake build; LiDryer / Fuego / Simple
./PeleC2d.<comp>.ex chamber-vent.inp   prob.pmf_datafile=$PWD/LiDryer_H2_p1_phi0_4000tu0300.dat
./PeleC2d.<comp>.ex chamber-plenum.inp prob.pmf_datafile=$PWD/LiDryer_H2_p1_phi0_4000tu0300.dat
```

`prob.pmf_datafile` is CWD-relative — always pass an absolute path. For MPI
builds set `USE_MPI = TRUE` (GNUmake) and run under `mpiexec`; the case is
bit-identical across decompositions like every case in this suite.

Cost model (measured at the shipped scale, serial Apple M-class core):
dt ≈ 6.7×10⁻⁸ s, ~120k steps to `stop_time = 8 ms`, ≈ 2×10⁵ cell-steps/s —
about an hour serial for the vent variant, ~2.5× that for the plenum; minutes
on a many-core node. Front-travel estimate for sizing `stop_time`: the front
leaves the kernel at ~σ_exp·S_L ≈ 100 cm/s and accelerates with its area
ratio; allow (L_chamber − r_kernel)/(150 cm/s) plus ~1 ms of ring-down.

## Measuring

```sh
./chamber_metrics.py runs/vent runs/plenum --xvent 1.2
```

prints, per plotfile: closed-end overpressure, chamber-mean overpressure,
burnt volume fraction, and the volumetric flux through the vent plane. The T7
QoIs are built from the **difference of the two runs' traces**: peak-aligned
overpressure (peaks aligned, not clocks — ignition transients differ),
V̇_comb − V̇_vent whose zero crossing marks the peak, and the post-peak
ringing frequency and decay rate, which is where the boundary's reflection
coefficient is directly visible (Lesson 9's punchline: repeat the vent
variant at σ = 16 and watch the ring-down change while the peak barely
moves).

## Boundary settings

The vent carries the flame-crossing recipe measured in `NSCBC-FlameOutflow`:
`extrap_temperature = 1`, `beta_s = 0`, `beta = 0.5`, σ = 0.25 —
acoustically open, so the ring-down shows the chamber rather than the
boundary. `extrap_material` stays off (a front crosses). The plenum variant
uses identical settings so placement is the only difference; with the
boundary two chamber lengths away, its settings are nearly irrelevant, which
is the point.

## Status

Production-complete, and the case has now killed THREE defective reversal
closures on its way here — it is the only configuration in the suite whose
boundary genuinely reverses (the vent breathes), so it finds what nothing
else can. The history, in order: the original hard ambient pin NaN'd the
first production run at 4.08 ms (backflow −144 → −2628 cm/s, a 138 K cell);
the soft pressure-only replacement survived but manufactured a spurious
0.32 atm chamber spike at 4.0–4.6 ms by drawing inflow against a
+0.3 atm adverse gradient (dropped `S_p`/`dR₊`/`T_in`, frozen ghost
velocity); and the first unification of the closure NaN'd again at 4.05 ms
by extrapolating material content into the backflow (the cold runaway,
385 → 89 K in 42 µs). Full forensics:
`Docs/NSCBC-reversal-branch-defect.md`; static gate: driver C13 a/b/c. The
shipped closure — unified acoustics, `w_mat`-upwinded material slopes —
completes all three variants to 8 ms.

The T7 comparison, measured with `chamber_qoi.py` on the 2026-08-25
production set (identical boundary settings, placement the only difference):

| variant | V̇ plateau [cm²/s] | end-of-burn peak [dyn/cm²] | f_ring [kHz] | decay λ [1/s] |
|---|---|---|---|---|
| vent, σ = 0.25 | 72.2 | **+37,473** | 3.33 | 674 |
| vent, σ = 16 | 71.6 | **+6,875** | 14.16 | 504 |
| plenum (same window) | 72.1 | +1,951 | 0.51 | 184 |

Three phases, three answers. **Buildup**: placement costs nothing — both
variants drive identical venting (V̇ 72.2 vs 72.1 cm²/s); the plenum's
standing +~3×10³ dyn/cm² offset is the inertia of the 2.4 cm duct column it
pushes, physics any boundary would see. **End of burn** (t ≈ 4.0–4.6 ms):
the last ~3% of charge burns while the vent flow chokes — a real
constant-volume pressurisation — and how much of it the boundary lets
through is where placement bites: +37k at σ = 0.25 against the plenum's
+2k, halved and halved again by σ = 16. This is the honest remaining
boundary cost, down 8.5× from the defective closure's 321k, and with NO
counter-gradient inflow (V̇ dips to −9 cm²/s where the defective closure
drew −374 against the gradient). **Ring-down** is where σ is directly
legible: at σ = 0.25 the "mode" is the relaxation itself
(3.33 kHz ≈ K/2π = σc/2πL), while σ = 16 stiffens the end toward a real
acoustic termination and the chamber rings near its burnt-gas quarter-wave
(14.2 measured vs ~18.8 kHz ideal) — the Lesson-9 punchline, measured. The
σ = 16 run also takes brief genuine backflow breaths (V̇ to −110 cm²/s)
through the fixed closure without incident.

Not registered in `Tests/CMakeLists.txt`: at ~1 CPU-hour per variant it is
an application capstone, not a CI gate.

## The EB follow-ups: the lateral plenum and the baffle

The two deliberate omissions, now built (2026-08-25):

* **`chamber-box.inp`** — the same 1.2 × 0.3 chamber as embedded walls
  INSIDE a 2.4 × 0.9 box: the vent at x = 1.4 opens into a plenum three
  chamber heights tall, and the characteristic outflow sits on the far x-hi
  face a chamber length downstream. This is the duct-plenum comparison with
  a *lateral* expansion — the geometry T7 actually poses.
* **`chamber-baffle.inp`** — `chamber-box.inp` plus the Sydney baffle: two
  plates across the chamber at x = 0.8 (thickness 0.05) leaving a central
  0.1 gap — 67% blockage. The A/B against `chamber-box.inp` is the baffle's
  contribution with everything else identical.

Both use one registered geometry (`eb2.geom_type = "nscbc-chamber-box"`,
built in `prob.cpp` from `prob.chamber_*`/`prob.baffle_*`; the ignition
kernel follows the EB closed end via `prob.kernel_x`). The whole chamber,
closed end included, is interior EB — nothing touches a domain face, which
is the documented EB-at-domain-boundary NaN limitation — and
`pelec.eb_zero_body_state = 1` is mandatory as in every EB+NSCBC
configuration. Boundary settings are the flame-crossing recipe, unchanged.
Both variants run `pelec.do_mol = 1` (the EB-supported hydro, as every EB
RegTest in the tree attests — the Godunov path died first), with the EB
surfaces placed OFF the grid lines (x.4·dx via `prob.chamber_yc`; a surface
on a cell face makes degenerate cut cells), and — the one that actually
mattered — **`pelec.eb_isothermal = 0`**. PeleC's EB walls default to
isothermal, and `eb_boundary_T` DEFAULTED TO 1 KELVIN when these variants
were built (since changed to 300 K because of this case); the first three
production attempts all died of it: the quadratic EB gradient faithfully conducts heat
into a 1 K wall, the cut cells refrigerate first (vfrac = 0.6 row at 268 K
by 25 µs, marching to 77 K by 1.6 ms), the corner where two cold walls meet
is coldest, and the arriving flame NaNs there — under either hydro scheme,
at any grid alignment. A 500-step probe pinned it: min fluid T 268.1 K with
default walls (identically with `bc_nscbc` off and either body state — not
the boundary, not zero-state poisoning), 298.0 K exactly with adiabatic
walls. Adiabatic EB matches the adiabatic domain walls, so the box is
thermally consistent with the duct variants. The scheme caveat stands: the
box-vs-duct-plenum comparison differs in hydro (MOL vs Godunov) as well as
geometry; the box-vs-baffle A/B is internally consistent.

Building these taught two counter lessons on day one. The reversal counter
read 314k in the first 100 steps, before the ignition wave was within a
centimetre of the boundary — measured at the face, that is not noise but a
coherent −0.28 cm/s micro-inflow across the whole outflow (Mach 8×10⁻⁶):
the EB small-cell fixup leaves the box a hair under-pressurised and the
σ-relaxation breathes it back in. Real backflow, honestly counted — so
read the counter with the magnitude of the flow in mind; equilibration
counts big and means nothing. Separately, the counter now carries a 10⁻⁹c
deadband (counting only — the closure is branch-free), so a face holding a
charge at *exactly* u = 0 cannot count pure sign-noise, and a zero count
means no physical backflow. Also expected: with `eb_zero_body_state = 1`
the covered cells hold zeroed states, so derived velocities in plotfile
dumps read NaN *inside the walls* — an artifact of the flag, not of the
solve.

Cost at dx = 0.0125: 192 × 72 = 13,824 cells, ~168k steps to 8 ms, 2.65 h
at 6 ranks. Metrics: `chamber_qoi.py --xvent 1.405 --xlo 0.205
--ymask 0.305,0.605` — the mask restricts every chamber reduction to the
cavity interior and drops zeroed EB body cells, which the box needs and the
duct variants don't.

**The lateral-expansion answer (first production set, 2026-08-25).** The
box against the straight-duct plenum, phase by phase (closed-end
overpressure in dyn/cm², vent flux in cm²/s; the ignition IC is clipped to
the chamber interior — an unclipped kernel pokes through the 0.05-thick
closed end and lights the plenum behind the chamber, which is how the
first attempt failed):

| | box (lateral plenum) | duct plenum |
|---|---|---|
| dp_closed at 0.6 ms | **1,239** | 4,251 |
| buildup character | ~900–1,400, flat | 1,600–4,300, humped |
| V̇_vent plateau | 46–53 | 67–73 |
| chamber burnout (bf → 0.99) | ~6 ms | ~4.3 ms |
| dp_closed at 8 ms | 3,705, rising | 2,679, past its 3,437 peak |

The lateral expansion is a ~3× better vent: the chamber sees a third of
the duct-plenum's overpressure throughout its burn, because the vented gas
escapes sideways instead of having to push a chamber-length column of duct
gas. The price is a slower burn — the weaker chamber–plenum pressure
difference drives ~30% less vent flow and the chamber takes ~40% longer to
empty. At late times both variants climb as the *extension* charge burns
(the whole domain carries the premixed charge, by the duct convention), so
neither shows a clean ring-down inside 8 ms and that QoI is honestly n/a
here. The remaining caveat from the build sequence stands: box-vs-duct
differs in hydro scheme (MOL vs Godunov) as well as geometry.

**The baffle answer (same production set).** Baffle against box, everything
identical but the plates (dp_closed in dyn/cm², V̇ in cm²/s):

| | box | baffle |
|---|---|---|
| buildup to 2 ms | ~900–1,400 | ~900–2,400, climbing |
| gap jet (2.5 ms) | V̇ = 47, bf = 0.46 | **V̇ = 281, bf = 0.71** |
| downstream burnout | bf = 0.63 at 3.5 ms | **bf = 0.95 at 3.5 ms** |
| peak dp_closed | 3,733 (the 8 ms plenum-burn tail) | **31,725 at 3.7 ms** |
| vent backflow | V̇ min −0.0 | **V̇ min −67** |
| reversal fills (whole run) | 1.5M | **14.1M** |
| state at 8 ms | 3,705, still climbing | settled to ~200–700 |

The Sydney phenomenology in miniature, all of it: the flame squeezes
through the 0.1 gap and accelerates 6× (V̇ 47 → 281), torches the
downstream half of the chamber in under a millisecond, and drives the
closed-end overpressure to 8.5× the unbaffled maximum — the
baffle-enhanced explosion peak this geometry exists to produce. The jet's
violence puts the vent plane into genuine sustained breathing (V̇ to −67;
14.1M reversal fills against the box's 1.5M — ten times the reversal duty,
all through the unified closure without incident), and because the jet
also consumes the plenum charge early, the baffle run is the only
configuration in the whole suite that reaches a settled state by 8 ms.
Both runs: zero NaN, ~168–174k steps, 2.7/3.4 h at 6 ranks.
