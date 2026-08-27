# NSCBC-FlameOutflow — a flame sitting *on* a characteristic outflow

A wrinkled, unanchored premixed sheet (LiDryer, φ = 0.4 H₂/air, S_L = 22.8 cm/s)
in a uniform stream at U = 2 S_L, with its mean position exactly on the outflow
plane. Half the boundary is in burnt gas, half in fresh, and the reaction zone
passes through the boundary cells at the two crossings — which is where the
wrinkle slope is largest, so the flame normal there is tilted 51° off the
boundary normal.

```sh
./PeleC2d.gnu.ex nscbc-flameoutflow.inp
./PeleC2d.gnu.ex nscbc-flameoutflow.inp pelec.bc_nscbc_sigma=16.0   # see below
```

The stream carries the sheet downstream at U − S_L = S_L, so nothing anchors it
and nothing has to be held in place. It is the configuration
`NSCBC-PMF/README.md` says is needed and cannot provide: a live heat-release
term, strong tangential gradients, and a multi-species diffusive structure, all
inside the boundary cells, from step 0.

What the case does **not** do is advect out within a measurement. The drift Mach
number is S_L/c ≈ 6×10⁻⁴, so clearing the boundary takes some 400 acoustic
transits. The sheet is initialised already straddling and is quasi-frozen over
the few relaxation times the boundary needs. That is not a compromise: the
quantity under test is the boundary's response to a reaction zone in its own
cells, and that is present throughout.

## The measurement protocol

Earlier attempts at this measurement were worthless, for two reasons worth
recording because both are easy to repeat.

**There is no useful self-referencing metric here.** The configuration has a
genuine mass imbalance — the front drifting downstream converts light burnt gas
into dense fresh gas, so the domain really does gain mass — and the resulting
pressure ramp is several times larger than any difference between boundary
settings. Comparing a run against its own domain mean measures the ramp, not the
boundary.

**A longer-domain reference is only a reference if it differs in one thing.**
The obvious reference — same problem, outflow moved downstream — differs in two
others: the physical ramp is spread over a larger volume, and, because
`L_ref = ProbHi − ProbLo`, *every* relaxation rate `K = σ(1−M²)c/L_ref` changes
with it, including the one at the inlet.

The protocol that works:

* `prob.nscbc_inflow = 0` in **every** run, test and reference. The inlet is a
  hard Dirichlet and the characteristic treatment is on the outflow alone, so
  the length-dependent inflow rate cannot contaminate anything.
* The reference puts its outflow at x = 3.0, i.e. 2.4 cm downstream of the test
  outflow. In burnt gas at c ≈ 8.6×10⁴ cm/s nothing from it can reach x < 0.6
  before t = 2.8×10⁻⁵ s. Runs stop at 2.4×10⁻⁵ s.
* Up to that time the reference restricted to x < 0.6 **is** the exact solution
  of the test problem: identical grid, identical scheme, identical inlet, and
  its outflow outside the domain of dependence. The difference field is the test
  run's outflow error.
* Because of that shielding the mean level is part of the error and must **not**
  be subtracted. Getting this wrong was the earlier mistake that made a hard
  Dirichlet outflow look better than it is.

`measure.py` implements it.

## What it measures

Errors against the shielded reference at t = 2.4×10⁻⁵ s (≈ 3.4 τ_relax at
σ = 1), in dyn/cm²; p_amb = 1.013×10⁶, so 1000 is 10⁻³ of ambient. `R` is the
acoustic reflection at the same σ from `Verification/NSCBC1D`. All
characteristic rows: β = 0.5, order 2, hard-Dirichlet inlet. `eT` is
`bc_nscbc_extrap_temperature`. Measured after the Phase-0 reaction-source sign
fix (`NSCBC.H`, `L_in += (1−β_s) S_p`) and the unified reversal closure
(`Docs/NSCBC-reversal-branch-defect.md`; gate: driver C13); the earlier
tables are in the git history. The recipe rows reproduce to the digit
(−506.5/568.3/638.1 and −521.9/581.5/646.2 against the soft-closure era's
−507/568/638 and −522/582/646) — never-reversing histories are bit-identical
across the closure changes, measured end-to-end. The failure-mode rows all
read HIGHER than under the soft closure (+14–22%). The reversal counter
(run with `pelec.sum_interval` set — this inp defaults it to −1, so counters
are silent unless you ask) shows why: in the σ = 1, eT = 0, β_s = 1 row the
ENTIRE outflow face is in reversal essentially every fill (12,800 per
25-step report = 64 rows × 4 layers × 2 stages × 25) — the piled-up ramp,
~2400 dyn/cm² against a 45.6 cm/s stream, stalls and reverses the
throughflow wholesale, so these failure rows are measurements *of* the
reversal closure. Each earlier closure bled pressure through that reversal
(the pin pinned it to ambient, the soft branch relaxed it with a frozen
ghost velocity), flattering exactly the rows it was corrupting; the unified
closure holds the reversed face with frozen material content and the same
σ-rated relaxation as forward flow, and these are the honest numbers.

| outflow | σ | eT | β_s | mean Δp | L2(Δp) | L2(Δp) bl | R [%] |
|---|---|---|---|---|---|---|---|
| hard `p = p_amb` | — | — | — | −527 | 586 | 631 | — |
| characteristic | 0.25 | 0 | 1 | +2752 | 2775 | 2896 | 0.76 |
| characteristic | 0.25 | 0 | 0 | +1592 | 1636 | 1873 | 0.76 |
| characteristic | 0.25 | 1 | 1 | +2429 | 2458 | 2655 | 0.76 |
| **characteristic** | **0.25** | **1** | **0** | **−507** | **568** | **638** | **0.76** |
| characteristic | 1 | 0 | 1 | +2377 | 2397 | 2455 | 2.56 |
| characteristic | 1 | 0 | 0 | +1167 | 1205 | 1299 | 2.56 |
| characteristic | 1 | 1 | 1 | +1981 | 2006 | 2113 | 2.56 |
| characteristic | 1 | 1 | 0 | −522 | 582 | 646 | 2.56 |
| characteristic | 16 | 0 | 1 | +39 | 255 | 142 | 28.1 |
| characteristic | 16 | 0 | 0 | −282 | 378 | 397 | 28.1 |
| **characteristic** | **16** | **1** | **1** | **−139** | **286** | **270** | **28.1** |
| characteristic | 16 | 1 | 0 | −511 | 570 | 618 | 28.1 |

(The pre-fix table also carried β = 1 rows, an order-1 row and `pin_farfield`;
none of those knobs changed in Phase 0 and their conclusions stand: β = 0.5 is
a free ~10%, order 2 is load-bearing, `pin_farfield` anchors to
p_target + ρc·u_out. See the git history for the numbers.)

What the corrected table says, in order of importance.

**β_s = 0 now earns its keep, and the old table was measuring the bug.** With
the sign corrected, switching the reaction source on cuts the error by 42%
at σ = 0.25 (2752 → 1592) and by half at σ = 1 (2377 → 1167). Under the old
sign the same switch did nothing (+2065 → +2074), because it was *doubling*
the heat-release push and the visible residue was only the difference between
2 S_p/K and the σ-dominated rest. The Sutherland–Kennedy cancellation this
term implements was always the right physics; the branch had its sign
inverted from the first commit to Phase 0.

**extrap_temperature and β_s = 0 together erase the model error, at every σ.**
The eT = 1, β_s = 0 rows read −507 / −522 / −511 for σ = 0.25 / 1 / 16: the
σ-dependence that dominated every other configuration is *gone*, and the
error lands at the hard-outflow level (−527 / 586) — but at σ = 0.25 it costs
0.76% reflection where the hard boundary reflects everything. This is the
consistent-closure result: extrap_temperature makes the ghosts carry the
diffusive dp/dt and β_s = 0 makes the incoming wave carry the chemical one,
and what remains no longer scales with the relaxation. Why the remaining
≈ −515 coincides with the hard boundary's error to 3% at all three σ is not
established; it smells like a shared floor, and nothing currently gates it.

**σ = 16 with β_s = 1 is still the best absolute error, and β_s = 0
overcorrects there.** +49 / 256 (eT = 0) and −139 / 286 (eT = 1) beat every
other row's L2, at the known price of 28% reflection. Adding β_s = 0 on top
pushes the mean through zero to −281 / −511: at high σ the relaxation is
already absorbing most of the chemical offset, and the source term then
pushes past the mark — the 1-D S_p accounting is not exact against this
tilted, 2-D front.

**σ is no longer the only control that matters.** In the old table only
σ = 16 beat the hard outflow. Now σ = 0.25 with the consistent closures
matches it while staying transparent, and the sentence "the default is the
worst entry in the table" applies only with the closures off.

## What the error is made of

The bisection below predates the sign fix but measures configurations the fix
does not touch (β_s never enters: chemistry off, or β_s = 1), so it stands.
At σ = 1 (t = 10⁻⁵ s, domain-mean pressure rise):

| configuration | mean Δp |
|---|---|
| uniform fresh gas, no flame in the domain | +25 |
| flame, `do_react = 0` (it just diffuses) | +417 |
| flame, reactions on, diffusion off | +793 |
| flame, reactions and diffusion | +1283 |
| flame, hard outflow | −23 |

**Half the error is there with chemistry switched off.** A purely non-reacting
density front already biases the boundary; reactions roughly double it. That
is why β_s = 0 alone halves the σ = 1 error rather than removing it (2377 →
1167 in the table above), and why extrap_temperature — the closure that lets
the ghosts carry the front's diffusive structure — is the other half of the
recipe.

The mechanism is in the ghost pressure. At an outflow

```
p_g = ½ ρc (R₊,g − R₋,g)
R₊,g = R₊,N + ℓ δR₊                        (extrapolated, limited)
R₋,g = R₋,N + ℓ Δn L_in / ((c − u_out) ρc)  (relaxed)
```

so, per ghost layer ℓ,

```
p_g − p_N  =  ½ ρc ℓ δR₊  −  (ℓ σ / 2 n_x) (p_N − p_∞)
                 ^^^^                ^^^^
              extrapolation        anchoring
```

The anchoring increment carries a factor 1/n_x — with n_x = 96 and σ = 1 it is
0.5% of the pressure offset it is trying to remove. The extrapolation term
carries no such factor. Setting the two equal gives the equilibrium offset

```
Δp  ≈  ρ c L_ref (du_out/dn)|_boundary / σ
```

which is the whole story: a **normal velocity gradient at the boundary** biases
the ghost pressure, and the relaxation can only fight it in proportion to σ. A
flame crossing the outflow is a large `du_out/dn` — the gas accelerates from
45.6 to 123.8 cm/s across 0.07 cm — and it is dilatational, not acoustic, so
there is nothing wrong with the invariant algebra. LODI simply has no way to
tell the two apart.

Two checks that the formula is the right one rather than a plausible story:

* It predicts a σ⁻¹ trend. Measured (β_s = 1, eT = 0): 2377 → 255 for
  σ = 1 → 16, and 2074 → 1450 → 218 for σ = 1 → 4 → 16 in the pre-fix table,
  whose β_s = 0 rows differ from β_s = 1 by under 1%.
* It predicts the offset is **grid-converged**, because δR₊ is a per-cell slope
  falling as 1/n_x while the anchoring carries 1/n_x explicitly. Measured
  (no chemistry, σ = 1, t = 10⁻⁵ s): n_x = 48 → +531, 96 → +409, 192 → +371.
  It converges to a finite value; it does not refine away.

## What to tell a user

This is the case that turns "place outflows away from flames" from received
wisdom into a number. If you must put one there:

* Turn on the consistent closures: `bc_nscbc_extrap_temperature = 1` and
  `bc_nscbc_beta_s = 0`. Together they hold the error at the hard-outflow
  level at ANY σ, so you can keep a transparent boundary (σ = 0.25, 0.76%
  reflection) instead of buying anchoring with reflection.
* If absolute pressure error matters more than transparency, σ = 16 with
  `beta_s = 1` is still the best row (L2 ≈ 260–290) at 28% reflection. Do
  not combine σ = 16 with β_s = 0 — the two corrections overshoot together.
* Or use `pin_farfield`, if a fixed offset of order ρc·u_out is acceptable and
  transparency is not.
* Set `bc_nscbc_beta = 0.5`. It is a free 10%.
* Leave `bc_nscbc_order = 2`.
* Run with `pelec.sum_interval > 0` at least once and read the fallback line.
  A `reaction source dropped` or `transverse dropped` count means a dial is
  silently doing nothing.

## Traps

**`pelec.allow_negative_energy = 0`** is harmless in an inert case and fatal
here — see `NSCBC-PMF/README.md`. The inputs file sets it explicitly to 1.

## Locating the reaction zone

Use a radical, not the temperature. The temperature midpoint sits in the preheat
zone, which is wide and diffusion-controlled; the radical peak sits in the
reaction zone, which is the thing that actually has to be inside the boundary
cells for `beta_s` to have anything to correct. `radicals.py` does both the
front-position tracking and the boundary-column profile.

Along the outflow column at t = 24 µs, H₂/air, σ = 1:

```
  y [cm]    T [K]         Y_H        Y_OH    q [erg/cm3/s]
  0.0031      299   1.494e-17   4.342e-14        9.387e+00
  0.1031      644   1.963e-07   4.007e-06        2.337e+08
  0.1281     1008   2.348e-05   2.532e-04        5.549e+09
  0.2531     1215   5.116e-05   9.520e-04        5.961e+09
  0.3781      301   6.946e-17   2.304e-13        8.336e+01
```

Twelve orders of magnitude in `Y_H` and nine in the heat release, along a single
column of boundary cells, with the two peaks at the two designed crossings. That
is the case doing what it claims.

The radical is also the better front tracker for the stability question below,
because the peak is 2–3 cells wide where the temperature midpoint is a ramp over
ten.

## Is the measurement contaminated by thermodiffusive instability?

A fair question of the H₂ case: at φ = 0.4, Le(H₂) ≈ 0.3, so an imposed wrinkle
is thermodiffusively unstable and would grow. It does not, because there is no
time for it to. Tracking the front by the H peak:

| t [s] | mean x_f [cm] | wrinkle half-amplitude [cm] |
|---|---|---|
| 0 | 0.65294 | 0.07990 |
| 9.78e−6 | 0.65262 | 0.07896 |
| 2.40e−5 | 0.65247 | 0.07809 |

The amplitude *decays* 2.3% over the run. The Darrieus–Landau e-folding time at
this wavenumber is ~9×10⁻³ s and the flame time δ/S_L is 3×10⁻³ s, against a
measurement window of 2.4×10⁻⁵ s — 130 to 390 times shorter than either. The
boundary equilibrates on the acoustic time, which is what makes the measurement
possible at all, and the flame is frozen on that time.

`NSCBC-FlameOutflow-DRM` runs the same problem with CH₄/air at φ = 0.75
(drm19, Le(CH₄) = 0.97, thermodiffusively neutral) and settles the question by
construction rather than by argument.

## A diffusive boundary condition (`bc_nscbc_extrap_temperature`)

The measurements above were made before the outflow had any diffusive treatment
at all. `Source/Diffusion.cpp` forms the conductive and species fluxes at a
physical boundary face from these ghost cells, and in the entropy closure the
ghost temperature is whatever the EOS returns from the extrapolated density and
the acoustically-set pressure — nothing chooses it with a heat flux in mind.

Check C8 in `Verification/NSCBC1D` puts a number on it. On a flame-like ramp
(2×10⁴ K/cm, 0.01 cm cells) the ghost overstates the face temperature gradient
by **40%**: T_ghost = 1680 K where the interior ramp continues to 1600 K.

`pelec.bc_nscbc_extrap_temperature = 1` closes the λ₀ family on temperature
instead — T is extrapolated on the same minmod slope as everything else and ρ
follows from the EOS — so the face gradient is the interior one, exactly. A
uniform state still comes back to 2×10⁻¹⁶, so nothing hyperbolic is given up.

Its measured effect is in the main matrix above: at β_s = 1 it buys 17% at
σ = 1 (2377 → 1981) and trades mean for near-boundary structure at σ = 16
(+49 → −139 mean, but L2(dT) in the boundary layer 4.6 → 2.7); combined with
the corrected β_s = 0 it is half of the closure pair that removes the
σ-dependence outright. (An earlier version of this section quoted −69% at
σ = 16 from β_s = 0 rows measured under the inverted source sign; those
numbers are in the git history and were measurements of the bug.)

The default is off, because it changes every outflow result. Turn it on when a
thermal or compositional structure is anywhere near the boundary.

## The material continuation (`bc_nscbc_extrap_material`), and what it settles

The ghost-pressure bias of the previous section has a targeted fix:
continue the *material* part of the incoming invariant's slope into the
ghost, bounded through the entropy family so the incoming model never feeds
on the waves it launches (`Source/NSCBC.H`; gated by checks C9(a) and C11 in
`Verification/NSCBC1D`, where it removes the bias exactly and holds a
manufactured sustained front that the entropy closure walks away from).

Measured here, σ = 1, β = 0.5, β_s = 0, t = 2.4×10⁻⁵ s, mean Δp:

| closure | mean Δp | change |
|---|---|---|
| entropy | +2074 | — |
| + material | +1972 | −5% |
| + temperature | +1894 | −9% |
| + **material and temperature** | **+1200** | **−42%** |

Two conclusions, and the first one closes an open question.

**The extrapolation bias is real but is not what dominates this case.** The
fix removes the C9(a) mechanism *exactly* in the driver and buys 5% here, so
the σ-suppressed error of this configuration is mostly something else — by
the bisection above and the C9(b) result, the unmodelled diffusive and
reactive enthalpy deposition in the boundary cells. That also answers the
question the C9/C10 commits left open, **against** the flux-form
reformulation's premise: a flux-form NSCBC removes the same extrapolation
bias, so it too would buy ~5% here, at many times the cost.

**The residual is not an unmodelled diffusive source either.** The natural next hypothesis — supply the missing
"viscous condition" as a dp/dt|_diffusion term in the modelled incoming wave — was built, verified exact on
quadratic profiles, and refuted by measurement: it moves this case from +1200 to +1771, and the 1-D conduction
test (C12) from +104 to −911. In the ghost-cell form the diffusion operator reads the ghost cells, so
`extrap_temperature` already carries the diffusive physics and an amplitude-side term counts it twice. What
remains of the σ = 1 error is therefore multi-dimensional (the tilted front's tangential structure) and/or the
2× under-continuation of the material bound at U = 2 S_L — not a missing 1-D source term.

**The two ghost closures compound.** Together they make the ghost a fully
consistent material continuation — p from the corrected acoustic pair, u
from the continued slope, T extrapolated, ρ from the EOS — and take 42%
where each alone takes single digits. Note two caveats: at U = 2 S_L the
quasi-steady bound in the material continuation underestimates the front's
true slope by 2× (the sheet drifts at S_L = U/2), so this configuration is
near its worst case; and the continuation is for *quasi-steady* structure —
see the exit test below before enabling it anywhere a front actually crosses.

## The exit test — `nscbc-flameexit.inp`

This case's original question was always "what happens when the flame
reaches the outflow", and the quasi-frozen configuration answers it only for
a front that *sits* there. `nscbc-flameexit.inp` makes the sheet actually
leave: U = 40 S_L, the sheet initialised fully inside, and the whole passage
— approach, crossing, exit, emptied domain — inside 6×10⁻⁴ s. The wrinkle
keeps the parent's slope (51° tilt at the crossings) and the grid keeps the
parent's dx. Physics fixes the truth during transit (the front must advect
at U − S_L ≈ 889 cm/s with its wrinkle frozen; the decay times are 100× the
window), and after the exit the exact solution is known outright: a uniform
fresh stream at (U, p_amb, T_in). `exit_metrics.py` tracks both.

All characteristic rows below: β = 0.5, order 2, `extrap_temperature = 1`
unless marked bare, `extrap_material = 0` — NOTE: the inp defaults
`extrap_material = 1` (the case's comment argues the bound is accurate at
40 S_L), so table rows must pass `pelec.bc_nscbc_extrap_material=0`
explicitly; with eM = 1 the transit over-vents catastrophically (front
expelled by 1.1×10⁻⁴ with a −0.37 atm excursion) in EVERY configuration,
which is the known eM-transit failure, not a boundary ranking. Re-measured
after the unified reversal closure (C13-gated), with the shipped (clamped)
tangential stencil. The recipe row reproduces **to the digit** (+2422 /
−8818 / −250, 6.7) with its reversal counter at zero through the whole
crossing and exit — a healthy transit never reverses the boundary — and the
failure rows' pre-reversal histories are also bit-preserved: the blocked and
bare-β_s = 0 transit peaks match to the digit (+24523, +23419; both predate
their reversals), while bare-β_s = 1's peak is late-time and moved with the
closure. What changed everywhere is the reversal-dominated late evolution
(see below).

| outflow | transit peak Δp | post-exit extreme | front | post-exit residual Δp, rms(u−U) |
|---|---|---|---|---|
| hard `p = p_amb` | −147 | −203 | on schedule | **+1.2, 0.01** |
| σ = 16, β_s = 1 | +853 | +515 | on schedule | **+8.5, 0.19** |
| **σ = 16, β_s = 0** | **−110** | −432 | on schedule | +8.5, 0.19 |
| σ = 1, β_s = 1 | +9119 | +15926 | slightly late | +134, 3.5 |
| σ = 1, β_s = 0 | −1224 | −4710 | on schedule | +133, 3.5 |
| σ = 0.25, β_s = 1 | +24523 | +54449 at 4.8e−4, then DECLINING | **BLOCKED: stalls, pushed back** | never exits |
| σ = 0.25, β_s = 0 | +2422 | −8818, recovering | on schedule | −250, 6.7 |
| σ = 0.25 bare (eT = 0), β_s = 1 | +58615 | — | pushed backwards, domain refills with burning | never recovers; maxT peaks 3444, easing |
| σ = 0.25 bare (eT = 0), β_s = 0 | +23419 | — | exits, slow, wrinkle destroyed; remnant re-anchors | +27500 oscillating offset |

Wrinkle amplitude holds at ≈ 0.034 in every completing run while the full 32
rows track the front, decaying only as rows leave the domain, so the wrinkle
column of the old table is subsumed by "front".

What the corrected table says, in order of importance.

**The corrected reaction source turns the transit from a σ = 16-only
manoeuvre into something any σ survives.** With eT = 1 and β_s = 0 the front
crosses on schedule at every σ measured, and the transit disturbance falls
7–10× against β_s = 1 at the same σ (+9119 → −1224 at σ = 1; +24523 → +2422
at σ = 0.25). At σ = 16 the characteristic boundary is now *quieter during
the crossing than the hard Dirichlet* (−110 against −147).

The physics, spelled out because the phrase "no model for dilatation" is
easy to over-read. The bare relaxation L_in = K(p − p_target) has its fixed
point at p = p_target with L_in = 0: the equilibrium it relaxes toward is a
uniform, source-free stream — correct for the LODI limit it comes from. A
reaction zone in the boundary cells breaks that assertion: the pressure
equation in characteristic form is ∂p/∂t = −½(L₊ + L₋) + S_p, so a steady
front on the face needs ½(L₊ + L₋) = S_p — sustained heat release is
sustained dilatation, and the expansion must leave through the waves. The
outgoing L₊, measured from the interior, carries its share automatically;
the bare incoming model supplies zero at p = p_target, so the unbalanced
source integrates into overpressure until the relaxation cancels it —
a standing offset ∝ S_p/K, i.e. ∝ σ⁻¹, which is exactly the measured
2377 → 255 column and why σ = 16 was the only working configuration before
Phase 0. β_s = 0 adds the missing term (+S_p on L_in, Sutherland–Kennedy):
the steady balance then closes AT p_target, and it is transit-safe because
S_p is an instantaneous pointwise measurement that asserts nothing about
persistence and vanishes as the reaction zone leaves.

This is not a formulation gap so much as a slot that was empty. Dilatation
is not a characteristic variable; it splits into a normal part (shared
between R₊ and R₋) and a transverse part, and the formulation has exactly
three places for it, all now populated: chemical normal dilatation on the
amplitude side (β_s, above); transverse dilatation on the amplitude side
(the γp ∇_t·u_t piece of the β term); and the front's mechanical structure
(the mass-conservation velocity rise) on the material side — order-2
slopes, `extrap_material`'s entropy-bounded continuation, the ghost T
closure. (The diffusive dp/dt deliberately stays OFF the amplitude side:
the diffusion operator reads the ghosts, so the ghost T closure already
carries it and an amplitude term double-counts — C12, 104 → −911.) What
genuinely cannot be added in this formulation is narrower: classifying the
measured du/dn at one instant as sustained structure (continue it),
transiting structure (release it — continuing over-vents, the measured
eM-during-transit failure), or acoustics (extrapolating feeds back, C5).
The fill is a pure function of the instantaneous interior state — the
determinism invariant — and those three have identical instantaneous
signatures with opposite correct responses. The measured ways out are the
C11x source-consistency bound (admit only the du/dn the local sources can
sustain, du/dn = S_p/ρc², per the same Sutherland–Kennedy relation) and
the queue-item-4 trend registers (history, kept outside the fill).

**The blocked state is the failure the closures half-on produce, and each
generation of the reversal closure reported it differently — only the
unified closure's number is gated.** σ = 0.25 with eT = 1 but β_s = 1
builds a ramp that stalls the front near x = 0.55 and pushes it back
upstream; the flame never leaves the domain. The row's early history is
bit-preserved across all three reversal closures (transit peak +24523 to
the digit — reversals only begin later), and then the closures diverge: the
hard pin read +0.04 atm (artificially clamping the ghosts to ambient during
backflow), the soft branch read +0.15 atm and climbing (its frozen ghost
velocity let the reversal feed the ramp), and the unified closure holds the
peak at +0.054 atm at 4.8×10⁻⁴ and DECLINING by 6×10⁻⁴ — the reversed face
(4.6 million reversal fills over the run) is finally being relaxed with the
same restoring σ-rated closure as forward flow. Survival of the *run* is
still not survival of the *physics*: the front remains blocked.

**The bare default still destroys the physics, with either β_s — it just no
longer NaNs.** Both bare rows used to die (1.9×10⁻⁴ and 4.1×10⁻⁴): the front
pushed back through the outflow reverses it, and the original closure
answered backflow with a hard ambient pin whose pressure discontinuity fed
NaNs — the same defect `NSCBC-Chamber` found in the vent ring-down. Under
the unified closure both runs reach 6×10⁻⁴ intact (15.2 and 2.1 million
reversal fills respectively) and show what the bare configuration actually
does: with β_s = 1 the front is driven backwards and the domain refills
with burning gas (maxT peaks at 3444 K at 5.6×10⁻⁴, easing by the end, with
the ramp topping out at +0.058 atm); with β_s = 0 the front exits slowly
with its wrinkle destroyed, then a remnant re-anchors at the boundary and
holds a +0.027 atm oscillating offset. Without the temperature closure the
boundary's diffusive fluxes are wrong through the whole transit and the
σ = 0.25 relaxation cannot pay that debt back. eT is the fidelity flag for
the transit; β_s is the fidelity flag for the source.

**Post-exit anchoring is σ's job alone.** After the flame leaves, β_s has
nothing to act on and both β_s rows land on identical residuals (+133 at
σ = 1, +8.5 at σ = 16, −250 recovering at σ = 0.25): the emptied domain
returns to ambient at the relaxation rate, so the residual ranking is the
anchoring ranking, exactly as before.

**`pin_farfield` and `extrap_material` conclusions are unchanged** — neither
knob changed in Phase 0. `pin_farfield` anchors a through-flow duct to
p + ρcu (stream 460 cm/s slow, +17500 persistent); `extrap_material` remains
for fronts that stay, not fronts that leave. The old rows are in the git
history.
