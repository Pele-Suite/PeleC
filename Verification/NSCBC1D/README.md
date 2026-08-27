# nscbc1d — standalone verification of `Source/NSCBC.H`

This driver compiles the production NSCBC kernel **unmodified** against a
minimal 1-D finite-volume Euler solver. Its purpose is to separate two failure
modes that are otherwise indistinguishable in a full PeleC run:

> *"the boundary condition is wrong"* vs *"the AMReX plumbing around it is wrong"*

It is cheap enough to sweep parameters, so it also produces the reference curves
that the AMReX-side regression tests are checked against.

## Build and run

```sh
cmake -S . -B build -DAMReX_DIR=<amrex-install>/lib/cmake/AMReX
cmake --build build
./build/nscbc1d          # run all checks
./build/nscbc1d sweep    # also dump the sigma sweep
```

Any AMReX build with `AMReX_SPACEDIM=1` works; MPI, OpenMP, EB and particles are
all unnecessary. The default mechanism is `air` (2 species, Fuego); override
with `-DPELE_MECHANISM=`. Build it a second time against a reacting mechanism
to exercise check C7, which is skipped otherwise:

```sh
cmake -S . -B build_lidryer -DAMReX_DIR=... -DPELE_MECHANISM=LiDryer
cmake --build build_lidryer && ./build_lidryer/nscbc1d      # 65/65 (air: 61/61)
```

Two further build axes, both run in CI (the `NSCBC-Driver` job):

```sh
cmake -S . -B build_adv -DAMReX_DIR=... -DPELE_NUM_ADV=2    # 61/61: pack_ghost's
cmake --build build_adv && ./build_adv/nscbc1d              # passive-scalar path
cmake -S . -B build_srk -DAMReX_DIR=... -DPELE_MECHANISM=LiDryer -DPELE_EOS=SRK
cmake --build build_srk && ./build_srk/nscbc1d              # 45/45 static checks
```

The SRK build is the EOS-portability claim made checkable: the kernel keeps
to density-carrying EOS entry points precisely so it runs under a real gas.
Under SRK the dynamic checks (C4/C5/C9b/C10/C11/C12) are skipped — they
integrate the mini solver for thousands of steps, every one paying several
Newton solves per cell, to re-verify algebra that is EOS-independent and
already gated under Fuego — and two gates change meaning: C1's tolerance sits
at the Newton round-trip floor (~1e-11) rather than machine epsilon, and C7
gates FD-vs-FD agreement instead of FD-vs-closed-form convergence, because
under SRK the kernel path *is* the finite difference.

Nothing in the driver may assume a particular mechanism. The first version of
`air_Y()` returned "0.233 for O2, else 0.767", which is right for the
two-species `air` mechanism and gives a composition summing to 6.6 for
LiDryer's nine — every check downstream then failed for reasons having nothing
to do with the boundary condition.

## Units

PeleC and PelePhysics work in **CGS**: cm, g, s, K, dyn/cm², erg/g.
`Constants::PATM = 1.01325e6`. The first version of this driver was written in
SI and every sound speed was 100× too large, which silently turned the
supersonic-outflow check into a subsonic one — that is exactly the class of
error this driver exists to catch, and it is why the checks assert on physical
relationships rather than on hard-coded numbers.

## What is checked

| | Check | What a failure means |
|---|---|---|
| **C1** | A uniform state is reproduced exactly in every ghost layer, at every σ, for both inflow and outflow | The kernel is manufacturing a gradient out of nothing; everything downstream is meaningless |
| **C2** | Every relaxation moves the boundary *toward* its target, on both the lo and hi face | A sign error. This is the assertion the legacy Fortran lacked, and its absence is why users had to be told "`relax_T` must be negative" |
| **C3** | Outflow *extrapolates* composition rather than imposing it; inflow imposes it exactly; `Σ Y = 1` and `UEDEN = UEINT + KE` hold to round-off | Species over-specification at outflow (the legacy defect), or a broken state identity |
| **C4** | Acoustic reflection of a pressure pulse: below 1% at σ=0.25, essentially zero at σ=0, and 2nd-order extrapolation no worse than 1st | The extrapolation or the invariant algebra is wrong |
| **C5** | The relaxation is a **rate**: grid-independent, and equal to `K = σ(1−M²)c/L` | The parameterisation has drifted to a value-blend, whose effective rate scales as `c/Δx` and therefore doubles when the mesh does |
| **C6** | Supersonic, reversed-flow and EB-body-state fallbacks each return a finite physical state and increment their counter; a supersonic INFLOW without `Target.p` is a counted substitution (`target_incomplete`) and with it an exact full-state imposition | A silent fallback — i.e. a bug that will not be found in production; or a supersonic inflow quietly borrowing pressure from a domain it should be causally upstream of |
| **C7** | The closed-form reaction source `dp/dt|_react` matches a directional finite difference along the reaction path, converging as τ→0; chemistry conserves mass; a cold state gives zero | The thermodynamics of `reaction_dpdt()` is wrong. Skipped automatically when the mechanism has no reactions |
| **C8** | On a flame-like temperature ramp the entropy closure's ghost overstates the face temperature gradient by 40%; `extrap_temperature` reproduces the ramp exactly and still returns a uniform state to round-off | The diffusion operator reads these ghosts, so this is the conductive heat flux leaving the domain being wrong |
| **C9** | (a) Extrapolating `R₊` across a normal velocity gradient manufactures ghost pressure `½ρc·ℓ·δu`, exactly, and order 1 gives exactly zero. With `extrap_material`, on the mass-conserving form of the same ramp, the bias vanishes while the ghost keeps the full `du/dn`. (b) A heat band at the boundary produces a σ-suppressed offset that the order control shows is *not* the extrapolation — reported, not gated | The ghost-pressure bias mechanism, isolated. (b) failing to isolate it dynamically is why C10 and C11 exist |
| **C10** | The source-free ramp: mass/momentum-consistent, no sustainer. Its own negative result — a source-free expansion cannot persist in a duct — plus a reported row showing `extrap_material` holds a *decaying* ramp alive at the face, which is its known cost | Nothing; the gated content moved to C11 |
| **C11** | The sustained ramp: C10's structure plus the manufactured energy source `S_E = ṁ dH/dx` that makes it an exact steady solution straddling the outflow — a flame's mechanical structure minus the chemistry. An *oracle* ghost fill (exact continuation) holds it, so the architecture is sound; the entropy closure drifts 15707 dyn/cm² in 0.7 relaxation times and distorts the face `du/dn` to 175% of exact; `extrap_material` holds those to 3175 and 80%, and cuts the static face-flux error 3.3× | The material-slope continuation is broken, or the late-time columns are being read without their caveat: the frozen source cannot follow a structure the boundary lets slip, so late-time drift is the MMS's artefact, not the boundary's |
| **C12** | With real conduction in the mini solver and a hot flank in the outflow cells, against a shielded reference: the entropy closure leaks 887 dyn/cm² of boundary error, `extrap_temperature` holds it to 104 | The diffusive boundary physics lives in the ghost **T closure**, not in the wave model — an amplitude-side diffusion source term was built, verified exact on quadratic profiles, measured to double-count (104 → −911 here; +1200 → +1771 in PeleC), and removed |
| **C13** | (a) The outflow closure is a *continuous* function of the interior normal velocity through `u_out = 0` — a swept flame-like stencil shows no outlier jump in ghost p, u or T at the crossing; (b) an over-pressured transient reversal still relaxes: ghost p toward target AND ghost `u_out` pushed outward; (c) a firm reversal does NOT extrapolate material content: with T falling toward the face, the ghost keeps the interior's own T | The reversal handling has drifted from either of its two proven properties. The 2026-08 NSCBC-Chamber production A/B (`Docs/NSCBC-reversal-branch-defect.md`) caught both failure modes in sequence: a dedicated reversal branch that dropped `dR₊`, `S_p` and `T_in` and froze the ghost velocity put a 4.2e3 dyn/cm² step at `u_out = 0` (vs 1.2e-2 sweep variation, fails a/b) and fed a growing dither and a spurious 0.3 atm spike; unifying the closure *without* upwinding the material slopes (fails c, ghost T 341.5 vs 400 K) instead fed the cold runaway — extrapolated outward-cooling ghosts advected back in, 385 → 89 K in 42 µs, NaN. The shipped closure passes all three: unified acoustics, `w_mat`-upwinded material slopes |
| **C14** | A duct held in steady INFLOW through an outflow face (lo face relaxes to a low-pressure sink, hi face is an outflow whose Target also carries the 300 K reservoir state; interior starts at 600 K): the frozen-material closure keeps the hi half at 597 K — the domain feeds on its own exhaust — while `backflow_material` flushes it to 299.5 K on the advective clock; a static probe confirms the Mach-[10⁻³,10⁻²] ramp is bit-inert at breathing amplitudes | Sustained recirculation is being fed recycled interior gas (the flag off is that, by design — but the FLUSH failing with the flag on means the reservoir ramp is broken), or the ramp has crept down into breathing amplitudes where the chamber's transient physics must stay frozen |

C4 also measures the **inflow** reflection curve: R = 2.3% / 4.8% / 19% / 57% at `relax_u` = 0.5 / 2 / 10 / 50 — soft inlets swallow acoustics, stiff ones are walls; the default reflects under 5%. And the kernel now carries a **transit guard**: an advisory counter (`material structure`) that fires when |dS| > 5% of ρ per cell sits in an outflow boundary cell — the configuration whose crossing the σ = 0.25 default does not survive.

The C11 oracle row is the load-bearing negative control: it separates "the ghost-cell *form* cannot do this" (false — the oracle holds the front indefinitely) from "this particular *closure* cannot" (true for the entropy closure, mostly fixed by `extrap_material`).

## C11x — the profile-fit experiment (reported, not gated)

C11 also carries an experiment on the question "if you own a 1-D profile of
the front, can the boundary use it?" The closure is given the profile *family*
(tanh between end states — the analogue of owning an unstretched flamelet) but
not its position or thickness; both are fitted per fill, statelessly, by
inverting T at the last two interior cells through the family (a closed-form
value-and-slope match). Ladder rows at σ = 1, ⟨p⟩ error at 0.7 relaxation
times / at t_end:

| ghost closure | 0.7 τ | t_end | what it supplies |
|---|---|---|---|
| entropy | 25007 | 143274 | nothing beyond the algebra |
| `extrap_temperature` | 26090 | 144414 | linear T |
| **fit** (profile T only, ρ from EOS at kernel p) | 25268 | 116569 | fitted-profile T |
| `extrap_material` | 3175 | 89132 | linear T *and* u |
| **fitU** (profile T and u; p stays relaxed) | **66.8** | **79.8** | fitted-profile structure |
| fitUX (same, family end-state 15% wrong) | 73.8 | 86.9 | robustness probe |
| oracle | −38.1 | −54.7 | the exact answer, placed exactly |

(The drifting rows' numbers grew when the outflow-reversal fallback became a
soft, σ-rated relaxation instead of a hard ambient pin — the old pin was
silently braking the frozen-source walk-off through the transient reversals
these badly-anchored runs experience. The equilibration-window values of
every closure that actually holds the front, and every gate, are unchanged
to the digit; the late-time columns sit in the artifact region either way.)

Three findings. **Material-only profile information buys nothing here** — fit
= `extrap_temperature` = entropy to 0.3%, because this MMS is inviscid and its
boundary error was never in the material content (C12 is where T-content
pays). **The fitted profile supplying T and u recovers 97% of the
`extrap_material` → oracle gap** — and, unlike `extrap_material`, it does not
walk off with the frozen source at late time (79.8 vs 27210 at t_end): the
per-fill re-fit re-locks the structure to the family and cuts the mismatch
feedback loop. **The fit is insensitive to a wrong family**: a 15% error in
the assumed end state costs 9%, because the value-and-slope match absorbs the
leading-order deformation — the stretch/curvature argument in miniature. The
pressure never comes from the profile in any row; it stays the relaxation's.

The release side is measured on C10's decaying ramp (the structure a correct
boundary must let die; `extrap_material`'s known failure). At t_end, ⟨p⟩
error: plain kernel **+1477**, `extrap_material` **−20406** (holds the ramp
alive), unbounded profile-fitU **−3499**. The stateless re-fit releases only
*partially* — the two-parameter family can shift and widen but cannot
represent a shrinking amplitude, so during the decay it keeps imposing
full-amplitude structure through weakened data.

The repair is the **source-consistency bound** (`fitB`): a steady front obeys
du/dn = (dp/dt)|src / ρc² (the Sutherland–Kennedy relation behind β_s), so
the continuation's amplitude is blended toward the plain kernel by
w = min(1, du/dn_sustainable / du/dn_measured), with the sustainable
dilatation computed from the measured local source. Measured, both horns:

| | sustained front (C11), 0.7 τ / t_end | decaying ramp (C10), t_end |
|---|---|---|
| plain kernel | 15707 / 27129 | **+1477** |
| `extrap_material` | 3175 / 27210 | −20406 |
| fitU, unbounded | **66.8 / 79.8** | −3499 |
| **fitU + source bound** | **66.8 / 79.8** | **+1477** |
| oracle | −38.1 / −54.7 | — |

On the sustained front the bound is inert to the printed digit — the
manufactured source sustains du/dn = 448 against the actual 449, so w ≈
0.998 — and on the source-free ramp it refuses the continuation outright and
reproduces the plain kernel exactly. One stateless closure now passes both
qualifications. (Here the closure is handed the manufactured source exactly;
a PeleC version would assemble dp/dt|src from `reaction_dpdt` plus the
diffusive term, and inherits their coverage and their gaps.)

The shape axis is also measured: a second sustained-front block replaces the
truth with a Richards curve (k = 3 — asymmetric, outside any tanh; the
manufactured source and the oracle follow it automatically) while the fit
still assumes tanh. At 0.7 τ / t_end: entropy 14632 / 86673, **fitU 54.9 /
71.7**, oracle −15.0 / −21.7 — the tanh fit through a non-tanh truth retains
**100%** of the recovery, and the source bound stays inert. The reason is
geometric: the ghosts extend 4 cells past the boundary while the front is ~20
cells wide, so any smooth monotone saturating family matched locally in value
and slope agrees with the truth to second order over the overhang. The
library's global shape barely matters; what carries the closure is (i)
monotone saturating structure with bounded end states, (ii) the local
value-and-slope match, refreshed statelessly, and (iii) the source gate.
That is the design statement for a 2-D/3-D version: it does not need the PMF
profile per se — it needs a one-parameter monotone family with measurable end
states.

Still untested: multi-species fitting on a progress variable, and the
interaction with β_s = 0, which feeds the same S_p into the incoming wave.

## Reference results

Measured with `air`, `n = 400`, `L = 10 cm`, `order = 2`, a 0.1% Gaussian
pressure pulse, integrated for 1.6 acoustic transit times. `R` is the peak
residual wave amplitude in the upstream half after the pulse has left, measured
against the *instantaneous domain mean* so that the σ-driven anchoring
adjustment is not miscounted as a reflected wave.

| σ | R [%] | mean p drift [dyn/cm²] | τ_relax [s] |
|---|---|---|---|
| 0.00 | 0.00078 | −0.0076 | ∞ |
| 0.05 | 0.160 | −1.68 | 5.75e−3 |
| 0.10 | 0.315 | −3.30 | 2.87e−3 |
| 0.15 | 0.466 | −4.88 | 1.92e−3 |
| 0.20 | 0.614 | −6.44 | 1.44e−3 |
| **0.25** | **0.758** | **−7.96** | **1.15e−3** |
| 0.30 | 0.899 | −9.44 | 9.58e−4 |
| 0.50 | 1.43 | −15.0 | 5.75e−4 |
| 1.00 | 2.56 | −26.8 | 2.87e−4 |
| 2.00 | 4.18 | −43.5 | 1.44e−4 |
| 4.00 | 7.20 | −60.3 | 7.19e−5 |
| 8.00 | 14.97 | −69.3 | 3.59e−5 |
| 16.00 | 28.14 | −70.8 | 1.80e−5 |
| — (`pin_farfield`) | 0.015 | +0.010 | n/a (value pin) |

The sweep runs past σ = 2 because `Exec/RegTests/NSCBC-FlameOutflow` needs
σ ≈ 10 to anchor an outflow with a flame crossing it, and the price of that has
to be quotable: 20-30% reflection, and a mean drift that has saturated — beyond
σ ≈ 8 the anchoring stops improving while the reflection keeps growing, so
σ > 16 buys nothing at all.

Three things to read out of that table.

**The reflection/anchoring trade-off is real and monotone.** `R` grows very
nearly linearly in σ while the pressure anchoring strengthens in step. σ = 0 is
perfectly non-reflecting and completely unanchored. This is the curve that makes
"σ = 0.25 is a good default" a measured statement rather than a received one.

**The relaxation is a genuine rate.** Doubling and quadrupling the resolution
changes the measured `K` by 3% (1925 → 1861 s⁻¹ from n=200 to n=800), and the
measured value sits within 7% of `σc/L`. A boundary condition parameterised as a
*value blend* instead has an effective rate of
`c/Δx`, which doubles when the mesh does, so its σ is not transferable and its
behaviour is not grid-converged. Keeping this check green is what keeps
literature σ values meaningful in PeleC.

**`pin_farfield` is not simply "σ = ∞".** The hard far-field pin is *both*
nearly non-reflecting (R = 0.015%) *and* anchored (drift +0.01 vs −7.96 at
σ=0.25) — because it constrains the incoming invariant's *value* rather than its
gradient, and a value constraint on the incoming characteristic alone reflects
nothing. Its cost is that it anchors to `p_target + ρc·u_out` rather than to
`p_target`, and that it does not converge under refinement. It is the right
choice for an open boundary onto a large quiescent reservoir and the wrong one
for a duct exhausting into a plenum whose true mean pressure is not `p_target`.

## Duct modes — the injection-fidelity bed (`t1` / `t3` / `t5`)

```sh
./nscbc1d t1              # step-change convergence vs relax_u (Dupuy)
./nscbc1d t3              # forced-inlet deterioration index (Daviller)
./nscbc1d t5 [relax_u=2]  # standing-wave P_RMS(x) profiles (Daviller)
```

One duct, run only when asked: NSCBC value-relaxation inflow at the lo face,
hard pressure ghost (FOExtrap, p pinned — fully reflecting, R ≈ −1) at the
hi face, mean inflow 2×10³ cm/s, forcing amplitude 10⁻³ c. These are Part
III of the design doc made runnable — the first *dynamic* measurements of
what the value-relaxation inlet does to injected signals, against the
Dupuy/Daviller phenomenology. Each mode carries its own relationship gates
and exits nonzero on failure; none is part of the default suite.

**t1 — step response (Dupuy Figs 5a/8a).** Inlet target stepped at t = 0
against the reflecting far end; t_conv = last time the domain-mean velocity
sits outside 0.1% of the step, in units of t_a = 2L/c, 40 t_a window:

| relax_u | 0.1 | 0.2 | **0.3** | 0.5 | 1 | 2 | 5 | 10–100 |
|---|---|---|---|---|---|---|---|---|
| t_conv/t_a | 31.3 | 13.5 | **6.06** | 6.84 | 11.2 | 20.3 | 39.9 | >40 (never, in-window) |

The CLR curve exactly as their DDE analysis predicts: an interior optimum
at 0.3, drift-limited convergence below it, reflection-delayed convergence
above, unusable by 5. Gated: interior minimum, argmin ∈ [0.1, 1].

**t3 — forced-inlet deterioration (Daviller §4–5).** Harmonic target
u^t = u₀ + A sin(2πft) at 0.8/1.0/1.2 × the duct's quarter-wave f₀
(Doppler-corrected, 866.7 Hz here). I_u = achieved u′ amplitude at the
boundary cell over A (1 = faithful injection); I_in = incoming-invariant
amplitude over the ideal injector's; the stiff (velocity-Dirichlet) duct
gain 1/|1+e^{iθ}| is the reference:

| relax_u | I_u at 0.8 f₀ | I_u at f₀ | I_u at 1.2 f₀ |
|---|---|---|---|
| 0 | 0.011 | 0.013 | 0.016 |
| 0.5 | 0.137 | 0.013 | 0.073 |
| 2 | 0.949 | 0.010 | 0.244 |
| 5 | 2.709 | 0.146 | 0.455 |

Three findings, all design-doc I.3e made quantitative. **relax_u = 0
injects nothing** (gated): the GC inflow has no amplitude slot — a
time-varying target enters only through the relaxation, so the injection
vanishes with the stiffness. This is the exact opposite of NRI's I ≡ 1 at
any K. **Off-resonance injection is monotone in relax_u but not faithful
anywhere** (gated on the monotonicity): 14% of the target at relax_u = 0.5,
95% at 2 — and 271% at 5, overshooting, because the duct's return wave and
the relaxation's response interfere and the closed loop has no notion of
the difference between "target not met" and "reflection arriving".
**On resonance the injection collapses instead of over-driving** (I_u ≈
0.01 at every stiffness while I_in stays order-1): the quarter-wave mode
has a velocity node at the driven end, and a velocity relaxation cannot
impose a signal at a velocity node — the incoming wave it launches returns
inverted and cancels it. An inlet that must both hold a mean and inject a
signal through one relaxation coefficient does neither faithfully near a
resonance; that is the measured cost of statelessness at the inlet, and the
NRI/register discussion (design doc II.2, queue item 4) is the literature's
answer to exactly this table.

**t3 also carries the queue-item-4 prototypes** (register + feed-forward
GATED at I_in ∈ [0.85, 1.15] across the matrix): the same matrix under
three inlet models, with the register held by the *driver loop* — updated
once per step outside the stages, read frozen — which is the
boundary-registers architecture with the kernel untouched. I_in
(Daviller's index, incoming amplitude over the ideal injector's):

| inlet model | I_in over the 9-point matrix |
|---|---|
| classical | 0.07 – 4.7 |
| NRI register alone | 0.22 – 0.95 |
| `Target::dudt` feed-forward alone | 0.85 – 5.24 |
| **register + feed-forward** | **0.89 – 1.08, incl. on resonance** |

The register-only NRI transplant fails structurally (the value form has no
amplitude slot to protect — I.3e measured); `dudt` is that slot, stateless
and now in the kernel; the register's real job is the *learned reference
mean* that lets the relaxation ignore the returning wave. Together they
put I_in within ~10% of unity at every stiffness and frequency, and the
residual velocity pattern is the duct's own `|1+e^{iθ}|/2` to three
digits. Full design: `Docs/NSCBC-boundary-registers-design.md`.

Two register footnotes. `./nscbc1d ndnr` is the outlet-side twin (gated):
relaxing toward the register's EMA mean instead of the fixed far-field
collapses the σ = 16 planar-pulse reflection 28.14% → 3.56% while keeping
the anchor. And the register invariant is referenced to the *frozen
ambient* ρc on purpose: `NSCBC1D_LOCAL_RC=1` re-runs t3 with the
instantaneous local impedance instead, degrading register+FF to
I_in = 1.61 at relax_u = 2 — the mean p/(ρc) times the coherent impedance
oscillation aliases into R₋ at signal amplitude. Kept as a runnable
forensic; PeleC carries the same reference as a slow-EMA'd `ema_rhoc`
(`Exec/RegTests/NSCBC-Acoustic/README.md` for its half of both tables).

**t5 — standing-wave pattern (Daviller §7).** Same runs; P_RMS(x) at 17
stations against the analytic envelope |sin(k_eff(L−x))|. Shape correlation
is 1.000 at all three frequencies (gated off-resonance at > 0.95): nodes
and antinodes sit exactly where the duct puts them, at every stiffness —
the *geometry* is never the problem. The amplitudes carry t3's story:
antinode P_RMS = 2819 / 1260 / 812 dyn/cm² at 0.8/1.0/1.2 f₀ (relax_u = 2),
i.e. the strongest response is NEAR resonance and the weakest ON it,
because the inlet cannot couple into the mode it is nominally driving
hardest. (Ideal injected amplitude for scale: ρcA ≈ 1.4×10³, doubled at a
standing-wave antinode.)

## Adding a check

Checks are plain functions that call `check(bool, name, detail)`. Prefer
assertions on *physical relationships* (monotonicity, grid-independence,
conservation identities, direction of relaxation) over assertions on numbers:
the numbers move when the mechanism, the mesh or the interior scheme changes,
and the relationships do not.
