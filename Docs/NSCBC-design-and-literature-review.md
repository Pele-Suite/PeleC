# The ghost-cell NSCBC: formulation from the ground up, its relation to the Poinsot–Lele flux form, and a test program drawn from the PAPERS folder

This document does three things. Part I describes the ghost-cell (GC) characteristic
boundary treatment adopted in `Source/NSCBC.H` from first principles — what is
computed, from which cells, with which differences — and sets it against the original
Poinsot & Lele (1992) flux-based specification, including the weaknesses this choice
creates, both the ones measured in this repository and the ones the literature
identifies. Parts II–IV review the papers in `PAPERS/` (Dupuy, Meziat Ramirez,
Douasbin & Poinsot, *JCP* 2026 — the NDNR method; Daviller, Oztarlik & Poinsot,
*Comput. Fluids* 2019 — the NRI-NSCBC inlet; Quillâtre *et al.* MCS7 2011,
Quillâtre *et al.* *IECR* 2013 and Vermorel, Quillâtre & Poinsot *Comb. Flame* 2017 —
the Sydney/SydGex vented-explosion program) for settings and concerns that bear on
this implementation, a catalogue of test problems we can set up, and a catalogue of
tests in those papers that we deliberately cannot set up, with the reasons.

Throughout, "our σ" is `pelec.bc_nscbc_sigma` with rate `K = σ(1−M²)c/L`; the CERFACS
papers use the reduced coefficient σ = KL/c. The two conventions coincide up to the
(1−M²) factor, which is 1 to a few parts in 10⁴ in every case discussed here, so
literature σ values transfer directly. This was a design goal (`Params::L_ref` exists
only so σ keeps its literature meaning) and the papers confirm it was the right one.

---

## Part I — The ghost-cell formulation, from the ground up

### I.1 What Poinsot & Lele actually specify

The original NSCBC (Poinsot & Lele, *JCP* 101:104, 1992, built on Thompson 1987) is a
**flux/amplitude-form, time-integrated** boundary condition. The Navier–Stokes
equations at the boundary *node* are recast in characteristic form. Along the
boundary normal the convective terms are written as wave *amplitude variations*

```
L1 = (u−c)(∂p/∂n − ρc ∂u/∂n)      incoming at a subsonic outflow
L2 = u (c²∂ρ/∂n − ∂p/∂n)          entropy
L3 = u ∂v/∂n,  L4 = u ∂w/∂n       tangential velocities
L5 = (u+c)(∂p/∂n + ρc ∂u/∂n)      outgoing
```

and the primitive-variable evolution at the boundary node reads (LODI form)

```
∂p/∂t + ½(L5 + L1) = 0
∂u/∂t + (1/2ρc)(L5 − L1) = 0        (+ entropy/tangential relations)
```

The method's contract has four parts:

1. **Outgoing amplitudes are measured.** Every `L` whose characteristic speed points
   out of the domain is evaluated with **one-sided normal differences of the resolved
   interior field at the boundary node**, once per Runge–Kutta stage, inside the PDE
   update.
2. **Incoming amplitudes are modeled.** At a subsonic outflow only `L1` enters, and
   the Rudy–Strikwerda relaxation `L1 = K(p − p_t)` (the "Classical Linear
   Relaxation", CLR, in Dupuy et al.'s nomenclature) replaces the unknowable exterior
   derivative. Later refinements add transverse terms (Yoo–Im; Lodato), source terms,
   and low-Mach-consistent forcing amplitudes for injection (Daviller et al.).
3. **The boundary state is prognostic.** The modeled and measured amplitudes are
   assembled into `∂U/∂t` at the boundary node, which is *integrated in time* with the
   interior scheme. The boundary therefore has *memory through the state itself*: what
   it is now depends on the entire history of amplitudes it was fed.
4. **Viscous conditions are a separate, explicit specification.** Because the
   amplitudes are inviscid constructs, Poinsot & Lele require additional conditions on
   the normal derivatives of the viscous fluxes at an open boundary (e.g. ∂q_n/∂n = 0,
   ∂τ_tn/∂n = 0), imposed on the diffusion operator directly.

Where the derivatives live, in summary: **all differences are formed at the boundary
node** — one-sided along the normal for outgoing amplitudes, centred/one-sided along
tangents for transverse terms — and there is a **time derivative at the heart of the
method**: the BC is an ODE on the boundary state.

### I.2 What our ghost-cell form does instead

`Source/NSCBC.H` never touches the boundary node's evolution equation. It fills
**ghost cells** during every `FillPatch`, as a **pure algebraic function of the
instantaneous interior state**, and then steps aside: PeleC's ordinary machinery — the
Godunov/PLM reconstruction and Riemann solve for the hyperbolic flux, and the diffusion
operator's face-gradient stencils — reads those ghost cells exactly as it reads
interior data. The boundary condition *is* the ghost state; the boundary flux is
whatever the interior scheme computes from it.

The fill works in an outward-normal frame (`u_out = n_sgn·u_n`) with linearised
characteristic **invariants** rather than amplitude variations:

```
R± = u_out ± p/(ρc)          acoustic pair          (ρc FROZEN at boundary cell N)
S  = ρ − p/c²                linearised entropy      (c frozen likewise)
```

Freezing the impedance is not an optimisation: `R+` is only an invariant of the system
linearised about a fixed state, so slopes of `R+` are meaningless unless every stencil
cell uses the same ρc.

**Where every difference is formed** (this is the complete list; there are no others):

* **Outgoing slopes.** `δR+`, `δS` (or `δT` under `extrap_temperature`), `δY_k`,
  `δu_t` are minmod-limited one-sided differences on the three interior cells
  `N, N−1, N−2` walking inward from the boundary, with tangential indices clamped into
  the domain ∩ FAB so only valid data is ever read. The ghost value at layer ℓ is
  `φ_N + ℓ·δφ` (`order=2`) or `φ_N` (`order=1`).
* **The incoming invariant.** `R−` receives the Poinsot–Lele relaxation *converted to
  a spatial gradient*: `dR−/dn = L_in/((c−u_out)ρc)` with
  `L_in = K(p_N − p_t) − (1−β)T_in − (1−β_s)S_react`, applied per ghost layer as
  `R−_g = R−_N + ℓ·Δn·dR−/dn`. Under `extrap_material`, `R−` additionally carries the
  material part of its own measured minmod slope, bounded through the entropy family
  (`du_mat = −u·δS/(ρ(1−M²))`) so nothing the incoming model consumes rides an
  incoming family.
* **Transverse terms.** `T_in` uses centred differences of `p`, `u_out`, `u_t` across
  the two tangential neighbours of cell `N` (spacing `2Δy`, degrading to one-sided or
  vanishing at clamped corners) — the ghost-cell transplantation of Motheau's T1/T4.
* **Point values, no derivatives.** The reaction source `S_react = dp/dt|_chem` is an
  EOS/chemistry evaluation at cell `N` alone. The EOS closes the fill
  (`(p_g, S_g or T_g, Y_g) → (ρ_g, e_g)`).
* **No time differences anywhere.** The fill has no state, no history, no ∂/∂t. It is
  idempotent under the repeated `FillPatch` calls of an SDC iteration — and, measured
  in this repository, bit-identical under box re-decomposition, MPI rank count, and
  checkpoint/restart.

### I.3 The fundamental differences, and what each one buys and costs

**(a) An ODE on the boundary state versus an algebraic map to ghost states.**
Poinsot–Lele integrate modeled dynamics *at* the boundary; we paint a modeled *spatial
structure* next to it and let the Riemann solver decide what crosses the face. The
relaxation "rate" K survives the translation — the per-layer increment is consumed
once per Δt by the scheme, so τ = 1/K is mesh-independent (check C5) — but its
*authority* does not: the anchoring increment per layer carries an explicit 1/n_x
that the outgoing slopes do not. Balancing the two gives the equilibrium offset

```
Δp ≈ ρ c L_ref (∂u_out/∂n)|_b / σ
```

which is the ghost-pressure bias of check C9(a), measured exactly in the driver and
confirmed (σ⁻¹ trend, grid-converged) in `NSCBC-FlameOutflow`. **This bias has no
analogue in the flux form** — there, a dilatational normal gradient enters the
*outgoing* amplitudes, which are measured, not extrapolated. It is the one structural
error the GC form adds, and this repository's measurements bound its importance: it is
real, exactly reproducible, and *not* the dominant term in a reacting front-crossing
outflow (removing it exactly, via `extrap_material`, buys 5% at σ=1).

**(b) Who computes the boundary flux.** In the flux form the modeler assembles it; in
the GC form the interior scheme does, from reconstructed face states. We get nonlinear
consistency, positivity handling, and species conservation for free, and the same code
path serves 1/2/3-D, all faces, EB cut cells and any mechanism. The cost is that the
model's intent is filtered through the reconstruction: the fill's carefully built
slopes pass through the scheme's own limiting before they become a flux. `order = 2`
is load-bearing for exactly this reason (first order clamps the front's structure and
flips the sign of the front-crossing error), and any future change to the hydro
reconstruction near `ext_dir` faces silently changes the boundary condition.

**(c) Statelessness.** No memory means: idempotence under SDC re-fills; restart
determinism (measured bit-identical); no per-point storage to allocate, communicate,
or checkpoint; trivially correct under AMR regrid and GPU relaunch. But it forecloses,
*by construction*, every boundary condition in the modern CERFACS line that needs a
time integral or filter at the boundary point: NRI's outgoing-velocity integral
`u− = (1/2ρc)∫L1 dt` (Daviller 2019, Eq. 23/37) and NDNR's exponential-moving-average
separation of mean-flow transients from acoustics (Dupuy 2026, Eq. 12/A.1) are both
*states*. This is the single most consequential architectural difference exposed by
the papers in `PAPERS/`, and Part II returns to it, because NDNR solves precisely the
σ trade-off this repository has measured most painfully.

**(d) Viscous conditions become ghost closures.** The flux form must specify viscous
conditions separately; the GC form cannot avoid specifying them — PeleC's diffusion
operator forms boundary-face conductive and species fluxes *from these same ghost
cells*, so whatever normal gradients the ghost carries **are** the diffusive fluxes.
This cuts both ways. Done wrong, it is a silent error no amplitude analysis will find:
the historical entropy closure overstates the face temperature gradient by 40% on a
flame-like ramp (C8) and leaks 887 dyn/cm² of boundary error under real conduction
(C12). Done right (`extrap_temperature`), the diffusive physics is carried correctly
with no amplitude-side term at all — and an amplitude-side "viscous condition"
transplanted into the GC form *double-counts*, measured: +104 → −911 in the driver,
+1200 → +1771 in PeleC. The flux form's viscous-condition machinery is therefore not a
missing feature here; it is the solution to a problem the GC form does not have,
replaced by a problem it does have (choosing the ghost closure), which is now gated.

**(e) Relaxation semantics for unsteady targets.** In the flux form, *injecting* a
signal and *holding* a mean are different operations with different amplitudes — an
acoustic injection enters as `−2ρc ∂u_a^t/∂t`, a vortical one as `−ρc ∂u_v^t/∂t`
(the factor-2 distinction of Prosser and Guézennec–Poinsot that Daviller's NRI-NSCBC
formalises), while the relaxation handles only the mean. Our GC inflow has no
amplitude slot at all: a time-varying target is imposed as a *value* through the same
relaxation that holds the mean. For slowly varying targets this is fine; for forcing
at acoustic frequencies (flame transfer functions, acoustic response studies) the
injected amplitude will be wrong by a frequency-dependent factor set by `relax_u` —
the measured inflow curve (R = 2.3/4.8/19/57% at relax_u = 0.5/2/10/50) is the
zero-injection limit of this transfer function, but the injection side has never been
measured. Tests T3 and T5 below are designed to measure it.

**(f) Multi-dimensional composition.** The GC form composes across faces, edges and
corners by per-point ownership plus clamped stencils, and the 3-D radial test shows
the composition is sound. Flux-form implementations need explicit corner-compatible
amplitude assembly (Lodato et al. 2008). Advantage: ghost cells.

### I.4 Known weaknesses of the GC form, consolidated

Measured here: the C9(a) extrapolation bias (I.3a; bounded, partially removable);
anchoring authority ∝ σ/n_x against slope terms, hence σ = O(10) at fronts with its
28% acoustic reflection price; the reconstruction filter (I.3b); the frozen-ρc
linearisation, weakest exactly where the boundary crosses an impedance jump (a flame);
`extrap_material`'s quasi-steady bound (worst case at U = 2 S_L, off by 2×; unstable
in a fast transit at small σ, where the transit guard now fires); level-local stencils
under AMR (warned, measured small at normal incidence); EB body-state detectability
(now mandatory `eb_zero_body_state`). From the literature: the CLR σ trade-off itself
(Dupuy Fig. 8: viable band σ ∈ [0.1, 0.3] with drift immediately below it — our sweep
reproduces both ends), and everything in I.3c/I.3e that statelessness forecloses.

---

## Part II — Settings, issues and concerns from the papers

### II.1 Dupuy, Meziat Ramirez, Douasbin, Poinsot (JCP 562:114987, 2026) — NDNR

The paper's diagnosis of CLR is a theory for what this repository measured
empirically. Their DDE analysis shows CLR's viable σ band is narrow (optimum ≈ 0.27–
0.3 for a step response, drift immediately below, reflection-delayed convergence
above), that convergence times scale on the acoustic delay t_a = 2L/c, and that in
complex cases users fall back to trial and error. Our σ sweep (driver + flame tables)
is the same curve seen from the reflection side; our flame-exit finding — σ ≈ 16 or
nothing — sits squarely inside their recommended NDNR band σ ∈ [10, L/(cΔt)], which
they can afford *because the EMA filter removes the acoustics from the relaxation's
deviation before σ acts on it*. In other words: **the literature has converged on the
anchoring strength our exit test demanded, and has a mechanism that removes the 28%
reflection we pay for it.** That mechanism is one EMA register per boundary point
(τ∞ ≈ 3 t_a; f_c ≈ f₀/10) applied to `p+` at outlets and `u−` at inlets — with
`p+ = −½∫L+dt`, i.e. two pieces of per-point time-integrated state.

Concerns to carry: (i) the flame-exit event in their Appendix B channel was an
*artifact of CLR's acoustic activity* — the flame left at 50 ms under CLR and was
still inside at 140 ms under NDNR. Boundary-generated acoustics can *move flames*;
our exit test measured the front being pushed *backwards* by the σ = 0.25 default,
the same coin's other face. Any PeleC vented-flame result obtained with strong
low-frequency boundary activity should be suspected of the same contamination.
(ii) Their t_a-based prescriptions assume the user can estimate the acoustic delay;
in a chamber whose sound speed field evolves (their B, our flame cases) σ-tuning by
formula fails — the argument for a broad-viable-band method, not for better tuning.

Implication for us: a **boundary-registers architecture** (a per-face MultiFab of
time-integrated `p+`/`u−` and their EMA, updated exactly once per level advance,
*read* by the stateless fill) would preserve FillPatch idempotence — the fill stays
algebraic given the registers; only the once-per-step register update carries state.
Restart requires checkpointing the registers. This is the single highest-leverage
capability the papers point at, and it is *additive*: CLR remains the registers-off
limit.

### II.2 Daviller, Oztarlik, Poinsot (Comput. Fluids 190:503, 2019) — NRI-NSCBC

Three results matter here. First, the injection amplitudes: acoustic forcing requires
the factor 2 (`−2ρc ∂u_a^t/∂t`), vortical injection the factor 1, and using either for
the other is measurably wrong — a distinction our value-relaxation inlet does not
currently express (I.3e). Second, the relaxation deviation should exclude the
outgoing-wave velocity `u−`, evaluated per point as `(1/2ρc)∫L1 dt` — with it, the
inlet is non-reflecting at *any* K (their deterioration index I ≡ 1), without it,
I reaches >20 at σ = 5 near duct resonances. Third, their Fig. 3/8 "operating zone"
diagrams formalise the same CLR compromise as Dupuy.

Concerns to carry: our inlet is the classical NSCBC column of their Table 1, so every
deterioration they quantify applies to us verbatim wherever an inlet faces a resonant
domain (their Fig. 18: NSCBC excited the duct's 1/4-, 3/4-, 5/4-wave modes hard enough
to pollute a flame-transfer-function measurement). PeleC users doing forced-response
or thermoacoustic work with `bc_nscbc` inlets will hit this. The TurbInflow path is
the vortical-injection case: `relax_u` low-pass filters the injected spectrum (already
documented), and their factor-1 analysis says the *correct* injection amplitude for
turbulence is gentler than for acoustics — worth measuring (T6) before anyone tunes
`relax_u` upward to "pass more turbulence" and turns the inlet into a wall.

### II.3 Quillâtre et al. (MCS7 2011; IECR 2013) and Vermorel, Quillâtre, Poinsot (C&F 183:207, 2017) — the Sydney/SydGex vented-explosion program

These are the application the boundary condition ultimately serves (the accompanying
press piece in `PAPERS/` makes the program's stakes explicit). The settled numerical
practice, in all three papers: **the NSCBC never sits at the vent.** The domain is
extended by a plenum "mimicking the atmosphere," the mesh coarsens toward its far
side, and NSCBC is applied on the plenum borders — so the flame, and the violent
venting flow, cross an *interior* plane, and the characteristic boundary only ever
sees the plenum's mild far field. Our flame-exit measurements are precisely the
quantification of what that practice avoids: the shipped σ = 0.25 dies on a transit,
and even the working recipe (σ = 16 + `extrap_temperature`) trades a 0.18%-of-ambient
transient for it.

Settings and sensitivities worth recording: the safety QoI (peak overpressure and its
timing) is governed by the *competition between combustion rate and venting rate* —
i.e. by the mass flux through the vent plane, which sits directly downstream of
whatever boundary treatment terminates the domain. Their post-peak pressure traces
ring at the first longitudinal duct mode (≈650 Hz at small scale) "considering a
non-zero impedance at its exit": the boundary's impedance is *in the measured
observable*. And the laminar-phase findings (two-step chemistry for T_ad; realistic
Lewis numbers for curvature response; errors up to 220% from Le = 1) are not boundary
matters but calibrate how accurate the rest of a PeleC explosion setup must be before
boundary errors of a few hundred dyn/cm² are even visible.

---

## Part III — Test problems we can set up

Ordered roughly by value per unit cost. None of these needs to be a CI case; T1–T3
are driver-level, T4–T6 reuse existing PeleC cases or utilities, T7–T8 are new small
cases in the existing NSCBC-* mould.

**T1 — Step-change convergence time versus σ (Dupuy §3, Figs 5a/8a).** 1-D driver: a
duct with a *fully reflecting* far end (hard pressure ghost), inlet target stepped at
t₀; measure t_conv (last time |u−u^t| > 0.1% of the step) versus σ over [0.1, 100].
Gates: a minimum near σ ≈ 0.3 and growth on both sides, matching the CLR curve their
DDE predicts. This validates our relaxation *dynamics* (everything so far gates
statics, rates, or reflection), and it is the baseline any future NDNR work must beat.

**T2 — Wave-packet recovery (Dupuy §3.3).** Same duct, boundary stabilised at its
target, then struck by a broadband packet built around the quarter-wave mode f₀ =
c/4L. Measure recovery time versus σ: CLR's t_conv explodes for σ ≳ 3 (their Fig. 8b).
Confirms — dynamically — why "raise σ" is not free even before a flame is involved.

**T3 — Forced-inlet deterioration index (Daviller §4–5).** Duct with harmonic inlet
forcing u_a^t(t) (our hook already receives `time`; the target is simply made
oscillatory) and a *reflecting* outlet (hard p). The analytic index
I = |L̂5/L̂5^t| = 1/(1 + R₁e^{2ikL}) is known in closed form for the classical
relaxation inlet; measure I at three frequencies straddling a duct resonance for
σ ∈ {0, 2, 5}. This measures, for the first time, what our *value-relaxation* inlet
does to an injected signal — the GC-form answer to the factor-2 question of I.3e —
against an exact reference. Driver-level first; PeleC 1-D after.

**T4 — Nozzle start-up convergence (Daviller §6).** PeleC already has
`EB-ConvergingNozzle` with a characteristic-outflow precedent. Add the `bcnormal_nscbc`
hook, initialise at rest, and measure time-to-steady-state and trapped-mode ringing
versus (σ, relax_u), reproducing their Figs 9–11 phenomenology (slow convergence at
small K; trapped 150 Hz oscillations at large K). This is the "compressible start-up
transient" scenario every PeleC user actually faces.

**T5 — Standing-wave pattern under forcing (Daviller §7, laminar core).** Duct,
harmonic inlet, reflecting outlet: the analytic p′(x) (their Eq. 33) gives a
quantitative P_RMS(x) profile. Classical NSCBC over-predicts the anti-node amplitudes
(their Fig. 15); measure ours. Cheap (1-D/2-D, no turbulence needed) and the sharpest
available gate on injection-amplitude fidelity.

**T6 — Turbulence through the characteristic inlet (Daviller §7, full).** 3-D channel
with PeleC's TurbInflow through the `bc_nscbc` inlet plus optional harmonic forcing;
gates: injected velocity spectrum matches the target signal (their Fig. 14 right),
pressure spectrum free of spurious cavity-mode peaks (their Fig. 14 left / Fig. 18).
Directly exercises the `relax_u`-filters-the-spectrum warning already in our docs.

**T7 — Mini-SydGex: vent-with-plenum versus vent-with-boundary.** A 2-D Sydney-like
chamber (closed end, ignition kernel, one baffle optional, open end) run two ways:
(a) vent opening into a plenum whose *far* side carries the characteristic boundary —
the CERFACS practice; (b) the boundary placed directly at the vent plane with our exit
recipe (σ = 16 + `extrap_temperature`; transit guard expected to fire). The QoI is the
overpressure trace: (a) is the reference, (b)−(a) is the price of dropping the plenum.
This is the application-shaped version of `nscbc-flameexit.inp`, it reuses all of that
machinery, and it converts "use a plenum" from imported practice into a measured
in-house number — including whether the post-peak ringing (their 650 Hz observation)
is reproduced with the right decay for each treatment.

**T8 — Outwardly-propagating laminar flame (Quillâtre MCS7, DNS sub-case).** Quarter
circle, symmetry sides, characteristic far-field arc; spherical flame from a 1 cm
kernel. Two boundary-relevant gates: the consumption-speed–curvature line (Markstein
slope vs Clavin–Joulin) must be *independent of the arc radius* (boundary
independence of a slowly arriving spherical dilatation), and the domain must not drift
while the front approaches all of the boundary at once. Also the natural PeleC
verification of the laminar phase all three explosion papers identify as controlling
the overpressure.

---

## Part IV — Tests in these papers that do **not** apply to our implementation

| Test (paper) | What it tests there | Why it cannot apply here |
|---|---|---|
| NDNR parameter maps and convergence surfaces over (σ, τ∞/t_a) (Dupuy Figs 7–10, Table 2) | The EMA-filtered relaxation: convergence time and reflection of a boundary condition with a per-point low-pass state | Our fill is stateless by design (I.3c): there is no τ∞, no EMA register, no `p+`/`u−` integral to filter. Applicable only after a boundary-registers architecture exists; until then only the CLR halves of their figures apply (= T1/T2). |
| NDNR reacting channel A/B comparison (Dupuy App. B) | NDNR vs CLR at both ends of a turbulence-forced flame channel, incl. the flame-exit artifact | Same statelessness barrier for the NDNR half; the CLR half is a heavier cousin of T6 and adds nothing T6/T7 don't measure more cheaply. |
| NRI non-reflection at arbitrary K (Daviller Figs 9–11, I ≡ 1 property) | The `u− = (1/2ρc)∫L1dt` subtraction making the inlet transparent at any relaxation strength | The integral is per-point time-integrated state; no slot exists in an algebraic ghost fill. The *classical-NSCBC* baselines of the same figures are exactly T3/T4. |
| Isolated-vortex injection amplitude check (Prosser/Guézennec line, cited and built on in Daviller §3) | That vortical injection needs the factor-1 amplitude in `L5`, not the acoustic factor 2 | The check is on a wave-amplitude *forcing slot* our inflow does not have — we impose values through relaxation, not `∂u^t/∂t` amplitudes. The GC-form question it corresponds to ("does a relaxed value-inlet distort injected vorticity?") is T6's spectrum gate. |
| Sub-grid combustion-model discrimination on SydGex (Vermorel §4–5: Colin vs Charlette efficiency functions, β sensitivity, cross-scale constant fitting; IECR/MCS7 chemistry & Lewis studies) | TFLES efficiency closures and reduced-chemistry/Le choices against the overpressure database across three scales | Combustion-model validation, orthogonal to the boundary condition; PeleC has no TFLES efficiency-function pair to discriminate. The *geometry and QoI* are reusable (T7), but matching their peak-overpressure tables tests a turbulent-combustion model, not a BC. |
| Scale-up study SS→MS→LS (Vermorel §5, 1:6:24.4 replicas) | A priori transferability of model constants across a 10⁴ volume ratio | Requires the GexCon medium/large rigs' data and billion-cell-class LES; nothing boundary-specific is isolated at the larger scales that T7 does not already probe at small scale. |
| Impedance-imposing boundaries (D-TDIBC and complex-impedance work cited in Dupuy's introduction) | Prescribing a frequency-dependent reflection coefficient via time-domain convolution at the boundary | Convolution kernels are the heaviest form of per-point boundary state; excluded by the same architecture that excludes NDNR/NRI, and further out on the same road. Catalogued so the door is marked, not shut. |

A closing remark on the catalogue: every "cannot apply" above traces to one root —
statelessness (I.3c) — except the combustion-model items, which are out of scope for a
boundary condition under any formulation. That is a clean situation: one architectural
decision (boundary registers, updated once per advance, checkpointed, read by the
otherwise-unchanged algebraic fill) converts the entire first column from
"inapplicable" to "implementable", with CLR as its registers-off limit.

---

## Part V — Recommended order of work

1. **T5 + T3** (injection fidelity, exact references, cheap): they measure the one
   semantic gap (I.3e) no current check covers, and they price the classical inlet for
   the forced-response use cases Daviller shows it fails.
2. **T7** (mini-SydGex): converts the plenum practice and our exit recipe into one
   comparable overpressure number; the application-shaped capstone for the existing
   exit work.
3. **Boundary-registers design note** (enables NDNR at outlets first — our measured σ
   trade-off is worst there — then NRI at inlets): one MultiFab of per-point
   `∫L_out dt` and its EMA per characteristic face; update in `PeleC::advance` at the
   new-time fill only; checkpoint alongside the state; fill reads registers as
   additional inputs and remains algebraic. Dupuy's Table 2 then supplies the defaults
   (σ > 10, τ∞ = 3 t_a) that our flame-exit table independently corroborates.
4. **T1/T2/T4/T8** as the verification bed fills out; **T6** when TurbInflow-through-
   NSCBC first meets a real user.
