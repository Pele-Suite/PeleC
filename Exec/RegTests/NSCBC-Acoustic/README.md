# NSCBC-Acoustic

A right-running isentropic acoustic pulse in a quiescent gas, launched toward a
subsonic outflow at `x-hi`. `x-lo` is a wall and `y` is periodic, so the only
thing that can absorb the pulse is the outflow boundary. This is the AMReX-side
counterpart of check C4 in `Verification/NSCBC1D`: the same physical problem,
solved by a completely different code path, as a cross-check that the boundary
condition survives contact with the framework.

## Running it

```sh
./PeleC-NSCBC-Acoustic nscbc-acoustic.inp                      # characteristic
./PeleC-NSCBC-Acoustic nscbc-acoustic.inp pelec.bc_nscbc=0     # hard boundary
```

With `bc_nscbc = 0` the problem's `bcnormal` imposes the ambient pressure
directly in the ghost cells — deliberately crude, so the comparison is stark.
With `bc_nscbc = 1` the `bcnormal_nscbc` hook in `prob.H` returns an outflow
target and the characteristic treatment takes over. Note how little the hook
has to say: for a pure non-reflecting outflow the target pressure is the only
quantity that may be specified, so it is the only quantity it sets.

## Results

`L = 10 cm`, `c = 34719 cm/s`, so one acoustic transit is `2.88e-4 s`. Run to
`t = 4.6e-4 s` (1.6 transit times), by which point the pulse has crossed the
outflow. `R` is the peak residual pressure disturbance in the upstream half,
measured against the instantaneous domain mean so the σ-driven anchoring
adjustment is not miscounted as a reflected wave, divided by the incident
amplitude.

| boundary treatment | R [%] | mean p at end [dyn/cm²] |
|---|---|---|
| hard: ambient p imposed in the ghost cells | **97.2** | 1013178.3 |
| NSCBC, σ = 0.25 | **0.81** | 1013241.8 |

Target pressure is 1013250.0. The hard boundary reflects essentially the whole
pulse; the characteristic treatment reduces the reflection by a factor of 120
and holds the mean pressure to 8 parts per million of the target.

**Cross-check.** The standalone 1-D driver predicts `R = 0.758%` at σ = 0.25
(`Verification/NSCBC1D/README.md`). PeleC's 2-D Godunov solve gives 0.810%.
Two independent solvers, two independent implementations of the surrounding
machinery, agreeing to 7% on a quantity that spans two orders of magnitude
between the good and bad boundary conditions. That agreement is the point of
having both.

## Things worth trying

* `pelec.bc_nscbc_sigma = 0` — perfectly non-reflecting; R drops further, and
  the mean pressure is then unanchored and free to drift over a long run.
* `pelec.bc_nscbc_sigma = 2.0` — over-relaxed; the reflection grows roughly
  linearly in σ.
* `pelec.bc_nscbc_order = 1` — zeroth-order extrapolation of the outgoing
  invariant instead of minmod-limited linear.
* `pelec.bc_nscbc_pin_farfield = 1` — the value-pin formulation, which is both
  non-reflecting and anchored but does not converge under mesh refinement.
* Refine `amr.n_cell` and confirm R does not grow: the relaxation is a rate,
  not a per-cell value blend, so it is grid-converged.

## Three dimensions

The same sources build in 1-D, 2-D and 3-D — `AMREX_SPACEDIM`, `AMREX_D_DECL`,
`AMREX_D_TERM` and `AMREX_D_EXPR` throughout — so only the inputs file changes.
The 2-D planar result is unmoved: mean pressure at the end is 1013178.3 with a
hard boundary and 1013241.8 with the characteristic one, the same to the last
printed digit as before the conversion.

| inputs | what it is |
|---|---|
| `nscbc-acoustic.inp` | 2-D planar; the historical regression |
| `nscbc-acoustic-corner.inp` | 2-D planar with y WALLS: the wall/NSCBC corner seam |
| `nscbc-acoustic-3d.inp` | 3-D planar; confirms the kernel builds and runs in 3-D |
| `nscbc-acoustic-3d-radial.inp` | 3-D radial, **all six faces characteristic** |
| `nscbc-acoustic-duct.inp` | forced duct: the PeleC half of driver t3/t5 (below) |

### The wall/NSCBC corner (`nscbc-acoustic-corner.inp`)

The planar pulse again, but with `NoSlipWall` above and below instead of
periodic — so the characteristic outflow at x-hi meets a wall at two
corners, and the fill's tangential stencils clamp at the wall-adjacent rows.
The pulse is uniform in y and carries no v, so the exact solution stays
y-uniform through the crossing: any y-structure in the residual, and any v
anywhere, is corner-manufactured, with no reference solution needed.
Measured (2026-08-25, unified reversal closure): upstream-half residual
R = 0.752% — identical to the periodic-y baseline to the third digit — with
y-structure 0.0000% of the incident amplitude and max|v| 0.000% of the
incident velocity. The corner is invisible.

### Why the radial case exists

A planar pulse loads one face at normal incidence, which is what the 2-D run
already tests. It says nothing about a ghost cell that lies outside the domain
in two or three directions at once — and that is the part of `BCfill.cpp` with
no 1-D analogue, so the standalone driver cannot reach it either.

A radial pulse reaches the six faces at normal incidence, the twelve edges at
45°, and the eight corners along the body diagonal, in one run. Because the
initial condition is exactly isotropic about the box centre, **any** departure
from sphericity in the departing front is boundary-generated; no reference
solution is needed to say so.

`metrics.py sphericity` reports the spread in the front radius over ~4000
directions, binned by χ, the angle to the nearest face normal: χ = 0 is a face
centre, 45° an edge, 54.7° a corner. Rays whose front has already left are
dropped, so as time goes on the surviving rays are exactly the oblique ones.

### What it measures

`c = 34719 cm/s`, box 5 cm, 96³. The front reaches the faces at 7.2×10⁻⁵ s, the
edges at 1.02×10⁻⁴ and the corners at 1.25×10⁻⁴.

| t [s] | rays | χ range | radius spread % | amplitude spread % |
|---|---|---|---|---|
| | | | NSCBC / hard | NSCBC / hard |
| 3.0×10⁻⁵ | 4000 | all | 7.355 / 7.355 | 1.782 / 1.782 |
| 6.0×10⁻⁵ | 4000 | all | 2.891 / 2.891 | 5.158 / 5.157 |
| 7.5×10⁻⁵ | 2260 | edges + corners | 1.964 / 2.452 | **4.70 / 7.34** |
| 9.0×10⁻⁵ | 555 | corners | 1.291 / 1.612 | **7.69 / 40.29** |
| 1.05×10⁻⁴ | 25 | corners | 0.600 / 0.523 | **2.25 / 22.14** |

Read the first two rows first: before the front reaches any face the two runs
are **identical to six figures**, as they must be, since nothing has touched the
boundary yet. That is the metric's own sanity check, and it is why the later
rows can be believed.

After the face crossing they separate, and the discriminator is the **amplitude**
spread, not the radius: the front arrives at the right time either way, but a
hard boundary corrupts its strength. At 9×10⁻⁵ s — corner-bound rays only, the
face and edge parts of the wave already gone — the amplitude around the
surviving arc varies by 40% with a hard boundary and 7.7% with the
characteristic one.

So corner and edge ownership works. The `apply()` algebra is dimension-agnostic
by construction, and this says the plumbing around it is too.

## AMR and EB variants

`nscbc-acoustic-amr.inp` runs the planar pulse from inside a 2× refined patch
kept away from the outflow: mean pressure and upstream residual match the
single-level run (1013241.8 / 0.751% vs 0.752%). PeleC now warns — once per
face — when a refined level touches a Hard/UserBC face with `bc_nscbc = 1`,
because the fill's stencil is level-local and a fine patch on a
characteristic face makes the boundary condition level-dependent. Measured
here the on-face artefact is small (0.744%, +0.1 dyn/cm²); nothing
guarantees that at higher ratios or oblique incidence.

`nscbc-acoustic-eb.inp` seats an EB solid inside the fill's stencil at the
outflow. It found two things: PeleC's *default* body state is a sampled
fluid state the fill cannot detect (the counter read zero with a solid in
the stencil — the silent fallback the counters exist to prevent), so
`bc_nscbc` with EB geometry now *requires* `pelec.eb_zero_body_state = 1`
and aborts otherwise; and a body cutting the domain face itself NaNs under
the characteristic *and* the hard boundary — a pre-existing PeleC
EB-at-domain-boundary limitation. With the flag set, the run counts ~56000
body-state stencil degradations and completes cleanly.

### Two things this case turned up

**The flow-reversal path is not hypothetical.** Running with
`pelec.sum_interval > 0`, the counters read zero everywhere until the pulse
leaves and then report `flow reversal 11456` in a 40-step window. That is
correct physics, not a bug: the rarefaction behind an outgoing spherical wave
pulls the pressure below ambient and draws gas back in through faces configured
as outflows. Under the unified reversal closure (no dedicated branch:
counted, then the same restoring relaxation as forward flow with the
material slopes upwinded off — `Docs/NSCBC-reversal-branch-defect.md`,
driver gate C13) the run is stable through it. Before this case, that path
had only ever been exercised by a synthetic state in check C6.

**Under the unified closure the reversal no longer caps the benefit.** The
pin-era closure made every reversed cell effectively a pressure Dirichlet,
and the residual advantage over a hard boundary collapsed to ~1.5×.
Re-measured 2026-08-25 with the unified closure (`metrics.py residual`,
t = 1.35×10⁻⁴): NSCBC max|dp| 0.0118 / L2 0.00060 of the incident amplitude
against the hard boundary's 0.0627 / 0.0147 — 5× and 24×, with the mean
pressure held to 0.3 ppm of target (1013249.7 vs 1013250). The boundary
keeps venting *through* the reversal instead of walling it, and it does so
while sustaining the full physical re-entry: 15.8 million reversal fills
over the run (the pin-era report of ~11k per window was the pin choking the
breathing off early, not less reversal happening). The sphericity table
above is bit-identical under the closure change — re-measured, every row to
the digit, both variants — because the front rays it tracks are causally
ahead of anything the reversed faces emit. The residual is still not purely
boundary error here (the imploding half of the split pulse is in the box at
the final time), and for flows that hold an outflow in *sustained*
recirculation the local-inflow treatment of work-queue item 2 remains the
right escalation.

## Duct mode — the PeleC half of driver t3/t5 (`nscbc-acoustic-duct.inp`)

```sh
./PeleC2d.<comp>.ex nscbc-acoustic-duct.inp \
    pelec.bc_nscbc_relax_u=<K> prob.force_freq=<f> \
    stop_time=<6 t_a + 4/f> amr.plot_per=<1/(24 f)>
./duct_metrics.py <rundir> --freq <f>
```

`pulse_type = 2`: a uniform stream enters through a characteristic INFLOW at
x-lo whose velocity target carries a harmonic forcing (the `time` argument of
`bcnormal_nscbc`, finally consumed), and x-hi falls back to the case's
hard-pressure `bcnormal` — the fully reflecting far end. This is the same
duct as the driver's `t1`/`t3`/`t5` modes (`Verification/NSCBC1D/README.md`:
what the value-relaxation inlet does to injected signals), solved by PeleC's
Godunov machinery instead of the mini solver. The banner prints the
Doppler-corrected quarter-wave `f0` (865.1 Hz here); the matrix is relax_u ∈
{0, 0.5, 2, 5} × f ∈ {0.8, 1.0, 1.2} f0. `duct_metrics.py` linearly detrends
the inlet series before projecting — with relax_u small the mean is weakly
anchored, and its wander otherwise leaks into a finite-window Fourier
projection as a spurious amplitude (measured: the leak read as I_u up to 0.73
before detrending, on runs whose P_RMS field is zero).

Measured I_u (achieved u′ amplitude at the inlet over the target amplitude),
PeleC against driver:

| relax_u | f/f0 | PeleC | driver |
|---|---|---|---|
| 0 | 0.8 / 1.0 / 1.2 | 0.000 / 0.000 / 0.000 | 0.011 / 0.013 / 0.016 |
| 0.5 | 0.8 / 1.0 / 1.2 | 0.135 / 0.001 / 0.073 | 0.137 / 0.013 / 0.073 |
| 2 | 0.8 / 1.0 / 1.2 | 0.912 / 0.008 / 0.241 | 0.949 / 0.010 / 0.244 |
| 5 | 0.8 / 1.0 / 1.2 | 2.722 / 0.141 / 0.441 | 2.709 / 0.146 / 0.455 |

And the t5 antinode amplitudes at relax_u = 2: P_RMS = 2823.5 / 1264.3 /
813.4 dyn/cm² against the driver's 2819.4 / 1260.5 / 812.5 — 0.1–0.3%
agreement, with the standing-wave envelope correlation at 1.000 in every
forced run. Two independent solvers and two independent implementations of
the surrounding machinery agree to 1–4% on the deterioration index (several
entries to the third digit) and to a few tenths of a percent on the
standing-wave amplitudes. Every driver conclusion survives the plumbing:
relax_u = 0 injects nothing, off-resonance injection is monotone in relax_u
but faithful nowhere (2.7× overshoot at relax_u = 5), and on resonance the
injection collapses at any stiffness because a velocity relaxation cannot
drive the velocity node it is aimed at. The deterioration is a property of
the formulation, not of either code.

## Boundary registers (NRI and NDNR) — the PeleC half of driver t3-register/ndnr

The stateful repairs from `Docs/NSCBC-boundary-registers-design.md`, wired
into PeleC (`pelec.bc_nscbc_nri=1` for inlets, `pelec.bc_nscbc_ndnr=1` for
outlets; registers update once per level-0 advance in
`PeleC::nscbc_update_registers`, and `BCfill` composes the effective target —
the kernel stays a pure function).

**NRI** (inlet reflection cancellation), duct point relax_u = 2, f = 0.8 f0,
with the hook's feed-forward (`prob.force_feedforward=1`):

| treatment | I_in | I_u | driver |
|---|---|---|---|
| feed-forward alone | 2.419 | 1.493 | 2.40 |
| feed-forward + register (`bc_nscbc_nri=1`) | **0.970** | 0.593 | 1.031 / 0.660 |

**NDNR** (EMA-mean relaxation at outlets), the planar σ = 16 pulse of the
main table:

| treatment | R [%] | mean p at end |
|---|---|---|
| σ = 16, classical | 30.10 | 1013178.2 |
| σ = 16 + `bc_nscbc_ndnr=1` | **3.66** | 1013212.8 |

The driver's collapse is 28.14% → 3.56%; PeleC reproduces the 8× reduction
while holding the mean to 37 ppm of the 1013250 target (the classical row
holds 71 ppm — NDNR anchors *better*, because the relaxation no longer
spends its authority fighting the acoustic wave).

**A measured implementation lesson.** The register invariant must be
referenced to a *time-frozen* impedance (the store carries a slow-EMA
`rho c` per point). Computing `R₋ = u − p/(ρc)` with the instantaneous
local impedance aliases the mean `p/(ρc)` (~10⁴ cm/s) times the coherent
0.1% impedance oscillation into the invariant — tens of cm/s, the same
order as the signal. Measured identically in both codes: local impedance
degrades I_in from 1.03 to 1.61 (driver, `NSCBC1D_LOCAL_RC=1`) and to 1.59
(PeleC). The kernel's per-fill frozen `rho_c` never had this problem
because it differences across the stencil at one instant; the register
differences across time.

**Determinism.** Straight-through vs checkpoint+restart (σ = 16 NDNR, split
at step 400 of 800) is bit-identical; the registers ride in the checkpoint
as `NSCBCRegisters`. Decomposition: n1 vs n6 of the same MPI binary is
bit-identical on every state field, including a 400×64 grid whose
`max_grid_size=32` splits the register faces tangentially across ranks.
(Note the case GNUmakefile defaults to `USE_MPI=FALSE`; build with
`USE_MPI=TRUE` before believing any `mpiexec` result, and compare
plotfiles field-by-field — a multi-rank run shards `Cell_D`.)
