# NSCBC boundary registers — design note (2026-08-25)

**Status: design with a measured prototype; the stateless half is already
in the kernel.** This is work-queue item 4: the design that must exist
before any per-point boundary state enters PeleC. It was written the way
this branch writes everything — the central claims were measured first, in
the driver, and one of them failed in an instructive way.

## Why registers, ranked by measured need

1. **Inlet injection fidelity.** Driver t3 (and its PeleC twin in
   `NSCBC-Acoustic`, agreeing to 1–4%) measured what the value-relaxation
   inlet does to an injected signal: nothing at `relax_u = 0`, an
   unfaithful K-dependent fraction off resonance (0.14 / 0.95 / 2.7 of
   target at relax_u 0.5 / 2 / 5), and collapse on a duct resonance at any
   stiffness — a velocity relaxation cannot drive a velocity node. Anyone
   doing forced-response or thermoacoustic work through `bc_nscbc` hits
   this (design doc II.2: classical-NSCBC excited duct modes hard enough
   to pollute a flame-transfer-function measurement).
2. **Outlet reflection where σ = 16 is still chosen.** The flame closures
   removed the σ-ramp for fronts (eT = 1, β_s = 0 at any σ), but where
   strong anchoring is genuinely wanted the 28% reflection price stands.
   NDNR (Dupuy 2026) removes it with one EMA register per boundary point
   (τ∞ ≈ 3 t_a) that strips acoustics from the relaxation's deviation
   before σ acts.
3. **Trend gates.** C11x's source-consistency closure wants a transit/decay
   gate; C14's transient-vs-sustained backflow classification is provably
   unavailable from an instantaneous snapshot (the sit-vs-transit
   discussion in the FlameOutflow README). A slow per-face trend register
   is the stateful escape both name.

## The architectural invariants (non-negotiable)

* The fill remains a **pure function of (interior state, Params, Target,
  registers)**. No state lives in or is written by the fill. FillPatch and
  SDC re-fills stay idempotent because the registers are frozen during an
  advance.
* Registers update **exactly once per level advance**, after the new-time
  state exists, outside any fill path. The update is per-point local, so
  it is MPI-decomposition-independent by construction.
* Registers are **checkpointed**; restart is bit-identical, gated.
* **CLR is the registers-off limit**: every register consumer is additive
  and off by default; the entire measured suite is unchanged with
  registers absent.

## The prototype, and what it measured

The driver's t3 duct (NSCBC inflow, hard-p reflecting far end, harmonic
forcing at 0.8 / 1.0 / 1.2 × the quarter-wave f₀) ran four inlet models.
The register lives in the *driver loop* — updated once per step from the
completed state, read frozen by the stages — which is precisely the
architecture proposed for PeleC; the kernel stayed pure throughout.
`I_in` is the incoming-invariant amplitude over the ideal injector's
(Daviller's index); `I_u` is the total velocity response over the target.

| inlet model | I_in over the 9-point matrix (3 K × 3 f) |
|---|---|
| classical value relaxation | 0.07 – 4.7, K- and f-dependent, resonance-wild |
| NRI register alone | 0.22 – 0.95 |
| `dudt` feed-forward alone | 0.85 – 5.24 (overshoots with K) |
| **register + feed-forward** | **0.89 – 1.08, including on resonance** |

Three findings:

* **A register-only NRI transplant fails, structurally.** Daviller's
  I ≡ 1 belongs to the flux form, where injection has its own amplitude
  slot (−2ρc ∂u_a^t/∂t) and NRI merely stops the relaxation from fighting
  the outgoing wave. Our value form has no slot (design doc I.3e);
  subtracting the outgoing wave from the deviation only exposes the bare
  injection bandwidth limit (measured: I_u *worse* than classical at
  K = 2–5 off resonance).
* **The amplitude slot is stateless, and it is now in the kernel.**
  `Target::dudt` — the analytic time derivative of the normal-velocity
  target, supplied by the problem hook, entering the incoming amplitude as
  `L_in += 2 ρc n_sgn dudt` (the 2 is u = (R₊+R₋)/2: imposing a velocity
  rate through one invariant needs twice it). Zero recovers the classical
  inlet identically.
* **The register's true job is the reference mean.** The outgoing-wave
  velocity u₋ = ½(R₋ − R̄₋) is algebraic *given a reference* R̄₋ — but an
  inlet has no pressure target, so the reference must be **learned**: an
  EMA of R₋ = u − p/ρc at the boundary cell, τ = 3 t_a. With it, the
  relaxation's deviation excludes the returning wave; with `dudt` beside
  it, injection is exact. Together: I_in within ~10% of unity at every
  stiffness and every frequency measured, and the residual I_u pattern
  equals the duct's own standing response `|1 + e^{iθ}|/2` to three digits
  — the boundary has stopped being part of the acoustics.

## Proposed PeleC design

**Storage.** One face-band `MultiFab` per characteristic face (the
boundary-adjacent cell layer), components:
`{ema_Rm}` now; `{ema_pplus, trend_dS}` reserved for phases B/C. Defined
only on levels whose grids touch the face (the AMR warning already
constrains characteristic faces to be level-consistent).

**Update.** At the end of `PeleC::advance` on the owning level, after the
new-time state: `R̄ ← R̄ + (Δt/(τ+Δt))(R − R̄)` per boundary-adjacent
cell, τ = 3 t_a with t_a = 2 L_ref/c estimated from the register row
itself. One kernel launch per face per step; GPU-clean; no
communication.

**Read.** `BCfill` — not the kernel — composes the effective target:
`u_tgt_eff = u_tgt + ½(R₋ − R̄₋)` for a registered inlet, and passes
`Target` to the unchanged pure fill. The problem hook opts in per face
(`Target` flag or a `pelec.bc_nscbc_nri` list); `dudt` is orthogonal and
already available to any hook today.

**Checkpoint.** The register MultiFab is written and read with the
checkpoint (PeleC's derived-MF machinery); a restart-bit-identity gate
joins the determinism suite.

**Gates.** Driver: the t3 mode rows become gated (I_in ∈ [0.85, 1.15]
across the matrix for register+FF; the classical rows stay as the
baseline). PeleC: the `NSCBC-Acoustic` duct mode re-runs with `dudt` and
the register — the cross-check standard is the 1–4% driver agreement the
classical table already achieved. C1 (uniform state) and the full
existing suite must be untouched with registers off — structurally
guaranteed, still gated.

**Phasing.**
* **A (inlet NRI)** — the register MF + update + BCfill composition +
  gates above. `dudt` is already usable without any of it.
* **B (NDNR outlets)** — same storage, `ema_pplus`; the relaxation acts on
  the EMA-high-passed deviation; measured against the σ sweep (target:
  σ = 16-class anchoring at ≪ 28% reflection; Dupuy's broad-band claim).
* **C (trend gates)** — a slow `|dS|`-trend component for C11x's
  transit/decay gate and C14's sustained-backflow classification;
  advisory first, closure-coupled only after measurement.

## Implementation record (phases A–C landed)

What shipped deviates from the proposal in one structural way: the store is
a flat `Gpu::DeviceVector<Real>` indexed by `(face, flattened tangential
index, component)` — not a face-band MultiFab. The band is a 2-D skin of
the domain; a MultiFab would buy distributed ownership at the price of
FillPatch-adjacent machinery, and the skin is small enough (≲1 MB at
chamber resolution) to replicate. Components:
`{ema_Rout, ema_p, u_minus, dp_ac, trend_dS, ema_rhoc, init}`.

Two lessons were paid for in measurement:

* **The invariant needs a time-frozen impedance** (`ema_rhoc`). With the
  instantaneous local `rho c`, the mean `p/(ρc)` (~10⁴ cm/s) times the
  coherent impedance oscillation aliases into `R_out` at signal amplitude.
  Driver with `NSCBC1D_LOCAL_RC=1`: I_in degrades 1.03 → 1.61; PeleC
  measured 1.59 before the fix and 0.970 after. The kernel's per-fill
  frozen `rho_c` is immune (stencil difference at one instant); a register
  differences across time and is not.
* **`amrex::DefaultGeometry().Domain()` is not a reliable level test** in
  the fill path: the comparison silently failed, the fill never took the
  register pointer, and the composition was a no-op (I_in bit-identical to
  feed-forward-only). The store now records the domain it was sized for
  and `nscbc_registers_peek(dom)` returns null for any other box — the
  level test and the staleness test are the same comparison.

Gates, all green: driver t3 register+FF I_in ∈ [0.85, 1.15] across the
matrix; driver ndnr σ = 16 reflection 28.14% → 3.56%; PeleC duct NRI point
I_in = 0.970 (driver 1.031); PeleC planar σ = 16 NDNR 30.10% → 3.66% with
the mean held to 37 ppm; restart bit-identity with `NSCBCRegisters` in the
checkpoint. Measured tables: `Exec/RegTests/NSCBC-Acoustic/README.md`.

Determinism, measured not assumed: n1 vs n6 of the same MPI binary is
bit-identical on every state field with NDNR on — including a 400×64 grid
with `max_grid_size=32`, which splits the register faces *tangentially*
across ranks. That case was run because the zero-on-unowned-ranks ownership
invariant looked vulnerable to corner-strip ghost fills reading unowned
(zeroed) register entries; the vulnerability does not materialize, so the
simple ownership design stands. (Two earlier "failures" of this test were
artifacts worth remembering: the case's GNUmakefile builds `USE_MPI=FALSE`
by default, so `mpiexec -n 6` launched six serial clones racing over one
plotfile; and a multi-rank plotfile shards `Cell_D` across files, so
file-level `cmp` against a serial run is meaningless — compare fields, not
files.)

## Not proposed

Raw time integrals without EMA (unbounded state, drift); any state inside
the fill or written during FillPatch; register-driven changes to the
outflow flame closures (eT/β_s solved that problem without state, and the
chamber killed three reversal branches that tried to be clever —
statefulness enters only where statelessness measurably cannot go).
