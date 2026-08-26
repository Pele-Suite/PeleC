# NSCBC-FlameOutflow-DRM — the same problem, CH₄/air chemistry

Identical to `../NSCBC-FlameOutflow` in every respect except the mechanism and
the flame: CH₄/air at φ = 0.75 under `drm19` instead of H₂/air at φ = 0.4 under
`LiDryer`. The problem definition is *the parent case's* — `prob.H`, `prob.cpp`
and `prob_parm.H` here are two-line forwarding headers — so the two cannot drift
apart.

```sh
./PeleC2d.gnu.ex nscbc-flameoutflow-drm.inp prob.nscbc_inflow=0
```

| | H₂/air, φ = 0.4 | CH₄/air, φ = 0.75 |
|---|---|---|
| mechanism | LiDryer (9 species) | drm19 (21 species) |
| S_L | 22.80 cm/s | 25.37 cm/s |
| T_ad | 1426 K | 1925 K |
| thermal thickness δ_T | 0.069 cm | 0.053 cm |
| expansion ratio ρ_u/ρ_b | 4.4 | 6.42 |
| ∂u/∂n through the front | 1130 s⁻¹ | 2672 s⁻¹ |
| **Le of the deficient reactant** | **≈ 0.3** | **0.970** |

The flame was computed with Cantera from `drm19/mechanism.yaml`
(freely-propagating, mixture-averaged transport, 1 atm, 300 K) and written in
the PMF format the parent case reads.

## Why this exists

Two reasons, and the second is the interesting one.

### It settles the thermodiffusive objection

At φ = 0.4, H₂/air runs at Le ≈ 0.3, where an imposed wrinkle is
thermodiffusively unstable. A reader is entitled to ask whether the parent
case's boundary numbers are measuring the flame rather than the boundary. At
φ = 0.75, CH₄/air has Le(CH₄) = 0.970 — neutral — so here the question does not
arise at all.

For the record it does not arise in the H₂ case either, and the argument is a
timescale one. Tracking the front by the radical peak (H for H₂, CH₃ for CH₄):

| | H₂, half-amplitude | CH₄, half-amplitude |
|---|---|---|
| t = 0 | 0.07990 | 0.08007 |
| mid-run | 0.07896 | 0.07982 |
| end | 0.07809 | 0.07947 |

Both **decay**, and the unstable one decays *more*. The Darrieus–Landau
e-folding time at this wavenumber is ~9×10⁻³ s and the flame time δ/S_L is
~3×10⁻³ s, against measurement windows of 1.5–2.4×10⁻⁵ s — 130 to 390 times
shorter. The boundary equilibrates on the acoustic time and the flame is frozen
on that time, which is the same fact that makes the measurement possible at all.

### It is an independent test of the mechanism

The parent case attributes its error to a ghost-pressure bias driven by the
normal velocity gradient,

```
Δp  ≈  ρ c L_ref (∂u_out/∂n)|_b / σ
```

CH₄/air at φ = 0.75 is a different flame in exactly the way that formula cares
about: its expansion ratio is 6.42 against 4.4, so ∂u/∂n through the front is
2672 s⁻¹ against 1130, while its burnt-gas impedance ρc is 15.2 against 19.8.
The formula therefore **predicts a 1.82× larger error at the same σ**, before
the run is made.

Measured, at t = 1.5×10⁻⁵ s, mean Δp against each case's own shielded
reference (rows with the reaction source *off*, which the Phase-0 sign fix
does not touch; the β_s = 0 rows of the original table were measurements of
the inverted sign and are superseded below):

| outflow | σ | β | β_s | H₂ | CH₄ | ratio |
|---|---|---|---|---|---|---|
| hard `p = p_amb` | — | — | — | −507 | −776 | 1.53 |
| characteristic | 1 | 1 | 1 | +1917 | +3254 | **1.70** |
| characteristic | 1 | 0.5 | 1 | +1572 | +2428 | **1.54** |

**1.62 measured against 1.82 predicted**, with no fitted quantity — the
prediction is from the flame's expansion ratio and thickness alone. These
β_s = 1 rows are failure-mode configurations whose piled-up structure
reverses the boundary (in the H₂ sitting case, run with counters on, the
ENTIRE outflow face is in reversal essentially every fill), so their
absolute values have moved with each generation of the reversal closure —
pin, soft, and now the C13-gated unified closure
(`Docs/NSCBC-reversal-branch-defect.md`) — while the H₂/CH₄ RATIO, which is
the check, has stayed within 2%: 1.65 (pin), 1.64 (soft), 1.62 (unified).
The pre-fix CH₄ β = 0.5 value (+2253) had been reproduced to the digit
(+2253.3) across a different machine, toolchain and build system, so the
shifts are the closure, not the environment.

## What it says about β_s, more strongly than the H₂ case did

CH₄/air puts **3.5× more heat release in the boundary cells** than H₂/air does —
peak `q` along the outflow column is 2.1×10¹⁰ against 6.0×10⁹ erg/cm³/s. If the
reaction-source correction works, this is where it should show hardest.

It does. At σ = 1, β = 0.5, t = 1.5×10⁻⁵ s, mean Δp against the shielded
reference:

| configuration | CH₄ |
|---|---|
| β_s = 1 (source off) | +2428 |
| β_s = 0 (source on) | **+919** |
| β_s = 0 + `extrap_temperature` | **−546.8** |
| hard `p = p_amb` | −776 |

Turning the corrected source on cuts the error 2.6× (the H₂ case measures
2.0×), and the closure pair — the ghosts carrying the diffusive dp/dt, the
incoming wave carrying the chemical one — lands at the hard-outflow level,
exactly as in the parent case; the recipe row reproduces the soft-closure
era's −546.8 to the digit (it never reverses, and forward-flow fills are
bit-identical across the closure change). The stronger the heat release,
the more the term is worth, which is what a correctly-signed
Sutherland–Kennedy cancellation must do.

`β = 0.5` is worth more in CH₄ than in H₂ — 20% against 10% — which is
consistent: the stronger dilatation makes the tangential structure at the
boundary stronger too.

## Build note

The GNUmakefile in this directory shipped with the parent case's
`Chemistry_Model := LiDryer` until the 2026-08 survey — every measurement
above was made through CMake, which always set `drm19`, and the GNUmake path
had never been exercised. A LiDryer build of this case reads the 21-species
datafile into 9-species chemistry and produces a garbage initial state
(mean molecular weight ≈ 1) that runs stably enough to look like physics.
If a DRM number ever looks alien, check the plotfile's species list first.

## The boundary column

`radicals.py` (in the parent case) at t = 1.5×10⁻⁵ s, σ = 1:

```
  y [cm]    T [K]       Y_CH3         Y_H       q [erg/cm3/s]
  0.0031      301   2.233e-09   5.031e-17          1.014e+02
  0.1031     1178   2.621e-04   3.268e-06          2.963e+09
  0.1281     1641   1.572e-04   6.399e-05          1.574e+10
  0.2031     1792   4.659e-10   6.199e-05          9.584e+08
  0.2781     1576   3.858e-04   4.367e-05          2.067e+10
  0.3781      303   5.851e-09   5.356e-16          7.672e+02
```

Six orders of magnitude in `Y_CH3` and eight in the heat release, along a single
column of boundary cells, with the CH₃ peaks at the two designed crossings and
the H pool filling the burnt half between them. Use CH₃ rather than temperature
to locate the front: the temperature midpoint sits in the preheat zone, which is
wide and diffusion-controlled, while CH₃ marks the fuel-consumption layer, which
is 2–3 cells wide and is the thing that actually has to be inside the boundary
cells.

## Cost and resolution

drm19 is 21 species against LiDryer's 9, and the case runs about 3.4× slower per
step. The grid is the parent's, so δ_T/Δx is 8.5 here against 11 there. The
parent case's resolution study (n_x = 48, 96, 192 giving +531, +409, +371)
establishes that the boundary offset is grid-converging by n_x = 96, so this is
adequate for the boundary measurement, which is what the case is for. It is not
adequate for a claim about the flame itself.
