# NSCBC-COVO — two multi-dimensional tests of the characteristic boundary

One executable, two problems, selected with `prob.probtype`. Both are here to
answer questions the 1-D tests cannot: what happens when a wave meets the
boundary *obliquely*, and what happens when the thing crossing the boundary is
not a wave at all.

```sh
./PeleC-NSCBC-COVO nscbc-covo.inp        # probtype 0, convected vortex
./PeleC-NSCBC-COVO nscbc-pulse2d.inp     # probtype 1, circular pulse
# add pelec.bc_nscbc=0 to either for the hard-boundary reference
./PeleC-NSCBC-COVO nscbc-wrapgate.inp    # the periodic-wrap gate; zero steps
```

A third input file lives here that is not a physics run. `nscbc-wrapgate.inp`
takes zero steps and exists only to give `PeleC::nscbc_check_periodic_wrap()`
a configuration it can discriminate on: a box spanning the periodic direction,
structure sitting on the periodic seam, and `bc_nscbc_beta < 1` so the
transverse stencil is assembled too. Expect a line like

```
  NSCBC periodic-seam check: dir 0 lo, periodic tangential dir 1 -- N image
  pairs agree to 2e-4 (boundary-row density spread S)
```

with a small mismatch — the clamped tangential stencil's measured residual
under amrex's corner-strip protocol, O(1e-4) here and up to 1.6e-3 with a
flame front on the seam corner — and a spread that is *not* zero. The gate
aborts only above 1e-2: a broken stencil (a naive wrap measures 2.2e-2 here)
fails outright, while the protocol residual passes and
is *reported*. A zero spread means the check passed vacuously and gated
nothing, which is what `nscbc-covo.inp` and `nscbc-acoustic.inp` do: their
states are uniform along the periodic direction near the seam. The shipped
file gates the inflow face; a second run with `prob.centre = 1.0 0.0` puts
the vortex on the outflow corner and gates that face. The stencil-design
history — clamp vs wrap vs strip-consistent band, and why simplicity chose
the clamp — is in the comment above the stencil in `Source/BCfill.cpp`.

---

## probtype 1 — circular Gaussian pulse

After Motheau et al. (2017). A quiescent square, all four faces non-reflecting,
a Gaussian pressure bump at the centre. The outgoing circular wave reaches every
face and every corner at once. The exact wavefront is a circle of radius `ct`,
so **any departure from circularity is boundary error** — and unlike an
amplitude metric it cannot be hidden by a uniform pressure offset.

Two numbers, both measured over the polar rays whose front is still inside the
domain (i.e. the corner-directed ones, once the face-normal parts have exited):
the spread in front **radius**, and the spread in front **amplitude**, around
the arc.

| | radius spread | amplitude spread |
|---|---|---|
| hard boundary | 0.98% | **21.4%** |
| NSCBC, β = 1 (transverse off) | 0.47% | 0.90% |
| NSCBC, β = 0.8 | 0.40% | **0.066%** |
| NSCBC, β = 0.5 | 0.40% | 1.08% |
| NSCBC, β = 0 (full transverse) | 0.36% | 2.95% |

The amplitude column is the striking one. The hard boundary does not so much
bend the front as *chew* it: by the time the corner-directed arc is all that is
left, its strength varies by a factor of 1.2 around the arc, because the parts
that passed near the faces have already been corrupted by reflections. The
characteristic treatment holds it to 0.9% without transverse terms and to
**0.066% with them at β = 0.8** — a further factor of 14.

The radius column saturates: 0.4% is roughly one cell, so all the β < 1 rows
are "as circular as this mesh can represent". Amplitude uniformity is the
discriminating measure here, not radius.

The picture says it faster than the table. With the hard boundary the
rarefaction ring squares off, kinks at the corners, and grows four bright
reflection lobes; with the NSCBC it stays round and simply leaves.

`W` is the domain half-width. The 1.28% "before contact" figure is the
Cartesian-mesh discretisation error for a circular front and is the floor for
this metric, not a boundary effect — note that it is measured over all 720 rays
while the crossing-time rows use only the ~52 corner-directed ones, so the rows
are comparable to each other but not to the first row.

---

## probtype 0 — convected isentropic vortex (COVO)

A Yee/Shu vortex carried at M = 0.2 from a subsonic inflow at x-lo to a subsonic
outflow at x-hi, y periodic. A vortex is a pure vorticity/entropy structure: it
carries essentially no acoustic content, so **an ideal boundary would let it
leave in silence** and any pressure signal left in the domain afterwards was
manufactured by the boundary. It is also the only test in this suite that
exercises the **inflow** path end to end.

This is the test that shows the current implementation's limit honestly.

| configuration | rms/δp | rms vs floor |
|---|---|---|
| **no boundary at all** (periodic; discretisation only) | **0.0047** | 1× |
| hard inflow + hard outflow | 0.0652 | 14.0× |
| NSCBC, β = 1 (transverse off) | 0.0468 | 10.0× |
| NSCBC, β = 0.8 | 0.0341 | 7.3× |
| NSCBC, **β = 0.5** | **0.0176** | **3.8×** |
| NSCBC, β = 0.35 | 0.0187 | 4.0× |
| NSCBC, β = 0.2 (constant) | 0.0325 | 6.9× |
| NSCBC, β < 0 (pointwise local Mach) | 0.0754 | 16.1× |
| NSCBC, β = 0 (full transverse) | 0.6184 | 132× |

(The Phase-0 re-sweep, run after the reaction-source and periodic-seam fixes,
reproduces every β row above to the digit — the inert 2-D baselines did not
move — and separates two things an earlier row conflated: a **constant**
β = 0.2 measures 6.9×, better than β = 1, while the legacy **pointwise
local-Mach** convention measures 16.1×. The convention's fault is not its
value but its wobble: the vortex's own swirl modulates M and therefore the
correction, point by point, as the structure crosses.)

Dial sensitivity at β = 1, for reference: σ = 0 gives 15×, `relax_u` = 0.2 gives
14×, and hard-inflow + NSCBC-outflow gives 14× — i.e. none of them matters
compared to β.

`δp` is the vortex's own pressure deficit (1.90e4 dyn/cm², ~1.9% of ambient).
The periodic row is obtained by comparing the final field with the initial field
shifted by `u·t`; it includes the shift-interpolation error, so the true floor
is lower still and the ratios below are conservative.

Three things follow, and they are the point of the test.

**Without transverse terms the treatment helps by only 1.4×**, and the residual
stays an order of magnitude above the no-boundary floor — against 120× for a
*planar* pulse at normal incidence in NSCBC-Acoustic. That gap is the whole
content of "1-D-normal LODI", and none of σ, `relax_u` or the choice of which
boundary is characteristic closes it. When knobs span a decade and the answer
does not move, the error is not in those knobs.

**Turning the transverse terms on closes most of it.** β = 0.5 brings the
residual to 3.8× the floor: 2.7× better than β = 1 and 3.7× better than a hard
boundary. That is the term the dial sweep was pointing at.

**But the response to β is sharply asymmetric.** Too little correction costs a
factor of a few; too much is catastrophic — β = 0 is close to unstable, at
132× the floor — and a *pointwise-varying* β (the local-Mach convention) is
worse than no transverse terms at all. That asymmetry is why the shipped
default is β = 1 rather than the optimum. A wrong β is far more dangerous than
an absent one.

**The optimum is predictable, not empirical.** For a plane wave meeting the
boundary at angle θ, the correction that exactly cancels the obliqueness error
in the 1-D characteristic decomposition is `(1−β) = cos θ / (1 + cos θ)`, so

    β_opt = 1 − cos θ / (1 + cos θ)

which is 0.5 at normal incidence and rises toward 1 at grazing incidence (0.67
at 60°, 0.79 at 75°). Both measured optima land on that curve: the vortex — a
broad, near-normal-incidence structure — optimises at **0.5**, and the circular
pulse, whose surviving arc crossed the faces near-tangentially, optimises at
**0.8**. Two independent problems, one formula, no fitting.

**The outflow was never the reflecting part.** Switching only the outflow to
the characteristic treatment cuts the upstream residual (mean |δp| over
x < 3 cm) from 1648 to 342 dyn/cm² — a factor of 4.8. It was not sending a
reflected wave back upstream; it was failing to absorb the vortex cleanly,
which is a different fault, and β is its fix.

Use `prob.nscbc_inflow = 0` to leave the inflow to the ordinary `bcnormal` and
isolate the outflow, which is how the third row above was obtained.

---

## Verification note

The control that makes all of this trustworthy: with `prob.vortex_strength = 0`
— uniform flow, nothing crossing anything — both boundary conditions produce
**exactly zero** pressure disturbance to machine precision, over 3000 steps.
That is the AMReX-side confirmation of check C1 in `Verification/NSCBC1D`, and
it means every number above is a response to the flow rather than noise the
boundary generates on its own.

## Reproducing the numbers

`Verification/NSCBCFields/fielddump` flattens one variable of a single-level 2-D
plotfile onto a regular array (`fextract` gives 1-D slices only, and the
quantity of interest here is a shape); `metrics.py` computes the metrics from
those dumps. Build it with whichever system you used for PeleC — set `DIM` and
`COMP` to match:

```sh
cd Verification/NSCBCFields
make DIM=2 COMP=llvm -j                                       # GNUmake
# or: cmake -S . -B build -DAMReX_DIR=<amrex>/lib/cmake/AMReX && cmake --build build
```

The runs. Each stops at step 3000, which is `t = 1.348e-3 s` — the vortex left
at `t ~ 1.08e-3`, so `plt03000` is the "after it has gone" state every number
in the vortex table refers to.

```sh
cd Exec/RegTests/NSCBC-COVO
EXE=./PeleC2d.llvm.ex          # or whatever your build produced

mkdir b05 && (cd b05 && $EXE ../nscbc-covo.inp pelec.bc_nscbc_beta=0.5)
mkdir b10 && (cd b10 && $EXE ../nscbc-covo.inp pelec.bc_nscbc_beta=1.0)
mkdir hard && (cd hard && $EXE ../nscbc-covo.inp pelec.bc_nscbc=0)
```

That trio already gives the headline comparison — 0.0176 / 0.0468 / 0.0652 —
without needing any reference:

```sh
F=../../Verification/NSCBCFields/fielddump.gnu.ex
for d in b05 b10 hard; do $F $d/plt03000 pressure $d.dat; done
python3 ../../Verification/NSCBCFields/metrics.py residual 1013250.0 b05/plt00000.dat b05.dat b10.dat hard.dat
```

The **3.8×** additionally needs the no-boundary floor, which is a fourth run
with the x faces made periodic so there is no boundary to get wrong:

```sh
mkdir per && (cd per && $EXE ../nscbc-covo.inp \
    geometry.is_periodic="1 1" pelec.lo_bc="Interior Interior" \
    pelec.hi_bc="Interior Interior" pelec.bc_nscbc=0)
```

The vortex never leaves in that run, so the floor is not the residual itself but
the difference between the final field and the initial field shifted by `u t`:

```python
import numpy as np, sys
sys.path.insert(0, "../../Verification/NSCBCFields")
from metrics import load
x, y, a0, _  = load("per0.dat")       # fielddump of per/plt00000
x, y, a1, t1 = load("per.dat")        # fielddump of per/plt03000
u  = 0.2 * 34719.02257                # M * c, printed by amrex_probinit
dx = x[1] - x[0]; L = x[-1] - x[0] + dx
n  = ((u * t1) % L) / dx; i0 = int(n); f = n - i0
ref = (1 - f) * np.roll(a0, i0, axis=1) + f * np.roll(a0, i0 + 1, axis=1)
print(np.sqrt(((a1 - ref) ** 2).mean()) / 1.89811e4)   # -> 0.0047
```

Then `3.8 = 0.0176 / 0.0047`. Note the shift interpolation contributes to that
floor, so the true floor is lower and every ratio quoted here is conservative.
