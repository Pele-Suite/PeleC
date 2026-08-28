# NSCBC-PMF — a lean H₂/air flame against a characteristic outflow

A 1-D premixed flame (LiDryer, φ = 0.4 H₂/air) initialised from the shipped PMF
profile, with a subsonic characteristic **inflow** of fresh gas at one end and a
characteristic **outflow** in burnt gas at the other. `prob.standoff` moves the
flame relative to the outflow; the domain is fixed so that `L_ref`, and hence the
relaxation rate `K = σ(1−M²)c/L`, is identical in every run.

```sh
./PeleC2d.llvm.ex nscbc-pmf.inp prob.standoff=-2.0    # flame 1.0 cm from the outflow
./PeleC2d.llvm.ex nscbc-pmf.inp prob.standoff=-1.2    # 0.2 cm
```

The flame is effectively stationary over the run: at S_L ≈ 23 cm/s it travels
2×10⁻³ cm in 10⁻⁴ s, while the acoustic transit is 2.5×10⁻⁵ s. It is a
quasi-frozen, actively-burning heat source, and what is being measured is how the
boundary copes.

## What this case is good for

It is the only case in the suite that exercises the **inflow** path end to end
against real chemistry and real transport — velocity, temperature and all nine
species imposed as incoming characteristics, with the outgoing acoustic taken
from the interior. It is a genuine reacting-flow regression for the boundary
condition as a whole.

## What it is *not* good for, and why

**It does not meaningfully test `bc_nscbc_beta_s`.** The reaction source term
exists to correct for heat release *in the boundary cell*, and in this
configuration there is none. The temperature profile at the outflow reads

```
1426  1426  1426  1426  1426  1426   K
```

— equilibrium burnt gas. Measured at `standoff = -2.0`, after 2.4 relaxation
times:

| `beta_s` | p at outflow − p∞ | domain mean − p∞ |
|---|---|---|
| 1.0 (source off) | 65.1 | 80.5 |
| 0.0 (source on) | 79.8 | 97.3 |

Those offsets are 6×10⁻⁵ of ambient and the difference between them is 15
dyn/cm² — noise. Pushing the flame closer (`standoff = -1.2`, 0.2 cm) only moves
the boundary into the slowly-recombining tail; it does not put the reaction zone
on the boundary.

A 1-D flame cannot put its reaction sheet *on* an outflow without the flame
leaving the domain. Testing the source term needs a configuration where the
reaction zone crosses the boundary — see `NSCBC-FlameOutflow`, where a wrinkled,
unanchored sheet is parked on the outflow plane so that half the boundary is
burnt and half is fresh, with the reaction zone in the boundary cells at the two
crossings and the flame normal tilted 51° there.

That case answers the question this one raised, and the answer is not the
expected one: `beta_s` is still not what matters. At a front-crossing outflow
the dominant error is a ghost-pressure bias driven by the *normal velocity
gradient*, it is present with chemistry switched off, and `sigma` is the control
that moves it.

## Two traps this case walked into

**`pelec.allow_negative_energy = 0` is harmless in an inert case and fatal
here.** With reactions on, it destroyed the state on the first step —
temperature collapsing from 298/1426 K to ~11 K across the entire domain. The
inputs file therefore sets it explicitly to 1 with a comment, because the value
was copied in from the inert NSCBC cases and cost an hour to find. Note that the
symptom looks exactly like a boundary blow-up: global NaNs a couple of steps in.
The bisection that identified it (`do_react=0` stable, `do_react=1` with
diffusion off unstable) is the reason it was not blamed on the boundary
condition.

**The `bcnormal` copied from `Exec/RegTests/PMF` evaluates the PMF profile at
the raw domain coordinate**, without subtracting `standoff`, while
`pc_initdata` uses `position − standoff`. With `standoff = 0` — the only value
the original case uses — the two agree. With any other value the boundary
imposes the wrong end of the profile; here that meant feeding 298 K fresh gas
into the burnt side, and the run NaN'd in two steps. Fixed in this copy; the
original in `Exec/RegTests/PMF` still has it, and it is a latent trap for anyone
who sets a non-zero standoff there.

`pc_initdata` here is also dimension-agnostic (the flame normal is the last
dimension), where the original hardcodes `prob_lo[2]` and is 3-D only.
