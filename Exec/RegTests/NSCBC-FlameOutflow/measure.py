#!/usr/bin/env python3
"""Outflow error for NSCBC-FlameOutflow, against a causally shielded reference.

    ./measure.py --ref runs/ref runs/sigma1 runs/sigma16 ...

Why a reference run is needed at all, and why it has to be built this way, is in
README.md; the two-line version is that this configuration has a real mass
imbalance whose pressure ramp dwarfs any difference between boundary settings,
so nothing self-referencing measures the boundary, and a longer-domain reference
is only a valid reference if it differs in exactly one thing.

Requirements on the runs, all of which the reference argument silently assumes:

  * every run, reference included, uses prob.nscbc_inflow=0, so the inlet is a
    hard Dirichlet and the length-dependent inflow rate K = relax_u c / L_ref
    is not part of what is being compared;
  * the reference has the same dx and the same inlet, and its outflow far
    enough downstream that nothing from it reaches the test domain within the
    comparison time (2.4 cm of burnt gas buys 2.8e-5 s);
  * the comparison time is inside that window.

Given all that, the reference restricted to the test domain is the exact
solution there, and the difference -- INCLUDING its mean level, which must not
be subtracted -- is the test run's outflow error.

Needs `fielddump` from Verification/NSCBCFields on PATH or via --fielddump.
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np

PAMB = 1013250.0
VARS = ("pressure", "Temp", "x_velocity", "y_velocity")


def dump(exe, plt, var, scratch):
    out = os.path.join(scratch, "f.dat")
    subprocess.run([exe, plt, var, out], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        fh.readline()
        nx, ny, xlo, ylo, dx, dy, t = fh.readline().split()
        nx, ny = int(nx), int(ny)
        a = np.loadtxt(fh)
    return (a.reshape(ny, nx),
            dict(nx=nx, ny=ny, xlo=float(xlo), ylo=float(ylo), dx=float(dx),
                 dy=float(dy), t=float(t)))


def plt_at(exe, d, target, scratch):
    """The plotfile in d closest in time to target."""
    best, bt, bd = None, None, np.inf
    for p in sorted(os.listdir(d)):
        f = os.path.join(d, p)
        if not (p.startswith("plt") and os.path.isdir(f)):
            continue
        _, m = dump(exe, f, "Temp", scratch)
        if abs(m["t"] - target) < bd:
            best, bt, bd = f, m["t"], abs(m["t"] - target)
    return best, bt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help="test run directories")
    ap.add_argument("--ref", required=True, help="shielded reference directory")
    ap.add_argument("--time", type=float, default=2.4e-5,
                    help="comparison time [s]; must be inside the causal "
                         "window of the reference")
    ap.add_argument("--fielddump", default="fielddump",
                    help="path to Verification/NSCBCFields/fielddump")
    ap.add_argument("--layer", type=int, default=16,
                    help="near-boundary width in cells for the local norms")
    a = ap.parse_args()

    with tempfile.TemporaryDirectory() as scratch:
        rp, tr = plt_at(a.fielddump, a.ref, a.time, scratch)
        if rp is None:
            sys.exit(f"no plotfiles in {a.ref}")
        print(f"reference {rp}   t = {tr:.5g} s")
        R = {v: dump(a.fielddump, rp, v, scratch)[0] for v in VARS}

        hdr = (f"{'case':<18} {'<dp>':>9} {'L2(dp)':>9} {'L2(dp) bl':>10} "
               f"{'max|dp|':>9} {'L2(dT) bl':>10} {'L2(du) bl':>10}")
        print()
        print(hdr)
        print("-" * len(hdr))
        prof = {}
        for d in a.runs:
            fp, tt = plt_at(a.fielddump, d, a.time, scratch)
            if fp is None:
                continue
            if abs(tt - tr) > 0.02 * max(tr, 1e-300):
                print(f"{os.path.basename(d):<18} SKIPPED: nearest plotfile is "
                      f"at t = {tt:.4g}, reference is at {tr:.4g}")
                continue
            Q = {}
            for v in VARS:
                Q[v], m = dump(a.fielddump, fp, v, scratch)
            nx = m["nx"]
            if R["pressure"].shape[1] < nx:
                sys.exit("reference domain is smaller than the test domain")
            bl = np.arange(nx) >= nx - a.layer
            dp = Q["pressure"] - R["pressure"][:, :nx]
            dT = Q["Temp"] - R["Temp"][:, :nx]
            du = Q["x_velocity"] - R["x_velocity"][:, :nx]
            prof[os.path.basename(d)] = np.sqrt((dp ** 2).mean(axis=0))
            print(f"{os.path.basename(d):<18} {dp.mean():9.1f} "
                  f"{np.sqrt((dp**2).mean()):9.1f} "
                  f"{np.sqrt((dp[:, bl]**2).mean()):10.1f} "
                  f"{np.abs(dp).max():9.1f} "
                  f"{np.sqrt((dT[:, bl]**2).mean()):10.3f} "
                  f"{np.sqrt((du[:, bl]**2).mean()):10.3f}")

        if prof:
            cells = [0, 1, 2, 4, 8, 16, 32, 64]
            print("\nrms_y(dp) vs cells in from the outflow:")
            print(f"{'':<18}" + "".join(f"{c:>9}" for c in cells))
            for n, pr in prof.items():
                print(f"{n:<18}" + "".join(
                    f"{pr[len(pr) - 1 - c]:9.1f}" for c in cells
                    if c < len(pr)))
        print(f"\ndyn/cm^2; p_amb = {PAMB:.4g}, so 1000 is 1e-3 of ambient. "
              f"The mean is part of the error and is not subtracted.")


if __name__ == "__main__":
    main()
