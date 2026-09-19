#!/usr/bin/env python3
"""Lesson-9 metrics for NSCBC-Chamber: the overpressure story in one table.

    ./chamber_metrics.py <rundir> [<rundir> ...] [--xvent 1.2]

Per plotfile: closed-end overpressure (mean p over the first 10% of the
CHAMBER length, minus p_amb), domain-mean overpressure, burnt volume fraction
of the chamber (progress variable from H2 depletion), and the volumetric flux
through the vent plane, per unit depth.  The T7 quantity of interest is the
DIFFERENCE of these traces between the vent variant and the plenum variant:
peak-aligned overpressure, V_comb' - V_vent' (zero crossing = the peak), and
the post-peak ringing.  Needs fielddump (Verification/NSCBCFields) on PATH or
via FIELDDUMP.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

FIELDDUMP = os.environ.get("FIELDDUMP", "fielddump")
PAMB = 1013250.0


def dump(plt, var, scratch):
    out = os.path.join(scratch, "f.dat")
    subprocess.run([FIELDDUMP, plt, var, out], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        fh.readline()
        nx, ny, xlo, ylo, dx, dy, t = fh.readline().split()
        nx, ny = int(nx), int(ny)
        a = np.loadtxt(fh).reshape(ny, nx)
    return a, float(dx), float(dy), float(t), float(xlo)


def main():
    argv = sys.argv[1:]
    xvent = 1.2
    if "--xvent" in argv:
        i = argv.index("--xvent")
        xvent = float(argv[i + 1])
        del argv[i:i + 2]
    args = [a for a in argv if not a.startswith("--")]
    with tempfile.TemporaryDirectory() as scratch:
        for d in args:
            plts = sorted(p for p in os.listdir(d)
                          if p.startswith("plt")
                          and os.path.isdir(os.path.join(d, p)))
            print(f"\n=== {d}  (vent plane at x = {xvent})")
            print(f"{'t [s]':>10} {'dp_closed':>10} {'dp_mean':>9} "
                  f"{'burnt frac':>10} {'Vdot_vent [cm^2/s]':>18}")
            for p in plts:
                f = os.path.join(d, p)
                pr, dx, dy, t, xlo = dump(f, "pressure", scratch)
                rho, *_ = dump(f, "density", scratch)
                rH2, *_ = dump(f, "rho_H2", scratch)
                u, *_ = dump(f, "x_velocity", scratch)
                nx_ch = int(round((xvent - xlo) / dx))   # chamber columns
                iv = min(nx_ch, pr.shape[1]) - 1         # vent-plane column
                y_h2 = rH2 / rho
                y_fresh = y_h2[:, :nx_ch].max()
                c = 1.0 - y_h2[:, :nx_ch] / max(y_fresh, 1e-30)
                print(f"{t:10.4e} {pr[:, :nx_ch//10 + 1].mean() - PAMB:10.1f} "
                      f"{pr[:, :nx_ch].mean() - PAMB:9.1f} "
                      f"{c.mean():10.4f} "
                      f"{u[:, iv].sum() * dy:18.3f}")


if __name__ == "__main__":
    main()
