#!/usr/bin/env python3
"""Metrics for the flame-exit test (NSCBC-FlameOutflow, fast-stream variant).

The exit test has no shielded reference -- at u_ratio = 40 the exit takes
~5e-4 s and nothing stays causally shielded that long.  It does not need one:

  * DURING transit, the physics fixes what the front must do: advect at
    U - S_L, keep its wrinkle amplitude (decay times are 100x the window),
    keep its radical peak.  Front kinematics and structure are tracked by the
    Y_H peak per transverse row, parabola-refined, exactly as radicals.py
    does for the parent case.
  * AFTER the exit the exact solution is known outright: uniform fresh stream
    at (U, p_amb, T_in).  Any residual field is boundary-made.  This is the
    strongest metric in the file and needs no reference at all.

Usage: exit_metrics.py rundir [rundir ...]
Needs fielddump (Verification/NSCBCFields) on PATH or in FIELDDUMP.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

FIELDDUMP = os.environ.get("FIELDDUMP", "fielddump")
PAMB = 1013250.0
UIN = 912.1643066  # 40 S_L, printed by amrex_probinit
TIN = 298.0


def dump(plt, var, scratch):
    out = os.path.join(scratch, "f.dat")
    subprocess.run([FIELDDUMP, plt, var, out], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        fh.readline()
        nx, ny, xlo, ylo, dx, dy, t = fh.readline().split()
        nx, ny = int(nx), int(ny)
        a = np.loadtxt(fh)
    return a.reshape(ny, nx), (nx, ny, float(xlo), float(ylo), float(dx),
                               float(dy), float(t))


def front(plt, scratch):
    """(mean x_f, wrinkle half-amp, resolved-row count, max Y_H)."""
    rho, m = dump(plt, "density", scratch)
    rhoY, _ = dump(plt, "rho_H", scratch)
    Y = rhoY / rho
    nx, ny, xlo, ylo, dx, dy, t = m
    x = xlo + (np.arange(nx) + 0.5) * dx
    xf = np.full(ny, np.nan)
    for j in range(ny):
        r = Y[j]
        i = int(np.argmax(r))
        if i in (0, nx - 1) or r[i] < 1e-8:  # peak gone or on the edge
            continue
        a, b, c = r[i - 1], r[i], r[i + 1]
        den = a - 2 * b + c
        xf[j] = x[i] + (0.5 * (a - c) / den * dx if den != 0 else 0.0)
    good = np.isfinite(xf)
    amp = 0.5 * (np.nanmax(xf) - np.nanmin(xf)) if good.sum() > 3 else np.nan
    mean = np.nanmean(xf) if good.any() else np.nan
    return mean, amp, int(good.sum()), float(Y.max()), m


def fields(plt, scratch):
    p, m = dump(plt, "pressure", scratch)
    T, _ = dump(plt, "Temp", scratch)
    u, _ = dump(plt, "x_velocity", scratch)
    return p, T, u, m


def main():
    for d in sys.argv[1:]:
        plts = sorted(p for p in os.listdir(d)
                      if p.startswith("plt") and os.path.isdir(f"{d}/{p}"))
        name = os.path.basename(os.path.normpath(d))
        print(f"\n=== {name}")
        print(f"{'t [s]':>10} {'<p>-pamb':>9} {'max|p-<p>|':>10} "
              f"{'x_f':>8} {'w-amp':>8} {'rows':>5} {'maxY_H':>9} "
              f"{'rms(u-U)':>9} {'rms(T-Tin)':>10} {'maxT':>7}")
        with tempfile.TemporaryDirectory() as sc:
            for pl in plts:
                plt = f"{d}/{pl}"
                xf, amp, rows, ymax, _ = front(plt, sc)
                p, T, u, m = fields(plt, sc)
                pm = p.mean()
                pac = np.abs(p - pm).max()
                # Post-exit cleanliness: once no row resolves a front, the
                # exact solution is the uniform fresh stream.
                rms_u = np.sqrt(((u - UIN) ** 2).mean())
                rms_T = np.sqrt(((T - TIN) ** 2).mean())
                print(f"{m[6]:10.3e} {pm - PAMB:9.1f} {pac:10.1f} "
                      f"{xf:8.4f} {amp:8.5f} {rows:5d} {ymax:9.2e} "
                      f"{rms_u:9.2f} {rms_T:10.2f} {T.max():7.0f}")


if __name__ == "__main__":
    main()
