#!/usr/bin/env python3
"""Duct-mode metrics for NSCBC-Acoustic: the PeleC half of driver t3/t5.

    ./duct_metrics.py <rundir> --freq F [--t0 T] [--u0 U] [--amp A]

Fourier-projects the inlet-column velocity and pressure at the forcing
frequency over t > t0 (default 6 acoustic round trips), exactly as the
driver's duct_forced() does, and accumulates P_RMS(x):

  I_u  = achieved u' amplitude at the inlet / target amplitude
  I_in = incoming-invariant (u + p/rho c) amplitude / ideal injector's 2A

plus the P_RMS(x) profile and its correlation against the analytic
|sin(k_eff (L-x))| envelope.  Sampling comes from the plotfiles, so
amr.plot_per must resolve the forcing period (the duct inp uses period/24).
Needs fielddump (Verification/NSCBCFields) on PATH or via FIELDDUMP.
"""
import argparse
import os
import re
import subprocess
import tempfile

import numpy as np

FIELDDUMP = os.environ.get("FIELDDUMP", "fielddump")
RHO0 = 1.1339884e-3  # printed by amrex_probinit (Null mechanism, 300 K, 1 atm)
C0 = 34719.0


def dump(plt, var, scratch):
    out = os.path.join(scratch, "f.dat")
    subprocess.run([FIELDDUMP, plt, var, out], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        fh.readline()
        nx, ny, xlo, ylo, dx, dy, t = fh.readline().split()
        nx, ny = int(nx), int(ny)
        a = np.loadtxt(fh).reshape(ny, nx)
    return a, float(dx), float(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--freq", type=float, required=True)
    ap.add_argument("--t0", type=float, default=None,
                    help="start of the window (default 6 t_a)")
    ap.add_argument("--u0", type=float, default=2000.0)
    ap.add_argument("--amp", type=float, default=34.7)
    a = ap.parse_args()

    plts = sorted(p for p in os.listdir(a.rundir)
                  if re.fullmatch(r"plt\d+", p)
                  and os.path.isdir(os.path.join(a.rundir, p)))
    om = 2.0 * np.pi * a.freq
    with tempfile.TemporaryDirectory() as scratch:
        rows = []
        for p in plts:
            f = os.path.join(a.rundir, p)
            pr, dx, t = dump(f, "pressure", scratch)
            u, _, _ = dump(f, "x_velocity", scratch)
            rows.append((t, u.mean(axis=0), pr.mean(axis=0)))
    rows.sort(key=lambda r: r[0])
    t = np.array([r[0] for r in rows])
    U = np.array([r[1] for r in rows])   # (nt, nx), y-averaged
    P = np.array([r[2] for r in rows])
    nx = U.shape[1]
    L = nx * dx

    t_a = 2.0 * L / C0
    t0 = a.t0 if a.t0 is not None else 6.0 * t_a
    m = t >= t0
    if m.sum() < 8:
        raise SystemExit(f"only {m.sum()} plotfiles after t0 = {t0:g}")
    tm, Um, Pm = t[m], U[m], P[m]
    w = np.gradient(tm)  # trapezoid-ish weights on the (nearly uniform) grid
    W = w.sum()

    def proj(sig):
        # Linear detrend before projecting: with relax_u small the mean is
        # weakly anchored and its slow wander leaks into a finite-window
        # Fourier projection as a spurious amplitude at every frequency.
        fit = np.polyfit(tm, sig, 1)
        s0 = sig - np.polyval(fit, tm)
        s = (2.0 / W) * np.sum(s0 * np.sin(om * tm) * w)
        c = (2.0 / W) * np.sum(s0 * np.cos(om * tm) * w)
        return np.hypot(s, c)

    amp_u = proj(Um[:, 0])
    amp_R = proj(Um[:, 0] + Pm[:, 0] / (RHO0 * C0))
    I_u = amp_u / a.amp
    I_in = amp_R / (2.0 * a.amp)

    pmean = np.sum(Pm * w[:, None], axis=0) / W
    prms = np.sqrt(np.maximum(
        np.sum((Pm - pmean) ** 2 * w[:, None], axis=0) / W, 0.0))
    keff = om * C0 / (C0 * C0 - a.u0 * a.u0)
    x = (np.arange(nx) + 0.5) * dx
    env = np.abs(np.sin(keff * (L - x)))
    corr = float(np.dot(prms, env) /
                 max(np.linalg.norm(prms) * np.linalg.norm(env), 1e-300))

    print(f"f = {a.freq:.1f} Hz   window t in [{t0:.4g}, {tm[-1]:.4g}] s "
          f"({m.sum()} samples)")
    print(f"I_u  = {I_u:.3f}")
    print(f"I_in = {I_in:.3f}")
    print(f"P_RMS antinode = {prms.max():.1f} dyn/cm^2, "
          f"envelope correlation = {corr:.3f}")
    st = np.linspace(0, nx - 1, 17).astype(int)
    print("x/L :", " ".join(f"{x[i]/L:6.2f}" for i in st))
    print("Prms:", " ".join(f"{prms[i]:6.1f}" for i in st))


if __name__ == "__main__":
    main()
