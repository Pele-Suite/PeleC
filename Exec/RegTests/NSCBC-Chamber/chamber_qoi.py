#!/usr/bin/env python3
"""Lesson-9 T7 quantities of interest: the two-run comparison table.

    ./chamber_qoi.py vent plenum [--sigma16 vent-sigma16] [--xvent 1.2]

Builds on chamber_metrics.py's per-plotfile trace (and shares its fielddump
protocol) but computes the DIFFERENCE quantities the README's Status table
needs:

  * peak-aligned overpressure difference (peaks aligned, not clocks --
    ignition transients differ between the two variants),
  * the volume balance V'_comb - V'_vent, whose zero crossing must mark the
    pressure peak (V'_comb = (sigma_exp - 1) dV_burnt/dt, with sigma_exp
    measured from the run's own burnt/fresh densities, not assumed),
  * the post-peak ring-down: dominant frequency and exponential decay rate
    of the closed-end trace, where the boundary's reflection coefficient is
    directly visible.

Traces are cached to <rundir>-traces.csv beside the run directory (NOT inside
it -- a live run may still be writing there); delete the cache to force
re-extraction.  Only directories matching plt<digits> exactly are read, so
plt*.old.* collision debris is ignored.  Needs numpy and fielddump
(Verification/NSCBCFields) on PATH or via FIELDDUMP.
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

FIELDDUMP = os.environ.get("FIELDDUMP", "fielddump")
PAMB = 1013250.0
GAMMA_EFF = 1.4          # only used for the compressibility cross-check line


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


def extract(rundir, xvent, scratch, xlo_ch=None, ymask=None):
    """One row per plotfile: t, dp_closed, dp_mean, burnt_frac, Vdot_vent,
    sigma_exp.  Cached beside (not inside) the run directory.

    For the EB box variants pass --xlo (the chamber's back wall) and
    --ymask ylo,yhi (the chamber's interior rows): the chamber region is
    then the interior cavity rather than "everything left of the vent",
    which in the box would sweep in plenum rows and zeroed EB wall cells.
    Cells with rho below 1e-8 (zeroed body state) are excluded everywhere.
    """
    cache = rundir.rstrip("/") + "-traces.csv"
    plts = sorted(p for p in os.listdir(rundir)
                  if re.fullmatch(r"plt\d+", p)
                  and os.path.isdir(os.path.join(rundir, p)))
    if os.path.exists(cache):
        rows = np.loadtxt(cache, delimiter=",", ndmin=2)
        if rows.shape[0] == len(plts):
            return rows
        print(f"# cache {cache} has {rows.shape[0]} rows, run has "
              f"{len(plts)} plotfiles -- re-extracting")
    rows = []
    for n, p in enumerate(plts):
        f = os.path.join(rundir, p)
        pr, dx, dy, t, xlo = dump(f, "pressure", scratch)
        rho, *_ = dump(f, "density", scratch)
        rH2, *_ = dump(f, "rho_H2", scratch)
        u, *_ = dump(f, "x_velocity", scratch)
        i0 = 0 if xlo_ch is None else int(round((xlo_ch - xlo) / dx))
        nx_ch = int(round((xvent - xlo) / dx))       # vent-plane column
        iv = min(nx_ch, pr.shape[1]) - 1
        ny = pr.shape[0]
        j0, j1 = 0, ny
        if ymask is not None:
            j0 = int(round(ymask[0] / dy))
            j1 = int(round(ymask[1] / dy))
        live = rho > 1.0e-8                          # excludes zeroed EB body
        ch = np.zeros_like(live)
        ch[j0:j1, i0:nx_ch] = True
        ch &= live
        with np.errstate(invalid="ignore", divide="ignore"):
            y_h2 = np.where(live, rH2 / np.maximum(rho, 1e-30), 0.0)
        y_fresh = y_h2[ch].max()
        c = np.where(ch, 1.0 - y_h2 / max(y_fresh, 1e-30), 0.0)
        burnt, fresh = ch & (c > 0.9), ch & (c < 0.1)
        sigma_exp = (rho[fresh].mean() / rho[burnt].mean()
                     if burnt.any() and fresh.any() else np.nan)
        closed = np.zeros_like(live)
        closed[j0:j1, i0:i0 + (nx_ch - i0) // 10 + 1] = True
        closed &= live
        vent = live[:, iv] if ymask is None else (live[:, iv] &
                                                  (np.arange(ny) >= j0) &
                                                  (np.arange(ny) < j1))
        rows.append([t,
                     pr[closed].mean() - PAMB,
                     pr[ch].mean() - PAMB,
                     c[ch].mean() if ch.any() else np.nan,
                     u[vent, iv].sum() * dy,
                     sigma_exp])
        if n % 40 == 0:
            print(f"#   {rundir}: {n}/{len(plts)}", file=sys.stderr)
    rows = np.array(rows)
    np.savetxt(cache, rows, delimiter=",", header=(
        "t,dp_closed,dp_mean,burnt_frac,Vdot_vent,sigma_exp"))
    return rows


def peak(tr):
    """Time and value of the overpressure peak (parabolic through the
    plotfile samples)."""
    t, dp = tr[:, 0], tr[:, 1]
    i = int(np.argmax(dp))
    if 0 < i < len(t) - 1:
        a, b, c = dp[i - 1], dp[i], dp[i + 1]
        den = a - 2 * b + c
        if den != 0.0:
            s = 0.5 * (a - c) / den
            return t[i] + s * (t[1] - t[0]), b - 0.25 * (a - c) * s
    return t[i], dp[i]


def vdot_balance(tr, area):
    """V'_comb - V'_vent on the trace grid, and its zero crossing nearest
    the pressure peak.  V'_comb = (sigma_exp-1) * area * d(burnt_frac)/dt."""
    t, bf, vv, se = tr[:, 0], tr[:, 3], tr[:, 4], tr[:, 5]
    se = np.where(np.isfinite(se), se, np.nanmedian(se))
    vcomb = (se - 1.0) * area * np.gradient(bf, t)
    bal = vcomb - vv
    tp, _ = peak(tr)
    sign = np.sign(bal)
    zc = np.nonzero(np.diff(sign) != 0)[0]
    if len(zc) == 0:
        return bal, np.nan, tp
    i = zc[np.argmin(np.abs(t[zc] - tp))]
    tz = t[i] - bal[i] * (t[i + 1] - t[i]) / (bal[i + 1] - bal[i])
    return bal, tz, tp


def ringdown(tr, t_from):
    """Dominant frequency (FFT, parabolic peak interpolation) and decay rate
    (log-linear fit through the rectified extrema) of dp_closed after t_from.
    Returns f [Hz], lambda [1/s], n_extrema."""
    m = tr[:, 0] >= t_from
    t, dp = tr[m, 0], tr[m, 1]
    if len(t) < 8:
        return np.nan, np.nan, 0
    # detrend with a one-sided moving mean wider than any acoustic period
    w = max(5, len(t) // 8) | 1
    base = np.convolve(dp, np.ones(w) / w, mode="same")
    s = dp - base
    dt = np.median(np.diff(t))
    sp = np.abs(np.fft.rfft(s * np.hanning(len(s))))
    fr = np.fft.rfftfreq(len(s), dt)
    i = 1 + int(np.argmax(sp[1:]))
    if 0 < i < len(sp) - 1:
        a, b, c = sp[i - 1], sp[i], sp[i + 1]
        den = a - 2 * b + c
        i = i + (0.5 * (a - c) / den if den != 0.0 else 0.0)
    f = float(i * fr[1])
    # extrema of |s|: local maxima of the rectified signal
    ex = [(t[j], abs(s[j])) for j in range(1, len(s) - 1)
          if abs(s[j]) >= abs(s[j - 1]) and abs(s[j]) >= abs(s[j + 1])
          and abs(s[j]) > 0]
    if len(ex) < 3:
        return f, np.nan, len(ex)
    te, ae = np.array(ex).T
    lam = -np.polyfit(te, np.log(ae), 1)[0]
    return f, float(lam), len(ex)


def main():
    argv = sys.argv[1:]
    xvent, s16dir, xlo_ch, ymask = 1.2, None, None, None
    if "--xvent" in argv:
        i = argv.index("--xvent")
        xvent = float(argv[i + 1])
        del argv[i:i + 2]
    if "--xlo" in argv:
        i = argv.index("--xlo")
        xlo_ch = float(argv[i + 1])
        del argv[i:i + 2]
    if "--ymask" in argv:
        i = argv.index("--ymask")
        ymask = tuple(float(v) for v in argv[i + 1].split(","))
        del argv[i:i + 2]
    if "--sigma16" in argv:
        i = argv.index("--sigma16")
        s16dir = argv[i + 1]
        del argv[i:i + 2]
    vent_dir, plenum_dir = argv[0], argv[1]
    x0 = 0.0 if xlo_ch is None else xlo_ch
    area = (xvent - x0) * 0.3  # chamber cross-section per unit depth [cm^2]

    with tempfile.TemporaryDirectory() as scratch:
        tv = extract(vent_dir, xvent, scratch, xlo_ch, ymask)
        tp_ = extract(plenum_dir, xvent, scratch, xlo_ch, ymask)
        ts = extract(s16dir, xvent, scratch, xlo_ch, ymask) if s16dir else None

    rows = {}
    for name, tr in [("vent", tv), ("plenum", tp_)] + (
            [("vent-sigma16", ts)] if ts is not None else []):
        t_pk, dp_pk = peak(tr)
        _, t_zc, _ = vdot_balance(tr, area)
        f, lam, nex = ringdown(tr, t_pk + 5e-4)
        rows[name] = (t_pk, dp_pk, t_zc, f, lam, nex)

    print("\nPer-variant:")
    print(f"{'variant':>14} {'t_peak [ms]':>12} {'dp_peak':>9} "
          f"{'t(Vdot=0) [ms]':>15} {'f_ring [kHz]':>13} "
          f"{'lambda [1/s]':>13} {'n_ext':>6}")
    for k, (t_pk, dp_pk, t_zc, f, lam, nex) in rows.items():
        print(f"{k:>14} {1e3*t_pk:12.4f} {dp_pk:9.1f} {1e3*t_zc:15.4f} "
              f"{1e-3*f:13.3f} {lam:13.1f} {nex:6d}")

    # peak-aligned difference, vent - plenum, on the vent trace's grid
    t_pv, _ = peak(tv)
    t_pp, _ = peak(tp_)
    shift = t_pv - t_pp
    tt = tv[:, 0]
    dpl = np.interp(tt - shift, tp_[:, 0], tp_[:, 1],
                    left=np.nan, right=np.nan)
    d = tv[:, 1] - dpl
    ok = np.isfinite(d)
    pre = ok & (tt < t_pv)
    post = ok & (tt >= t_pv)
    print("\nPeak-aligned difference (vent - plenum), closed-end trace:")
    print(f"  peak shift        : {1e3*shift:+.4f} ms "
          f"(vent {1e3*t_pv:.4f}, plenum {1e3*t_pp:.4f})")
    if pre.any():
        print(f"  buildup  max|d|   : {np.abs(d[pre]).max():9.1f} dyn/cm^2")
    if post.any():
        print(f"  ringdown max|d|   : {np.abs(d[post]).max():9.1f} dyn/cm^2")
        print(f"  ringdown rms d    : {np.sqrt(np.mean(d[post]**2)):9.1f}"
              " dyn/cm^2")

    print("\nREADME table (markdown):")
    print("| variant | t_peak [ms] | dp_peak [dyn/cm^2] | t(Vdot=0) [ms] |"
          " f_ring [kHz] | decay [1/s] |")
    print("|---|---|---|---|---|---|")
    for k, (t_pk, dp_pk, t_zc, f, lam, _) in rows.items():
        print(f"| {k} | {1e3*t_pk:.3f} | {dp_pk:.0f} | {1e3*t_zc:.3f} |"
              f" {1e-3*f:.2f} | {lam:.0f} |")


if __name__ == "__main__":
    main()
