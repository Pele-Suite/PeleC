#!/usr/bin/env python3
"""
Metrics for the multi-dimensional NSCBC tests.

Reads the regular arrays written by fielddump -- 2-D or 3-D, detected from the
header -- and reports:

  circularity  -- for the circular-pulse case, the spread in the radius of the
                  outgoing wavefront over polar angle, restricted to the rays
                  along which the front is still inside the domain.  The exact
                  solution is a circle, so any spread is boundary error.  This
                  is the diagnostic used by Motheau et al. (2017): a good
                  non-reflecting boundary lets the pressure contours stay
                  circular as the wave crosses; a poor one flattens and buckles
                  them, worst at the corners where the wave arrives obliquely.

  sphericity   -- the 3-D form, for Exec/RegTests/NSCBC-Acoustic with
                  prob.pulse_type = 1.  Same idea, but the rays are spread over
                  the sphere and the report is binned by CHI, the angle between
                  the ray and the nearest face normal: chi = 0 is a face centre,
                  chi = 45 deg an edge, chi = 54.7 deg a corner.  That single
                  coordinate covers all three, which is what makes this the test
                  of corner and edge ownership rather than of the 1-D algebra.

                  Read it as follows.  Early on every ray is still inside and
                  the spread measures nothing but discretisation.  Once the
                  front has crossed the faces only the oblique rays survive the
                  "still inside" filter -- and those are exactly the ones aimed
                  at edges and corners.  So a spread that grows after the face
                  crossing is error injected at the faces and carried into the
                  still-interior part of the front, which is the thing a corner
                  bug produces.

  residual     -- what is left in the domain after the wave has gone, relative
                  to the incident amplitude.  For the vortex case this is the
                  whole story: a vortex carries no acoustic content, so any
                  pressure signal left behind was manufactured by the boundary.

Usage:
    metrics.py circularity <p_amb> <c> <file> [<file> ...]
    metrics.py sphericity  <p_amb> <c> <file> [<file> ...]
    metrics.py residual    <p_amb> <ref_amp_file> <file> [<file> ...]
"""
import sys
import numpy as np


def load(fn):
    """2-D loader, kept exactly as it was so the 2-D callers are untouched."""
    with open(fn) as f:
        f.readline()
        nx, ny, xlo, ylo, dx, dy, t = f.readline().split()
        nx, ny = int(nx), int(ny)
        a = np.loadtxt(f)
    a = a.reshape(ny, nx)
    x = float(xlo) + (np.arange(nx) + 0.5) * float(dx)
    y = float(ylo) + (np.arange(ny) + 0.5) * float(dy)
    return x, y, a, float(t)


def load_nd(fn):
    """Dimension-agnostic loader: returns (axes, array, t) with axes a list of
    coordinate vectors, slowest axis first in the array as fielddump writes it.
    """
    with open(fn) as f:
        f.readline()
        tok = f.readline().split()
        a = np.loadtxt(f)
    if len(tok) == 7:
        nx, ny = int(tok[0]), int(tok[1])
        xlo, ylo, dx, dy, t = (float(v) for v in tok[2:])
        axes = [xlo + (np.arange(nx) + 0.5) * dx,
                ylo + (np.arange(ny) + 0.5) * dy]
        return axes, a.reshape(ny, nx), t
    if len(tok) == 10:
        nx, ny, nz = (int(v) for v in tok[:3])
        xlo, ylo, zlo, dx, dy, dz, t = (float(v) for v in tok[3:])
        axes = [xlo + (np.arange(nx) + 0.5) * dx,
                ylo + (np.arange(ny) + 0.5) * dy,
                zlo + (np.arange(nz) + 0.5) * dz]
        return axes, a.reshape(nz, ny, nx), t
    raise ValueError(f"{fn}: header has {len(tok)} fields, expected 7 or 10")


def trilinear(axes, a, q):
    """Sample a 3-D array on a regular grid at scattered points q[:, 3]."""
    x, y, z = axes
    d = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
    o = (x[0], y[0], z[0])
    n = (len(x), len(y), len(z))
    f = [(q[:, m] - o[m]) / d[m] for m in range(3)]
    i0 = [np.clip(np.floor(f[m]).astype(int), 0, n[m] - 2) for m in range(3)]
    tt = [np.clip(f[m] - i0[m], 0.0, 1.0) for m in range(3)]
    out = 0.0
    for bz in (0, 1):
        for by in (0, 1):
            for bx in (0, 1):
                wgt = ((tt[0] if bx else 1 - tt[0]) *
                       (tt[1] if by else 1 - tt[1]) *
                       (tt[2] if bz else 1 - tt[2]))
                out = out + wgt * a[i0[2] + bz, i0[1] + by, i0[0] + bx]
    return out


def fibonacci_directions(n):
    """n roughly equal-area directions on the unit sphere."""
    k = np.arange(n) + 0.5
    cz = 1.0 - 2.0 * k / n
    r = np.sqrt(np.maximum(0.0, 1.0 - cz * cz))
    phi = np.pi * (1.0 + 5.0 ** 0.5) * k
    return np.stack([r * np.cos(phi), r * np.sin(phi), cz], axis=1)


def sphericity(p_amb, c, files, ndir=4000):
    print()
    print("  Wavefront sphericity -- radius of the outgoing front vs direction")
    print("  (exact solution is a sphere; every deviation is boundary error)")
    print("  chi = angle to the nearest face normal: 0 deg face, 45 edge, "
          "54.7 corner")
    print()
    hdr = ("  %-18s %9s %7s %10s %9s %9s   %s" %
           ("file", "t [s]", "rays", "r_mean", "spread%", "amp_sprd%",
            "spread% by chi [0-20|20-40|40-55]"))
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for fn in files:
        axes, a, t = load_nd(fn)
        if a.ndim != 3:
            print(f"  {fn}: not a 3-D dump")
            continue
        x, y, z = axes
        ctr = np.array([0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1]),
                        0.5 * (z[0] + z[-1])])
        half = np.array([0.5 * (x[-1] - x[0]), 0.5 * (y[-1] - y[0]),
                         0.5 * (z[-1] - z[0])])
        dp = a - p_amb
        r_th = c * t
        nh = fibonacci_directions(ndir)
        # Distance from the centre to the boundary along each ray.
        with np.errstate(divide="ignore"):
            dbnd = np.min(half / np.maximum(np.abs(nh), 1e-300), axis=1)
        keep = r_th < 0.90 * dbnd
        if (keep.sum() < 16) or (r_th <= 0):
            print("  %-18s %9.3e %7s" % (fn.split("/")[-1], t, "front gone"))
            continue
        nk = nh[keep]
        # chi: angle to the nearest face normal, i.e. to the largest |n|
        chi = np.degrees(np.arccos(np.max(np.abs(nk), axis=1)))
        rr = np.linspace(0.55 * r_th, 1.45 * r_th, 2400)
        r_peak = np.empty(len(nk))
        amp = np.empty(len(nk))
        # Chunked so the sample array stays a sane size.
        step = max(1, 2_000_000 // len(rr))
        for s0 in range(0, len(nk), step):
            sub = nk[s0:s0 + step]
            q = (ctr[None, None, :] +
                 rr[:, None, None] * sub[None, :, :]).reshape(-1, 3)
            vals = trilinear(axes, a * 0 + dp, q).reshape(len(rr), len(sub))
            kmax = np.argmax(np.abs(vals), axis=0)
            r_peak[s0:s0 + len(sub)] = rr[kmax]
            amp[s0:s0 + len(sub)] = np.abs(vals[kmax, np.arange(len(sub))])
        spread = (r_peak.max() - r_peak.min()) / r_peak.mean()
        amp_spread = (amp.max() - amp.min()) / amp.mean()
        byc = []
        for lo_, hi_ in ((0, 20), (20, 40), (40, 55)):
            m = (chi >= lo_) & (chi < hi_)
            byc.append("%.3f" % (100 * (r_peak[m].max() - r_peak[m].min()) /
                                 r_peak[m].mean()) if m.sum() > 8 else "  --  ")
        print("  %-18s %9.3e %7d %10.4f %9.3f %9.3f   %s" %
              (fn.split("/")[-1], t, keep.sum(), r_peak.mean(), 100 * spread,
               100 * amp_spread, " | ".join(byc)))


def bilinear(x, y, a, xq, yq):
    """Sample a on a regular grid at scattered (xq, yq)."""
    dx, dy = x[1] - x[0], y[1] - y[0]
    fi = (xq - x[0]) / dx
    fj = (yq - y[0]) / dy
    i0 = np.clip(np.floor(fi).astype(int), 0, len(x) - 2)
    j0 = np.clip(np.floor(fj).astype(int), 0, len(y) - 2)
    tx = np.clip(fi - i0, 0.0, 1.0)
    ty = np.clip(fj - j0, 0.0, 1.0)
    return ((1 - tx) * (1 - ty) * a[j0, i0] + tx * (1 - ty) * a[j0, i0 + 1] +
            (1 - tx) * ty * a[j0 + 1, i0] + tx * ty * a[j0 + 1, i0 + 1])


def dist_to_boundary(theta, x, y):
    """Distance from the origin to the domain boundary along each ray."""
    xmax, ymax = x[-1], y[-1]
    xmin, ymin = x[0], y[0]
    ct, st = np.cos(theta), np.sin(theta)
    big = 1e30
    with np.errstate(divide="ignore", invalid="ignore"):
        dx1 = np.where(ct > 0, xmax / ct, np.where(ct < 0, xmin / ct, big))
        dy1 = np.where(st > 0, ymax / st, np.where(st < 0, ymin / st, big))
    return np.minimum(dx1, dy1)


def circularity(p_amb, c, files, ntheta=720):
    print()
    print("  Wavefront circularity -- radius of the outgoing front vs polar angle")
    print("  (exact solution is a circle; every deviation is boundary error)")
    print()
    print("  %-16s %8s %8s %10s %10s %10s %10s" %
          ("file", "t/tau", "rays", "r_mean", "spread%", "amp_sprd%", "peak|dp|"))
    print("  " + "-" * 82)
    for fn in files:
        x, y, a, t = load(fn)
        dp = a - p_amb
        r_th = c * t
        theta = np.linspace(0.0, 2 * np.pi, ntheta, endpoint=False)
        dbnd = dist_to_boundary(theta, x, y)
        # Only rays whose front is still comfortably inside the domain.
        keep = r_th < 0.90 * dbnd
        if keep.sum() < 8 or r_th <= 0:
            print("  %-16s %8.3f %8s" % (fn.split("/")[-1], t, "front gone"))
            continue
        th = theta[keep]
        # Scan each retained ray for the front.
        # The radial scan must be finer than the effect being measured: at 400
        # samples one increment is ~0.2 % of r_th, which quantises the radius
        # spread and makes two genuinely different boundaries report the same
        # number.  2400 puts the quantisation an order of magnitude below the
        # grid spacing.
        rr = np.linspace(0.55 * r_th, 1.45 * r_th, 2400)
        R, TH = np.meshgrid(rr, th, indexing="ij")
        vals = bilinear(x, y, dp, R * np.cos(TH), R * np.sin(TH))
        k = np.argmax(np.abs(vals), axis=0)
        r_peak = rr[k]
        amp = np.abs(vals[k, np.arange(len(th))])
        spread = (r_peak.max() - r_peak.min()) / r_peak.mean()
        amp_spread = (amp.max() - amp.min()) / amp.mean()
        print("  %-16s %8.3f %8d %10.4f %10.3f %10.3f %10.4e" %
              (fn.split("/")[-1], r_th / max(x[-1], 1e-30), keep.sum(),
               r_peak.mean(), 100 * spread, 100 * amp_spread, amp.mean()))


def residual(p_amb, ref_file, files):
    _, a0, _ = load_nd(ref_file)
    inc = np.abs(a0 - p_amb).max()
    print()
    print("  Residual pressure disturbance, relative to the incident amplitude")
    print("  incident |dp|_max = %.5e dyn/cm^2" % inc)
    print()
    print("  %-16s %10s %12s %12s %14s" %
          ("file", "t [s]", "max|dp|/inc", "L2|dp|/inc", "mean p"))
    print("  " + "-" * 68)
    for fn in files:
        _, a, t = load_nd(fn)
        dp = a - p_amb
        l2 = np.sqrt((dp ** 2).mean())
        print("  %-16s %10.4e %12.5f %12.6f %14.4f" %
              (fn.split("/")[-1], t, np.abs(dp).max() / inc, l2 / inc, a.mean()))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "circularity":
        circularity(float(sys.argv[2]), float(sys.argv[3]), sys.argv[4:])
    elif mode == "sphericity":
        sphericity(float(sys.argv[2]), float(sys.argv[3]), sys.argv[4:])
    elif mode == "residual":
        residual(float(sys.argv[2]), sys.argv[3], sys.argv[4:])
    else:
        print(__doc__)
        sys.exit(1)
