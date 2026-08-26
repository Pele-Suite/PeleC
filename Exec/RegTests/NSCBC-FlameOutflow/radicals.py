"""Reaction-zone diagnostics from a radical rather than from the temperature.

The temperature midpoint sits in the preheat zone, which is wide and diffusion-
controlled; the radical peak sits in the reaction zone, which is what actually
has to cross the boundary for the beta_s term to have anything to correct.  For
the front POSITION the radical is also sharper: the T-midpoint contour is a
smooth ramp over ~10 cells, the CH3 peak is 2-3.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure import dump as _dump   # noqa: E402


FIELDDUMP = os.environ.get('FIELDDUMP', 'fielddump')


def var(pltf, name):
    import tempfile
    with tempfile.TemporaryDirectory() as sc:
        return _dump(FIELDDUMP, pltf, name, sc)


def massfrac(pltf, sp):
    rho, m = var(pltf, 'density')
    rhoY, _ = var(pltf, f'rho_{sp}')
    return rhoY / rho, m


def peak_locus(pltf, sp):
    """x of the radical peak in each transverse row, parabola-refined."""
    Y, m = massfrac(pltf, sp)
    nx, ny, xlo, ylo, dx, dy, t = m
    x = xlo + (np.arange(nx) + 0.5) * dx
    out = np.full(ny, np.nan)
    for j in range(ny):
        r = Y[j]
        i = int(np.argmax(r))
        if i in (0, nx - 1) or r[i] <= 0:
            continue
        a, b, c = r[i - 1], r[i], r[i + 1]
        den = a - 2 * b + c
        out[j] = x[i] + (0.5 * (a - c) / den * dx if den != 0 else 0.0)
    return out, m


def report(pltfs, sp, tag):
    print(f'--- {tag}: front located by the {sp} peak')
    for p in pltfs:
        xf, m = peak_locus(p, sp)
        good = np.isfinite(xf)
        if good.sum() < 4:
            print(f'   {os.path.basename(p)}: {sp} peak not resolved '
                  f'({good.sum()} rows)')
            continue
        amp = 0.5 * (np.nanmax(xf) - np.nanmin(xf))
        print(f'   t = {m[6]:9.3e}   mean x_f = {np.nanmean(xf):.5f} cm   '
              f'wrinkle half-amp = {amp:.5f} cm')


def boundary_profile(pltf, species):
    """Radicals and heat release along the last interior column."""
    out = {}
    for sp in species:
        Y, m = massfrac(pltf, sp)
        out[sp] = Y[:, -1]
    q, m = var(pltf, 'heatRelease')
    out['q'] = q[:, -1]
    T, _ = var(pltf, 'Temp')
    out['T'] = T[:, -1]
    out['y'] = m[3] + (np.arange(m[1]) + 0.5) * m[5]
    out['t'] = m[6]
    return out


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'ch4'
    if which == 'h2':
        d = '/tmp/fo/final/ref'
        ps = [f'{d}/plt00000', f'{d}/plt00400', f'{d}/plt00974']
        report(ps, 'H', 'H2/air phi=0.40, LiDryer, Le~0.3')
        b = boundary_profile(ps[-1], ['H', 'OH'])
    else:
        d = sys.argv[2] if len(sys.argv) > 2 else '/tmp/ch4/ic'
        ps = sorted(f'{d}/{p}' for p in os.listdir(d)
                    if p.startswith('plt') and os.path.isdir(f'{d}/{p}'))
        report(ps, 'CH3', 'CH4/air phi=0.75, drm19, Le=0.97')
        b = boundary_profile(ps[-1], ['CH3', 'H'])
    print(f'\n--- boundary column at t = {b["t"]:.3e} s')
    keys = [k for k in b if k not in ('y', 't', 'q', 'T')]
    hdr = f'{"y [cm]":>8} {"T [K]":>8} ' + ' '.join(f'{"Y_"+k:>10}' for k in keys)
    print(hdr + f' {"q [erg/cm3/s]":>14}')
    for j in range(0, len(b['y']), 4):
        print(f'{b["y"][j]:8.4f} {b["T"][j]:8.0f} ' +
              ' '.join(f'{b[k][j]:10.3e}' for k in keys) +
              f' {b["q"][j]:14.3e}')
