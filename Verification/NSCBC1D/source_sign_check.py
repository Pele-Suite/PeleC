"""
Convention-free check of the sign of the reaction source term in the modelled
incoming wave of a ghost-cell NSCBC.

1-D linear acoustics with a volumetric pressure source S (= dp/dt|_{rho,e} of
heat release) confined to the cells next to a subsonic outflow:

    dp/dt + rho c^2 du/dx = S(x)
    du/dt + (1/rho) dp/dx = 0

Left end: rigid wall (u = 0).  Right end: the branch's ghost-cell fill, written
exactly as Source/NSCBC.H does it in its own variables --

    R+ = u + p/(rho c)   extrapolated with its minmod slope (order 2)
    R- = u - p/(rho c)   R-_g = R-_N + dx * L_in / (c rho c)
    L_in = K (p_N - p_t) + s * S_p(N)        s = +1, 0, -1

Godunov (exact linear Riemann) fluxes everywhere, including the boundary face.
The exact steady state of the PDE with a far-field pressure p_t has p == p_t
and du/dx = S/(rho c^2).  The equilibrium offset (p_N - p_t) K / S_p tells
which sign is right: 0 means the boundary reproduces the exact solution.
"""
import numpy as np

rho, c = 1.17e-3, 3.48e4        # cgs air
L, N = 40.0, 400
dx = L / N
x = (np.arange(N) + 0.5) * dx
p_t = 1.0e6
sigma = 1.0
K = sigma * c / L

S = np.zeros(N)
S[-1] = 2.0e7                   # dp/dt|chem in the boundary cell, dyn/cm2/s
Sp_N = S[-1]

def minmod(a, b):
    return np.where(a * b > 0, np.sign(a) * np.minimum(abs(a), abs(b)), 0.0)

def run(s, nsteps=60000, cfl=0.5):
    p = np.full(N, p_t)
    u = np.zeros(N)
    dt = cfl * dx / c
    rc = rho * c
    for _ in range(nsteps):
        # ghost at right (the kernel's fill, one layer)
        Rp = u + p / rc
        Rm = u - p / rc
        dRp = minmod(Rp[-1] - Rp[-2], Rp[-2] - Rp[-3])
        Rp_g = Rp[-1] + dRp
        L_in = K * (p[-1] - p_t) + s * Sp_N
        Rm_g = Rm[-1] + dx * L_in / (c * rc)
        u_g = 0.5 * (Rp_g + Rm_g)
        p_g = 0.5 * rc * (Rp_g - Rm_g)
        # left ghost: rigid wall (reflect u)
        pe = np.concatenate(([p[0]], p, [p_g]))
        ue = np.concatenate(([-u[0]], u, [u_g]))
        pL, pR = pe[:-1], pe[1:]
        uL, uR = ue[:-1], ue[1:]
        pf = 0.5 * (pL + pR) + 0.5 * rc * (uL - uR)
        uf = 0.5 * (uL + uR) + 0.5 * (pL - pR) / rc
        p = p - dt * rho * c * c * (uf[1:] - uf[:-1]) / dx + dt * S
        u = u - dt * (pf[1:] - pf[:-1]) / (rho * dx)
    return p, u

print(f"K = {K:.3g} 1/s, S_p = {Sp_N:.3g}, predicted offsets S_p/K = {Sp_N/K:.4g}")
for s, name in [(+1, "+S_p (NSCBC.H:1029, Phase 0 onward)"), (0, "term off (beta_s=1)"),
                (-1, "-S_p (as written before Phase 0)")]:
    p, u = run(s)
    off = p[-1] - p_t
    print(f"s={s:+d}  {name:32s}  p_N - p_t = {off:10.2f}   (p_N-p_t)K/S_p = {off*K/Sp_N:6.3f}"
          f"   du/dx at face = {(u[-1]-u[-2])/dx:.4g}  exact {Sp_N/(rho*c*c):.4g}")
