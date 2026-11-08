"""C161 (Claude Opus agent, 2026-11-08). Interaction of tilted solitons (soliton_valley.py functional) on a ring.
Ring of length L = m*T with m solitons (each winding +2 pi in phi), tilts all equal (same valley) or alternating.
U(T) = E_ring/m - E_sol, with E_sol the isolated tilted soliton energy (both relative to the uniform +x state).
"""
import sys
import numpy as np
from scipy.optimize import minimize

def ring_energy(v, dx, h, kappa, mu):
    V = v.reshape(-1, 3); nr = np.linalg.norm(V, axis=1); n = V / nr[:, None]
    n1 = np.roll(n, -1, 0); d = (n1 - n) / dx
    cz = n[:, 0] * n1[:, 1] - n[:, 1] * n1[:, 0]
    W = -h * n[:, 0] + 0.5 * kappa * n[:, 2] ** 2 + 0.5 * mu * n[:, 1] ** 2 + h
    E = np.sum(0.5 * np.sum(d ** 2, 1)) * dx - np.sum(cz) + np.sum(W) * dx
    G = -d + np.roll(d, 1, 0)
    n0 = np.roll(n, 1, 0)
    G[:, 0] += -n1[:, 1] + n0[:, 1]
    G[:, 1] += n1[:, 0] - n0[:, 0]
    G[:, 0] += -h * dx; G[:, 1] += mu * n[:, 1] * dx; G[:, 2] += kappa * n[:, 2] * dx
    G = (G - np.sum(G * n, 1)[:, None] * n) / nr[:, None]
    return E, G.ravel()

def ring_min(h, kappa, mu, T, m, tilts, dx=0.08):
    L = m * T; N = int(round(L / dx)); dx = L / N
    x = np.arange(N) * dx
    ph = np.zeros(N); nz = np.zeros(N)
    for j in range(m):
        c = (j + 0.5) * T
        ph += np.pi * (1 + np.tanh(1.5 * (x - c)))
        nz += tilts[j] * np.exp(-(x - c) ** 2 / 4)
    n0 = np.stack([np.cos(ph), np.sin(ph), nz], 1); n0 /= np.linalg.norm(n0, axis=1)[:, None]
    r = minimize(ring_energy, n0.ravel(), args=(dx, h, kappa, mu), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=60000, gtol=1e-11, ftol=1e-15))
    V = r.x.reshape(-1, 3); n = V / np.linalg.norm(V, axis=1)[:, None]
    return r.fun, n

if __name__ == '__main__':
    kappa, mu, h = 0.05, 2.3, float(sys.argv[1]) if len(sys.argv) > 1 else 0.0265
    Einf, _ = ring_min(h, kappa, mu, 60.0, 1, [0.5])
    print(f'kappa={kappa} mu={mu} h={h}  E_sol(ring T=60) = {Einf:+.6f}')
    print('T      U_same        U_alt        maxnz_same  signs_alt')
    for T in (8.0, 10.0, 12.0, 15.0, 18.0, 22.0):
        Es, ns = ring_min(h, kappa, mu, T, 2, [0.5, 0.5])
        Ea, na = ring_min(h, kappa, mu, T, 2, [0.5, -0.5])
        print(f'{T:5.1f}  {Es/2-Einf:+.6e}  {Ea/2-Einf:+.6e}  {np.abs(ns[:,2]).max():.3f}  {np.sign(na[:,2].max())*1:+.0f}/{np.sign(na[:,2].min()):+.0f}')
        sys.stdout.flush()
