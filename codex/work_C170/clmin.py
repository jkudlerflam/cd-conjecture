#!/usr/bin/env python3
"""C170 (independent code): classical ring minimisation for the C165 family, to test the referee's
claim that charge-1 DM harmonics give a concave N^4 term and shrink or remove the coexistence window.
usage: clmin.py g Dx G  (a = 2 pi 13/100, L = 100, B = 0)"""
import sys, numpy as np
from scipy.optimize import minimize
from lr2 import parent, harm

def run(J, L, rng, a):
    j = np.arange(L)
    def f(x):
        v = x.reshape(L, 3); r = np.linalg.norm(v, axis=1); n = v / r[:, None]
        nr = np.roll(n, -1, 0); nl = np.roll(n, 1, 0)
        E = np.einsum('ia,ab,ib->', n, J, nr)
        g = nr @ J.T + nl @ J
        gv = (g - n * np.sum(g * n, 1)[:, None]) / r[:, None]
        return E, gv.ravel()
    best = None
    starts = []
    for p in (11, 12, 13, 14):
        q = -2 * np.pi * p / L
        for th in (np.pi / 2, 1.3, 1.0, 0.7, 0.4, 0.15, 0.0):
            starts.append(np.stack([np.sin(th) * np.cos(q * j), np.sin(th) * np.sin(q * j), (-1) ** j * np.cos(th)], 1).ravel())
    for _ in range(4):
        starts.append(rng.standard_normal(3 * L))
    for x0 in starts:
        x0 = x0 + 1e-3 * rng.standard_normal(3 * L)
        r = minimize(f, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-11, ftol=1e-15))
        if best is None or r.fun < best.fun - 1e-12:
            best = r
    v = best.x.reshape(L, 3); n = v / np.linalg.norm(v, axis=1)[:, None]
    N = abs(np.mean((-1) ** j * n[:, 2]))
    Fp = np.abs(np.fft.fft(n[:, 0] + 1j * n[:, 1])) / L
    k = np.argmax(Fp)
    return best.fun / L, N, k, Fp[k]

if __name__ == '__main__':
    g, Dx, G = map(float, sys.argv[1:4])
    L = 100; a = 2 * np.pi * 13 / L
    rng = np.random.default_rng(170)
    dJ, _ = harm(g=g, Dx=Dx, G=G)
    for Dz in np.round(np.arange(0.995, 1.0451, 0.0025), 4):
        J = parent(a, Dz) + dJ
        e, N, k, w = run(J, L, rng, a)
        print(f'g={g} Dx={Dx} G={G} Dz={Dz:.4f} e={e:.8f} N^2={N*N:.4f} spiral k={k} weight={w:.3f}', flush=True)
