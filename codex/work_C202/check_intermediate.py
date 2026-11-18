"""Sections 5-6 intermediate quantities on actual Bethe roots, with my own restricted continuum:
Q_L, L sup|Cum R|, L|nu - m(Q_L)|, L sup|sigma_L - r_Q| on [-Q,Q], compared with s, s/2, D."""
import numpy as np, json
from scipy.optimize import brentq
from ref148 import *
rows = []
for D_ in [1.01, 1.2, 2.0, float(np.cosh(1.36)), 5.0, 20.0]:
    lam = float(np.arccosh(D_)); s = 1.5 if lam < 1.36 else 0.5; g = 13.5 * s
    Dc = g * 4 / (pi * np.sinh(2 * lam)) + s / (2 * pi * np.tanh(lam))
    for L in [8, 9, 32, 33, 128, 129]:
        for n in sorted({2, L // 4, L // 2 - 2, L // 2 - 1, L // 2}):
            if n < 2: continue
            E, x = bethe_energy(L, n, lam); nu = n / L
            f = lambda y: z_L(np.array([y]), x, L, lam)[0] - nu / 2
            if f(pi / 2) <= 1e-14:
                Q = pi / 2
            else:
                # first crossing beyond x_n: scan then refine
                ys = np.linspace(x[-1], pi / 2, 4001); fs = z_L(ys, x, L, lam) - nu / 2
                i = np.argmax(fs >= 0); Q = brentq(f, ys[i - 1], ys[i], xtol=1e-15)
            N = int(min(4000, 40 * Q / lam + 200)); t, w = gl(N); y = Q * t; wy = Q * w
            r = np.linalg.solve(np.eye(N) + Kk(y[:, None] - y[None, :], lam) * wy[None, :], a1(y, lam))
            m = float(wy @ r)
            sg = sigma_L(y, x, L, lam)
            # Cum R on a fine grid: (# roots <= x)/L - (z_L(x) + nu/2)
            xs = np.linspace(-Q, Q, 20001)
            cum = np.searchsorted(x, xs, side='right') / L - (z_L(xs, x, L, lam) + nu / 2)
            rows.append(dict(Delta=D_, L=L, n=n, Q=Q, x_n=float(x[-1]), L_cumR=float(L * np.abs(cum).max()), s=s,
                             L_mismatch=float(L * abs(nu - m)), L_dens=float(L * np.abs(sg - r).max()), D=float(Dc),
                             min_sigma_on_in_Q=float(sg.min())))
bad = [r for r in rows if r['L_cumR'] > r['s'] + 1e-9 or r['L_mismatch'] > r['s'] / 2 or r['L_dens'] > r['D']]
summ = dict(cases=len(rows), violations=len(bad), max_L_cumR=max(r['L_cumR'] for r in rows), max_L_mismatch=max(r['L_mismatch'] for r in rows),
            max_dens_over_D=max(r['L_dens'] / r['D'] for r in rows), min_sigma_on_minusQ_Q=min(r['min_sigma_on_in_Q'] for r in rows))
json.dump(dict(summary=summ, rows=rows), open('check_intermediate.json', 'w'), indent=1)
print(summ)
