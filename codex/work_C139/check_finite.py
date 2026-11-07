"""Finite-L checks of the C139 intermediate inequalities for 0<Delta<1. J=1, M0=1/2.
Usage: python3 check_finite.py DELTA OUT.json [Lmax]
Diagnostics in floating point, not certified."""
import sys, json, time
import numpy as np
from scipy.integrate import cumulative_trapezoid
from bethe_scan import solve
from continuum_rep import Model, Restricted, Q_of_nu, e_bulk

delta = float(sys.argv[1]); out = sys.argv[2]
Lmax = int(sys.argv[3]) if len(sys.argv) > 3 else 1024
M = Model(delta); g = M.g
cache = {}
def bethe(L, n):
    if (L, n) not in cache:
        if n == 0:
            cache[(L, n)] = (np.longdouble(L*delta/4), np.array([]))
        else:
            E, d, x = solve(L, n, delta, return_roots=True)
            cache[(L, n)] = (E, x)
    return cache[(L, n)]

def counting_density(x, roots, L):
    res = np.empty_like(x)
    for i in range(0, len(x), 4000):
        xs = x[i:i+4000]
        res[i:i+4000] = M.a1(xs) - M.K(xs[:, None]-roots[None, :]).sum(1)/L
    return res

rows = []
Ls = [L for L in [8, 9, 16, 17, 32, 33, 64, 65, 128, 129, 256, 257, 512, 513, 1024, 1025, 2048] if L <= Lmax]
for L in Ls:
    ns = sorted(set([1, 2, int(np.sqrt(L)), L//4, int(0.4*L), L//2-2, L//2-1, L//2]))
    for n in ns:
        if n < 1 or n > L//2: continue
        t0 = time.time()
        E, x = bethe(L, n); nu = n/L
        # continuum at mass nu
        if nu < 0.5:
            Q = Q_of_nu(M, nu); R = Restricted(M, Q); eb = R.E
            s_fun = R.s_full
        else:
            Q = np.inf; eb = e_bulk(M, 0.5); s_fun = M.r0
        # z_nu on a fine grid
        X = max(np.max(np.abs(x)), 0 if not np.isfinite(Q) else Q) + 6.0
        grid = np.linspace(0, X, 40001)
        sg = s_fun(grid)
        zg = cumulative_trapezoid(sg, grid, initial=0.0)
        znu = np.sign(x)*np.interp(np.abs(x), grid, zg)
        Cnu = np.clip(znu, -nu/2, nu/2)
        j = np.arange(1, n+1)
        Dplus = j/L - nu/2 - Cnu; Dminus = Dplus - 1/L
        LD = float(L*max(np.max(np.abs(Dplus)), np.max(np.abs(Dminus))))
        # counting density: global min and sup deviation
        xx = np.concatenate([np.linspace(-X, X, 20001), x])
        sig = counting_density(xx, x, L)
        dev = float(L*np.max(np.abs(sig - s_fun(np.abs(xx)))))
        sig_roots = sig[-n:]
        # energy
        LdF = float(L*(E - L*eb))
        row = dict(delta=delta, L=L, n=n, nu=nu, Q=float(Q) if np.isfinite(Q) else None,
                   L_supD=LD, B_over_2=M.B/2, uniform_ref=(1+M.qJ)/2,
                   L_sup_sigma_dev=dev, D0=M.B/np.tan(g)/np.pi,
                   min_sigma=float(sig.min()), min_sigma_roots=float(sig_roots.min()),
                   min_L_sigma_roots=float(L*sig_roots.min()),
                   L_F_minus_Le=LdF, xmax=float(np.max(x)))
        if 2*L <= 2*Lmax:
            E2, _ = bethe(2*L, 2*n)
            row['defect_over_M0'] = float(L*(E2-2*E)/0.5)
        rows.append(row)
        print(json.dumps(row), f'{time.time()-t0:.1f}s', flush=True)
    json.dump(rows, open(out, 'w'), indent=1)
