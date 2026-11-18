"""Bethe scan on Delta > 1: doubling defect, sector error against my own continuum, and the Section 4 diagnostics."""
import numpy as np, json, sys, time
from multiprocessing import Pool
from ref148 import *

def kKe(lam):
    q = np.exp(-lam); n = np.arange(0, 200)
    th2 = 2 * np.sum(q ** ((n + 0.5) ** 2))
    return pi / 2 * th2 ** 2

def fillings(L):
    c = {1, 2, max(1, L // 8), max(1, L // 4), max(1, L // 3), L // 2 - 3, L // 2 - 2, L // 2 - 1, L // 2}
    return sorted(n for n in c if 1 <= n <= L // 2)

def diag(L, n, lam, x):
    out = {}
    if n >= 2:
        xg = np.concatenate([np.linspace(x[j], x[j + 1], 24) for j in range(n - 1)]) if n <= 300 else np.linspace(x[0], x[-1], 8000)
        sg = sigma_L(xg, x, L, lam)
        r0n = float(r0(x[-1], lam)[0])
        out['min_sigma_over_r0xn_in'] = float(sg.min() / r0n)
        out['min_sigma_over_r0_in'] = float((sg / r0(xg, lam)).min())
        out['L_r0_xn'] = L * r0n
        out['Loc_bound'] = 4 * np.sin(pi / 8) * kKe(lam) / pi ** 2
        yo = np.linspace(x[-1], pi - x[-1], 4001)
        so = sigma_L(yo, x, L, lam)
        out['min_sigma_out'] = float(so.min())
        out['L_negmass_out_half'] = float(L * np.trapezoid(np.minimum(so, 0), yo) / 2)
        # sawtooth: z_L at roots must equal I_j/L
        out['z_err'] = float(np.abs(z_L(x, x, L, lam) - (np.arange(n) - (n - 1) / 2) / L).max())
    return out

def run(D):
    lam = float(np.arccosh(D)); m0 = M0(lam); rows = []
    Ls = [2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 16, 17, 32, 33, 64, 65, 128, 129, 256, 257, 512, 513, 1024]
    do_bulk = lam >= 0.044
    cache = {}
    for L in Ls:
        for n in fillings(L):
            try:
                E1, x1 = bethe_energy(L, n, lam)
                E2, x2 = bethe_energy(2 * L, 2 * n, lam)
            except Exception as ex:
                rows.append(dict(Delta=D, L=L, n=n, fail=str(ex))); continue
            row = dict(Delta=D, lam=lam, L=L, n=n, defect=L * (E2 - 2 * E1) / m0)
            nu = n / L
            if do_bulk or 2 * n == L:
                if nu not in cache:
                    cache[nu] = ebulk_half(lam) if 2 * n == L else ebulk(nu, lam)
                eb = cache[nu]
                row['secerr_L'] = L * (E1 - L * eb) / m0
                row['secerr_2L'] = 2 * L * (E2 - 2 * L * eb) / m0
            row.update(diag(L, n, lam, x1))
            d2 = diag(2 * L, 2 * n, lam, x2)
            row.update({k + '_2L': v for k, v in d2.items()})
            rows.append(row)
    return rows

if __name__ == '__main__':
    Deltas = [1.0001, 1.001, 1.01, 1.05, 1.2, 1.5, 2.0, float(np.cosh(1.36)), 3.0, 5.0, 10.0, 20.0, 100.0]
    t0 = time.time()
    with Pool(7) as p:
        res = p.map(run, Deltas)
    rows = [r for rr in res for r in rr]
    with open('scan_bethe.jsonl', 'w') as f:
        for r in rows: f.write(json.dumps(r) + '\n')
    print('rows', len(rows), 'fails', sum('fail' in r for r in rows), 'time', time.time() - t0)
