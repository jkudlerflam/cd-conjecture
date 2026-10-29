#!/usr/bin/env python3
"""Evaluate A/(pi M log2 d) along parametrized bond families.

    python3 scan.py --family tj --par mu=0.5 --Ls 6 8 10 12 --out data/tj.jsonl

Each run appends one JSON line: family, parameters, M, E_0(L), fits for A, the gap times L,
the velocity estimate from the current state and c_half, plus the ratio.  Symmetric charge
columns (spin-like, listed per family) are restricted to the two smallest |q| values; filling
charges are scanned in full so the chemical potential picks the density.
"""
import argparse, json, time, sys
import numpy as np
import families as F
from symring import symring
from core import hnorm, fit_cas, LOG2

FAMILIES = {
    'heis': (lambda **p: F.heisenberg(), []),
    'xxz': (F.xxz, [0]),
    'tj': (F.tj, [1]),
    'osp12': (F.osp12, [0]),
    'blbq': (F.blbq, [0]),
    'spin1xxzd': (F.spin1_xxz_d, [0]),
    'hubbard': (F.hubbard, []),
    'kk': (F.kugel_khomskii, [0, 1]),
    'twochain': (F.two_chain_coupled, [0, 1]),
    'gperm': (F.graded_perm, []),
}


def parse_par(items):
    out = {}
    for it in items:
        k, v = it.split('=')
        try:
            out[k] = json.loads(v)
        except Exception:
            out[k] = v
    return out


def run(family, par, Ls, nlev=3, entropy=True, tol=1e-10, symcols=None):
    make, defsym = FAMILIES[family]
    h, d, charges = make(**par)
    charges = np.asarray(charges)
    if charges.ndim == 1:
        charges = charges[:, None]
    symcols = defsym if symcols is None else symcols

    def filt(t, L):
        for c in symcols:
            if c < len(t) and abs(t[c]) > 1:
                return False
        return True

    t0 = time.time()
    M = float(hnorm(h, d))
    rows = []
    for L in Ls:
        r = symring(h, d, L, charges, nlev=nlev, tol=tol, entropy=entropy,
                    target_filter=filt if charges.shape[1] else None)
        rows.append(r)
    E0 = [r['E0'] for r in rows]
    A, e0, b = fit_cas(Ls, E0)
    A2 = fit_cas(Ls[-2:], E0[-2:], quartic=False)[0]
    res = {
        'family': family, 'par': par, 'd': d, 'Ls': list(Ls), 'E0': E0, 'k0': [r['k0'] for r in rows],
        'target': [list(r['target']) for r in rows], 'gapL': [r['gap'] * r['L'] for r in rows],
        'v_est': [r['v_est'] for r in rows], 'M': M, 'A_fit': A, 'e0_fit': e0, 'A_last2': A2,
        'ratio': A / (np.pi * M * np.log(d) / LOG2), 'ratio_last2': A2 / (np.pi * M * np.log(d) / LOG2),
        'time': time.time() - t0,
    }
    if len(Ls) >= 3:
        res['A_last3'] = fit_cas(Ls[-3:], E0[-3:])[0]
    if entropy:
        S = [r['S_half'] for r in rows]
        res['S_half'] = S
        if len(Ls) >= 2:
            res['c_half'] = 3 * (S[-1] - S[-2]) / np.log(Ls[-1] / Ls[-2])
            res['A_cv'] = res['c_half'] * rows[-1]['v_est']
            res['ratio_cv'] = res['A_cv'] / (np.pi * M * np.log(d) / LOG2)
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True)
    ap.add_argument('--par', nargs='*', default=[])
    ap.add_argument('--Ls', nargs='+', type=int, required=True)
    ap.add_argument('--nlev', type=int, default=3)
    ap.add_argument('--no-entropy', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = run(a.family, parse_par(a.par), a.Ls, nlev=a.nlev, entropy=not a.no_entropy)
    with open(a.out, 'a') as f:
        f.write(json.dumps(res) + '\n')
    print(json.dumps({k: res[k] for k in ('family', 'par', 'M', 'A_fit', 'A_last2', 'ratio', 'ratio_last2', 'c_half', 'gapL', 'time') if k in res}))
