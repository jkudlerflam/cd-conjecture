#!/usr/bin/env python3
"""Both estimators of A for every (model, par) in the given JSON-lines files.

  python3 analyze.py data/calib.jsonl [--Lmod 3] [--Lmin 8]

A_cas    fit of E_0(L)/L = e0 - pi A/(6L^2) + b/L^4 over all sizes kept
A_cas3   the same fit on the three largest sizes (drift check)
v        L * (lowest level at k0 +- 2pi/L - E0)/(2 pi), extrapolated linearly in 1/L^2 over the
         two largest sizes (Ising-type chains without a current: subtract the k0 level below it)
c_half   3 [S(L2/2) - S(L1/2)] / ln(L2/L1) over the two largest sizes
A_cv     c_half * v
"""
import argparse, json
from collections import defaultdict
import numpy as np
import core


def load(paths):
    recs = defaultdict(dict)
    for p in paths:
        for line in open(p):
            r = json.loads(line)
            recs[(r['model'], tuple(r['par']))][r['L']] = r
    return recs


def level(r, k, n=0, target='gs'):
    t = r['target'] if target == 'gs' else target
    for tg, kk, e in r['spec']:
        if kk % r['L'] == k % r['L'] and tg == t and len(e) > n:
            return e[n]
    return np.nan


def vel(r, nocurrent=False):
    L, k0 = r['L'], r['k0']
    e1 = np.nanmin([level(r, k0 + 1), level(r, k0 - 1)])
    base = level(r, k0, 1) if nocurrent else r['E0']
    return L * (e1 - base) / (2 * np.pi)


def analyze(rs, Lmod=1, Lmin=0, Lpar=None, nocurrent=False, e0=None):
    Ls = sorted(L for L in rs if L % Lmod == 0 and L >= Lmin and (Lpar is None or L % 2 == Lpar))
    E = [rs[L]['E0'] for L in Ls]
    out = dict(Ls=Ls, hnorm=rs[Ls[-1]]['hnorm'])
    if len(Ls) >= 3:
        out['A_cas'], out['e0'], _ = core.fit_cas(Ls, E)
        out['A_cas3'] = core.fit_cas(Ls[-3:], E[-3:])[0]
        out['A_loc'] = core.cas_local(Ls, E)[0]
    if e0 is not None:
        out['A_e0'] = [-6 * L * (rs[L]['E0'] - e0 * L) / np.pi for L in Ls]
    v = [vel(rs[L], nocurrent) for L in Ls]
    out['v_L'] = v
    if len(Ls) >= 2:
        x = 1 / np.array(Ls[-2:], float) ** 2
        out['v'] = float(v[-1] - (v[-1] - v[-2]) / (x[1] - x[0]) * x[1])
        S = [rs[L]['S_half'] for L in Ls]
        if S[-1] is not None and S[-2] is not None:
            out['c_half'] = core.c_half(S[-2], Ls[-2], S[-1], Ls[-1])
            out['A_cv'] = out['c_half'] * out['v']
    for key in ('A_cas', 'A_cv'):
        if key in out:
            out['R_' + key[2:]] = out[key] / (2 * np.pi * out['hnorm'])
    return out


def fmt(d):
    s = []
    for k, v in d.items():
        if isinstance(v, float):
            s.append('%s=%.4f' % (k, v))
        elif isinstance(v, list) and v and isinstance(v[0], float):
            s.append('%s=[%s]' % (k, ' '.join('%.4f' % x for x in v)))
        else:
            s.append('%s=%s' % (k, v))
    return ' '.join(s)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--Lmod', type=int, default=1)
    ap.add_argument('--Lmin', type=int, default=0)
    ap.add_argument('--Lpar', type=int, default=None)
    ap.add_argument('--nocurrent', action='store_true')
    ap.add_argument('--model', default=None)
    a = ap.parse_args()
    for key, rs in sorted(load(a.files).items()):
        if a.model and key[0] != a.model:
            continue
        try:
            print(key, fmt(analyze(rs, a.Lmod, a.Lmin, a.Lpar, a.nocurrent)))
        except Exception as e:  # noqa
            print(key, 'failed', e)
