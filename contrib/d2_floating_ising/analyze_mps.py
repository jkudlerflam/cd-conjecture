#!/usr/bin/env python3
"""Fits for the two calibrated estimators.
FES (data/fes/*.json): c = 6 dS/dln xi, joint fit over the largest four chi, and consecutive pair slopes;
  S is the mean over the bonds of the cell, runs with gradient error above 1e-3 are dropped;
  kappa = dln xi/dln chi over the same window (z = 1 CFT: 1.34, 1.05 for c = 1, 3/2).
OBC (data/obc/*.json): S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + const, fitted separately on even and odd
  bonds over L/8 <= l <= 7L/8 and averaged (largest chi); c_dbl = 6 [Sbar(L2) - Sbar(L1)]/ln(L2/L1) from
  the window means at two sizes.
usage: analyze_mps.py [--json data/mps_summary.json]
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, glob, argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def fes(path, nfit=4, errmax=1e-3):
    d = json.load(open(path))
    R = [r for r in d['results'] if r['err'] < errmax and np.isfinite(r['xi']) and r['xi'] > 0]
    out = dict(name=d['name'], n=d['n'], chis=[r['chi'] for r in R], xi=[r['xi'] for r in R],
               S=[float(np.mean(r['S'])) for r in R], err=[r['err'] for r in R])
    if len(R) >= 2:
        lx = np.log(out['xi']); S = np.array(out['S']); lc = np.log(out['chis'])
        out['pairs'] = [float(6 * (S[i + 1] - S[i]) / (lx[i + 1] - lx[i])) if abs(lx[i + 1] - lx[i]) > 1e-6 else None
                        for i in range(len(R) - 1)]
        k = min(nfit, len(R))
        out['c_fit'] = float(6 * np.polyfit(lx[-k:], S[-k:], 1)[0])
        out['kappa'] = float(np.polyfit(lc[-k:], lx[-k:], 1)[0])
        out['xi_max'] = float(out['xi'][-1])
    return out


def obc_fit(S, L):
    S = np.asarray(S); l = np.arange(1, L)
    x = np.log(2 * L / np.pi * np.sin(np.pi * l / L)) / 6
    w = (l >= L / 8) & (l <= 7 * L / 8)
    cs = []
    for par in (0, 1):
        m = w & (l % 2 == par)
        cs.append(np.polyfit(x[m], S[m], 1)[0])
    return float(np.mean(cs)), float(np.mean(S[w & (np.abs(l - L / 2) <= L / 8)]))


def obc(paths):
    rows = {}
    for p in paths:
        d = json.load(open(p))
        r = d['results'][-1]
        c, Sbar = obc_fit(r['S'], d['L'])
        rows.setdefault(d['name'], {})[d['L']] = dict(c_obc=c, Sbar=Sbar, chi=r['chi'], trunc=r['trunc'])
    for name, R in rows.items():
        Ls = sorted(R)
        if len(Ls) >= 2:
            a, b = Ls[0], Ls[-1]
            R['c_dbl'] = float(6 * (R[b]['Sbar'] - R[a]['Sbar']) / np.log(b / a))
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', default=os.path.join(HERE, 'data', 'mps_summary.json'))
    a = ap.parse_args()
    F = {os.path.basename(p)[:-5]: fes(p) for p in sorted(glob.glob(os.path.join(HERE, 'data', 'fes', '*.json')))}
    O = obc(sorted(glob.glob(os.path.join(HERE, 'data', 'obc', '*.json'))))
    json.dump(dict(fes=F, obc=O), open(a.json, 'w'), indent=1, default=float)
    print('FES')
    for k, f in F.items():
        if 'c_fit' in f:
            print(f"  {k:60s} n={f['n']:2d} chi<={f['chis'][-1]:4d} xi={f['xi_max']:9.1f} c={f['c_fit']:.3f} "
                  f"kappa={f['kappa']:.2f} pairs={[round(x, 2) for x in f['pairs'][-4:] if x is not None]}")
    print('OBC')
    for k, R in O.items():
        print(f'  {k:60s}', {L: round(v['c_obc'], 3) for L, v in R.items() if isinstance(L, int)},
              'c_dbl', round(R.get('c_dbl', np.nan), 3))
