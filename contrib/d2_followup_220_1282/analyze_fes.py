#!/usr/bin/env python3
"""Finite-entanglement scaling fits from the data/fes_*.json files written by run_fes.py.

For every run: S = mean entropy over the bonds of the cell, xi = -n / ln|lambda_1/lambda_0| of the
cell transfer matrix.  Reported:
  c_pair(chi)  = 6 dS / d ln xi between consecutive chi,
  c_fit        = 6 x slope of S against ln xi over the last `nfit` chi values (and its spread
                 over windows shifted by one point),
  kappa        = slope of ln xi against ln chi over the same points (FES theory, z = 1 CFT:
                 kappa = 6 / (c (sqrt(12/c) + 1))).
usage: analyze_fes.py [files...] [--nfit 4] [--json out.json]
"""
import sys, json, glob, os
import numpy as np


def kappa_th(c):
    return 6 / (c * (np.sqrt(12 / c) + 1))


def analyse(path, nfit=4, errmax=1e-3):
    d = json.load(open(path))
    byc = {}
    for r in d['results']:                      # a continued run repeats its first chi: keep the last
        byc[r['chi']] = r
    rs = [byc[c] for c in sorted(byc)]
    rs = [r for r in rs if r['err'] < errmax and np.isfinite(r['xi']) and r['xi'] < 1e7]
    chi = np.array([r['chi'] for r in rs], float)
    S = np.array([np.mean(r['S']) for r in rs])
    xi = np.array([r['xi'] for r in rs])
    e = np.array([r['e'] for r in rs])
    out = dict(file=os.path.basename(path), model=d['model'], n=d['n'], mode=d['mode'], c_exact=d.get('c_exact'),
               chi=chi.tolist(), S=S.tolist(), xi=xi.tolist(), e=e.tolist(), err=[r['err'] for r in rs])
    if len(rs) >= 2:
        out['c_pair'] = (6 * np.diff(S) / np.diff(np.log(xi))).tolist()
    if len(rs) >= 3:
        k = min(nfit, len(rs))
        fits = []
        for off in range(0, len(rs) - k + 1):
            sl = slice(off, off + k)
            fits.append(6 * np.polyfit(np.log(xi[sl]), S[sl], 1)[0])
        out['c_fit'] = fits[-1]
        out['c_fit_windows'] = fits
        out['kappa'] = float(np.polyfit(np.log(chi[-k:]), np.log(xi[-k:]), 1)[0])
    return out


if __name__ == '__main__':
    args = sys.argv[1:]
    nfit = 4
    jout = None
    if '--nfit' in args:
        i = args.index('--nfit'); nfit = int(args[i + 1]); del args[i:i + 2]
    if '--json' in args:
        i = args.index('--json'); jout = args[i + 1]; del args[i:i + 2]
    files = args or sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'fes_*.json')))
    allres = []
    for f in files:
        o = analyse(f, nfit)
        allres.append(o)
        print(f"{o['model']:14s} n={o['n']:2d} {o['mode']:6s} c_exact={o['c_exact']}  chi_max={o['chi'][-1] if o['chi'] else '-'}")
        for i in range(len(o['chi'])):
            cp = o['c_pair'][i - 1] if i > 0 else float('nan')
            print(f"    chi={o['chi'][i]:5.0f} e={o['e'][i]:.10f} err={o['err'][i]:.0e} S={o['S'][i]:.5f} "
                  f"xi={o['xi'][i]:9.3f} c_pair={cp:6.3f}")
        if 'c_fit' in o:
            print(f"    c_fit(last {nfit}) = {o['c_fit']:.3f}   windows {np.round(o['c_fit_windows'], 3).tolist()}"
                  f"   kappa = {o['kappa']:.3f}")
    if jout:
        json.dump(allres, open(jout, 'w'), indent=1)
