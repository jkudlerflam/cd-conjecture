#!/usr/bin/env python3
"""Fits for the refinement stage.
  fes:  data/fes/NAME_nN.json -> c_pair = 6 dS/dln xi between consecutive chi (S = cell mean),
        c_fit over the last 4 chi (spread over shifted windows of 4), kappa = dln xi/dln chi.
  obc:  data/obc/NAME_L*_s*[_edge].json (largest chi) -> window fits of d2_followup_220_1282/analyze_obc.py
        (plain, bond-energy, oscillating term, smoothed over one period).
usage: analyze.py fes [files]   |   analyze.py obc [files]
"""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np


def fes(path, nfit=4, errmax=1e-3):
    d = json.load(open(path))
    byc = {}
    for r in d['results']:
        byc[r['chi']] = r
    rs = [byc[c] for c in sorted(byc) if byc[c]['err'] < errmax and np.isfinite(byc[c]['xi'])]
    S = np.array([np.mean(r['S']) for r in rs]); xi = np.array([r['xi'] for r in rs]); chi = np.array([r['chi'] for r in rs])
    o = dict(name=d['name'], n=d['n'], chi=chi.tolist(), xi=xi.tolist(), S=S.tolist())
    if len(rs) >= 2:
        o['c_pair'] = (6 * np.diff(S) / np.diff(np.log(xi))).tolist()
    if len(rs) >= 3:
        k = min(nfit, len(rs))
        fits = [6 * np.polyfit(np.log(xi[i:i + k]), S[i:i + k], 1)[0] for i in range(len(rs) - k + 1)]
        o['c_fit'] = float(fits[-1]); o['c_spread'] = float(np.ptp(fits[-3:])) if len(fits) > 1 else 0.0
        o['kappa'] = float(np.polyfit(np.log(chi[-k:]), np.log(xi[-k:]), 1)[0])
    return o


def obc(path):
    sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
    from analyze_obc import fits          # plain, bond-energy, oscillating-term and one-period-smoothed fits
    d = json.load(open(path))
    r = d['results'][-1]
    L = d['L']
    f = fits(np.array(r['S']), np.array(r['Eb']), L)
    return dict(name=d['name'], L=L, chi=r['chi'], edge=d['edge'], trunc=r['trunc'], c_plain=float(f['c_plain']),
                c_eb=float(f['c_Eb']), c_osc=float(f['c_osc']), k_osc=float(f['k_osc']), c_smooth=float(f['c_smooth']),
                Sbar=f['Sbar'])


if __name__ == '__main__':
    kind = sys.argv[1]
    files = sys.argv[2:] or sorted(glob.glob(os.path.join(HERE, 'data', kind, '*.json')))
    for f in files:
        if f.endswith('_cells.json'):
            continue
        o = fes(f) if kind == "fes" else obc(f)
        if kind == "fes" and not o["chi"]:
            continue
        if kind == 'fes':
            print(f"{o['name']:22s} n={o['n']:2d} chi={o['chi'][-1]:4.0f} xi={o['xi'][-1]:8.1f} c_fit={o.get('c_fit', np.nan):.3f}"
                  f" spread={o.get('c_spread', np.nan):.3f} kappa={o.get('kappa', np.nan):.2f} c_pair={np.round(o.get('c_pair', []), 3).tolist()}")
        else:
            print(f"{o['name']:22s} L={o['L']:4d} chi={o['chi']:4d} edge={o['edge']!s:5s} trunc={o['trunc']:.1e} c_plain={o['c_plain']:.3f} "
                  f"c_eb={o['c_eb']:.3f} c_osc={o['c_osc']:.3f} (k={o['k_osc']:.3f}) c_smooth={o['c_smooth']:.3f}")
