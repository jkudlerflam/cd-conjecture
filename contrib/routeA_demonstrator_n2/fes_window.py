#!/usr/bin/env python3
"""FES and block-entropy central charges in a common correlation-length window, for the microscopic chain and
the calibration chains (c = 2: two XX chains, effective chain; c = 1: one XX chain).
usage: fes_window.py XIMIN XIMAX.  Claude (Opus 5.5), 2026-10-28."""
import os, sys, json, glob, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def points(names):
    P = []
    for nm in names:
        f = os.path.join(HERE, 'data', 'fes', nm + '.json')
        if os.path.exists(f):
            P += [(x['chi'], x['xi_cells'], x['S'][0]) for x in json.load(open(f))['results'] if x['err'] < 1e-4]
    return sorted(set(P))


def fit(P, lo, hi):
    Q = [p for p in P if lo <= p[1] <= hi]
    if len(Q) < 2:
        return None, len(Q)
    x = np.log([q[1] for q in Q]); y = np.array([q[2] for q in Q])
    return float(6 * np.polyfit(x, y, 1)[0]), len(Q)


def block(prefix):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'block', prefix + '_chi*.json'))):
        d = json.load(open(f)); S = dict(d['S'])
        out.append((d['chi'], round(d['xi_cells'], 2), round(3 * (S[2] - S[1]) / math.log(2), 3)))
    return sorted(out)


if __name__ == '__main__':
    lo, hi = float(sys.argv[1]), float(sys.argv[2])
    sets = {'XX (c=1)': ['xx1s'], 'two XX (c=2)': ['xx2s', 'xx2'], 'effective 0.15 (c=2)': ['effgs_l0.15', 'effg_l0.15'],
            'effective 0.20 (c=2)': ['effgs_l0.20', 'effg_l0.20'], 'microscopic 0.15': ['micg_l0.15'], 'microscopic 0.20': ['micg_l0.20']}
    res = {}
    for k, v in sets.items():
        c, n = fit(points(v), lo, hi)
        res[k] = dict(c_fes=c, npts=n)
        print(f'{k:22s} FES c in xi [{lo}, {hi}] = {c if c is None else round(c, 3)} ({n} points)')
    for k, p in [('XX (c=1)', 'xx1s'), ('two XX (c=2)', 'xx2s'), ('effective 0.15 (c=2)', 'effgs_l0.15'),
                 ('effective 0.20 (c=2)', 'effgs_l0.20'), ('microscopic 0.15', 'micg_l0.15'), ('microscopic 0.20', 'micg_l0.20')]:
        b = block(p)
        res[k]['block'] = b
        print(f'{k:22s} block c(1,2) (chi, xi, c):', b)
    json.dump(res, open(os.path.join(HERE, 'data', f'window_{lo}_{hi}.json'), 'w'), indent=1)
