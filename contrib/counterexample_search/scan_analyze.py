#!/usr/bin/env python3
"""Rank scan points by entanglement slope and attach CFT diagnostics.

For each point and size: c_ent (fit l >= 3), lowest gap times L/(2 pi), ground
momentum, ground-state degeneracy across momenta. A point is FLAGGED when, at the
largest size, c_ent exceeds the SU(d)_1 calibration value at the same size
(c_ref), c_ent does not decrease between the two largest sizes by more than 0.1,
and z_eff of the lowest gap between the two largest sizes lies in [0.7, 1.4].
Usage: python3 scan_analyze.py --d 4 --cref 3.47 --top 25 files...
"""
import argparse, json, glob
from collections import defaultdict
import numpy as np
from est import c_ent, lowest_gap, v_adj

p = argparse.ArgumentParser()
p.add_argument('files', nargs='+')
p.add_argument('--d', type=int, required=True)
p.add_argument('--cref', type=float, required=True, help='c_ent of SU(d)_1 at the largest size')
p.add_argument('--top', type=int, default=25)
p.add_argument('--lmin', type=int, default=3)
p.add_argument('--json', default=None)
p.add_argument('--pair', type=int, nargs=2, default=None)
a = p.parse_args()

pts = defaultdict(dict)
coeffs = {}
for f in a.files:
    for line in open(f):
        r = json.loads(line)
        pts[r['point']][r['L']] = r
        coeffs[r['point']] = r['coeff']

rows = []
for i, byL in pts.items():
    Ls = sorted(byL)
    if len(Ls) < 2:
        continue
    L1, L2 = (Ls[-2], Ls[-1]) if a.pair is None else tuple(a.pair)
    if L1 not in byL or L2 not in byL:
        continue
    r1, r2 = byL[L1], byL[L2]
    ce = {L: c_ent(byL[L], a.lmin) for L in Ls}
    g1, g2 = lowest_gap(r1), lowest_gap(r2)
    z = float(-np.log(g2 / g1) / np.log(L2 / L1)) if g1 > 0 and g2 > 0 else np.nan
    ndeg = len(r2['ground_degenerate_momenta'])
    chalf = float(3 * (r2['S'][-1] - r1['S'][-1]) / np.log(L2 / L1))
    flag = (ce[L2] > a.cref) and (ce[L2] >= ce[L1] - 0.1) and (0.7 <= z <= 1.4)
    # ferro-multiplet artefact: entropy identical across sizes up to the polarized-multiplet value
    ferro = abs(r2['S'][-1] - r1['S'][-1]) < 0.15 and z > 1.6
    rows.append(dict(point=i, coeff=coeffs[i], ce=ce, z=z, gapL=g2 * L2 / (2 * np.pi), k0=r2['k0'], ndeg=ndeg,
                     flag=bool(flag), E0L=r2['E0'] / L2, chalf=chalf, ferro=bool(ferro)))
rows.sort(key=lambda x: -x['chalf'])
nflag = sum(x['flag'] for x in rows)
print('points %d, flagged %d (c_ref %.3f)' % (len(rows), nflag, a.cref))
for x in rows[:a.top]:
    Ls = sorted(x['ce'])
    print('%4d %-28s ' % (x['point'], ' '.join('%6.3f' % c for c in x['coeff'])) +
          ' '.join('c%d=%5.2f' % (L, x['ce'][L]) for L in Ls) +
          '  chalf=%5.2f z=%5.2f gapL=%6.3f k0=%d ndeg=%d %s %s' % (x['chalf'], x['z'], x['gapL'], x['k0'], x['ndeg'], 'FLAG' if x['flag'] else '', 'ferro?' if x['ferro'] else ''))
if a.json:
    json.dump(rows, open(a.json, 'w'))
