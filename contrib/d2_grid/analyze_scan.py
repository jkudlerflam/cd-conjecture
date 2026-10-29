#!/usr/bin/env python3
"""Summaries of the sphere scans: distribution of c_S(16), gapless fraction, local maxima.

Per point: c_S(L) (ring chord fit), the scaled gaps g_L = L (E_n - E_0) for the first three excitations
among all momentum sectors (two levels per sector), and the ratio g_16/g_12 of the first gap that is
not a near-degeneracy (E_n - E_0 > 1e-3 at L = 16).  Classes:
  gapped     c_S(16) < 0.25 or the scaled gap grows by more than 1.25 from L = 12 to 16
  ssb        ground-state near-degeneracy and flat entanglement
  gapless    otherwise (c_S(16) >= 0.25 and scaled gap not growing)
Local maxima of c_S(16) on each sphere use the 8 nearest Fibonacci neighbours.
usage: analyze_scan.py [--top 40] [--json data/scan_summary.json]
"""
import os, sys, json, glob, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import strata

ap = argparse.ArgumentParser()
ap.add_argument('--top', type=int, default=40)
ap.add_argument('--json', default=os.path.join(HERE, 'data', 'scan_summary.json'))
ap.add_argument('--thr', type=float, default=1.1)
a = ap.parse_args()


def gaps(E, L):
    e = np.sort(np.concatenate([np.asarray(x) for x in E]))
    d = e - e[0]
    return L * d[1:4], d


def classify(r):
    c16 = r['L16']['cS']
    g12, d12 = gaps(r['L12']['E'], 12)
    g16, d16 = gaps(r['L16']['E'], 16)
    nd = np.flatnonzero(d16[1:] > 1e-3)
    if len(nd) == 0:
        return 'ssb', np.nan
    n = nd[0] + 1
    ratio = (16 * d16[n]) / max(12 * d12[n], 1e-12)
    if c16 < 0.25 or ratio > 1.25:
        return 'gapped', ratio
    if n > 1 and c16 < 0.5:
        return 'ssb', ratio
    return 'gapless', ratio


recs = {}
for f in sorted(glob.glob(os.path.join(HERE, 'data', 'scan', '*.jsonl'))):
    key = os.path.basename(f)[:-6]
    rows = [json.loads(l) for l in open(f) if l.strip()]
    recs[key] = rows

summary = dict(spheres={}, top=[], maxima=[])
allpts = []
for key, rows in recs.items():
    if key == 'stag' or not rows or 'L16' not in rows[0]:
        continue
    st, seed = key.split('_')[0], key.split('_')[1]
    if st not in strata.COORDS:
        continue
    N = max(int(r['tag'].split('_')[-1]) for r in rows) + 1
    X = strata.fibonacci({'NS': 600, 'Z2': 400, 'TR': 500, 'RA': 500}[st])
    c16 = {int(r['tag'].split('_')[-1]): r['L16']['cS'] for r in rows}
    c12 = {int(r['tag'].split('_')[-1]): r['L12']['cS'] for r in rows}
    cls = {}
    for r in rows:
        i = int(r['tag'].split('_')[-1])
        cl, ratio = classify(r)
        cls[i] = cl
        allpts.append(dict(tag=r['tag'], st=st, c16=r['L16']['cS'], c12=r['L12']['cS'], cls=cl, ratio=ratio,
                           v=r['v'], m0=r['L16']['m0']))
    idx = np.array(sorted(c16))
    Xi = X[idx]
    D = Xi @ Xi.T
    nb = np.argsort(-D, 1)[:, 1:9]
    vals = np.array([c16[i] for i in idx])
    for k, i in enumerate(idx):
        if vals[k] >= vals[nb[k]].max() and vals[k] > 0.8:
            summary['maxima'].append(dict(tag=f'{st}_{seed}_{i}', c16=vals[k], c12=c12[i], cls=cls[i]))
    nc = {c: sum(1 for x in cls.values() if x == c) for c in ('gapped', 'ssb', 'gapless')}
    summary['spheres'][key] = dict(n=len(rows), max_c16=float(vals.max()), **nc)

byst = {}
for p in allpts:
    byst.setdefault(p['st'], []).append(p)
print('stratum  points  gapless  ssb  gapped  max c_S(16)  #c_S(16)>1.1  #>1.2  #>1.3  max over gapless with c12<c16')
for st, ps in byst.items():
    c = np.array([p['c16'] for p in ps])
    gl = [p for p in ps if p['cls'] == 'gapless']
    up = [p['c16'] for p in gl if p['c16'] > p['c12']]
    print(f"{st:6s} {len(ps):7d} {len(gl):8d} {sum(p['cls']=='ssb' for p in ps):5d} {sum(p['cls']=='gapped' for p in ps):7d}"
          f"  {c.max():10.3f}  {np.sum(c>1.1):12d} {np.sum(c>1.2):6d} {np.sum(c>1.3):6d}  {max(up) if up else float('nan'):.3f}")
top = sorted(allpts, key=lambda p: -p['c16'])[:a.top]
print('\ntop points by c_S(16):')
for p in top:
    print(f"{p['tag']:12s} c12={p['c12']:.3f} c16={p['c16']:.3f} {p['cls']:8s} ratio={p['ratio']:.2f} m0={p['m0']}")
summary['top'] = top
summary['n_points'] = len(allpts)
summary['hist'] = {st: np.histogram([p['c16'] for p in ps], bins=np.arange(-0.1, 2.01, 0.1))[0].tolist() for st, ps in byst.items()}
json.dump(summary, open(a.json, 'w'), indent=1, default=float)
