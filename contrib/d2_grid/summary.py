#!/usr/bin/env python3
"""Tables for README.md: survey counts, screening statistics, FES verdicts for every flagged point.

FES verdict for a point (run with the most converged chi values among its cells):
  c1       c_fit within 0.06 of 1 over at least five converged chi, kappa in [1.15, 1.45]
  gapped   the last two pair slopes are below 0.5 and falling (S saturates while xi grows)
  c1*      not c1, but the median of the last five pair slopes is within 0.06 of 1 and only one of those
           five deviates from 1 by more than 0.3 (a single jump when the spiral unlocks in the cell)
  open     anything else (listed individually)
Writes data/summary.json.
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from analyze import fes

cand = json.load(open(os.path.join(HERE, 'candidates.json')))
best = {}
for f in glob.glob(os.path.join(HERE, 'data', 'fes', '*_n*.json')):
    o = fes(f)
    if o['name'] not in cand:
        continue
    if o['name'] not in best or len(o['chi']) > len(best[o['name']]['chi']) or \
            (len(o['chi']) == len(best[o['name']]['chi']) and o['chi'] and o['chi'][-1] > best[o['name']]['chi'][-1]):
        best[o['name']] = o


def verdict(o):
    cp = o.get('c_pair', [])
    if len(o['chi']) >= 5 and abs(o.get('c_fit', 9) - 1) < 0.06 and 1.15 <= o.get('kappa', 0) <= 1.45:
        return 'c1'
    last = np.array(cp[-5:])
    if len(last) == 5 and abs(np.median(last) - 1) < 0.06 and np.sum(np.abs(last - 1) > 0.3) <= 1:
        return 'c1*'
    if len(cp) >= 3 and cp[-1] < 0.5 and cp[-2] < 0.7 and cp[-1] <= cp[-2] + 0.05:
        return 'gapped'
    return 'open'


rows = []
for k, o in sorted(best.items()):
    rows.append(dict(name=k, n=o['n'], chi=o['chi'][-1] if o['chi'] else 0, xi=o['xi'][-1] if o['xi'] else 0,
                     c_fit=o.get('c_fit'), spread=o.get('c_spread'), kappa=o.get('kappa'),
                     c_pair=o.get('c_pair', []), verdict=verdict(o), nchi=len(o['chi']),
                     c_med=float(np.median(o['c_pair'][-5:])) if len(o.get('c_pair', [])) >= 5 else None))
missing = [k for k in cand if k not in best]
cnt = {v: sum(r['verdict'] == v for r in rows) for v in ('c1', 'c1*', 'gapped', 'open')}
c1 = [r['c_fit'] for r in rows if r['verdict'] == 'c1']
print('candidates', len(cand), 'with FES', len(rows), 'missing', missing)
print('verdicts', cnt)
if c1:
    print('c = 1 group: c_fit mean %.3f, std %.3f, min %.3f, max %.3f' % (np.mean(c1), np.std(c1), min(c1), max(c1)))
for r in rows:
    if r['verdict'] not in ('c1', 'gapped'):
        print(f"{r['name']:24s} {r['verdict']:7s} n={r['n']:2d} nchi={r['nchi']} chi={r['chi']:4.0f} xi={r['xi']:8.1f} "
              f"c_fit={r['c_fit'] if r['c_fit'] is None else round(r['c_fit'],3)} kappa={r['kappa'] if r['kappa'] is None else round(r['kappa'],2)} "
              f"c_med={r['c_med'] if r['c_med'] is None else round(r['c_med'],3)} c_pair={np.round(r['c_pair'],2).tolist()}")
json.dump(dict(rows=rows, missing=missing, counts=cnt), open(os.path.join(HERE, 'data', 'summary.json'), 'w'), indent=1)
