#!/usr/bin/env python3
"""Search every computed point for the tricritical-singlet fingerprint (units 2 pi / L):
  phi   = lowest spin-0 level at k0 + pi               (x = 1/5 at the glued point)
  F     = lowest excited spin-0 level at k0            (x = 6/5; ratio F/phi = 6)
  dphi  = lowest spin-0 level at k0 + pi -+ 2pi/L       (descendant of phi, ratio 6 as well)
Control (O'Brien-Fendley TCI, controls/tci_obf.py): eps'/eps = 8.1, 6.9, 6.5 and descendant
ratio 6.8, 6.5, 6.3 at L = 8, 12, 16.  A candidate needs phi scale invariant between L = 8 and 12
(|phi12/phi8 - 1| < 0.2), phi soft compared with the SU(2) current, and both ratios above 4.5.
Prints all soft scale-invariant singlets, sorted by F/phi at L = 12."""
import glob, sys
import numpy as np
from analyze import load, levels

files = sys.argv[1:] or glob.glob('data/grid/*.jsonl') + glob.glob('data/grid_local/*.jsonl') + glob.glob('data/fine/*.jsonl')
recs = load(files)


def fp(r):
    L = r['L']
    lv = [x for x in levels(r) if x[0] > 1e-7]
    g = lambda dk, s: [x[0] for x in lv if x[1] == dk and x[2] == s]
    first = lambda l: l[0] if l else np.nan
    return dict(phi=first(g(L // 2, 0)), F=first(g(0, 0)), dphi=first(g(L // 2 - 1, 0)), cur=first(g(1, 1)),
                S=r['S'][-1], Sgs=levels(r)[0][2])


rows = []
for key, b in recs.items():
    if 8 not in b or 12 not in b:
        continue
    a, c = fp(b[8]), fp(b[12])
    if a['Sgs'] != 0 or c['Sgs'] != 0 or not np.isfinite(a['phi']) or not np.isfinite(c['phi']):
        continue
    inv = c['phi'] / a['phi']
    soft = c['phi'] / c['cur'] if np.isfinite(c['cur']) else np.nan
    if abs(inv - 1) < 0.2 and (not np.isfinite(soft) or soft < 0.6):
        rows.append((c['F'] / c['phi'], c['dphi'] / c['phi'], a['F'] / a['phi'], a['dphi'] / a['phi'], key[0] / np.pi, key[1],
                     inv, soft, 3 * (c['S'] - a['S']) / np.log(1.5)))
rows.sort(reverse=True)
print('%6s %6s | %6s %6s | %6s %6s | %5s %5s | %5s' % ('th/pi', 'eta', 'F/p12', 'd/p12', 'F/p8', 'd/p8', 'p12/8', 'p/cur', 'chalf'))
for r in rows:
    print('%6.3f %6.2f | %6.2f %6.2f | %6.2f %6.2f | %5.2f %5.2f | %5.2f' % (r[4], r[5], r[0], r[1], r[2], r[3], r[6], r[7], r[8]))
print(len(rows), 'soft scale-invariant singlets;', sum(1 for r in rows if r[0] > 4.5 and r[1] > 4.5), 'with both ratios > 4.5')
