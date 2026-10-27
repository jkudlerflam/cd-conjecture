#!/usr/bin/env python3
"""Size dependence of the TCI-relevant levels (units 2 pi / L) at L = 8 and 12 for chosen points.
usage: detail.py 'th_lo th_hi eta_lo eta_hi' files..."""
import sys, glob
import numpy as np
from analyze import load, summary
lo, hi, elo, ehi = map(float, sys.argv[1].split())
files = sys.argv[2:] or glob.glob('data/grid/*.jsonl') + glob.glob('data/grid_local/*.jsonl')
recs = load(files)
print('%6s %6s | %5s | %5s %5s %5s | %5s %5s %5s | %5s %5s | %5s %5s | z' % ('th', 'eta', 'chalf', 'cur8', 'pi08', 'k008', 'cur12', 'pi012', 'k0012', 'r8', 'r12', 'k/p8', 'k/p12'))
for key in sorted(recs):
    t, e = key[0] / np.pi, key[1]
    if not (lo <= t <= hi and elo <= e <= ehi):
        continue
    b = recs[key]
    if 8 not in b or 12 not in b:
        continue
    s8, s12 = summary(b[8]), summary(b[12])
    ch = 3 * (s12['S_half'] - s8['S_half']) / np.log(1.5)
    z = -np.log((s12['gap'] / 12) / (s8['gap'] / 8)) / np.log(1.5)
    print('%6.3f %6.2f | %5.2f | %5.3f %5.3f %5.3f | %5.3f %5.3f %5.3f | %5.2f %5.2f | %5.2f %5.2f | %5.2f' % (
        t, e, ch, s8['cur'], s8['pi0'], s8['k00'], s12['cur'], s12['pi0'], s12['k00'],
        s8['pi0'] / s8['cur'], s12['pi0'] / s12['cur'], s8['k00'] / s8['pi0'], s12['k00'] / s12['pi0'], z))
