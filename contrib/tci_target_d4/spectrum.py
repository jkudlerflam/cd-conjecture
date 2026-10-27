#!/usr/bin/env python3
"""Print the labelled low spectrum (units 2 pi / L, momentum offset from k0, spin) of one point.
usage: spectrum.py theta_over_pi eta [emax] [files...]"""
import sys, glob
import numpy as np
from analyze import load, levels
t, e = float(sys.argv[1]), float(sys.argv[2])
emax = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
files = sys.argv[4:] or glob.glob('data/grid/*.jsonl') + glob.glob('data/grid_local/*.jsonl') + glob.glob('data/fine/*.jsonl')
recs = load(files)
key = min(recs, key=lambda k: abs(k[0] / np.pi - t) + abs(k[1] - e))
print('point', key[0] / np.pi, key[1])
for L in sorted(recs[key]):
    r = recs[key][L]
    print('L=%d k0=%d E0/L=%.6f S_half=%.4f' % (L, r['k0'], r['E0'] / L, r['S'][-1]))
    lv = [x for x in levels(r) if x[0] < emax]
    # group by momentum offset
    for dk in sorted(set(x[1] for x in lv)):
        print('   dk=%2d: ' % dk + '  '.join('%.3f(S%s)' % (x[0], ('%g' % x[2]) if x[2] >= 0 else '?') for x in lv if x[1] == dk))
