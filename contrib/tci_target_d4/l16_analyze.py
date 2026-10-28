#!/usr/bin/env python3
"""Read the L = 16 output of ed_lean.py (data/l16_k0.jsonl with --entropy, l16_k8, l16_k7) and the
L = 8, 12 chart data, and print for each point (units 2 pi / L): phi = E(k=8) - E0, F = E1(k=0) - E0,
desc = E(k=7) - E0, the ratios F/phi and desc/phi, x_eff = 1/(desc/phi - 1), and c_half(12 -> 16).
Assumes, as at L = 8 and 12 for these points, that the ground state is at k = 0 and the lowest
levels at k = 8 and 7 and the first excited level at k = 0 are singlets (check with --sz 2 if needed).
Tricritical control at L = 16: F/phi = 6.49, desc/phi = 6.34 (x_eff = 0.19).
usage: l16_analyze.py [dir]"""
import json, sys, glob, os
import numpy as np
from analyze import load, levels
d = sys.argv[1] if len(sys.argv) > 1 else 'data'
recs = {}
for k in (0, 8, 7):
    f = os.path.join(d, 'l16_k%d.jsonl' % k)
    if os.path.exists(f):
        for line in open(f):
            r = json.loads(line)
            recs.setdefault((round(r['theta'], 6), round(r['eta'], 6)), {})[k] = r
small = load(glob.glob('data/grid/*.jsonl') + glob.glob('data/fine/*.jsonl') + glob.glob('data/left/*.jsonl'))
f = 16 / (2 * np.pi)
for key, b in recs.items():
    if 0 not in b:
        continue
    e0 = b[0]['E'][0]
    phi = (b[8]['E'][0] - e0) * f if 8 in b else np.nan
    F = (b[0]['E'][1] - e0) * f
    de = (b[7]['E'][0] - e0) * f if 7 in b else np.nan
    line = '(%.3f pi, %.2f): phi %.4f F/phi %.2f desc/phi %.2f x_eff %.2f' % (key[0] / np.pi, key[1], phi, F / phi, de / phi, 1 / (de / phi - 1))
    s = small.get(key, {})
    if 12 in s and 'S_half' in b[0]:
        line += '  c_half(12->16) %.3f' % (3 * (b[0]['S_half'] - s[12]['S'][-1]) / np.log(16 / 12))
    print(line)
