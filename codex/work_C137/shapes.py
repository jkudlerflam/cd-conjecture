"""Describe optimal windows: relative sup distance of c_opt to the sine-square and C135 sine windows,
and the most negative weight."""
import json, glob
import numpy as np
from common import HERE, ssd_window, sine_window
from report import load_survey
rows = []
names = ['Heisenberg', 'XX', 'XXZ_0.5', 'XXZ_2', 'XXZ_-0.9', 'XXZ_-0.99', 'TFIM_1', 'TFIMc-2', 'TFIMc0.5', 'TFIMc5'] + [f'random_{i:02d}' for i in range(30)]
for name in names:
    for r in load_survey(name):
        n = r['n']
        if n < 4: continue
        c = np.array(r['opt']['c']); s1 = ssd_window(n); s2 = sine_window(n)
        rows.append(dict(name=name, n=n, d_ssd=float(np.max(np.abs(c - s1)) / s1.max()), d_sine=float(np.max(np.abs(c - s2)) / s2.max()), cmin=float(c.min())))
json.dump(rows, open(HERE / 'results' / 'shapes.json', 'w'), indent=1)
import collections
agg = collections.defaultdict(list)
for r in rows: agg[r['name']].append(r)
for k, v in agg.items():
    ev = [x for x in v if x['n'] % 2 == 0]; od = [x for x in v if x['n'] % 2 == 1]
    print(k, 'even: max d_ssd %.3g max d_sine %.3g | odd: max d_ssd %.3g max d_sine %.3g | min c %.3g' % (
        max([x['d_ssd'] for x in ev] or [np.nan]), max([x['d_sine'] for x in ev] or [np.nan]),
        max([x['d_ssd'] for x in od] or [np.nan]), max([x['d_sine'] for x in od] or [np.nan]), min(x['cmin'] for x in v)))
