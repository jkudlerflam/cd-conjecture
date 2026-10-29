#!/usr/bin/env python3
"""From each junction patch take the point with the largest min(c_S(12), c_S(16)); keep it if that
value exceeds 1.1, and write the list to patch_screen_points.json for obcscreen.py."""
import json, glob, os
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
G = defaultdict(list)
for f in glob.glob(os.path.join(HERE, 'data', 'scan', 'patch_points_*.jsonl')):
    for l in open(f):
        r = json.loads(l)
        G[r['tag'].rsplit('_p', 1)[0]].append(r)
out = []
for k, rows in sorted(G.items()):
    r = max(rows, key=lambda r: min(r['L12']['cS'], r['L16']['cS']))
    m = min(r['L12']['cS'], r['L16']['cS'])
    if m > 1.1:
        out.append(dict(name=r['tag'], v=r['v'], c12=r['L12']['cS'], c16=r['L16']['cS']))
json.dump(out, open(os.path.join(HERE, 'patch_screen_points.json'), 'w'))
print(len(G), 'patches,', len(out), 'selected')
