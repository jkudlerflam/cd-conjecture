#!/usr/bin/env python3
"""Merge ED result files from several directories (laptop and Typhon) into results/.
Sizes present in more than one source must agree in E0 to 1e-8 (checked)."""
import json, glob, os, sys
srcs = sys.argv[1:] or ['results_local', 'results_cluster']
merged = {}
for s in srcs:
    for f in glob.glob(os.path.join(s, '*.json')):
        n = os.path.basename(f)
        if n.startswith('dmrg_') or n in ('summary.json', 'junction.json'):
            continue
        d = json.load(open(f))
        if 'sizes' not in d:
            continue
        m = merged.setdefault(n, {'name': d['name'], 'meta': d['meta'], 'sizes': {}})
        for L, rec in d['sizes'].items():
            if L in m['sizes']:
                assert abs(m['sizes'][L]['E0'] - rec['E0']) < 1e-8 * max(1, abs(rec['E0'])), (n, L)
            else:
                m['sizes'][L] = rec
for n, m in merged.items():
    json.dump(m, open(os.path.join('results', n), 'w'))
print('merged', len(merged), 'files')
