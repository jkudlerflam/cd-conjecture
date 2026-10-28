#!/usr/bin/env python3
"""Print the summary rows used in README.md from data/report_*.json.
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import json, sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for fn in ('report_spatial.json', 'report_compile.json'):
    p = os.path.join(HERE, 'data', fn)
    if not os.path.exists(p):
        continue
    for rec in json.load(open(p)):
        t = rec['true']; info = rec['info']
        fps = t['fermi_points']
        vs = [abs(f['v']) for f in fps] or [np.nan]
        cr = t.get('closest_resonances') or []
        top = cr[0] if cr else {}
        pt = [c for c in cr if c['kind'] == 'pair transfer']
        w = t.get('worst_channel') or {}
        print('%-38s FP=%d  v %.1e..%.1e  fwd %.1e  K-1 %s  closest %s %.1e (g %.1e)  pairT %s  worst %s %.1e g %.1e flow %.2f  rep %s' % (
            rec['spec'], t['n_fermi_points'], min(vs), max(vs), t.get('forward_max_offdiag', np.nan),
            None if not t.get('K_bands_estimate') else '%.1e' % max(abs(k - 1) for k in t['K_bands_estimate'] if k is not None),
            top.get('kind', '-')[:14], top.get('mismatch', np.nan), top.get('gamma', np.nan),
            ('%.1e' % pt[0]['mismatch']) if pt else '-',
            w.get('kind', '-')[:14], w.get('mismatch', np.nan), w.get('gamma', np.nan), w.get('flow', np.nan),
            None if not t.get('closest_repeated_point_channel') else '%.1e' % t['closest_repeated_point_channel']['mismatch']))
        print('      bands', [(b['band'], None if b.get('centre') is None else round(b['centre'], 4), round(b.get('filling', -1), 4)) for b in t['bands']],
              ' k', [round(f['k'], 4) for f in fps], ' body', ['%.0e' % x for x in t['body_norms'][1:4]])
        if 'proxy' in rec:
            q = rec['proxy']; wq = q.get('worst_channel') or {}
            print('      proxy: fwd %.1e worst g %.1e flow %.2f  v %.2f..%.2f' % (q.get('forward_max_offdiag', np.nan), wq.get('gamma', np.nan), wq.get('flow', np.nan),
                  min(abs(f['v']) for f in q['fermi_points']) if q['fermi_points'] else np.nan, max(abs(f['v']) for f in q['fermi_points']) if q['fermi_points'] else np.nan))
