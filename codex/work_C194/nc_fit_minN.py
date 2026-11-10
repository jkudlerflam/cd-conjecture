"""C194 (unregistered variant): drop sizes with N <= 2 (few-body, pair size ~ 21 sites comparable to L/N) and
redo the registered L-fit and the KT n_c fit.  Reads gap_table.json, writes nothing."""
import json, numpy as np
from fitgap import fit_gap
import nc_fit
g = json.load(open('gap_table.json'))
pts = []
for nom, v in sorted(g.items(), key=lambda kv: kv[1]['n']):
    rows = [r for r in v['rows'] if r['N'] > 2]
    if len(rows) >= 3:
        f = fit_gap([r['L'] for r in rows], [r['D'] for r in rows])
        print(f"n={nom:5s} sizes {[r['L'] for r in rows]}: D_inf = {f['D']:+.5f} +- {f['err']:.5f} ({f['D']/f['err']:+.1f} err)")
        pts.append((v['n'], f['D'], f['err']))
nc_fit.run('KT fit on the N > 2 variant', *zip(*pts))
