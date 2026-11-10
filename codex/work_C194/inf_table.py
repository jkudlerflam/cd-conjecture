"""C194 infinite-chain table (registered readings, RESULT.md definition 3).
Usage: python3 -E -s inf_table.py <dir> [<dir> ...]   (runs with the same tag in several dirs: largest chi list wins)
Per chi: iterations and whether the run stopped on tolerance or on the cap, S_bar, xi_0, xi_1 (T1f), xi_2 (T2).
Per run: growth of xi_1, xi_2, xi_0 over the last chi doubling, c_fit over the largest four chi (C189 definition),
and the G1 test (xi_1 growth < 1.15 and xi_2 growth > 1.5)."""
import json, sys, os
import numpy as np
from summarize_inf import load, xis, cfit

np.seterr(all='ignore')
runs = {}
for d in sys.argv[1:]:
    for tag, recs in load(d, 'inf_').items():
        if tag not in runs or len(recs) > len(runs[tag]):
            runs[tag] = recs
out = {}
for tag in sorted(runs, key=lambda t: runs[t][0]['cfg']['p'] / runs[t][0]['cfg']['q']):
    recs = sorted(runs[tag], key=lambda r: r['chi'])
    c = recs[0]['cfg']; n = c['p'] / c['q']; cap = c.get('max_it', 400)
    print(f"=== {tag}  n = {c['p']}/{c['q']} = {n:.4f}  (cap {cap})")
    rows = []
    for r in recs:
        x = xis(r); Sb = 0.5 * (r['S'] + r['S_prev'])
        it = r['conv']['iters']; how = 'CAP' if it >= cap else 'tol'
        rows.append(dict(chi=r['chi'], S=Sb, xi0=x.get('xi0'), xi1=x.get('xi1f'), xi2=x.get('xi2'), xi=x['xi'],
                         it=it, how=how, dS=r['conv']['dS'], e=r['e_site'], trunc=r['trunc']))
        print(f"  chi={r['chi']:4d} it={it:4d}({how}) dS={r['conv']['dS']:.1e} e={r['e_site']:.10f} trunc={r['trunc']:.1e} "
              f"S={Sb:.4f} xi0={x.get('xi0', np.nan):8.1f} xi1={x.get('xi1f', np.nan):8.1f} xi2={x.get('xi2', np.nan):8.1f} "
              f"xi1/xi2={x.get('xi1f', np.nan) / x.get('xi2', np.nan):.3f}")
    res = dict(n=n, rows=rows)
    if len(rows) >= 4:
        last = rows[-4:]
        res['c_fit'] = float(cfit([r['xi'] for r in last], [r['S'] for r in last]))
        res['c_fit_xi2'] = float(cfit([r['xi2'] for r in last], [r['S'] for r in last]))
    byc = {r['chi']: r for r in rows}
    top = rows[-1]['chi']
    if top // 2 in byc:
        a, b = byc[top // 2], rows[-1]
        g = {k: b[k] / a[k] for k in ('xi0', 'xi1', 'xi2')}
        res['growth'] = g; res['doubling'] = (top // 2, top)
        res['G1'] = bool(g['xi1'] < 1.15 and g['xi2'] > 1.5)
        print(f"  last doubling {top // 2}->{top}: xi1 x{g['xi1']:.3f}, xi2 x{g['xi2']:.3f}, xi0 x{g['xi0']:.3f}; "
              f"c_fit(xi=max) = {res.get('c_fit', np.nan):.3f}, c_fit(xi2) = {res.get('c_fit_xi2', np.nan):.3f}; "
              f"G1 reading (xi1 < 1.15 and xi2 > 1.5): {res['G1']}")
    out[tag] = res
json.dump(out, open('inf_table.json', 'w'), indent=1, default=float)
