"""Optimal, SSD and sine windows by exact diagonalization. Usage: python survey.py NAME nmin nmax"""
import sys, time, warnings, json
warnings.filterwarnings('ignore')
from common import *
import models

name = sys.argv[1]; nmin, nmax = int(sys.argv[2]), int(sys.argv[3])
import os
out_path = HERE / 'results' / (f'{name}_n{nmin}.json' if os.environ.get('C137_SPLIT') else f'{name}.json')
recs = json.load(open(out_path)) if out_path.exists() else []
done = {r['n'] for r in recs}
for n in range(nmin, nmax + 1):
    if n in done:
        continue
    m = models.get(name, n)
    h = m['h']; t0 = time.time()
    rec = dict(name=name, n=n, e0_exact=m['e0'], g=m.get('g'))
    fam = WindowFamily(h, n, udirs=m['udirs'], reflect=m['reflect'])
    r = optimize_window(fam, maxit=300)
    rec['opt'] = dict(value=r['value'], c=r['c'].tolist(), y={str(a): v.tolist() for a, v in r['y'].items()},
                      iters=r['iters'], pred_gap=r['pred_gap'], evals=r['evals'])
    for wname, c in [('ssd', ssd_window(n)), ('sine', sine_window(n))]:
        famc = WindowFamily(h, n, udirs=m['udirs'], reflect=m['reflect'], cfix=c)
        rc = optimize_window(famc, maxit=300)
        rec[wname + '_u'] = dict(value=rc['value'], pred_gap=rc.get('pred_gap', 0.0))
        ev, _ = lowest(WindowFamily(h, n, udirs=(), cfix=c).matrix(np.zeros(0)), 2)
        rec[wname + '_0'] = dict(value=float(ev[0]))
    rec['seconds'] = time.time() - t0
    recs.append(rec)
    recs.sort(key=lambda r: r['n'])
    write_json(out_path, recs)
    print(name, n, rec['opt']['value'], rec['ssd_u']['value'], rec['sine_u']['value'], round(rec['seconds'], 1), flush=True)
