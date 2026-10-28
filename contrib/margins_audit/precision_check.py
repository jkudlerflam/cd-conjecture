"""Precision check: recompute the main schedule point (and C33's exact point at
n = 23) with 60 more decimal digits and compare roots, centres and margins with
results/n{n}.json.  Usage: python3 -B precision_check.py n [n ...]
Appends to precision_check.json."""
import sys, os, json, time
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hpmat                                    # noqa: E402
from leading import leading, default_params, beta   # noqa: E402
from bands import analyse                       # noqa: E402
from run_schedule import prec_for, summarize    # noqa: E402

EXTRA = 60
out = []
fn = os.path.join(HERE, 'precision_check.json')
if os.path.exists(fn):
    out = json.load(open(fn))
for n in [int(a) for a in sys.argv[1:]]:
    R = json.load(open(os.path.join(HERE, 'results', 'n%d.json' % n)))
    labels = ['C33_kappa1e-3'] + (['C33_exact_n23'] if n == 23 else [])
    for lab in labels:
        t0 = time.time()
        r = R['runs'][lab]
        hpmat.setprec(prec_for(n) + EXTRA)
        P = default_params(n)
        b = beta(n, P['a'], P['Dx'])
        P['td'] = b * Decimal(2) ** (-r['p']) if r['p'] is not None else Decimal(0)
        eps = 2 * b * b * Decimal(2) ** (-r['f'])
        L = leading(n, P)
        o, wh = analyse(L['T'], L['mu'], eps)
        s = summarize(o, b)
        dr = max(abs(Decimal(a) - Decimal(c)) for a, c in zip(s['roots_k'], r['roots_k']))
        dc = max(abs(Decimal(a['centre_k']) - Decimal(c['centre_k'])) for a, c in zip(s['bands'], r['bands']))
        rec = dict(n=n, run=lab, prec=prec_for(n), prec_check=prec_for(n) + EXTRA,
                   max_root_diff=float(dr), max_centre_diff=float(dc),
                   delta4=[r['delta4'], s['delta4']], dPT=[r['dPT'], s['dPT']],
                   log10_lmin_F=[r['log10_lambda_min_F'], s['log10_lambda_min_F']],
                   seconds=time.time() - t0)
        print(json.dumps(rec), flush=True)
        out = [x for x in out if not (x['n'] == n and x['run'] == lab)] + [rec]
        json.dump(out, open(fn, 'w'), indent=1)
