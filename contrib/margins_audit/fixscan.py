"""Parameter fixes at one n: scan the field power f (cheap, T and mu fixed)
and the loop strength p (one coefficient per value), keeping everything else
at the main schedule.  Usage: python3 -B fixscan.py n [p offsets ...]
Writes fixscan_n{n}.json."""
import sys, os, json, time
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hpmat                                    # noqa: E402
from leading import leading, default_params, beta   # noqa: E402
from bands import analyse                       # noqa: E402
from run_schedule import prec_for               # noqa: E402

n = int(sys.argv[1])
poffs = [int(x) for x in sys.argv[2:]] or [-3, -1, 1]
R = json.load(open(os.path.join(HERE, 'results', 'n%d.json' % n)))
r = R['runs']['C33_kappa1e-3']
p0, f0 = r['p'], r['f']
hpmat.setprec(prec_for(n))
P = default_params(n)
b = beta(n, P['a'], P['Dx'])
out = []


def rec(o, p, f):
    x = dict(n=n, p=p, f=f, delta2=float(o['delta2']), delta4=float(o['delta4']), witness4=o['w4'],
             dPT=float(o['dPT']), dU1=float(o['dU1']), dU2=float(o['dU2']), dF=float(o['dF']),
             unimodal=o['unimodal'], loop_whitened=o['loop_whitened'], field_whitened=o['field_whitened'])
    print(json.dumps(x), flush=True)
    out.append(x)
    json.dump(out, open(os.path.join(HERE, 'fixscan_n%d.json' % n), 'w'), indent=1)


for dp in [0] + poffs:
    t0 = time.time()
    p = p0 + dp
    P['td'] = b * Decimal(2) ** (-p)
    L = leading(n, P)
    fs = [f0 - 3, f0 - 2, f0 - 1, f0, f0 + 1, f0 + 3] if dp == 0 else [f0]
    for f in fs:
        o, wh = analyse(L['T'], L['mu'], 2 * b * b * Decimal(2) ** (-f))
        rec(o, p, f)
