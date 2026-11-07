"""Lesson 2 on the free TFIM: windows at g = 1 + c/n for large n (Majorana method)."""
import sys, json, warnings, time
warnings.filterwarnings('ignore')
from common import *
from free import *
from e0 import e0_tfim

cs = [float(a) for a in sys.argv[4].split(",")] if len(sys.argv) > 4 else [-40, -20, -10, -5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20, 40]
ns_opt = [int(a) for a in sys.argv[1].split(',')]
ns_fix = [int(a) for a in sys.argv[2].split(',')]
out = []
for c in cs:
    for n in sorted(set(ns_opt + ns_fix)):
        g = 1 + c / n; h = tfim(g); e0 = e0_tfim(g); M = np.sqrt(1 + g * g)
        rec = dict(c=c, n=n, g=g, e0=e0, M=M)
        ssd = ssd_window(n); sn = sine_window(n)
        rec['ssd_0'] = window_energy(h, ssd)
        rec['ssd_site'] = window_energy(h, ssd, onsite_weights=ssd_site_weights(n))
        rec['sine_0'] = window_energy(h, sn)
        rec['sine_v2'] = window_energy(h, sn, onsite_weights=sine_site_weights(n))
        if n in ns_opt:
            t = time.time()
            rec['opt'] = optimize_free_fast(h, n)['value']
            rec['ssd_u'] = optimize_free_fast(h, n, cfix=ssd)['value']
            rec['sine_u'] = optimize_free_fast(h, n, cfix=sn)['value']
            rec['sec'] = time.time() - t
        for k in list(rec):
            if k in ('ssd_0', 'ssd_site', 'sine_0', 'sine_v2', 'opt', 'ssd_u', 'sine_u'):
                rec['K_' + k] = n * n * (e0 - rec[k]) / M
        out.append(rec)
        print(c, n, {k: round(v, 4) for k, v in rec.items() if k.startswith('K_')}, flush=True)
    write_json(HERE / 'results' / f'lesson2_free_{sys.argv[3]}.json', out)
