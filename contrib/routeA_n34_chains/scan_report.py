#!/usr/bin/env python3
"""Weak-coupling report (Hartree-Fock Fermi points, velocities, centres, resonances, dimensionless
forward and non-forward couplings) for a list of chain specs (models34.py, without the level), and
for the velocity-equalized proxy of each.  usage: scan_report.py OUT.json SPEC6 [SPEC6 ...]
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import models34 as Mo
from chain_report import report, jsonable

out = []
for spec6 in sys.argv[2:]:
    K, info = Mo.load(spec6 + ':full')
    print('#####', spec6, 'mu', info['mu'], flush=True)
    r = report(K, info['n'], mu=0.0, top=8)
    r.pop('T_hf'); r.pop('M_hf')
    rec = dict(spec=spec6, info={k: v for k, v in info.items()}, true=r)
    try:
        Kp, pinfo = Mo.load(spec6 + ':pfull')
        rp = report(Kp, info['n'], mu=0.0, top=8, verbose=False)
        rp.pop('T_hf'); rp.pop('M_hf')
        rec['proxy'] = rp
        rec['proxy_info'] = dict(phi=pinfo['proxy_phi'], scales=pinfo['proxy_scales'])
        print('   proxy: #FP %d forward min eig %.4f max offdiag %.2e worst gamma %.2e flow %.2e' % (
            rp['n_fermi_points'], rp.get('forward_matrix_min_eig', np.nan), rp.get('forward_max_offdiag', np.nan),
            (rp.get('worst_channel') or {}).get('gamma', np.nan), (rp.get('worst_channel') or {}).get('flow', np.nan)), flush=True)
    except Exception as e:
        print('   proxy failed', e)
    out.append(jsonable(rec))
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
