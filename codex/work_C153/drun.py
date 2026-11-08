#!/usr/bin/env python3
"""One open-chain point of the staggered chain with udmrg: q=0 ground state (S, F profiles) and nexc excited
states, q=1 and q=2 ground states. Output JSON to --outdir. Written for codex/work_C153 (Claude Opus agent, 2026-11-07).
usage: drun.py --u U --b B --L L --chis 64 128 256 --nexc 2 --outdir DIR"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from udmrg import model, run_state
ap = argparse.ArgumentParser()
ap.add_argument('--u', type=float, required=True); ap.add_argument('--b', type=float, required=True)
ap.add_argument('--L', type=int, required=True); ap.add_argument('--chis', type=int, nargs='+', default=[64, 128])
ap.add_argument('--nexc', type=int, default=2); ap.add_argument('--qs', type=int, nargs='+', default=[0, 1, 2])
ap.add_argument('--etol', type=float, default=1e-10); ap.add_argument('--outdir', default='data')
ap.add_argument('--tag', default='')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
out = os.path.join(a.outdir, f'u{a.u:.6f}_b{a.b:.4f}_L{a.L}{a.tag}.json')
Ws, d, n = model('stag', a.L, a.u, a.b, False)
N0 = a.L // 2
res = dict(u=a.u, b=a.b, L=a.L, chis=a.chis, sectors={}); t0 = time.time()
for q in a.qs:
    states, rec = [], []
    for k in range(1 + (a.nexc if q == 0 else 0)):
        dm, hist = run_state(Ws, d, n, N0 + q, a.chis, ortho=states, etol=a.etol, seed=1 + k)
        E = dm.energy() / dm.norm2(); r = dict(E=float(E), hist=hist)
        if q == 0 and k == 0:
            S, F, spec = dm.measure(); r.update(S=S, F=F)
        elif k < a.nexc:
            dm.measure()
        states.append(dm); rec.append(r)
        res['sectors'][str(q)] = rec; res['time_s'] = time.time() - t0
        json.dump(res, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
print('done', out, f'{time.time()-t0:.0f}s')
