#!/usr/bin/env python3
"""Finite-entanglement scaling (damped VUMPS of vumps.py, copied unchanged from contrib/d2_followup_220_1282)
for H = sum [XX - YY + u ZZ + b Z] (original frame, two-site cell) and controls xx, isxx.
Adapted from contrib/d2_staggered/fes_point.py for codex/work_C153 (Claude Opus agent, 2026-11-07).
usage: fes.py stag --u U --b B --chis 8 12 ... --outdir DIR"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
ap = argparse.ArgumentParser()
ap.add_argument('name'); ap.add_argument('--u', type=float, default=0.0); ap.add_argument('--b', type=float, default=0.0)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64, 96, 128])
ap.add_argument('--maxit', type=int, default=600); ap.add_argument('--tol', type=float, default=1e-9)
ap.add_argument('--outdir', default='fes_data')
a = ap.parse_args()
d = 2
if a.name == 'stag':
    h = V.bond(np.diag([1., -1., a.u]), [0, 0, a.b]); tag = f'stag_u{a.u:.6f}_b{a.b:.4f}'
elif a.name == 'xx':
    h = V.bond(np.diag([1., 1., 0.]), [0, 0, 0]); tag = 'xx'
elif a.name == 'isxx':
    hi = V.bond(np.diag([-1., 0, 0]), [0, 0, -1.]); hx = V.bond(np.diag([1., 1., 0.]), [0, 0, 0])
    one = np.eye(4).reshape(2, 2, 2, 2)
    T = (np.einsum('abcd,efgh->aebfcgdh', hi.reshape(2, 2, 2, 2), one)
         + np.einsum('abcd,efgh->aebfcgdh', one, hx.reshape(2, 2, 2, 2)))
    h = T.reshape(16, 16); d = 4; tag = 'isxx'
W = V.mpo_from_term(h, 2, d=d)
os.makedirs(a.outdir, exist_ok=True)
out = os.path.join(a.outdir, f'{tag}.json')
data = dict(name=a.name, u=a.u, b=a.b, results=[]); t0 = time.time()
def on_chi(psi, rec):
    rec['time_s'] = time.time() - t0; data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w'), default=float); os.replace(out + '.tmp', out)
V.run([W] * 2, 2, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=0, mode='damped', on_chi=on_chi)
print('done', out, time.time() - t0)
