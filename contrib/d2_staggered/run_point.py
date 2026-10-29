#!/usr/bin/env python3
"""Charge-resolved DMRG data for one (model, u, b, L, boundary condition).

usage: run_point.py stag --u 1.35 --b 2 --L 64 [--pbc] --chis 64 128 192 [--nexc 3] [--qs 0 1 2]
       run_point.py tfi --L 64 ...        run_point.py isxx --L 64 ...
For each charge offset q in --qs (sector N = N0 + q, N0 = L/2 for stag and isxx, 0 for tfi) the
ground state is computed; in the q = 0 sector also --nexc excited states (and --nexc1 in q = 1).
Output data/<bc>/<tag>.json: energies, chi histories, entanglement S and charge variance F on every
bond of the q = 0 ground state, and the charge-resolved Schmidt spectrum at the central bond.
Written for contrib/d2_staggered (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from udmrg import model, run_state

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--u', type=float, default=0.0)
ap.add_argument('--b', type=float, default=0.0)
ap.add_argument('--L', type=int, required=True)
ap.add_argument('--pbc', action='store_true')
ap.add_argument('--chis', type=int, nargs='+', default=[32, 64, 128])
ap.add_argument('--nexc', type=int, default=3)
ap.add_argument('--nexc1', type=int, default=1)
ap.add_argument('--qs', type=int, nargs='+', default=[0, 1, 2])
ap.add_argument('--etol', type=float, default=1e-10)
ap.add_argument('--out', default=None)
a = ap.parse_args()

bc = 'pbc' if a.pbc else 'obc'
tag = a.out or (f'{a.name}_u{a.u:.4f}_b{a.b:.4f}_L{a.L}' if a.name == 'stag' else f'{a.name}_L{a.L}')
out = os.path.join(HERE, 'data', bc, tag + '.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
Ws, d, n = model(a.name, a.L, a.u, a.b, a.pbc)
N0 = 0 if a.name == 'tfi' else a.L // 2
qs = [0] if a.name == 'tfi' else a.qs
res = dict(name=a.name, u=a.u, b=a.b, L=a.L, bc=bc, chis=a.chis, N0=N0, sectors={})
t0 = time.time()


def save():
    res['time_s'] = time.time() - t0
    json.dump(res, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)


log = lambda s: print(s, flush=True)
for q in qs:
    N = N0 + q
    states, rec = [], []
    nst = 1 + (a.nexc if q == 0 else (a.nexc1 if q == 1 else 0))
    for k in range(nst):
        dm, hist = run_state(Ws, d, n, N, a.chis, ortho=states, etol=a.etol, seed=1 + k, log=log)
        E = dm.energy() / dm.norm2()
        r = dict(E=E, hist=hist, overlaps=[dm.overlap(s) for s in states])
        if q == 0 and k == 0:
            S, F, spec = dm.measure()
            r.update(S=S, F=F, spec_mid=spec[a.L // 2], spec_quarter=spec[a.L // 4])
        else:
            dm.measure()
        states.append(dm)
        rec.append(r)
        log(f'q={q} k={k} E={E:.12f} t={time.time()-t0:.0f}s')
        res['sectors'][str(q)] = rec
        save()
print('done', out)
