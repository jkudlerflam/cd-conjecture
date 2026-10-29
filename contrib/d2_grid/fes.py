#!/usr/bin/env python3
"""Finite-entanglement scaling with damped VUMPS (vumps.py of contrib/d2_followup_220_1282, unchanged).

usage: fes.py NAME --cellscan 1 2 3 ... [--chi 16]      (energy per site for each cell size n)
       fes.py NAME --n N --chis 8 12 16 24 32 48 64 96 128
Output data/fes/NAME_cells.json or data/fes/NAME_nN.json (S on every bond of the cell, transfer-matrix
spectrum, xi = -n/ln|l1/l0|, energy, gradient error).  c_fit = 6 dS/dln xi (analyze.py).
Written for contrib/d2_grid (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "d2_followup_220_1282"))
import numpy as np
from models_grid import model
import vumps as V

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--cellscan', type=int, nargs='*', default=None)
ap.add_argument('--chi', type=int, default=16)
ap.add_argument('--n', type=int, default=1)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64, 96, 128])
ap.add_argument('--mode', default='damped')
ap.add_argument('--maxit', type=int, default=400)
ap.add_argument('--tol', type=float, default=1e-8)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--tag', default='', help='suffix for the output file name')
a = ap.parse_args()
h, r, c, note, _ = model(a.name)
W = V.mpo_from_term(h, r)
os.makedirs(os.path.join(HERE, 'data', 'fes'), exist_ok=True)
t0 = time.time()
if a.cellscan is not None:
    out = os.path.join(HERE, 'data', 'fes', f'{a.name}_cells.json')
    res = []
    for n in a.cellscan:
        try:
            psi, rs = V.run([W] * n, n, [8, a.chi], tol=1e-7, maxit=150, seed=a.seed, mode=a.mode, log=lambda r: None)
            rr = rs[-1]
            res.append(dict(n=n, e=rr['e'], err=rr['err'], S=float(np.mean(rr['S'])), xi=rr['xi'],
                            lam1=rr['lam'][1]))
        except Exception as ex:
            res.append(dict(n=n, error=str(ex)))
        print(a.name, res[-1], f't={time.time()-t0:.0f}', flush=True)
        json.dump(dict(name=a.name, note=note, chi=a.chi, res=res), open(out, 'w'))
    sys.exit()
out = os.path.join(HERE, 'data', 'fes', f'{a.name}_n{a.n}{a.tag}.json')
data = dict(name=a.name, note=note, c_exact=c, n=a.n, mode=a.mode, results=[])


def save(psi, rec):
    rec['time_s'] = time.time() - t0
    rec['S'] = [float(x) for x in rec['S']]
    rec['mag'] = V.local_expect(psi, list(V.PAULI)).tolist()
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
    print(f"{a.name} n={a.n} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
          f"S={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)


V.run([W] * a.n, a.n, a.chis, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode, log=lambda r: None, on_chi=save)
