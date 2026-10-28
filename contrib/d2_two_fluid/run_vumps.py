#!/usr/bin/env python3
"""Finite-entanglement scaling with VUMPS (vumps.py, damped mode) in an n-site cell, for a
model spec of models.py.  Adapted from ../d2_followup_220_1282/run_fes.py.

usage: run_vumps.py SPEC --n N --chis 8 12 16 24 32 48 64 --out data/vu_NAME.json
Per chi: energy per site, gradient error, entropies on every bond of the cell, the leading
transfer-matrix eigenvalues (xi = -n / ln|l1/l0|), and the cell expectation values of the
Pauli matrices.
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
from models import get

ap = argparse.ArgumentParser()
ap.add_argument('spec')
ap.add_argument('--n', type=int, default=4)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64])
ap.add_argument('--maxit', type=int, default=1500)
ap.add_argument('--tol', type=float, default=1e-8)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--out', required=True)
a = ap.parse_args()
h, d = get(a.spec)
W = V.mpo_from_term(h, 2, d=d)
data = dict(spec=a.spec, n=a.n, d=d, results=[])
t0 = time.time()


def save(psi, rec):
    rec['time_s'] = time.time() - t0
    if d == 2:
        rec['mag'] = V.local_expect(psi, list(V.PAULI)).tolist()
    data['results'].append(rec)
    json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
    print(f"{a.spec} n={a.n} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
          f"Smean={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)


V.run([W] * a.n, a.n, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode='damped', krylov=16,
      log=lambda rec: None, on_chi=save, stol=1e-7)
