#!/usr/bin/env python3
"""Finite-entanglement scaling run: VUMPS with an n-site cell over a list of bond dimensions.

usage: run_fes.py MODEL --n N --chis 8 16 32 ... [--mode damped|vumps] [--maxit 400] [--tol 1e-8]
                  [--seed 0] [--out data/fes_MODEL_nN.json]
Writes after every chi: energy per site, gradient error, entropies on every bond of the cell,
the leading eigenvalues of the cell transfer matrix (xi = -n / ln|l1/l0|), and the run time.
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
from models import model

ap = argparse.ArgumentParser()
ap.add_argument('model')
ap.add_argument('--n', type=int, default=1)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 16, 32])
ap.add_argument('--mode', default='damped')
ap.add_argument('--maxit', type=int, default=400)
ap.add_argument('--tol', type=float, default=1e-8)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--krylov', type=int, default=16)
ap.add_argument('--stol', type=float, default=1e-7, help='stop when S changes by less than this over 20 steps')
ap.add_argument('--out', default=None)
ap.add_argument('--save_state', action='store_true')
ap.add_argument('--load', default=None, help='npz state to start from (its chi must equal the first entry of --chis)')
a = ap.parse_args()
h, r, c, note = model(a.model)
W = V.mpo_from_term(h, r)
out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', f'fes_{a.model}_n{a.n}_{a.mode}.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
data = dict(model=a.model, note=note, c_exact=c, n=a.n, mode=a.mode, seed=a.seed, results=[])
t0 = time.time()


def save(psi, rec):
    rec['time_s'] = time.time() - t0
    rec['mag'] = V.local_expect(psi, list(V.PAULI)).tolist()
    data['results'].append(rec)
    if a.save_state:
        np.savez_compressed(out.replace('.json', f"_chi{rec['chi']}.npz"),
                            **{f'AL{k}': psi.AL[k] for k in range(psi.n)}, **{f'AR{k}': psi.AR[k] for k in range(psi.n)},
                            **{f'C{k}': psi.C[k] for k in range(psi.n)})
    json.dump(data, open(out + '.tmp', 'w'))
    os.replace(out + '.tmp', out)
    print(f"{a.model} n={a.n} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
          f"Smean={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)


psi0 = None
if a.load:
    z = np.load(a.load)
    psi0 = V.UMPS([z[f'AL{k}'] for k in range(a.n)], [z[f'AR{k}'] for k in range(a.n)], [z[f'C{k}'] for k in range(a.n)])
    if os.path.exists(out):
        data = json.load(open(out))
V.run([W] * a.n, a.n, a.chis, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode, krylov=a.krylov,
      log=lambda rec: None, on_chi=save, psi=psi0, stol=a.stol)
