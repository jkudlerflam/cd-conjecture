#!/usr/bin/env python3
"""Finite-entanglement scaling of a logical chain with VUMPS (vumps.py, copied unchanged from
contrib/d2_followup_220_1282).  usage:
  run_fes34.py SPEC --nc 2 --chis 8 16 32 ... [--mode damped|vumps] [--out file.json]
SPEC as in models34.py.  The bond is a two-site term on local dimension 2^n (one logical cell
per site); nc is the VUMPS unit cell in logical cells.  Written by Claude (Opus 5.5) for
contrib/routeA_n34_chains, 2026-10-28."""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
import models34 as Mo

ap = argparse.ArgumentParser()
ap.add_argument('spec')
ap.add_argument('--nc', type=int, default=2)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 16, 32])
ap.add_argument('--mode', default='damped')
ap.add_argument('--maxit', type=int, default=400)
ap.add_argument('--tol', type=float, default=1e-8)
ap.add_argument('--stol', type=float, default=1e-7)
ap.add_argument('--tau0', type=float, default=10.0)
ap.add_argument('--taumax', type=float, default=1e5)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--mpotol', type=float, default=1e-13)
ap.add_argument('--out', default=None)
a = ap.parse_args()
K, info = Mo.load(a.spec)
n = info['n']; d = 2 ** n
W = V.mpo_from_term(K, 2, d=d, tol=a.mpotol)
info['mpo_D'] = int(W.shape[0])
tag = a.spec.replace(':', '_').replace(',', '-')
out = a.out or os.path.join(Mo.HERE, 'data', 'fes', f'fes_{tag}_nc{a.nc}_{a.mode}.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
data = dict(info=info, nc=a.nc, mode=a.mode, results=[])
t0 = time.time()
print('spec', a.spec, 'MPO D', W.shape[0], flush=True)


def save(psi, rec):
    rec['time_s'] = time.time() - t0
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w'))
    os.replace(out + '.tmp', out)
    print(f"chi={rec['chi']} e={rec['e']:.14f} err={rec['err']:.1e} it={rec['iters']} "
          f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)


V.run([W] * a.nc, a.nc, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode,
      log=lambda rec: None, on_chi=save, stol=a.stol, tau0=a.tau0, tau_max=a.taumax)
