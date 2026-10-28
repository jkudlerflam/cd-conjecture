#!/usr/bin/env python3
"""Bulk finite-entanglement scaling on a long open chain (fdmrg.FDMRG), for a model spec of
models.py.  At bond dimension chi the chi-induced correlation length xi is much shorter than
the distance to the ends, so the central region is the infinite-chain finite-chi state and
needs no commensurate unit cell (two incommensurate modes are allowed).

usage: run_fdmrg.py SPEC --L 400 --chis 8 12 16 24 32 48 64 --win 40 80 --out data/x.json
Per chi: energy, S on every bond, the mean S over the central 40 bonds, the leading
eigenvalues of the transfer matrix of central windows of m sites (xi = -m / ln|l1/l0|), and,
for d = 2, the profiles of <n>, <sigma^x>, <sigma^y>.
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fdmrg import FDMRG
from vumps import mpo_from_term
from models import get

ap = argparse.ArgumentParser()
ap.add_argument('spec')
ap.add_argument('--L', type=int, default=400)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64])
ap.add_argument('--sweeps', type=int, default=4)
ap.add_argument('--maxsweeps', type=int, default=20)
ap.add_argument('--etol', type=float, default=1e-9)
ap.add_argument('--win', type=int, nargs='+', default=[40, 80])
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--out', required=True)
a = ap.parse_args()
h, d = get(a.spec)
W = mpo_from_term(h, 2, d=d)
dm = FDMRG(W, a.L, d=d, seed=a.seed)
data = dict(spec=a.spec, L=a.L, d=d, results=[])
t0 = time.time()
n_op = np.diag([0., 1.])
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]])
for chi in a.chis:
    Eprev = None
    for sw in range(a.maxsweeps):
        E, tr = dm.sweep(chi)
        if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < a.etol * a.L:
            break
        Eprev = E
    mid = a.L // 2
    rec = dict(chi=chi, E=E, e_site=E / a.L, sweeps=sw + 1, trunc=tr, S=dm.S.tolist(),
               S_mid=float(np.mean(dm.S[mid - 20:mid + 20])))
    if d == 2:
        loc = dm.local_ops([n_op, X, Y])
        rec['n'] = loc[0].tolist(); rec['sx'] = loc[1].tolist(); rec['sy'] = loc[2].tolist()
        rec['n_mid'] = float(np.mean(loc[0][mid - 20:mid + 20]))
    spec = {}
    for m in a.win:
        j0 = mid - m // 2
        try:
            w = dm.window_spectrum(j0, m)
            spec[str(m)] = [[float(z.real), float(z.imag)] for z in w]
            rec[f'xi_{m}'] = float(-m / np.log(abs(w[1]) / abs(w[0])))
        except Exception as ex:
            spec[str(m)] = str(ex)
    rec['win'] = spec
    rec['time_s'] = time.time() - t0
    data['results'].append(rec)
    json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
    print(f"{a.spec} chi={chi} e={E/a.L:.10f} sweeps={sw+1} trunc={tr:.1e} S_mid={rec['S_mid']:.5f} "
          + ' '.join(f"xi{m}={rec.get(f'xi_{m}', float('nan')):.2f}" for m in a.win)
          + (f" n_mid={rec['n_mid']:.4f}" if d == 2 else '') + f" t={time.time()-t0:.0f}s", flush=True)
