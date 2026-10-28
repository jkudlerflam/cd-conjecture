#!/usr/bin/env python3
"""Finite-entanglement scaling with two-site iDMRG (idmrg.py) for a model spec of models.py.

usage: run_idmrg.py SPEC --chis 16 24 32 48 64 96 128 [--steps 400] --out data/x.jsonl [--label L]
Each line of the output: chi, energy per site, S on both cell bonds, xi from the two-site
transfer matrix, and (d = 2) the densities <n> on both sites, <sigma^x>, <sigma^y>, and the
adjacent-pair density <n n>.
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import idmrg
from models import get

_orig = idmrg.measure


def measure(A, B, lam, lam_prev, h, d):
    res = _orig(A, B, lam, lam_prev, h, d)
    if d == 2:
        inv = np.where(lam_prev > 1e-10, 1 / np.maximum(lam_prev, 1e-300), 0.0)
        Ga = inv[:, None, None] * A
        Gb = B * inv[None, None, :]
        th = np.tensordot(lam_prev[:, None, None] * Ga * lam[None, None, :], Gb * lam_prev[None, None, :], axes=([2], [0]))
        nrm = np.vdot(th, th).real
        n = np.diag([0., 1.]); sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]])
        def one(O, site):
            t = np.tensordot(O, th, axes=([1], [1 + site]))
            t = np.moveaxis(t, 0, 1 + site)
            return float(np.vdot(th, t).real / nrm)
        nn = np.kron(n, n).reshape(2, 2, 2, 2)
        t = np.tensordot(nn, th, axes=([2, 3], [1, 2])).transpose(2, 0, 1, 3)
        res.update(n=[one(n, 0), one(n, 1)], sx=[one(sx, 0), one(sx, 1)], sy=[one(sy, 0), one(sy, 1)],
                   nn=float(np.vdot(th, t).real / nrm))
    return res


idmrg.measure = measure

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('spec')
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64, 96, 128])
    ap.add_argument('--steps', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--label', default='')
    a = ap.parse_args()
    h, d = get(a.spec)
    t0 = time.time()

    class Log:
        pass
    res = idmrg.run(h, d, a.chis, steps_per_chi=a.steps, tol_e=a.tol, seed=a.seed, log=False)
    with open(a.out, 'a') as f:
        for r in res:
            r.update(spec=a.spec, label=a.label, time_s=time.time() - t0)
            f.write(json.dumps(r) + '\n')
            print(f"{a.spec} chi={r['chi']} e={r['e_site']:.10f} S={r['S_AB']:.5f},{r['S_BA']:.5f} xi={r['xi']:.2f} "
                  f"n={r.get('n')} steps={r['steps']}", flush=True)
