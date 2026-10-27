#!/usr/bin/env python3
"""Spectrum in every momentum sector plus ground-state entanglement for one chain.

usage: run_point.py NAME --sizes 8 9 ... 26 [--k 4] [--procs P] [--out results/NAME.json]
NAME is a key of candidates_params.json, controls.CONTROLS, cluster_ising[_lam],
or extra_points.json.  Momentum sectors are distributed over P worker processes.
Resumes: sizes already in the output file are skipped.
"""
import os
for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "1")
import argparse, json, sys, time
import multiprocessing as mp
import numpy as np
import ed
from controls import local_term

_H = None
_L = None


def _work(m):
    H = ed.sector_hamiltonian(_L, _H, m)
    e, u = ed.lowest(H, k=_K)
    return m, e, u[:, 0]


def main():
    global _H, _L, _K
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--sizes', type=int, nargs='+', required=True)
    ap.add_argument('--k', type=int, default=4)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--out', default=None)
    ap.add_argument('--no-ent', action='store_true')
    a = ap.parse_args()
    h, meta = local_term(a.name)
    out = a.out or f'results/{a.name}.json'
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    data = json.load(open(out)) if os.path.exists(out) else {'name': a.name, 'meta': meta, 'sizes': {}}
    _H, _K = h, a.k
    for L in a.sizes:
        if str(L) in data['sizes']:
            continue
        t0 = time.time()
        _L = L
        ed.orbits(L)  # computed before fork, shared copy-on-write
        if a.procs > 1:
            with mp.get_context('fork').Pool(min(a.procs, L)) as pool:
                res = pool.map(_work, range(L), chunksize=1)
        else:
            res = [_work(m) for m in range(L)]
        spec = {m: e for m, e, _ in res}
        E0 = min(v[0] for v in spec.values())
        tol = 1e-8 * max(1.0, abs(E0))
        gms = sorted(m for m, v in spec.items() if v[0] < E0 + tol)
        rec = {'E0': float(E0), 'ground_momenta': gms,
               'spectrum': {str(m): list(map(float, v)) for m, v in spec.items()}, 'time_s': None}
        if not a.no_ent:
            u0 = next(u for m, e, u in res if m == gms[0])
            del res
            S = ed.entropies(ed.expand(L, gms[0], u0), L)
            rec['S'] = S.tolist()
        rec['time_s'] = time.time() - t0
        data['sizes'][str(L)] = rec
        tmp = out + '.tmp'
        json.dump(data, open(tmp, 'w'))
        os.replace(tmp, out)
        print(f"{a.name} L={L} E0={E0:.10f} m0={gms} t={rec['time_s']:.1f}s", flush=True)


if __name__ == '__main__':
    main()
