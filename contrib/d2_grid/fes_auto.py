#!/usr/bin/env python3
"""Cell scan followed by finite-entanglement scaling, for one survey point.

1. For each cell size n in --cells, damped VUMPS at chi = 8 then 16 (vumps.py of
   contrib/d2_followup_220_1282, unchanged); record the energy per site.
2. Among converged injective runs, take the cell with the lowest energy (smaller n wins within --etol per site), and grow chi through
   --chis in that cell, writing data/fes/NAME_nN.json after every chi (same format as fes.py).
usage: fes_auto.py NAME [--cells 1 2 ... 12] [--chis 8 12 16 24 32 48 64 96]
NAME is a key of candidates.json (9-vector convention of strata.py) or a control of models_grid.py.
Written for contrib/d2_grid (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from models_grid import model
import vumps as V

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--cells', type=int, nargs='+', default=list(range(1, 13)))
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64, 96])
ap.add_argument('--etol', type=float, default=1e-6)
ap.add_argument('--maxit', type=int, default=400)
a = ap.parse_args()
h, r, c, note, _ = model(a.name)
W = V.mpo_from_term(h, r)
os.makedirs(os.path.join(HERE, 'data', 'fes'), exist_ok=True)
t0 = time.time()
cells = []
for n in a.cells:
    try:
        _, rs = V.run([W] * n, n, [8, 16], tol=1e-7, maxit=150, seed=0, mode='damped', log=lambda r: None)
        rr = rs[-1]
        cells.append(dict(n=n, e=rr['e'], err=rr['err'], S=float(np.mean(rr['S'])), xi=rr['xi']))
    except Exception as ex:
        cells.append(dict(n=n, e=np.inf, error=str(ex)))
    print(a.name, 'cell', cells[-1], f't={time.time()-t0:.0f}', flush=True)
json.dump(dict(name=a.name, chi=16, res=cells), open(os.path.join(HERE, 'data', 'fes', f'{a.name}_cells.json'), 'w'),
          default=float)
# a cell run is usable only if it converged to an injective state: small gradient, finite xi, and an
# energy above -|h| (non-injective fixed points give meaningless energies)
hn = np.linalg.norm(h, 2)
good = [x for x in cells if np.isfinite(x['e']) and x.get('err', 1) < 1e-4 and x.get('xi', np.inf) < 1e5 and x['e'] > -hn]
if not good:
    good = [x for x in cells if np.isfinite(x['e']) and x.get('err', 1) < 1e-2 and x['e'] > -hn]
emin = min(x['e'] for x in good)
best = min(x['n'] for x in good if x['e'] <= emin + a.etol)
out = os.path.join(HERE, 'data', 'fes', f'{a.name}_n{best}.json')
data = dict(name=a.name, note=note, c_exact=c, n=best, mode='damped', results=[])


def save(psi, rec):
    rec['time_s'] = time.time() - t0
    rec['S'] = [float(x) for x in rec['S']]
    rec['mag'] = V.local_expect(psi, list(V.PAULI)).tolist()
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
    print(f"{a.name} n={best} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
          f"S={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)


V.run([W] * best, best, a.chis, tol=1e-8, maxit=a.maxit, seed=0, mode='damped', log=lambda r: None, on_chi=save)
