#!/usr/bin/env python3
"""Finite-entanglement scaling with the damped VUMPS of contrib/d2_followup_220_1282/vumps.py (unchanged).

NAME is a control (ising, xx, isxx) or a family point 'FAM:k=v,...' (models_fi.parse_name).
With --cells, a chi = 8, 16 cell scan picks the lowest-energy converged injective cell (smaller n wins
within --etol per site), as in contrib/d2_grid/fes_auto.py; otherwise --n is used.  Then chi grows
through --chis and data/fes/TAG_nN.json is written after every chi.  c = 6 dS/dln xi (analyze.py).
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
import vumps as V
import models_fi as MF

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--tag', default=None)
ap.add_argument('--n', type=int, default=2)
ap.add_argument('--cells', type=int, nargs='*', default=None)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64, 96, 128])
ap.add_argument('--etol', type=float, default=1e-6)
ap.add_argument('--maxit', type=int, default=500)
a = ap.parse_args()
if ':' in a.name:
    fam, p = MF.parse_name(a.name)
    J, f = MF.family(fam, p)
    h, d, c = V.bond(J, f), 2, None
else:
    h, d, c = MF.control_term(a.name)
tag = a.tag or a.name.replace(':', '_').replace(',', '_').replace('=', '')
W = V.mpo_from_term(h, 2, d=d)
os.makedirs(os.path.join(HERE, 'data', 'fes'), exist_ok=True)
t0 = time.time()
n = a.n
cells = []
if a.cells:
    hn = np.linalg.norm(h, 2)
    for m in a.cells:
        try:
            _, rs = V.run([W] * m, m, [8, 16], d=d, tol=1e-7, maxit=150, seed=0, mode='damped', log=lambda r: None)
            rr = rs[-1]
            cells.append(dict(n=m, e=rr['e'], err=rr['err'], xi=rr['xi'], S=float(np.mean(rr['S']))))
        except Exception as ex:
            cells.append(dict(n=m, e=np.inf, err=1.0, error=str(ex)))
        print(tag, 'cell', cells[-1], f't={time.time()-t0:.0f}', flush=True)
    good = [x for x in cells if np.isfinite(x['e']) and x['err'] < 1e-4 and x.get('xi', np.inf) < 1e5 and x['e'] > -hn]
    if not good:
        good = [x for x in cells if np.isfinite(x['e']) and x['err'] < 1e-2 and x['e'] > -hn]
    emin = min(x['e'] for x in good)
    n = min(x['n'] for x in good if x['e'] <= emin + a.etol)
out = os.path.join(HERE, 'data', 'fes', f'{tag}_n{n}.json')
data = dict(name=a.name, c_exact=c, n=n, cells=cells, results=[])


def save(psi, rec):
    rec['time_s'] = time.time() - t0
    rec['S'] = [float(x) for x in rec['S']]
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w'), default=float); os.replace(out + '.tmp', out)
    print(f"{tag} n={n} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
          f"S={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)


V.run([W] * n, n, a.chis, d=d, tol=1e-9, maxit=a.maxit, seed=0, mode='damped', log=lambda r: None, on_chi=save)
