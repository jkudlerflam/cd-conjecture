#!/usr/bin/env python3
"""Open-chain DMRG (FDMRG class of contrib/d2_followup_220_1282/fdmrg.py, unchanged) for survey points.
The end sites carry half the uniform field (sum of bond terms) unless --edgefield (full field).
usage: obc.py NAME --L 96 --chis 32 64 96 128 [--seed 1] [--edgefield]
Output data/obc/NAME_L{L}_s{seed}[_edge].json with S on every bond and bond energies.
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
from vumps import mpo_from_term, PAULI
from fdmrg import FDMRG

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--L', type=int, required=True)
ap.add_argument('--chis', type=int, nargs='+', default=[32, 64, 96, 128])
ap.add_argument('--sweeps', type=int, default=4)
ap.add_argument('--etol', type=float, default=1e-9)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--edgefield', action='store_true')
a = ap.parse_args()
h, r, c, note, f = model(a.name)
W = mpo_from_term(h, r)
edge = 0.5 * np.einsum('a,aij->ij', np.asarray(f), PAULI) if (a.edgefield and f is not None) else None
tag = '_edge' if edge is not None else ''
out = os.path.join(HERE, 'data', 'obc', f'{a.name}_L{a.L}_s{a.seed}{tag}.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
data = dict(name=a.name, note=note, c_exact=c, L=a.L, seed=a.seed, edge=edge is not None, results=[])
dm = FDMRG(W, a.L, seed=a.seed, edge=edge)
t0 = time.time()
for chi in a.chis:
    Eprev = None
    for sw in range(a.sweeps * 5):
        E, tr = dm.sweep(chi)
        if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < a.etol * a.L:
            break
        Eprev = E
    rec = dict(chi=chi, E=E, sweeps=sw + 1, trunc=tr, S=dm.S.tolist(), time_s=time.time() - t0)
    if r == 2:
        rec['Eb'] = dm.bond_energies(h).tolist()
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
    print(f"{a.name} L={a.L} chi={chi} E={E:.10f} sweeps={sw+1} trunc={tr:.1e} Smid={dm.S[a.L//2-1]:.5f} "
          f"t={time.time()-t0:.0f}s", flush=True)
