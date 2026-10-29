#!/usr/bin/env python3
"""Open-chain DMRG (FDMRG of contrib/d2_followup_220_1282/fdmrg.py, unchanged), entanglement on every bond.

NAME as in fes_fi.py.  For family points the two end sites carry the full uniform field (edge term f.s/2
added to the bond sum), as with --edgefield in contrib/d2_grid/obc.py.
usage: obc_fi.py NAME --L 128 --chis 32 64 128 [--tag T]
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
from fdmrg import FDMRG
import models_fi as MF

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--L', type=int, required=True)
ap.add_argument('--chis', type=int, nargs='+', default=[32, 64, 128])
ap.add_argument('--sweeps', type=int, default=4)
ap.add_argument('--etol', type=float, default=1e-9)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--tag', default=None)
a = ap.parse_args()
edge = None
if ':' in a.name:
    fam, p = MF.parse_name(a.name)
    J, f = MF.family(fam, p)
    h, d, c = V.bond(J, f), 2, None
    edge = 0.5 * np.einsum('a,aij->ij', np.asarray(f, float), V.PAULI)
else:
    h, d, c = MF.control_term(a.name)
    if a.name in ('ising', 'isxx'):          # full transverse field on the end sites as well
        hz = 0.5 * V.PAULI[2] * -1.0
        edge = hz if d == 2 else np.kron(hz, np.eye(2))
tag = a.tag or a.name.replace(':', '_').replace(',', '_').replace('=', '')
W = V.mpo_from_term(h, 2, d=d)
out = os.path.join(HERE, 'data', 'obc', f'{tag}_L{a.L}.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
data = dict(name=a.name, c_exact=c, L=a.L, results=[])
dm = FDMRG(W, a.L, d=d, seed=a.seed, edge=edge)
t0 = time.time()
for chi in a.chis:
    Eprev = None
    for sw in range(a.sweeps * 5):
        E, tr = dm.sweep(chi)
        if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < a.etol * a.L:
            break
        Eprev = E
    rec = dict(chi=chi, E=float(E), sweeps=sw + 1, trunc=float(tr), S=[float(x) for x in dm.S], time_s=time.time() - t0)
    if d == 2:
        rec['Eb'] = [float(np.real(x)) for x in dm.bond_energies(h)]
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
    print(f"{tag} L={a.L} chi={chi} E={E:.10f} sweeps={sw+1} trunc={tr:.1e} Smid={dm.S[a.L//2-1]:.5f} "
          f"t={time.time()-t0:.0f}s", flush=True)
