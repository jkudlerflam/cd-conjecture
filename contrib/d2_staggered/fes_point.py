#!/usr/bin/env python3
"""Finite-entanglement scaling in the thermodynamic limit with the damped VUMPS of
contrib/d2_followup_220_1282/vumps.py (imported unchanged; dense, no symmetry imposed).

Model in the ORIGINAL frame, H = sum_j [X X - Y Y + u Z Z + b Z] (one-site translation invariant),
two-site unit cell.  Controls: xx (c=1), ising (c=1/2), isxx (Ising (x) XX on a d=4 site, c=3/2).
usage: fes_point.py stag --u 1.35 --b 2 --chis 8 12 16 24 32 48 64 96 128
Output data/fes/<tag>.json with S, xi = -n/ln|l1/l0| and the transfer spectrum at each chi.
c is fitted by analyze.py from S = (c/6) ln xi + const.
Written for contrib/d2_staggered (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
import vumps as V

X, Y, Z = V.PAULI
I2 = np.eye(2)
ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--u', type=float, default=0.0)
ap.add_argument('--b', type=float, default=0.0)
ap.add_argument('--n', type=int, default=2)
ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64, 96, 128])
ap.add_argument('--maxit', type=int, default=600)
ap.add_argument('--tol', type=float, default=1e-9)
ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args()

d = 2
if a.name == 'stag':
    h = V.bond(np.diag([1., -1., a.u]), [0, 0, a.b]); tag = f'stag_u{a.u:.4f}_b{a.b:.4f}'
elif a.name == 'xx':
    h = V.bond(np.diag([1., 1., 0.]), [0, 0, 0]); tag = 'xx'
elif a.name == 'ising':
    h = V.bond(np.diag([-1., 0, 0]), [0, 0, -1.]); tag = 'ising'
elif a.name == 'isxx':
    hi = V.bond(np.diag([-1., 0, 0]), [0, 0, -1.]); hx = V.bond(np.diag([1., 1., 0.]), [0, 0, 0])
    # sites (i1 x1)(i2 x2): rearrange kron(hi on i1 i2, hx on x1 x2)
    # H_bond = hi (x) 1 + 1 (x) hx, both bond terms already carry half of each on-site field
    one = np.eye(4).reshape(2, 2, 2, 2)
    T = (np.einsum('abcd,efgh->aebfcgdh', hi.reshape(2, 2, 2, 2), one)
         + np.einsum('abcd,efgh->aebfcgdh', one, hx.reshape(2, 2, 2, 2)))
    h = T.reshape(16, 16); d = 4; tag = 'isxx'
else:
    raise KeyError(a.name)
W = V.mpo_from_term(h, 2, d=d)
out = os.path.join(HERE, 'data', 'fes', f'{tag}_n{a.n}.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
data = dict(name=a.name, u=a.u, b=a.b, n=a.n, results=[])
t0 = time.time()


def on_chi(psi, rec):
    rec['time_s'] = time.time() - t0
    data['results'].append(rec)
    json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)


V.run([W] * a.n, a.n, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode='damped', on_chi=on_chi)
print('done', out, time.time() - t0)
