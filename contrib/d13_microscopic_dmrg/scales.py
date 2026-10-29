#!/usr/bin/env python3
"""Logical energy scales of small members of the C88 material bond, by exact diagonalization.

Frozen material ring (full occupancy, so markers and cell lengths are fixed) of M cells of
p + 1 sites each, projected onto the legal (zero-penalty, K = infinity) words reachable from a
seed.  In the sector with one payload fermion (one bit set, charge 2) the M lowest levels form
the logical band mu + 2|T| cos(k + phi); in the charge-0 sector the gap above the unique ground
state is the cell gap.  Imports codex/work_C88/bond13.py read-only.

usage: python3 -B scales.py --p 2 --M 3 [--D 1 --g .2 --lam .2 --td .1 --Delta x]
Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-29.
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'codex', 'work_C88'))
sys.dont_write_bytecode = True
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import bond13 as B


def legal_closure(seed, h, d, allowed_idx):
    trans = {i: [(int(j), h[j, i]) for j in np.flatnonzero(np.abs(h[:, i]) > 0) if j != i] for i in range(d * d)}
    diag = np.real(np.diag(h))
    seed = tuple(seed); states = [seed]; index = {seed: 0}
    rows, cols, vals = [], [], []
    L = len(seed); cur = 0
    while cur < len(states):
        w = states[cur]
        e = 0.0
        for i in range(L):
            j = (i + 1) % L
            e += diag[w[i] * d + w[j]]
            for q, v in trans[w[i] * d + w[j]]:
                z = list(w); z[i], z[j] = divmod(q, d); z = tuple(z)
                if not all(allowed_idx[z[k], z[(k + 1) % L]] for k in range(L)):
                    continue
                if z not in index:
                    index[z] = len(states); states.append(z)
                rows.append(index[z]); cols.append(cur); vals.append(v)
        rows.append(cur); cols.append(cur); vals.append(e)
        cur += 1
    H = sp.coo_matrix((vals, (rows, cols)), shape=(cur, cur)).tocsr()
    return H, states


def lowest(H, k):
    if H.shape[0] <= 2000:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(eigsh(H, k=k, which='SA', tol=1e-12)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--M', type=int, default=3)
    ap.add_argument('--D', type=float, default=1.0)
    ap.add_argument('--Delta', type=float, default=None)
    ap.add_argument('--g', type=float, default=0.2)
    ap.add_argument('--lam', type=float, default=0.2)
    ap.add_argument('--td', type=float, default=0.1)
    ap.add_argument('--flux', type=float, default=None)
    a = ap.parse_args()
    kw = dict(D=a.D, g=a.g, td=a.td, K=0.0, musel=0.0)
    if a.Delta is not None: kw['Delta'] = a.Delta
    if a.flux is not None: kw['flux'] = a.flux
    P = B.params(n=a.p, **kw)
    h = B.material_bond(P, a.lam)
    d = 11
    allowed = np.array([[B.allowed(x, y) for y in B.MAT] for x in B.MAT])
    ix = {x: i for i, x in enumerate(B.MAT)}
    vac = [ix[('R', 0)]] * a.p + [ix[('e', 0)]]       # all bits zero, walker at the left end
    out = dict(params={k: (float(v) if np.isscalar(v) else v) for k, v in P.items()}, lam=a.lam, p=a.p, M=a.M)
    H0, st0 = legal_closure(vac * a.M, h, d, allowed)
    e0 = lowest(H0, 3)
    one = [ix[('R', 0)], ix[('R', 1)]] + [ix[('R', 0)]] * (a.p - 2) + [ix[('e', 0)]]
    H1, st1 = legal_closure(one + vac * (a.M - 1), h, d, allowed)
    e1 = lowest(H1, a.M + 2)
    band = e1[:a.M] - e0[0]
    out.update(dim0=H0.shape[0], dim1=H1.shape[0], E0=e0.tolist(), cell_gap=float(e0[1] - e0[0]),
               band=band.tolist(), bandwidth=float(band.max() - band.min()), next_gap=float(e1[a.M] - e1[a.M - 1]))
    print(json.dumps(out))


if __name__ == '__main__':
    main()
