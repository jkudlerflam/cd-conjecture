#!/usr/bin/env python3
"""Check cqdmrg.py on the d = 13 bond against exact diagonalization of a short open chain.

The exact space is the closure of the seed under all moves of the bond (finite K, illegal words
included), so it is the full charge sector connected to the seed.  Writes data/check_ed.json.
Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-29.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import run_mobile as RM
import cqdmrg as Q
B = RM.B


def ed(h, seed, d=13):
    trans = {i: [(int(j), h[j, i]) for j in np.flatnonzero(np.abs(h[:, i]) > 0)] for i in range(d * d)}
    states = [tuple(seed)]; index = {states[0]: 0}; rows, cols, vals = [], [], []
    L = len(seed); cur = 0
    while cur < len(states):
        w = states[cur]
        for i in range(L - 1):
            for q, v in trans[w[i] * d + w[i + 1]]:
                z = list(w); z[i], z[i + 1] = divmod(q, d); z = tuple(z)
                if z not in index:
                    index[z] = len(states); states.append(z)
                rows.append(index[z]); cols.append(cur); vals.append(v)
        cur += 1
    H = sp.coo_matrix((vals, (rows, cols)), shape=(cur, cur)).tocsr()
    e, v = eigsh(H, k=1, which='SA', tol=1e-13)
    return float(e[0]), cur


def main():
    out = []
    for (p, L, N, Nlog, K) in [(2, 6, 6, 1, 30.0), (2, 7, 6, 1, 30.0), (2, 8, 6, 1, 30.0)]:
        P = B.params(n=p, D=0.0, Delta=0.5, g=0.3, td=0.2, K=K, musel=0.0)
        h = B.bond(P, lam=0.3, J=0.5, t=1.0)
        cfg, Mk = RM.config(p, L, N, Nlog)
        E, dim = ed(h, cfg)
        dm, hist = Q.run(h, RM.charges(), L, cfg, [64, 256], sweeps=3, maxsweeps=10, etol=1e-13, log=lambda s: None)
        rec = dict(p=p, L=L, N=N, Nlog=Nlog, K=K, dim=dim, E_ed=E, E_dmrg=hist[-1]['E'], diff=hist[-1]['E'] - E)
        print(json.dumps(rec), flush=True); out.append(rec)
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'data', 'check_ed.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
