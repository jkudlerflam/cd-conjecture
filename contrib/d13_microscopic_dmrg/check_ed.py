#!/usr/bin/env python3
"""Check cqdmrg.py on the d = 13 bond against exact diagonalization of a short open chain.

The exact space is the full (Q, N, Mk) sector of the seed (finite K, illegal words included),
the same space the DMRG explores.  Writes data/check_ed.json.
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
    """Ground energy in the full (Q, N, Mk) sector of the seed on an open chain (sparse kron build)."""
    L = len(seed)
    q = RM.charges()
    hs = sp.csr_matrix(np.where(np.abs(h) > 0, h, 0))
    H = None
    for i in range(L - 1):
        t = sp.kron(sp.kron(sp.identity(d ** i, format='csr'), hs, format='csr'), sp.identity(d ** (L - i - 2), format='csr'), format='csr')
        H = t if H is None else H + t
    idx = np.indices([d] * L).reshape(L, -1)
    tot = q[idx].sum(axis=0)
    want = q[list(seed)].sum(axis=0)
    sel = np.flatnonzero(np.all(tot == want, axis=1))
    Hs = H[sel][:, sel]
    e, v = eigsh(Hs, k=1, which='SA', tol=1e-12)
    return float(e[0]), len(sel)


def main():
    out = []
    for (p, L, N, Nlog, K) in [(2, 5, 3, 1, 30.0), (2, 6, 3, 1, 30.0), (2, 6, 6, 1, 30.0)]:
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
