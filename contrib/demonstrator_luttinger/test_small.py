#!/usr/bin/env python3
"""Small-size checks of dmrg.py against exact diagonalization (logical M = 12, micro M = 3)."""
import numpy as np, scipy.sparse as sps
from scipy.sparse.linalg import eigsh
import chain as CH, dmrg as DM

def dense_H(ch, Q0):
    L = ch.L; dims = ch.dims
    def emb(ops):
        out = sps.identity(1, format='csr')
        for s in range(L):
            out = sps.kron(out, sps.csr_matrix(ops.get(s, np.eye(dims[s]))), format='csr')
        return out
    H = sum(emb({s: ch.onsite[s]}) for s in range(L))
    for s in range(L - 1):
        for A, B in ch.bonds[s]:
            H = H + emb({s: A, s + 1: B})
    q = np.zeros(1)
    for s in range(L):
        q = (q[:, None] + ch.charge[s][None, :]).ravel()
    idx = np.where(np.abs(q - Q0) < 1e-9)[0]
    return H.tocsr()[idx][:, idx]

for name, ch, Q0, chis in [('logical', CH.logical(12, -0.28125, 0.0208, lam=0.1), 2.0, [16, 32, 64]),
                           ('micro', CH.micro(3, 0.15, 2.0, -1.0), 0.5, [16, 32, 64, 128])]:
    Hs = dense_H(ch, Q0)
    w = eigsh(Hs, k=3, which='SA', tol=1e-14)[0]; w.sort()
    Ws = ch.mpo(Q0=Q0, w=0.05)
    dm, E, _ = DM.run_ground(Ws, ch.dims, chis, log=None, sweeps_min=4)
    dm1, E1, _ = DM.run_ground(Ws, ch.dims, chis, log=None, sweeps_min=4, ortho=[(1.0, dm.M)])
    print(name, 'ED', w[:2], 'DMRG', E, E1, 'diff', E - w[0], E1 - w[1], 'overlap', dm1.overlap(dm.M))
