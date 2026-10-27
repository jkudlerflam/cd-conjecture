#!/usr/bin/env python3
"""TCI spectroscopy control: O'Brien-Fendley chain
   H = sum_j [ -Z_j Z_{j+1} - X_j + lam (X_j Z_{j+1} Z_{j+2} + Z_j Z_{j+1} X_{j+2}) ],
periodic, Z2-resolved (parity prod X), tricritical Ising at lam_c ~ 0.428.
Prints (E - E0)/(E_eps - E0) for the lowest levels in each parity sector, where eps is the
lowest even excitation (x = 1/5): CFT values even: 1/5, 6/5 (eps'), 1 + 1/5 (desc.), 2 (T);
odd: 3/40 (sigma), 7/8 (sigma').  So ratios: sigma 0.375, eps 1, sigma' 4.375, eps' 6."""
import sys
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh


def ops(L):
    X = sparse.csr_matrix(np.array([[0, 1], [1, 0]], float))
    Z = sparse.csr_matrix(np.array([[1, 0], [0, -1]], float))
    I = sparse.identity(2, format='csr')

    def site(o, j):
        out = None
        for k in range(L):
            m = o if k == j else I
            out = m if out is None else sparse.kron(out, m, format='csr')
        return out
    return [site(X, j) for j in range(L)], [site(Z, j) for j in range(L)]


def run(L, lam, nlev=12):
    Xs, Zs = ops(L)
    H = sparse.csr_matrix((2 ** L, 2 ** L))
    for j in range(L):
        a, b, c = j, (j + 1) % L, (j + 2) % L
        H = H - Zs[a] @ Zs[b] - Xs[a] + lam * (Xs[a] @ Zs[b] @ Zs[c] + Zs[a] @ Zs[b] @ Xs[c])
    # parity prod X is diagonal in the X basis; use projector via basis rotation-free trick: P = prod X
    P = Xs[0]
    for j in range(1, L):
        P = P @ Xs[j]
    out = {}
    for s in (+1, -1):
        Hs = H + 50.0 * (sparse.identity(2 ** L) - s * P)   # push the other parity up
        e = np.sort(eigsh(Hs, k=nlev, which='SA', tol=1e-10)[0])
        out[s] = e
    return out


lam = float(sys.argv[1]) if len(sys.argv) > 1 else 0.428
for L in (8, 12, 16):
    o = run(L, lam)
    e0 = o[1][0]
    eps = o[1][1] - e0
    print('L=%d  even: %s' % (L, np.round((o[1][1:7] - e0) / eps, 3).tolist()))
    print('      odd : %s   L*eps/2pi=%.4f' % (np.round((o[-1][:5] - e0) / eps, 3).tolist(), L * eps / 2 / np.pi))
