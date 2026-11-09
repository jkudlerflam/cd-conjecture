#!/usr/bin/env python3
"""C170 (independent code): sparse ED of H = sum_j [s_j^T J s_{j+1} + B.s_j] (Pauli), open or periodic."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
P = [sp.csr_matrix(np.array(m, complex)) for m in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
I2 = sp.identity(2, format='csr', dtype=complex)

def op(o, j, L):
    out = sp.identity(1, format='csr', dtype=complex)
    for k in range(L):
        out = sp.kron(out, o if k == j else I2, format='csr')
    return out

def ham(J, B, L, pbc=False, Js=None, Bs=None):
    S = [[op(P[a], j, L) for a in range(3)] for j in range(L)]
    H = sp.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    nb = L if pbc else L - 1
    for j in range(nb):
        k = (j + 1) % L
        Jj = J if Js is None else Js[j]
        for a in range(3):
            for b in range(3):
                if abs(Jj[a][b]) > 0:
                    H = H + Jj[a][b] * (S[j][a] @ S[k][b])
    for j in range(L):
        Bj = B if Bs is None else Bs[j]
        for a in range(3):
            if abs(Bj[a]) > 0:
                H = H + Bj[a] * S[j][a]
    return H, S

def low(H, k=1):
    w = sla.eigsh(H, k=k, which='SA', return_eigenvectors=False, tol=1e-12, maxiter=20000)
    return np.sort(w.real)
