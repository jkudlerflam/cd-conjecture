"""Sparse ED for translation-invariant nearest-neighbour qubit chains with PBC.
H = sum_j [ S_j^T J S_{j+1} + B . S_j ], S = sigma/2 (or a bond-dependent list of J's).
Used to check the Kitaev spin chain in a [110] field (Sorensen et al. 2022) in its uniform gauge.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
S = [sx, sy, sz]


def site_op(o, j, L):
    return sp.kron(sp.kron(sp.identity(2 ** j, format='csr'), sp.csr_matrix(o)), sp.identity(2 ** (L - j - 1), format='csr'), format='csr')


def build(Js, Bs, L, pbc=True):
    """Js: function bond index -> 3x3 J; Bs: function site -> B vector."""
    ops = [[site_op(S[a], j, L) for a in range(3)] for j in range(L)]
    H = sp.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    nb = L if pbc else L - 1
    for j in range(nb):
        J = Js(j); k = (j + 1) % L
        for a in range(3):
            for b in range(3):
                if abs(J[a, b]) > 1e-14:
                    H = H + J[a, b] * (ops[j][a] @ ops[k][b])
    for j in range(L):
        B = Bs(j)
        for a in range(3):
            if abs(B[a]) > 1e-14:
                H = H + B[a] * ops[j][a]
    return H


def ent(psi, L, l):
    m = psi.reshape(2 ** l, 2 ** (L - l))
    s = np.linalg.svd(m, compute_uv=False) ** 2
    s = s[s > 1e-14]
    return float(-(s * np.log(s)).sum())
