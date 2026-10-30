"""Fermionic (Majorana) decomposition of a two-cell n-flavour payload bond.

Qubit basis of one cell: index z = sum_f bit_f 2^f (flavour f = bit f).
Two-cell index i*2^n + j (left cell i, right cell j). Jordan-Wigner order:
(left, f=0), (left, f=1), ..., (right, f=0), ... ; this is the order in
which the rebuilt bond carries its fermion signs (checked by `jw_check`).
The decomposition is Hilbert-Schmidt orthogonal: K = sum_S k_S gamma_S over
even Majorana monomials gamma_S, so the quadratic part is the one-body
operator relative to infinite temperature and the rest is interaction.
"""
import itertools
import numpy as np


def ops(nm):
    """c_m (annihilation) on nm modes, JW strings, basis index = sum b_m 2^m
    with the mode order m = 0..nm-1 matching the qubit layout below."""
    dim = 2 ** nm
    C = []
    for m in range(nm):
        M = np.zeros((dim, dim))
        for s in range(dim):
            if (s >> m) & 1:
                sign = (-1) ** bin(s & ((1 << m) - 1)).count('1')
                M[s ^ (1 << m), s] = sign
        C.append(M)
    return C


def mode_layout(n):
    """Map two-cell qubit index (i*2^n + j) to the JW mode bit string.
    JW mode m: left flavour f -> m = f, right flavour f -> m = n + f."""
    d = 2 ** n
    P = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            s = i | (j << n)
            P[s, i * d + j] = 1
    return P          # P maps qubit-layout vectors to JW-mode vectors


def majoranas(nm):
    C = ops(nm)
    G = []
    for c in C:
        G.append(c + c.T)
        G.append(1j * (c.T - c))
    return G


def decompose(K, n):
    """Return (scalar, h, rest) with K = scalar + sum_ab h_ab c_a^dag c_b - tr(h)/2 + rest,
    h the 2n x 2n one-body matrix on JW modes (left flavours, right flavours)."""
    nm = 2 * n
    P = mode_layout(n)
    Kj = P @ K @ P.T
    dim = 2 ** nm
    C = ops(nm)
    scalar = np.trace(Kj).real / dim
    # one-body coefficient: for HS-orthogonal decomposition,
    # h_ab = Tr(K [c_a^dag c_b - delta_ab/2]) / (dim/4) ... use projection on the
    # operator space spanned by {c_a^dag c_b - delta/2}, which is HS-orthogonal to 1
    # and to all quartic Majorana monomials.
    basis, idx = [], []
    for a in range(nm):
        for b in range(nm):
            O = C[a].T @ C[b] - (0.5 if a == b else 0.0) * np.eye(dim)
            basis.append(O)
            idx.append((a, b))
    Gm = np.array([[np.trace(x.conj().T @ y) for y in basis] for x in basis])
    rhs = np.array([np.trace(x.conj().T @ Kj) for x in basis])
    coef = np.linalg.solve(Gm, rhs)
    h = np.zeros((nm, nm), complex)
    for (a, b), v in zip(idx, coef):
        h[a, b] = v
    quad = sum(v * O for v, O in zip(coef, basis))
    rest = Kj - scalar * np.eye(dim) - quad
    return scalar, h, P.T @ rest @ P


def quad_operator(h, n):
    nm = 2 * n
    C = ops(nm)
    dim = 2 ** nm
    Q = sum(h[a, b] * (C[a].T @ C[b] - (0.5 if a == b else 0) * np.eye(dim)) for a in range(nm) for b in range(nm))
    P = mode_layout(n)
    return P.T @ Q @ P
