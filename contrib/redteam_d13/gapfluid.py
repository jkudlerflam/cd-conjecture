"""Gap-space vacancy fluid of the projected d = 13 chain (C91 eq. for H0).

N material letters on a ring, N_v vacancies distributed in the N gaps.
A vacancy moves from gap j to gap j+1 (it crosses letter j+1) with amplitude -t,
unit norm (C91: ||E^+-|| = 1). Contact on material bond j is c_j = 1(g_j = 0).
The cyclic centre-of-mass gluing is dropped (periodic gap ring); local
quantities converge to the physical half-density values as N grows.
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def compositions(N, Nv):
    """All (g_0..g_{N-1}) with g_j >= 0 and sum = Nv, as an int8 array."""
    out = []
    for bars in itertools.combinations(range(Nv + N - 1), N - 1):
        prev = -1
        g = []
        for b in bars:
            g.append(b - prev - 1)
            prev = b
        g.append(Nv + N - 1 - prev - 1)
        out.append(g)
    return np.array(out, dtype=np.int16)


def encode(G, base):
    w = base ** np.arange(G.shape[1], dtype=np.int64)
    return G.astype(np.int64) @ w


def hop_matrix(G):
    """Sparse -sum_j (E_j^- E_{j+1}^+ + h.c.) on the gap ring (t = 1)."""
    M, N = G.shape
    base = int(G.max()) + 2
    codes = encode(G, base)
    order = np.argsort(codes)
    sc = codes[order]
    rows, cols = [], []
    w = base ** np.arange(N, dtype=np.int64)
    for j in range(N):
        jp = (j + 1) % N
        ok = G[:, j] > 0
        src = np.nonzero(ok)[0]
        newc = codes[src] - w[j] + w[jp]
        pos = np.searchsorted(sc, newc)
        dst = order[pos]
        assert np.all(sc[pos] == newc)
        rows.append(dst)
        cols.append(src)
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    A = sp.csr_matrix((-np.ones(len(r)), (r, c)), shape=(M, M))
    return (A + A.T).tocsr()


def ground(H):
    vals, vecs = spl.eigsh(H, k=2, which='SA', tol=1e-13)
    i = np.argsort(vals)
    return vals[i], vecs[:, i]


def static_susceptibility(H, E0, psi, ops):
    """S_ab = sum_{n != 0} <0|A_a|n><n|A_b|0>/(E_n-E0) for diagonal ops.

    Solved as (H-E0) x_b = Q A_b psi on the complement of psi with MINRES.
    """
    M = H.shape[0]
    X = []
    for a in ops:
        v = a * psi
        v = v - psi * (psi @ v)
        Hs = H - E0 * sp.identity(M, format='csr')
        x, info = spl.minres(Hs, v, rtol=1e-12, maxiter=20000)
        x = x - psi * (psi @ x)
        X.append((x, v))
    n = len(ops)
    S = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            S[a, b] = X[a][1] @ X[b][0]
    return S


def response_vector(H, E0, psi, a):
    """x = (H-E0)^-1 Q a psi on the complement of the ground state psi."""
    M = H.shape[0]
    v = a * psi
    v = v - psi * (psi @ v)
    op = spl.LinearOperator((M, M), matvec=lambda y: H @ y - E0 * y + psi * (psi @ y), dtype=float)
    x, info = spl.minres(op, v, rtol=1e-13, maxiter=50000)
    assert info == 0
    return x - psi * (psi @ x)


def free_fermion_alpha(L, N):
    """<n_i n_{i+1}>/rho for N free fermions on L sites, closed shell chosen by parity."""
    # periodic for odd N, antiperiodic for even N (Jordan-Wigner of hard-core bosons)
    shift = 0.0 if N % 2 == 1 else 0.5
    ks = 2 * np.pi * (np.arange(L) + shift) / L
    ks = np.where(ks > np.pi, ks - 2 * np.pi, ks)
    occ = np.argsort(-np.cos(ks))[:N]
    k = ks[occ]
    G0 = N / L
    G1 = np.sum(np.exp(1j * k)).real / L
    nn = G0 * G0 - abs(G1) ** 2
    return nn / G0
