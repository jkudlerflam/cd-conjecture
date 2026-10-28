#!/usr/bin/env python3
"""Logical XXZ rings H = Jxy sum (s+s- + h.c.) + Jzz sum ZZ (+ a0 per bond), charge and momentum resolved.

Used for the calibration (Jzz = 0, the XX chain) and for the effective logical chain (exact two-cell
couplings) at the same ring sizes as the microscopic ring spectra.  Returns, per (2Q, m), the lowest levels.
Claude (Opus 5.5), 2026-10-27.
"""
import itertools
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh, LinearOperator


def sector(M, nup):
    states = [sum(1 << i for i in c) for c in itertools.combinations(range(M), nup)]
    return np.array(states, dtype=np.int64), {s: i for i, s in enumerate(states)}


def build(M, nup, Jxy, Jzz, a0=0.0):
    st, ix = sector(M, nup)
    rows, cols, vals = [], [], []
    for i, s in enumerate(st):
        diag = a0 * M
        for j in range(M):
            jn = (j + 1) % M
            bj, bn = (s >> j) & 1, (s >> jn) & 1
            diag += Jzz * (1 if bj == bn else -1)
            if bj != bn:
                t = s ^ ((1 << j) | (1 << jn))
                rows.append(ix[t]); cols.append(i); vals.append(Jxy)
        rows.append(i); cols.append(i); vals.append(diag)
    H = sps.csr_matrix((vals, (rows, cols)), shape=(len(st), len(st)))
    # translation j -> j+1
    T_rows = []
    for i, s in enumerate(st):
        t = ((s << 1) | (s >> (M - 1))) & ((1 << M) - 1)
        T_rows.append(ix[t])
    T = sps.csr_matrix((np.ones(len(st)), (T_rows, np.arange(len(st)))), shape=(len(st), len(st)))
    return H, T


def levels(M, Q2, m, Jxy, Jzz, a0=0.0, nev=3, shift=10.0):
    nup = (M + Q2) // 2                       # up spins carry charge +1/2
    H, T = build(M, nup, Jxy, Jzz, a0)
    N = H.shape[0]
    k = 2 * np.pi * m / M
    fac = 1.0 if (m % M == 0 or 2 * m == M) else 2.0
    Ts = [sps.identity(N, format='csr')]
    for r in range(1, M):
        Ts.append(T @ Ts[-1])
    Pm = fac / M * sum(np.cos(k * r) * Ts[r] for r in range(M))
    if N <= 400:
        A = (Pm @ (H - shift * sps.identity(N)) @ Pm).toarray()
        w = np.linalg.eigvalsh((A + A.T) / 2)
        w = w[w < -1e-9][:nev] + shift
        return np.sort(w)
    op = LinearOperator((N, N), matvec=lambda v: Pm @ (H @ v - shift * v), dtype=float)
    v0 = Pm @ np.random.default_rng(0).normal(size=N)
    w = eigsh(op, k=nev, which='SA', v0=v0, tol=1e-12, return_eigenvectors=False)
    return np.sort(w) + shift


if __name__ == '__main__':
    for M in (4, 6, 8):
        print(M, 'XX Q=0 m=0', levels(M, 0, 0, -1.0, 0.0), 'm=1', levels(M, 0, 1, -1.0, 0.0), 'Q=1 m=0', levels(M, 2, 0, -1.0, 0.0))
        # free-fermion check of E0: Jxy = -1 hopping -> energies -2 cos k
