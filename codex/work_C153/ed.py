#!/usr/bin/env python3
"""Exact diagonalization of the staggered-charge chain on even rings, rotated frame.

H' = sum_j [X_j X_j+1 + Y_j Y_j+1 - u Z_j Z_j+1] + b sum_j (-1)^j Z_j      (Pauli matrices)
unitarily equivalent (pi rotation about x on odd sites) to H = sum [XX - YY + u ZZ + b Z].
Q = sum_j (-1)^j Z_j (original) = sum_j Z_j (rotated) = L - 2N, N = number of 1-bits (Z = -1).
The wrap bond (L-1, 0) carries the twist phase exp(i phi) on sigma^+_{L-1} sigma^-_0.
Written for codex/work_C153 (Claude Opus agent, 2026-11-07).
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from itertools import combinations


def basis(L, N):
    st = np.array([sum(1 << i for i in c) for c in combinations(range(L), N)], dtype=np.int64)
    st.sort()
    return st


def hamiltonian(L, N, u, b, phi=0.0):
    st = basis(L, N)
    D = len(st)
    bits = ((st[:, None] >> np.arange(L)[None, :]) & 1).astype(np.int8)   # 1 means Z = -1
    z = 1 - 2 * bits.astype(np.float64)
    diag = -u * np.sum(z * np.roll(z, -1, axis=1), axis=1) + b * (z @ ((-1.0) ** np.arange(L)))
    rows, cols, vals = [np.arange(D)], [np.arange(D)], [diag.astype(complex)]
    for j in range(L):
        k = (j + 1) % L
        diff = bits[:, j] != bits[:, k]
        src = np.nonzero(diff)[0]
        new = st[src] ^ ((1 << j) | (1 << k))
        dst = np.searchsorted(st, new)
        amp = np.full(len(src), 2.0, dtype=complex)
        if k == 0 and phi != 0.0:
            # sigma^+_{L-1} sigma^-_0 e^{i phi} + h.c.; sigma^+ raises Z (bit 1 -> 0)
            # bit_{L-1} = 1, bit_0 = 0 -> after: bit_{L-1} = 0, bit_0 = 1 : that is s+_{L-1} s-_0
            up = bits[src, j] == 1
            amp = np.where(up, 2.0 * np.exp(1j * phi), 2.0 * np.exp(-1j * phi))
        rows.append(dst); cols.append(src); vals.append(amp)
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
    if phi == 0.0:
        H = H.real.tocsr()
    return H, st


def lowest(L, N, u, b, phi=0.0, k=1, vecs=False):
    H, st = hamiltonian(L, N, u, b, phi)
    if H.shape[0] <= 400:
        w, v = np.linalg.eigh(H.toarray())
        w, v = w[:k], v[:, :k]
    else:
        w, v = eigsh(H, k=k, which='SA', tol=1e-12)
        o = np.argsort(w); w, v = w[o], v[:, o]
    return (w, v, st) if vecs else w


def luttinger(L, u, b, phi=0.05):
    """K, v from charge gap and twist curvature on a ring; returns dict."""
    N0 = L // 2
    E0 = lowest(L, N0, u, b)[0]
    E1 = lowest(L, N0 + 1, u, b)[0]
    Ep = lowest(L, N0, u, b, phi=phi)[0]
    dE1 = E1 - E0
    dEt = Ep - E0
    vK = 2 * np.pi * L * dEt / phi ** 2          # E(phi) - E0 = v K phi^2 / (2 pi L)
    v_over_K = 2 * L * dE1 / np.pi              # E(N0+1) - E0 = pi v / (2 K L)
    if vK <= 0 or v_over_K <= 0:
        return dict(L=L, E0=E0, dE1=dE1, dEt=dEt, K=np.nan, v=np.nan)
    return dict(L=L, E0=E0, dE1=dE1, dEt=dEt, K=float(np.sqrt(vK / v_over_K)), v=float(np.sqrt(vK * v_over_K)))


if __name__ == '__main__':
    import sys
    # calibration: XXZ at b = 0, exact K = pi/(2 arccos u), v = 2 pi sin(g)/g with cos g = -u
    for u in [0.0, 0.5, 0.9]:
        g = np.arccos(-u)
        for L in [12, 16, 20]:
            d = luttinger(L, u, 0.0)
            print(f'u={u} L={L} K={d["K"]:.4f} (exact {np.pi/(2*np.arccos(u)):.4f}) v={d["v"]:.4f} (exact {2*np.pi*np.sin(g)/g:.4f})', flush=True)
