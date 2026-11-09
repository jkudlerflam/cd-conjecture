#!/usr/bin/env python3
"""C170: Neel plateau edge of the C165 parent (AF XXZ, Pauli units: H = sum mu.mu_perp + Dz mu^z mu^z + hz sum mu^z).
Exact spin gap of the S=1/2 XXZ chain (J = 1, Delta = cosh gamma > 1), des Cloizeaux-Gaudin:
  gap = (pi sinh gamma / gamma) sum_n 1/cosh((2n+1) pi^2/(2 gamma)).
In Pauli units the coupling is J = 4 and the Zeeman term is 2 hz S^z, so hz_c1 = 2 gap(J=1).
Also the saturation field hz_sat = 2 (1 + Delta). ED check on rings in the mu frame with S^z sectors."""
import numpy as np
from itertools import combinations

def gap(D):
    g = np.arccosh(D)
    n = np.arange(-50, 50)
    return np.pi * np.sinh(g) / g * np.sum(1 / np.cosh((2 * n + 1) * np.pi ** 2 / (2 * g)))

def ed_sector(L, D, M):
    # mu frame, Pauli units, ring; sector with M up spins
    basis = [sum(1 << i for i in c) for c in combinations(range(L), M)]
    idx = {b: i for i, b in enumerate(basis)}
    import scipy.sparse as sp, scipy.sparse.linalg as sla
    rows, cols, vals = [], [], []
    for i, b in enumerate(basis):
        diag = 0.0
        for j in range(L):
            k = (j + 1) % L
            sj = 1 if b >> j & 1 else -1; sk = 1 if b >> k & 1 else -1
            diag += D * sj * sk
            if sj != sk:
                b2 = b ^ (1 << j) ^ (1 << k)
                rows.append(idx[b2]); cols.append(i); vals.append(2.0)   # XX + YY = 2 (s+s- + s-s+)
        rows.append(i); cols.append(i); vals.append(diag)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)))
    return sla.eigsh(H, k=1, which='SA', return_eigenvectors=False)[0]

if __name__ == '__main__':
    for D in (1.1, 1.25, 1.3, 1.4, 1.6, 2.0):
        print(f'Delta={D}: exact hz_c1 = {2*gap(D):.4f}   hz_sat = {2*(1+D):.2f}   classical spin-flop hz = {2*np.sqrt(D*D-1):.3f}')
    for L in (12, 16, 20):
        for D in (2.0, 1.6):
            e0 = ed_sector(L, D, L // 2); e1 = ed_sector(L, D, L // 2 + 1)
            print(f'ED ring L={L} Delta={D}: E(M=1)-E(0) = {e1-e0:.4f} (field term -2 hz per flipped spin -> hz_c1(L) = {(e1-e0)/2:.4f}), e0/L = {e0/L:.6f}')
