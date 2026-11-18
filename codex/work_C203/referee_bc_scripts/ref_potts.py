"""Referee C203: Theorem C checks with own code.
Potts ring H = -sum_j [E_j + F_{j,j+1}] on Q^L, E = sqrt(Q)|u><u|, F = sqrt(Q) sum_a |aa><aa|.
XXZ ring H_X(N; Delta, phi) on N = 2L spins, sector Sz = 0 (and others), as defined in C201 Section 4.1.
"""
import itertools, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import quad
from ref_common import digits, gauge_M


def potts_H(Q, L, sign=+1):
    D = digits(L, Q)
    N = Q ** L
    pw = Q ** (L - 1 - np.arange(L))
    sq = np.sqrt(Q)
    diag = np.zeros(N)
    rows, cols, vals = [], [], []
    for j in range(L):
        diag += -sign * sq * (D[:, j] == D[:, (j + 1) % L])
        diag += -sign / sq  # diagonal of E_j
        for dt in range(1, Q):
            new = (D[:, j] + dt) % Q
            rows.append(np.arange(N) + (new - D[:, j]) * pw[j])
            cols.append(np.arange(N))
            vals.append(np.full(N, -sign / sq))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return H + sp.diags(diag)


def xxz_sector(N, Delta, phi, m2):
    """H_X on N spins in sector with sum sigma^z = m2 (so S^z = m2/2)."""
    nup = (N + m2) // 2
    states = [s for s in itertools.combinations(range(N), nup)]
    pos = {s: k for k, s in enumerate(states)}
    dim = len(states)
    rows, cols, vals = [], [], []
    diag = np.zeros(dim)
    for k, s in enumerate(states):
        up = np.zeros(N, int); up[list(s)] = 1
        for i in range(N):
            j = (i + 1) % N
            if up[i] != up[j]:
                diag[k] += -Delta  # (Delta/2)(szsz - 1) = -Delta on antiparallel
                # hop: amplitude 1 in bulk; seam (N-1, 0): e^{i phi} sigma+_{N-1} sigma-_0
                new = up.copy(); new[i], new[j] = up[j], up[i]
                kk = pos[tuple(np.nonzero(new)[0])]
                if j != 0:
                    amp = 1.0
                else:
                    # i = N-1, j = 0. sigma+_{N-1} sigma-_0 raises N-1: from (N-1 down, 0 up) to (N-1 up, 0 down)
                    amp = np.exp(1j * phi) if (up[i] == 0 and up[j] == 1) else np.exp(-1j * phi)
                rows.append(kk); cols.append(k); vals.append(amp)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim), dtype=complex) + sp.diags(diag)
    return H


def low(H, k=1):
    if H.shape[0] <= 2500:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(spla.eigsh(H, k=k, which='SA', tol=1e-13)[0])


def eps0_xxz(Delta):
    """Bulk ground energy per site of H_X = 2 sum (SxSx+SySy+Delta SzSz) - N Delta/2."""
    if abs(Delta - 1) < 1e-14:
        e = 0.25 - np.log(2)
    else:
        g = np.arccos(Delta)
        a_ = np.pi - g
        f = lambda x: ((np.exp((a_ - np.pi) * x) - np.exp(-(a_ + np.pi) * x)) / (1 - np.exp(-2 * np.pi * x))
                       * 2 * np.exp(-g * x) / (1 + np.exp(-2 * g * x))) if x > 1e-9 else a_ / np.pi
        I = 2 * quad(f, 0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
        e = Delta / 4 - np.sin(g) / 2 * I
    return 2 * e - Delta / 2


def main():
    for Q, Lmax in [(2, 9), (3, 8), (4, 7)]:
        g = np.arccos(np.sqrt(Q) / 2)
        a = np.cos(g)
        sq = np.sqrt(Q)
        # bond and M
        u = np.ones(Q) / np.sqrt(Q)
        E = sq * np.outer(u, u)
        F = np.zeros((Q * Q, Q * Q))
        for x in range(Q):
            F[x * Q + x, x * Q + x] = sq
        I = np.eye(Q)
        C0 = F + 0.5 * (np.kron(E, I) + np.kron(I, E))
        ev = np.linalg.eigvalsh(C0)
        M = gauge_M(-C0, Q)
        e0 = 2 * eps0_xxz(a)
        print(f"\nQ={Q}: gamma={g:.6f} a={a:.6f}  spec(C0)=[{ev.min():.6f},{ev.max():.6f}]  M_SDP={M:.6f}  (sqrtQ+1)/2={(sq + 1) / 2:.6f}  e0=2eps0={e0:.10f}")
        for L in range(2, Lmax + 1):
            EP = low(potts_H(Q, L))[0]
            N = 2 * L
            row = {}
            for name, phi in [("2g", 2 * g), ("0", 0.0), ("pi", np.pi), ("g", g)]:
                row[name] = low(xxz_sector(N, a, phi, 0))[0]
            # lowest over other S^z sectors at twist 2g, for reference
            other = min(low(xxz_sector(N, a, 2 * g, 2 * mm))[0] for mm in range(1, L + 1))
            cas = L * (L * e0 - EP)
            print(f" L={L}: E_P={EP:.12f}  XXZ(Sz=0,phi=2g)={row['2g']:.12f}  diff={EP - row['2g']:.1e}  "
                  f"phi=0:{row['0'] - EP:+.4f} phi=pi:{row['pi'] - EP:+.4f} phi=g:{row['g'] - EP:+.4f}  "
                  f"min Sz!=0 @2g:{other - EP:+.4f}  L(Le0-E)={cas:.6f}  /M={cas / M:.6f}")


if __name__ == "__main__":
    main()
