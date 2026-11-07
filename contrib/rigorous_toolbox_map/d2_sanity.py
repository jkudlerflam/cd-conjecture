"""Small sanity check for the d = 2 direction (DRAFT_C131).

Computes Q(L) = L (E_{2L} - 2 E_L) / M for spin-1/2 nearest-neighbour rings,
even L = 2, 4, 6, 8 (so 2L <= 16), by sparse exact diagonalization.

Families:
  (a) the U(1) family h = Sx Sx + Sy Sy + Delta Sz Sz + D (Sx Sy - Sy Sx)
      - (B/2)(Sz 1 + 1 Sz), on a grid of (Delta, D, B);
  (b) random Hermitian bonds (Gaussian), normalized.
M is the gauge-optimized half spread, computed by convex minimization over
o (traceless Hermitian 2x2) of (lambda_max - lambda_min)/2 of h + o x 1 - 1 x o.

numpy and scipy only. Runtime: a few minutes on a laptop.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize

sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
id2 = np.eye(2, dtype=complex)
paulis = [sx, sy, sz]


def gauge_M(h):
    def spread(p):
        o = sum(p[i] * paulis[i] for i in range(3))
        hh = h + np.kron(o, id2) - np.kron(id2, o)
        ev = np.linalg.eigvalsh(hh)
        return (ev[-1] - ev[0]) / 2
    best = None
    for s in range(6):
        x0 = np.random.default_rng(s).normal(size=3) * (0 if s == 0 else 1)
        r = minimize(spread, x0, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 4000})
        if best is None or r.fun < best:
            best = r.fun
    return best


def ring_H(h, L):
    """Sparse ring Hamiltonian sum_x h_{x,x+1} on (C^2)^L."""
    dim = 2 ** L
    hs = sp.csr_matrix(h)
    H = sp.csr_matrix((dim, dim), dtype=complex)
    for x in range(L - 1):
        H = H + sp.kron(sp.kron(sp.identity(2 ** x), hs), sp.identity(2 ** (L - x - 2)))
    # closing bond (L-1, 0): build via permutation of the bond operator
    h4 = h.reshape(2, 2, 2, 2)  # out_a out_b in_a in_b with a = site L-1, b = site 0
    rows, cols, vals = [], [], []
    for s in range(dim):
        b0 = (s >> (L - 1)) & 1      # site 0 is the most significant bit
        a0 = s & 1                   # site L-1 is the least significant bit
        base = s & ~(1 << (L - 1)) & ~1
        for a1 in range(2):
            for b1 in range(2):
                v = h4[a1, b1, a0, b0]
                if v != 0:
                    t = base | (b1 << (L - 1)) | a1
                    rows.append(t); cols.append(s); vals.append(v)
    H = H + sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim))
    return H


def E0(h, L):
    H = ring_H(h, L)
    if 2 ** L <= 256:
        return np.linalg.eigvalsh(H.toarray())[0]
    return spla.eigsh(H, k=1, which="SA", tol=1e-10)[0][0]


def Q_profile(h, Ls=(2, 4, 6, 8)):
    M = gauge_M(h)
    E = {}
    for L in sorted(set(list(Ls) + [2 * l for l in Ls])):
        E[L] = E0(h, L)
    return M, [L * (E[2 * L] - 2 * E[L]) / M for L in Ls]


if __name__ == "__main__":
    out = []
    print("(a) U(1) family: Delta, D, B -> M, Q(2), Q(4), Q(6), Q(8)")
    for Delta in [-0.95, -0.5, 0.0, 0.5, 1.0, 1.5]:
        for D in [0.0, 0.3]:
            for B in [0.0, 0.4, 0.8]:
                h = (np.kron(sx, sx) + np.kron(sy, sy) + Delta * np.kron(sz, sz)
                     + D * (np.kron(sx, sy) - np.kron(sy, sx))
                     - (B / 2) * (np.kron(sz, id2) + np.kron(id2, sz)))
                M, Q = Q_profile(h)
                line = f"{Delta:6.2f} {D:4.1f} {B:4.1f}  M={M:.4f}  " + " ".join(f"{q:8.4f}" for q in Q)
                print(line); out.append(line)
    print("(b) random bonds: seed -> M, Q(2), Q(4), Q(6), Q(8)")
    for seed in range(12):
        g = np.random.default_rng(100 + seed).normal(size=(4, 4)) + 1j * np.random.default_rng(200 + seed).normal(size=(4, 4))
        h = (g + g.conj().T) / 2
        M, Q = Q_profile(h)
        line = f"seed {seed:2d}  M={M:.4f}  " + " ".join(f"{q:8.4f}" for q in Q)
        print(line); out.append(line)
    with open("d2_sanity_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
