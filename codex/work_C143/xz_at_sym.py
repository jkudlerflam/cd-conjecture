"""Symmetry-resolved ground states in the XZ family h = XZ + alpha ZX + beta YY.
For each ring length, the ground manifold is resolved by momentum and by the global parity P = prod_j Y_j;
the half-ring entanglement entropy is computed for every symmetry-resolved ground state (all states within
1e-8 of E0 that are simultaneous eigenvectors of T and P). We report min and max of S over them, and c
estimates from the minimum (a canonical choice that avoids cat superpositions of degenerate sectors)."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from ed import bond, ham_sparse, momentum_projector

def h_of(alpha, beta):
    J = np.zeros((3, 3)); J[0, 2] = 1.0; J[2, 0] = alpha; J[1, 1] = beta
    return bond(J)

def parityY(L):
    N = 1 << L; x = np.arange(N)
    y = (~x) & (N - 1)
    # Y|0> = i|1>, Y|1> = -i|0>  (|0> = up = bit 0): amplitude i^(#zeros) (-i)^(#ones)
    ones = np.array([bin(v).count("1") for v in range(N)])
    amp = (1j) ** (L - ones) * (-1j) ** ones
    return sp.csr_matrix((amp, (y, x)), shape=(N, N))

def entropy(psi, L):
    M = psi.reshape(2 ** (L // 2), -1)
    s = np.linalg.svd(M, compute_uv=False) ** 2; s = s[s > 1e-14]
    return float(-np.sum(s * np.log(s)))

def resolved(h, L, nev=4):
    H = ham_sparse(h, L); P = parityY(L)
    cands = []
    for m in range(L):
        Pk = momentum_projector(L, m)
        if Pk.shape[1] == 0: continue
        Hk = Pk.conj().T @ H @ Pk
        if Hk.shape[0] <= 400:
            w, v = np.linalg.eigh(Hk.toarray()); w = w[:nev]; v = v[:, :nev]
        else:
            w, v = spla.eigsh(Hk, k=nev, which="SA", tol=1e-12); o = np.argsort(w); w = w[o]; v = v[:, o]
        for i in range(len(w)):
            cands.append((w[i], m, Pk @ v[:, i]))
    E0 = min(c[0] for c in cands)
    G = [c for c in cands if c[0] < E0 + 1e-8]
    out = []
    # within each momentum, diagonalize parity on the degenerate vectors
    for m in sorted(set(c[1] for c in G)):
        V = np.array([c[2] for c in G if c[1] == m]).T
        V, _ = np.linalg.qr(V)
        Pm = V.conj().T @ (P @ V)
        pw, pv = np.linalg.eig(Pm)
        for j in range(len(pw)):
            psi = V @ pv[:, j]; psi /= np.linalg.norm(psi)
            out.append((m, complex(pw[j]), entropy(psi, L)))
    return E0, out

if __name__ == "__main__":
    pts = [tuple(map(float, p.split(","))) for p in sys.argv[1:]]
    Ls = [8, 10, 12, 14, 16]
    for (al, be) in pts:
        h = h_of(al, be); Smin = []; Smax = []; degs = []
        for L in Ls:
            E0, out = resolved(h, L)
            Ss = [o[2] for o in out]; Smin.append(min(Ss)); Smax.append(max(Ss)); degs.append(len(out))
        c8_12_16 = [3 * (Smin[i + 2] - Smin[i]) / np.log(Ls[i + 2] / Ls[i]) for i in range(len(Ls) - 2)]
        print(f"alpha={al:+.3f} beta={be:+.3f}: deg={degs} Smin={np.round(Smin,4).tolist()} Smax={np.round(Smax,4).tolist()} "
              f"c_est(L, L+4 from Smin)={np.round(c8_12_16,3).tolist()}", flush=True)
