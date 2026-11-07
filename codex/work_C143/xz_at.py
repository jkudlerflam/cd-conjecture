"""The two-parameter family through the XZ point,  h = X(x)Z + alpha Z(x)X + beta Y(x)Y.

Exact structure (PROOF section 4.2): with a_j = X_j Z_{j+1} and b_j = Z_j X_{j+1}, [a_i, b_j] = 0 for all i, j,
Y_j Y_{j+1} = a_j b_j, and each of {a_j}, {b_j} generates a Majorana chain (neighbours anticommute).
So H = sum_j [a_j + alpha b_j + beta a_j b_j]: two Majorana chains with velocity ratio alpha and an
energy-energy (Ashkin-Teller) coupling beta. XZ is alpha = beta = 0.

Diagnostics on rings L = 8..16: ground-state degeneracy, lowest gaps times L, soft count S_L,
and the half-ring entanglement entropy of the ground state (when unique), from which c is estimated by
the difference S(L/2) at successive L against (1/3) log L.
"""
import sys
import numpy as np
import scipy.sparse.linalg as spla
from ed import bond, ham_sparse, momentum_projector, spectrum, all_levels

def h_of(alpha, beta):
    J = np.zeros((3, 3)); J[0, 2] = 1.0; J[2, 0] = alpha; J[1, 1] = beta
    return bond(J)

def ground_entropy(h, L):
    """Lowest state over all momentum sectors; returns (E0, degeneracy within 1e-9, S(L/2)) for the lowest vector."""
    from ed import ham_sparse
    H = ham_sparse(h, L)
    best = None
    for m in range(L):
        P = momentum_projector(L, m)
        if P.shape[1] == 0: continue
        Hk = (P.conj().T @ H @ P)
        if Hk.shape[0] < 300:
            w, v = np.linalg.eigh(Hk.toarray()); w = w[:3]; v = v[:, :3]
        else:
            w, v = spla.eigsh(Hk, k=3, which="SA", tol=1e-12)
            o = np.argsort(w); w = w[o]; v = v[:, o]
        for i in range(len(w)):
            if best is None or w[i] < best[0] - 1e-12:
                best = (w[i], P @ v[:, i], [w])
            best[2].append(w)
    E0 = best[0]
    allw = np.concatenate(best[2])
    deg = int(np.sum(np.abs(allw - E0) < 1e-9))
    psi = best[1] / np.linalg.norm(best[1])
    M = psi.reshape(2 ** (L // 2), 2 ** (L - L // 2))
    s = np.linalg.svd(M, compute_uv=False) ** 2
    s = s[s > 1e-14]
    return E0, deg, float(-np.sum(s * np.log(s)))

if __name__ == "__main__":
    pts = [(1.0, 0.0), (1.0, 0.5), (0.3, 0.0), (0.0, 0.3), (0.0, 0.6), (0.0, 1.0), (-0.19, 0.3), (-0.38, 0.6), (-0.64, 1.0)]
    if len(sys.argv) > 1:
        pts = [tuple(map(float, p.split(","))) for p in sys.argv[1:]]
    Ls = [8, 10, 12, 14, 16]
    for (al, be) in pts:
        h = h_of(al, be)
        res = [ground_entropy(h, L) for L in Ls]
        Ss = [r[2] for r in res]
        cest = [3 * (Ss[i + 1] - Ss[i]) / np.log(Ls[i + 1] / Ls[i]) for i in range(len(Ls) - 1)]
        print(f"alpha={al:+.3f} beta={be:+.3f}: deg={[r[1] for r in res]}, S(L/2)={np.round(Ss,4).tolist()}, c_est(pairs)={np.round(cest,3).tolist()}", flush=True)
