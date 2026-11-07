"""Checks at the XZ (uniform compass) point h = X (x) Z.

1. Ground-state degeneracy and low spectrum on rings, against the Majorana prediction:
   odd Majoranas form a critical chain, even Majoranas are free zero modes.
2. Local conserved operators Q_j = Z_j X_{j+1}.
3. Rank-one family h = (a.s)(x)(b.s) with angle theta between a and b: gap and degeneracy.
4. Tangent space: which of the 11 directions (J 9 + B 3, minus scale) keep the degeneracy at first order.
"""
import numpy as np
from ed import bond, spectrum, all_levels, ham_sparse, SX, SY, SZ, I2

np.set_printoptions(precision=6, suppress=True, linewidth=140)


def op_on(L, ops):
    """ops: dict site -> 2x2."""
    M = np.array([[1.0 + 0j]])
    for i in range(L):
        M = np.kron(M, ops.get(i, I2))
    return M


def levels(h, L, tol=1e-8):
    w, k = all_levels(spectrum(h, L))
    o = np.argsort(w)
    return w[o], k[o]


print("== 1. XZ point degeneracies")
J = np.zeros((3, 3)); J[0, 2] = 1.0
h = bond(J)
for L in [4, 6, 8, 10, 12]:
    w, k = levels(h, L)
    e0 = w[0]
    deg = np.sum(np.abs(w - e0) < 1e-8)
    gaps = np.unique(np.round(w - e0, 8))[:5]
    # Majorana prediction: H = sum_j i g_{2j-1} g_{2j+1} (odd chain, L Majoranas, periodic or antiperiodic by parity)
    print(f"L={L}: E0={e0:.10f}, E0/L={e0/L:.6f}, ground degeneracy={deg}, 2^(L/2)={2**(L//2)}, first distinct gaps={gaps[1:]}")
print("expected E0/L -> -4/pi =", -4 / np.pi, "(critical Majorana chain with hopping 2: E0 = -sum_{k} |2 sin k| over half BZ...)")

print("== 2. Q_j = Z_j X_{j+1} commutes with H (L=6)")
L = 6
H = ham_sparse(h, L).toarray()
for j in range(2):
    Q = op_on(L, {j: SZ, (j + 1) % L: SX})
    print(" ||[H,Q_%d]|| =" % j, np.linalg.norm(H @ Q - Q @ H))

print("== 3. rank-one family (a.s)(x)(b.s), a = x, b = cos(t) x + sin(t) z")
for t in [0.0, 0.5, 1.0, 1.3, 1.5, 1.55, np.pi / 2]:
    J = np.zeros((3, 3)); J[0, 0] = np.cos(t); J[0, 2] = np.sin(t)
    h = bond(J)
    row = []
    for L in [8, 10, 12]:
        w, k = levels(h, L)
        e0 = w[0]
        deg = np.sum(np.abs(w - e0) < 1e-8)
        n_low = np.sum(w - e0 < 0.05)
        row.append((deg, n_low))
    w, k = levels(h, 12)
    print(f" theta={t:.3f}: (exact deg, #states within 0.05) at L=8,10,12: {row}")

print("== 4. tangent directions at XZ (first-order splitting of the ground manifold, L=8)")
L = 8
J0 = np.zeros((3, 3)); J0[0, 2] = 1.0
h0 = bond(J0)
H0 = ham_sparse(h0, L).toarray()
w, V = np.linalg.eigh(H0)
G = V[:, np.abs(w - w[0]) < 1e-8]
print(" ground degeneracy", G.shape[1])
names = []
dirs = []
lab = "xyz"
for a in range(3):
    for b in range(3):
        if (a, b) == (0, 2):
            continue
        Jd = np.zeros((3, 3)); Jd[a, b] = 1.0
        names.append(lab[a] + lab[b]); dirs.append(bond(Jd))
for a in range(3):
    Bd = np.zeros(3); Bd[a] = 1.0
    names.append("B" + lab[a]); dirs.append(bond(np.zeros((3, 3)), Bd))
for nm, hd in zip(names, dirs):
    Hd = ham_sparse(hd, L).toarray()
    P = G.conj().T @ Hd @ G
    e = np.linalg.eigvalsh(P)
    print(f"  {nm}: first-order splitting width {e[-1]-e[0]:.4f}, mean shift {e.mean():.4f}")
