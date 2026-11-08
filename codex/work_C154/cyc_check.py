"""C154: exact checks of the cyclic stratum h = a XZ + b YX + g ZY (Claude Opus agent, 2026-11-07).
(1) J = D P with D = diag(a,b,g), P the cyclic permutation (P e_z = e_x, P e_x = e_y, P e_y = e_z).
(2) sigma_j = P^{-j} tau_j maps the uniform ring (L = 0 mod 3) to the XYZ ring with bond tensor D_j = P^j D P^-j.
(3) Twisted Klein generators U_n = prod_j (n_j . sigma_j), n_j = P^{-j} e_n, commute with H (L = 0 mod 3).
(4) Kramers time reversal Theta = prod (i Y) K commutes with H (B = 0)."""
import numpy as np
from functools import reduce
I2 = np.eye(2); SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]]); SZ = np.diag([1., -1]).astype(complex)
PA = [SX, SY, SZ]
P = np.zeros((3, 3)); P[0, 2] = 1; P[1, 0] = 1; P[2, 1] = 1   # P_xz = P_yx = P_zy = 1

def site_op(op, j, L):
    return reduce(np.kron, [op if k == j else I2 for k in range(L)])

def ring(Jlist, L):
    H = 0
    for j in range(L):
        J = Jlist[j % len(Jlist)]
        for a in range(3):
            for b in range(3):
                if J[a, b] != 0:
                    H = H + J[a, b] * site_op(PA[a], j, L) @ site_op(PA[b], (j + 1) % L, L)
    return H

rng = np.random.default_rng(1)
for (a, b, g) in [(1, 0.1, 0), (1, 0, 0.1), (1, 0.37, -0.52), tuple(rng.normal(size=3))]:
    D = np.diag([a, b, g]); J = D @ P
    assert abs(J[0, 2] - a) < 1e-15 and abs(J[1, 0] - b) < 1e-15 and abs(J[2, 1] - g) < 1e-15
    Pinv = P.T
    Dj = [np.linalg.matrix_power(P, j) @ D @ np.linalg.matrix_power(Pinv, j) for j in range(3)]
    for L in (6, 9):
        Hs = ring([J], L); Ht = ring(Dj, L)
        es, et = np.linalg.eigvalsh(Hs), np.linalg.eigvalsh(Ht)
        # twisted Klein and Kramers
        res = []
        for n in range(3):
            U = reduce(lambda x, y: x @ y, [site_op(sum((np.linalg.matrix_power(Pinv, j) @ np.eye(3)[n])[c] * PA[c] for c in range(3)), j, L) for j in range(L)])
            res.append(np.abs(U @ Hs - Hs @ U).max())
        Th = reduce(lambda x, y: x @ y, [site_op(1j * SY, j, L) for j in range(L)])
        kr = np.abs(Th @ Hs.conj() - Hs @ Th).max()
        print(f"(a,b,g)=({a:+.3f},{b:+.3f},{g:+.3f}) L={L}: diag(D_j)={[np.round(np.diag(d),3).tolist() for d in Dj]} "
              f"spec diff={np.abs(es-et).max():.1e}  [U_n,H]={['%.0e'%r for r in res]}  [Theta,H]={kr:.0e}")
