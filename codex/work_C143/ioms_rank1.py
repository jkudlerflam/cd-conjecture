"""Local integrals of motion (strictly local operators Q with [H_L, Q] = 0) for rank-one bonds.

PROOF section 2.6 shows: a nontrivial strictly local IOM of any range forces rank J <= 1, and for range two the
only cases are the classical bonds and the XZ point. Here we test ranges r = 2, 3, 4 numerically on the
four-parameter rank-one family J = m n^T, m = x, n = cos(t) x + sin(t) z, field B = (bx, by, bz):
for each bond the dimension of the kernel of Q -> [H_L, Q] on operators supported on sites 0..r-1
(identity excluded) is the number of independent IOMs of range <= r there.
Usage: python3 ioms_rank1.py r L NSTARTS SEED
"""
import sys
import numpy as np
from scipy.optimize import minimize
from ed import bond, ham_sparse, PAULI, I2

P4 = [I2] + PAULI


def basis_ops(r, L):
    ops = []
    import itertools
    for idx in itertools.product(range(4), repeat=r):
        if all(i == 0 for i in idx):
            continue
        M = np.array([[1.0 + 0j]])
        for i in idx:
            M = np.kron(M, P4[i])
        ops.append((idx, np.kron(M, np.eye(2 ** (L - r)))))
    return ops


def bond_of(p):
    t = p[0]
    q = np.asarray(p[1:], float); nq = np.linalg.norm(q)
    bx, by, bz = (3.0 * np.tanh(nq) * q / nq) if nq > 0 else (0.0, 0.0, 0.0)   # |B| < 3 |J|
    J = np.zeros((3, 3)); J[0, 0] = np.cos(t); J[0, 2] = np.sin(t)
    return bond(J, (bx, by, bz))


def svals(p, ops, L):
    H = ham_sparse(bond_of(p), L).toarray()
    cols = []
    for _, Q in ops:
        C = H @ Q - Q @ H
        cols.append(np.concatenate([C.real.ravel(), C.imag.ravel()]))
    s = np.linalg.svd(np.array(cols).T, compute_uv=False)
    return s / np.linalg.norm(H)


if __name__ == "__main__":
    r, L, n, seed = map(int, sys.argv[1:5])
    ops = basis_ops(r, L)
    print("known points: XZ (t=pi/2,B=0) smallest s.v.:", svals([np.pi / 2, 0, 0, 0], ops, L)[-3:])
    print("classical (t=0, B along x):", svals([0, 0.7, 0, 0], ops, L)[-3:])
    rng = np.random.default_rng(seed)
    for i in range(n):
        p0 = np.array([rng.uniform(0, np.pi), *rng.normal(size=3)])
        r_ = minimize(lambda p: svals(p, ops, L)[-1], p0, method="Nelder-Mead",
                      options=dict(maxfev=800, xatol=1e-7, fatol=1e-10))
        t = r_.x[0] % np.pi
        q = r_.x[1:]; nq = np.linalg.norm(q); Bv = 3 * np.tanh(nq) * q / nq
        print(f"start {i}: min sigma {r_.fun:.2e} at theta/pi={t/np.pi:.4f}, B={np.round(Bv, 4)}", flush=True)
