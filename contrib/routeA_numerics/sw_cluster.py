#!/usr/bin/env python3
"""Numerically exact Schrieffer-Wolff (direct-rotation / des Cloizeaux) effective Hamiltonians of
open clusters of n identical route-A cells, and their linked-cluster decomposition.

Cluster Hamiltonian (open, n cells):  H_n = sum_j H0_j + lam sum_{j<n} sum_a FR_a(j) FL_a(j+1).
Low band: the N^n lowest eigenvectors Psi (energies E).  With B = (G^{(x)n})^dagger Psi and the polar
factor W = (B B^dagger)^{-1/2} B, the effective logical Hamiltonian is H_eff = W E W^dagger.  This is
the Bravyi-DiVincenzo-Loss Schrieffer-Wolff Hamiltonian, which is linked-cluster additive, so
    h2 = H_eff(2),     h3 = H_eff(3) - h2 (x) 1 - 1 (x) h2,     ...
are the bond and three-cell terms of the logical chain  H_log = sum_j h2_{j,j+1} + sum_j h3_{j,j+1,j+2}.

usage:  python3 sw_cluster.py [--cell c8] [--lams 0.01 0.02 ...] [--n3] [--out data/heff_c8.npz]
Written by Claude (Opus 5.5) for contrib/routeA_numerics, 2026-10-27.
"""
import os, sys, argparse, time, json
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cell import load_c8, swap, hermitian_basis  # noqa: E402


def cluster_op(cell, n, lam):
    e, FR, FL, D = cell['e'], cell['FR'], cell['FL'], cell['D']
    H0 = cell['H0']
    diag0 = np.allclose(H0, np.diag(np.diag(H0)))
    dim = D ** n

    def mv(x):
        x = x.reshape([D] * n + [-1])
        out = np.zeros_like(x)
        for j in range(n):
            if diag0:
                sh = [1] * (n + 1); sh[j] = D
                out += e.reshape(sh) * x
            else:
                out += np.moveaxis(np.tensordot(H0, x, ([1], [j])), 0, j)
        for j in range(n - 1):
            for a in range(3):
                t = np.moveaxis(np.tensordot(FR[a], x, ([1], [j])), 0, j)
                t = np.moveaxis(np.tensordot(FL[a], t, ([1], [j + 1])), 0, j + 1)
                out += lam * t
        return out.reshape(dim, -1)
    return LinearOperator((dim, dim), matvec=lambda v: mv(v)[:, 0], matmat=mv, dtype=complex), dim


def low_band(cell, n, lam, tol=1e-13, extra=4):
    N, D = cell['N'], cell['D']
    k = N ** n
    op, dim = cluster_op(cell, n, lam)
    if dim <= 3000:
        Hd = op.matmat(np.eye(dim, dtype=complex))
        w, U = linalg.eigh((Hd + Hd.conj().T) / 2)
        return w[:k], U[:, :k], w[k:k + extra]
    # start in the code space (index i < N on every cell in the eigenbasis)
    rng = np.random.default_rng(1)
    v0 = np.zeros([D] * n, complex)
    v0[tuple([slice(0, N)] * n)] = rng.normal(size=[N] * n)
    w, U = eigsh(op, k=k + extra, which='SA', tol=tol, v0=v0.ravel(), ncv=max(3 * (k + extra), 60), maxiter=20000)
    idx = np.argsort(w)
    w, U = w[idx], U[:, idx]
    # Rayleigh-Ritz inside the band to clean near-degenerate clusters
    Q, _ = np.linalg.qr(U[:, :k])
    Hq = Q.conj().T @ op.matmat(Q)
    wq, uq = linalg.eigh((Hq + Hq.conj().T) / 2)
    return wq, Q @ uq, w[k:]


def heff(cell, n, lam, **kw):
    N, D = cell['N'], cell['D']
    E, Psi, above = low_band(cell, n, lam, **kw)
    # code space: first N eigenbasis states on each cell
    Psi_t = Psi.reshape([D] * n + [-1])[tuple([slice(0, N)] * n)].reshape(N ** n, -1)
    B = Psi_t                                   # (logical, eigen) = (G^n)^dagger Psi
    u, s, vh = linalg.svd(B)
    W = u @ vh
    H = (W * E) @ W.conj().T
    H = (H + H.conj().T) / 2
    info = dict(band_top=float(E.max()), next=float(above.min()) if len(above) else None,
                overlap_min_sv=float(s.min()))
    return H, info


def decompose2(h2, N):
    """h2 = a 1 + b P + (on-site) + (connected non-swap).  Returns coefficients and norms."""
    P = swap(N)
    I = np.eye(N * N)
    T = hermitian_basis(N)                       # T[0] = 1/sqrt N
    coef = np.array([[np.trace(np.kron(A, B) @ h2).real for B in T] for A in T])   # h2 = sum c_mn T_m T_n
    # SU(N)-invariant part: a 1 + b P,  P = 1/N + sum_{mu>0} T_mu (x) T_mu^*  ... use projection
    Pc = P - np.trace(P) / N ** 2 * I
    b = np.trace(Pc @ h2).real / np.trace(Pc @ Pc).real
    a = (np.trace(h2).real - b * np.trace(P).real) / N ** 2
    R = h2 - a * I - b * P
    onsite_L = sum(coef[m, 0] * np.kron(T[m], T[0]) for m in range(1, N * N))
    onsite_R = sum(coef[0, m] * np.kron(T[0], T[m]) for m in range(1, N * N))
    conn = R - onsite_L - onsite_R
    return dict(a=a, b=b, R=R, coef=coef, onsite_L=onsite_L, onsite_R=onsite_R, conn=conn,
                nR=np.linalg.norm(R, 2), n_on=np.linalg.norm(onsite_L + onsite_R, 2), n_conn=np.linalg.norm(conn, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lams', type=float, nargs='+', default=[0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2])
    ap.add_argument('--n3', action='store_true')
    ap.add_argument('--out', default=os.path.join(HERE, 'data', 'heff_c8.npz'))
    a = ap.parse_args()
    cell = load_c8()
    N = cell['N']
    res = {}
    log = []
    for lam in a.lams:
        t0 = time.time()
        h2, i2 = heff(cell, 2, lam)
        d = decompose2(h2, N)
        rec = dict(lam=lam, a=d['a'], b=d['b'], b_over_lam2=d['b'] / lam ** 2, nR=d['nR'], n_on=d['n_on'],
                   n_conn=d['n_conn'], R_over_b=d['nR'] / d['b'], info2=i2)
        res[f'h2_{lam}'] = h2
        if a.n3:
            H3, i3 = heff(cell, 3, lam)
            h3 = H3 - np.kron(h2, np.eye(N)) - np.kron(np.eye(N), h2)
            res[f'h3_{lam}'] = h3
            rec.update(n_h3=float(np.linalg.norm(h3 - np.trace(h3) / N ** 3 * np.eye(N ** 3), 2)), info3=i3)
            rec['h3_over_b'] = rec['n_h3'] / d['b']
        rec['t'] = time.time() - t0
        log.append(rec)
        print(json.dumps(rec), flush=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    np.savez(a.out, lams=np.array(a.lams), **res)
    json.dump(log, open(a.out.replace('.npz', '.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
