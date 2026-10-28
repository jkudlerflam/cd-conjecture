"""Referee check of the leading hopping formula, PROOFS eq. (8).

(a) The exact g^2 hopping of the microscopic bond (ref_bond.py, physical basis)
    equals -2 g^2 (q^(2n+2)/S^2) [Gfull(w) + i eta Gfull(wv)] at finite D, where
    Gfull is the two-token Gram of the full (2n+2)-vertex scalar folded bus,
    left token entering at B, right token at A.
(b) D^(2n+4) Gfull(w) -> Gamma(w) as D -> infinity, in 300-digit decimals,
    with Gamma(w) from the independent eigenmode formula in ref_pencil.py.
"""
import json, sys
import numpy as np
from decimal import Decimal as Dm, getcontext
import ref_bond, ref_pencil as rp


def scalar_graph(n, w, D, a, num=float):
    m = 2 * n + 2  # order: A, d_0..d_{n-1}, x_0..x_{n-1}, B
    K = [[num(0)] * m for _ in range(m)]
    for i in range(m):
        K[i][i] = num(a) + (num(D) if (i == 0 or i > n) else num(0))
    edges = [(0, 1, 1), (1, n + 1, 1), (2 * n, 2 * n + 1, 1)]
    edges += [(1 + j, 2 + j, w) for j in range(n - 1)] + [(n + 1 + j, n + 2 + j, 1) for j in range(n - 1)]
    for i, j, t in edges:
        K[i][j] = K[j][i] = -num(t) if not isinstance(t, num) else -t
    return K


def gfull_np(n, w, D, a):
    K = np.array(scalar_graph(n, w, D, a))
    m = len(K); I = np.eye(m)
    KK = np.kron(K, I) + np.kron(I, K)
    R = np.linalg.inv(KK)
    e = np.zeros(m * m); e[(m - 1) * m + 0] = 1  # left token at B, right at A
    z = R @ e
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            Pi = np.kron(np.diag(np.eye(m)[1 + i]), I); Pj = np.kron(I, np.diag(np.eye(m)[1 + j]))
            G[i, j] = 2 * z @ Pi @ R @ Pj @ z
    return G


def gfull_dec(n, w, D, a):
    """Same quantity in Decimal, Gaussian elimination."""
    m = 2 * n + 2
    K = scalar_graph(n, Dm(w), Dm(D), Dm(a), num=Dm)
    M = m * m
    def kk(p, q):
        i1, i2 = divmod(p, m); j1, j2 = divmod(q, m)
        return (K[i1][j1] if i2 == j2 else 0) + (K[i2][j2] if i1 == j1 else 0)
    Amat = [[Dm(kk(p, q)) for q in range(M)] for p in range(M)]
    # LU without pivoting (diagonally dominant)
    L = [r[:] for r in Amat]
    for c in range(M):
        piv = L[c][c]
        for r in range(c + 1, M):
            if L[r][c] != 0:
                f = L[r][c] / piv; L[r][c] = f
                Lr = L[r]; Lc = L[c]
                for k in range(c + 1, M):
                    if Lc[k] != 0: Lr[k] -= f * Lc[k]
    def solve(b):
        y = b[:]
        for r in range(M):
            y[r] -= sum((L[r][k] * y[k] for k in range(r) if L[r][k] != 0), Dm(0))
        x = [Dm(0)] * M
        for r in range(M - 1, -1, -1):
            x[r] = (y[r] - sum((L[r][k] * x[k] for k in range(r + 1, M) if L[r][k] != 0), Dm(0))) / L[r][r]
        return x
    e = [Dm(0)] * M; e[(m - 1) * m] = Dm(1)
    z = solve(e)
    G = [[Dm(0)] * n for _ in range(n)]
    for j in range(n):
        b = [z[p] if (p % m) == 1 + j else Dm(0) for p in range(M)]
        y = solve(b)
        for i in range(n):
            G[i][j] = 2 * sum((z[p] * y[p] for p in range(M) if (p // m) == 1 + i), Dm(0))
    return G


if __name__ == '__main__':
    out = {}
    # (a) finite-D exactness of the g^2 formula
    a, q, w, v, eta = 3.5, .8, .8, .65, .6
    rows = []
    for n in (1, 2, 3):
        for D in (.2, 2.):
            g = 1e-4
            K = ref_bond.bond(n, w=w, v=v, D=D, g=g, q=q, a=a, eta=eta)
            N = 2 ** n
            T = K[np.ix_([(1 << j) * N for j in range(n)], [1 << j for j in range(n)])] / g ** 2
            S = sum(q ** (2 * l) for l in range(n + 2))
            pred = -2 * q ** (2 * n + 2) / S ** 2 * (gfull_np(n, w, D, a) + 1j * eta * gfull_np(n, w * v, D, a))
            rows.append(dict(n=n, D=D, rel_err=float(np.abs(T - pred).max() / np.abs(pred).max())))
            print(rows[-1], flush=True)
    out['finite_D_g2_formula'] = rows
    # (b) D -> infinity, exact-ish decimals
    getcontext().prec = 300
    rows = []
    for n in (1, 2, 3, 4):
        Gam = rp.gram(n, Dm(1) / 4)
        prev = None
        for D in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 6):
            G = gfull_dec(n, Dm(1) / 4, D, 4)
            sc = Dm(D) ** (2 * n + 4)
            err = max(abs(sc * G[i][j] - Gam[i][j]) / abs(Gam[i][j]) for i in range(n) for j in range(n))
            rows.append(dict(n=n, D=D, max_entry_rel_err=float(err), times_D=float(err * D)))
            print(rows[-1], flush=True)
    out['large_D_limit'] = rows
    json.dump(out, open('ref_largeD.json', 'w'), indent=1)
