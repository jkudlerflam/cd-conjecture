"""The effective B chain on the tuned Ashkin-Teller cone of XZ.

Second-order Schrieffer-Wolff on the degenerate manifold (PROOF.md section 3):
    H_B = t sum_j b_j  -  beta^2 sum_{j != l} K(l - j) b_j b_l  + const,     t = alpha + beta <a>_L,
with K the kernel of kernel.py (ring of the same length L). Nearest-neighbour pairs cancel exactly
(K real symmetric, b_j b_{j+1} = - b_{j+1} b_j), so in units of beta^2, with lam = t / beta^2,
    h_B(lam) = lam sum_j b_j - 2 sum_j sum_{r=2}^{L/2} w_r K(r) b_j b_{j+r},   w_r = 1/2 at r = L/2 else 1.
Representation on n = L/2 spins (after an on-site Hadamard): b_{2i-1} = Z_i, b_{2i} = X_i X_{i+1},
closing bond b_L = c X_n X_1 with c = +-1, and the parity P = prod Z_i conserved.
Pauli strings P(x, z) = prod X^x Z^z; P(x1,z1)P(x2,z2) = (-1)^{|z1 & x2|} P(x1^x2, z1^z2).
"""
import numpy as np
import numba as nb
import scipy.sparse.linalg as spla
from kernel import kernel_row


def b_string(j, n, c):
    """b_j (j = 0..L-1, b_0 = b_1 in 1-based notation) as (x, z, coeff)."""
    if j % 2 == 0:  # 1-based odd: Z_i, i = j//2
        return 0, 1 << (j // 2), 1.0
    i = j // 2
    k = (i + 1) % n
    return (1 << i) | (1 << k), 0, (c if k == 0 else 1.0)


def mult(p, q):
    x1, z1, c1 = p; x2, z2, c2 = q
    s = -1.0 if bin(z1 & x2).count("1") % 2 else 1.0
    return x1 ^ x2, z1 ^ z2, c1 * c2 * s


def terms(L, lam, K, c=1):
    n = L // 2
    out = {}

    def add(t, w):
        key = (t[0], t[1])
        out[key] = out.get(key, 0.0) + w * t[2]

    for j in range(L):
        add(b_string(j, n, c), lam)
    for j in range(L):
        for r in range(2, L // 2 + 1):
            w = 0.5 if (L % 2 == 0 and r == L // 2) else 1.0
            l = (j + r) % L
            p = mult(b_string(j, n, c), b_string(l, n, c))
            add(p, -2.0 * w * K[r])
    xs = np.array([k[0] for k in out], dtype=np.int64)
    zs = np.array([k[1] for k in out], dtype=np.int64)
    cs = np.array([out[k] for k in out], dtype=np.float64)
    keep = np.abs(cs) > 1e-14
    return xs[keep], zs[keep], cs[keep]


@nb.njit(cache=True)
def _popc(v):
    c = 0
    while v:
        v &= v - 1
        c += 1
    return c


@nb.njit(parallel=True, cache=True)
def _mv(states, index, xs, zs, cs, v, out):
    D = states.shape[0]
    for a in nb.prange(D):
        y = states[a]
        acc = 0.0
        for t in range(xs.shape[0]):
            s = y ^ xs[t]
            sg = 1.0 - 2.0 * (_popc(zs[t] & s) & 1)
            acc += cs[t] * sg * v[index[s]]
        out[a] = acc


def sector(n, parity):
    s = np.arange(1 << n, dtype=np.int64)
    pc = np.array([bin(v).count("1") & 1 for v in range(1 << n)])
    st = s[pc == parity]
    index = -np.ones(1 << n, dtype=np.int64)
    index[st] = np.arange(len(st))
    return st, index


def spectrum(L, lam, K, c, parity, k=6):
    n = L // 2
    xs, zs, cs = terms(L, lam, K, c)
    st, index = sector(n, parity)
    D = len(st)

    def mv(v):
        out = np.empty(D)
        _mv(st, index, xs, zs, cs, np.ascontiguousarray(v.real), out)
        return out

    if D <= 1024:
        H = np.zeros((D, D))
        for a in range(D):
            e = np.zeros(D); e[a] = 1
            H[:, a] = mv(e)
        assert np.allclose(H, H.T)
        return np.linalg.eigvalsh(H)[:k]
    op = spla.LinearOperator((D, D), matvec=mv, dtype=float)
    w = spla.eigsh(op, k=k, which="SA", tol=1e-9, return_eigenvectors=False, ncv=max(30, 3 * k))
    return np.sort(w)


def ground_vec(L, lam, K, c, parity):
    n = L // 2
    xs, zs, cs = terms(L, lam, K, c)
    st, index = sector(n, parity)
    D = len(st)

    def mv(v):
        out = np.empty(D)
        _mv(st, index, xs, zs, cs, np.ascontiguousarray(v.real), out)
        return out

    op = spla.LinearOperator((D, D), matvec=mv, dtype=float)
    w, V = spla.eigsh(op, k=2, which="SA", tol=1e-10, ncv=30)
    o = np.argsort(w)
    return w[o], V[:, o], st, index, (xs, zs, cs)


def expect_b(L, st, index, psi, c):
    """<b_j> for j = 0..L-1 in a sector eigenvector psi (real)."""
    n = L // 2
    vals = []
    for j in range(L):
        x, z, cc = b_string(j, n, c)
        s = st ^ x
        sg = 1.0 - 2.0 * (np.array([bin(int(v)).count("1") & 1 for v in (z & s)]))
        vals.append(cc * np.sum(psi[index[s]] * sg * psi))
    return np.array(vals)
