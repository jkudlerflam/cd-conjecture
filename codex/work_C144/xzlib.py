"""Shared utilities for C144 (XZ point and its tuned Ashkin-Teller cone).

Bond convention (as in C143 ed.py and specs/_d2_background.md):
    h = sum_{ab} J[a,b] s^a (x) s^b + sum_a B[a] (s^a (x) 1 + 1 (x) s^a)/2 ,  H_L = sum_i h_{i,i+1} on a ring.
Basis: site i is bit (L-1-i) of the integer.
"""
import numpy as np
import scipy.sparse.linalg as spla
from scipy.optimize import minimize
import numba as nb

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [SX, SY, SZ]
NAMES = "xyz"


def bond(J, B=(0, 0, 0)):
    J = np.asarray(J, float).reshape(3, 3)
    h = np.zeros((4, 4), complex)
    for a in range(3):
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
        h += B[a] * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a])) / 2
    return h


def term(name):
    """Bond matrix for a named direction: 'xz','zx','yy',..., 'Bx','By','Bz'."""
    if name[0] == "B":
        B = np.zeros(3); B[NAMES.index(name[1])] = 1.0
        return bond(np.zeros((3, 3)), B)
    J = np.zeros((3, 3)); J[NAMES.index(name[0]), NAMES.index(name[1])] = 1.0
    return bond(J)


def gauge_M(h):
    """M = min over gauge o of (lmax - lmin)/2 of h + o.(s(x)1 - 1(x)s), with the trace removed."""
    G = [np.kron(P, I2) - np.kron(I2, P) for P in PAULI]

    def f(o):
        w = np.linalg.eigvalsh(h + sum(o[a] * G[a] for a in range(3)))
        return (w[-1] - w[0]) / 2

    best = None
    for o0 in [np.zeros(3), 0.3 * np.ones(3), -0.3 * np.ones(3), np.array([0.2, -0.4, 0.1])]:
        r = minimize(f, o0, method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-10, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun


@nb.njit(parallel=True, cache=True)
def _matvec_ring(hr, hi, L, v_re, v_im, o_re, o_im):
    N = 1 << L
    for x in nb.prange(N):
        acc_re = 0.0
        acc_im = 0.0
        for i in range(L):
            a = L - 1 - i
            b = L - 1 - ((i + 1) % L)
            bi = (x >> a) & 1
            bj = (x >> b) & 1
            outidx = 2 * bi + bj
            base = x & ~((1 << a) | (1 << b))
            for inn in range(4):
                c_re = hr[outidx, inn]
                c_im = hi[outidx, inn]
                if c_re == 0.0 and c_im == 0.0:
                    continue
                y = base | ((inn >> 1) << a) | ((inn & 1) << b)
                acc_re += c_re * v_re[y] - c_im * v_im[y]
                acc_im += c_re * v_im[y] + c_im * v_re[y]
        o_re[x] = acc_re
        o_im[x] = acc_im


def ring_operator(h, L):
    """LinearOperator for H_L = sum_i h_{i,i+1} (gather form: out[x] = sum_y H[x,y] v[y])."""
    N = 1 << L
    hr = np.ascontiguousarray(h.real); hi = np.ascontiguousarray(h.imag)

    def mv(v):
        v = np.asarray(v).ravel()
        vr = np.ascontiguousarray(v.real); vi = np.ascontiguousarray(v.imag)
        o_re = np.empty(N); o_im = np.empty(N)
        _matvec_ring(hr, hi, L, vr, vi, o_re, o_im)
        return o_re + 1j * o_im

    return spla.LinearOperator((N, N), matvec=mv, dtype=complex)


def ground_energy(h, L, k=1, tol=1e-10):
    if L <= 10:
        from scipy.sparse import csr_matrix
        N = 1 << L
        Hd = np.zeros((N, N), complex)
        op = ring_operator(h, L)
        for x in range(N):
            e = np.zeros(N, complex); e[x] = 1.0
            Hd[:, x] = op.matvec(e)
        w = np.linalg.eigvalsh(Hd)
        return w[:k] if k > 1 else w[0]
    op = ring_operator(h, L)
    w = spla.eigsh(op, k=max(k, 2), which="SA", tol=tol, return_eigenvectors=False, ncv=40)
    w = np.sort(w)
    return w[:k] if k > 1 else w[0]


# ---------------- free Majorana utilities ----------------

def majorana_modes(A):
    """A real antisymmetric, H = (i/2) sum_{mn} A_mn g_m g_n = (1/4) g^T (2iA) g, {g_m,g_n} = 2 delta.
    Returns (eps>0 excitation energies, V columns) with g_m = sum_mu sqrt2 (V_m mu d_mu + conj d_mu^dag)
    up to the convention fixed by the self-test, and E0 = -(1/2) sum eps."""
    h = 2j * A
    w, U = np.linalg.eigh(h)
    pos = w > 1e-12
    return w[pos], U[:, pos], w


def chain_A(L, bc=-1, t=1.0):
    """Uniform Majorana ring sum_j i t g_j g_{j+1}, with g_{L+1} = bc g_1."""
    A = np.zeros((L, L))
    for j in range(L):
        k = (j + 1) % L
        s = bc if j == L - 1 else 1.0
        A[j, k] += t * s
        A[k, j] -= t * s
    return A
