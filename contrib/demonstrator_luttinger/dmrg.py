#!/usr/bin/env python3
"""Two-site finite DMRG on open chains with site-dependent dimensions (real arithmetic).

MPO W[s] has shape (Dl, Dr, d, d) = (left channel, right channel, out, in), upper triangular.
Excited states: H + w_ortho sum_k |psi_k><psi_k| through overlap environments.
Measurements: entropy on every bond, local charges, left-block charge moments, window transfer spectrum.
Claude (Opus 5.5), 2026-10-28.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import eigs, LinearOperator


def lt(L, A, W):
    t = np.tensordot(L, A, ([2], [0]))               # a', w, s, b
    t = np.tensordot(t, W, ([1, 2], [0, 3]))         # a', b, w', s'
    t = np.tensordot(A, t, ([0, 1], [0, 3]))         # b', b, w'   (real: conj = id)
    return t.transpose(0, 2, 1)


def rt(R, B, W):
    t = np.tensordot(B, R, ([2], [0]))               # a, s, w', b'
    t = np.tensordot(t, W, ([1, 2], [3, 1]))         # a, b', w, s'
    t = np.tensordot(t, B, ([1, 3], [2, 1]))         # a, w, a'
    return t


def hac2(L, W1, W2, R, x):
    t = np.tensordot(L, x, ([2], [0]))               # a', w, s1, s2, b
    t = np.tensordot(t, W1, ([1, 2], [0, 3]))        # a', s2, b, w1, s1'
    t = np.tensordot(t, W2, ([3, 1], [0, 3]))        # a', b, s1', w2, s2'
    t = np.tensordot(t, R, ([1, 3], [0, 1]))         # a', s1', s2', b'
    return t


def lanczos(apply, x0, kmax=12, tol=1e-13):
    """Lowest Ritz pair in a Krylov space of size <= kmax (full reorthogonalization)."""
    n0 = np.linalg.norm(x0)
    V = [x0 / n0]
    a, b = [], []
    w = apply(V[0])
    best = None
    for k in range(kmax):
        ak = float(np.vdot(V[k], w)); a.append(ak)
        w = w - ak * V[k] - (b[-1] * V[k - 1] if k > 0 else 0)
        for v in V:
            w = w - np.vdot(v, w) * v
        bk = np.linalg.norm(w)
        T = np.diag(a) + np.diag(b, 1) + np.diag(b, -1)
        e, U = np.linalg.eigh(T)
        best = (e[0], U[:, 0])
        if bk < tol or k == kmax - 1 or bk * abs(U[-1, 0]) < tol:
            break
        b.append(bk)
        V.append(w / bk)
        w = apply(V[-1])
    e0, u = best
    x = sum(c * v for c, v in zip(u, V[:len(u)]))
    return e0, x / np.linalg.norm(x)


class DMRG:
    def __init__(self, Ws, dims, mps=None, chi0=8, seed=1):
        self.Ws, self.dims, self.L = Ws, dims, len(dims)
        L = self.L
        if mps is None:
            rng = np.random.default_rng(seed)
            lg = np.log(np.array(dims, float))
            cum_l = np.concatenate([[0.0], np.cumsum(lg)])
            cum_r = np.concatenate([np.cumsum(lg[::-1])[::-1], [0.0]])
            bd = [int(round(min(chi0, np.exp(cum_l[i]), np.exp(cum_r[i])))) for i in range(L + 1)]
            mps = [rng.normal(size=(bd[i], dims[i], bd[i + 1])) for i in range(L)]
        self.M = [m.copy() for m in mps]
        self.right_canonicalize()
        self.ortho = []       # list of (weight, mps tensors) for excited states
        self.S = np.zeros(L - 1)
        self.build_envs()

    def right_canonicalize(self):
        M = self.M
        for i in range(self.L - 1, 0, -1):
            a, s, b = M[i].shape
            q, r = np.linalg.qr(M[i].reshape(a, s * b).T)
            M[i] = q.T.reshape(-1, s, b)
            M[i - 1] = np.tensordot(M[i - 1], r.T, ([2], [0]))
        M[0] /= np.linalg.norm(M[0])

    def build_envs(self):
        L = self.L
        self.Le = [None] * (L + 1); self.Re = [None] * (L + 1)
        self.Le[0] = np.ones((1, 1, 1)); self.Re[L] = np.ones((1, 1, 1))
        for i in range(L - 1, 0, -1):
            self.Re[i] = rt(self.Re[i + 1], self.M[i], self.Ws[i])
        self.Lo = [[None] * (L + 1) for _ in self.ortho]
        self.Ro = [[None] * (L + 1) for _ in self.ortho]
        for k, (_, P) in enumerate(self.ortho):
            self.Lo[k][0] = np.ones((1, 1)); self.Ro[k][L] = np.ones((1, 1))
            for i in range(L - 1, 0, -1):
                self.Ro[k][i] = self._ro(self.Ro[k][i + 1], self.M[i], P[i])

    @staticmethod
    def _ro(R, B, P):     # R (cur, other) -> through site
        t = np.tensordot(P, R, ([2], [1]))           # p, s, c
        return np.tensordot(B, t, ([1, 2], [1, 2]))  # c', p

    @staticmethod
    def _lo(L, A, P):
        t = np.tensordot(L, P, ([1], [0]))           # a, s, p'
        return np.tensordot(A, t, ([0, 1], [0, 1]))  # a', p'

    def add_ortho(self, mps, weight):
        self.ortho.append((weight, [m.copy() for m in mps]))
        self.build_envs()

    def sweep(self, chi, kmax=12, cutoff=1e-13):
        L, Ws, M = self.L, self.Ws, self.M
        maxtr = 0.0
        E = None
        for direction in (1, -1):
            rng = range(L - 1) if direction == 1 else range(L - 2, -1, -1)
            for i in rng:
                th = np.tensordot(M[i], M[i + 1], ([2], [0]))
                shp = th.shape
                Le, Re, W1, W2 = self.Le[i], self.Re[i + 2], Ws[i], Ws[i + 1]
                projs = []
                for k, (wk, P) in enumerate(self.ortho):
                    t = np.tensordot(self.Lo[k][i], P[i], ([1], [0]))
                    t = np.tensordot(t, P[i + 1], ([2], [0]))
                    t = np.tensordot(t, self.Ro[k][i + 2], ([3], [1]))
                    projs.append((wk, t.ravel()))

                def f(v):
                    out = hac2(Le, W1, W2, Re, v.reshape(shp)).ravel()
                    for wk, p in projs:
                        out = out + wk * np.dot(p, v) * p
                    return out
                E, x = lanczos(f, th.ravel(), kmax=kmax)
                a, s1, s2, b = shp
                try:
                    u, s, vh = linalg.svd(x.reshape(a * s1, s2 * b), full_matrices=False, lapack_driver='gesdd')
                except linalg.LinAlgError:
                    u, s, vh = linalg.svd(x.reshape(a * s1, s2 * b), full_matrices=False, lapack_driver='gesvd')
                keep = min(chi, max(1, int((s > cutoff * s[0]).sum())))
                maxtr = max(maxtr, float((s[keep:] ** 2).sum()))
                s = s[:keep] / np.linalg.norm(s[:keep])
                p = s ** 2; p = p[p > 1e-30]
                self.S[i] = float(-(p * np.log(p)).sum())
                if direction == 1:
                    M[i] = u[:, :keep].reshape(a, s1, keep)
                    M[i + 1] = (s[:, None] * vh[:keep]).reshape(keep, s2, b)
                    self.Le[i + 1] = lt(self.Le[i], M[i], Ws[i])
                    for k, (_, P) in enumerate(self.ortho):
                        self.Lo[k][i + 1] = self._lo(self.Lo[k][i], M[i], P[i])
                else:
                    M[i + 1] = vh[:keep].reshape(keep, s2, b)
                    M[i] = (u[:, :keep] * s[None, :]).reshape(a, s1, keep)
                    self.Re[i + 1] = rt(self.Re[i + 2], M[i + 1], Ws[i + 1])
                    for k, (_, P) in enumerate(self.ortho):
                        self.Ro[k][i + 1] = self._ro(self.Ro[k][i + 2], M[i + 1], P[i + 1])
        return float(E), maxtr

    # ---------------------------------------------------------------- measurements (centre at site 0)
    def overlap(self, P):
        X = np.ones((1, 1))
        for i in range(self.L):
            X = self._lo(X, self.M[i], P[i])
        return float(X[0, 0])

    def local(self, ops_per_site):
        """<O_s> for one diagonal operator (vector) per site."""
        out = np.zeros(self.L)
        X = np.ones((1, 1))
        for j in range(self.L):
            A = self.M[j]
            t = np.tensordot(X, A, ([1], [0]))
            out[j] = np.tensordot(A, t * ops_per_site[j][None, :, None], ([0, 1, 2], [0, 1, 2]))
            X = np.tensordot(A, t, ([0, 1], [0, 1]))
        return out

    def block_moments(self, q):
        """<Q_A>, <Q_A^2> for the left blocks A = sites [0, l), l = 1..L-1; q[s] = diagonal charge vector."""
        E0 = np.ones((1, 1)); E1 = np.zeros((1, 1)); E2 = np.zeros((1, 1))
        m1, m2 = [], []
        for j in range(self.L - 1):
            A = self.M[j]
            Aq = A * q[j][None, :, None]
            Aqq = A * (q[j] ** 2)[None, :, None]
            T = lambda X, B: np.tensordot(A, np.tensordot(X, B, ([1], [0])), ([0, 1], [0, 1]))
            E2 = T(E2, A) + 2 * T(E1, Aq) + T(E0, Aqq)
            E1 = T(E1, A) + T(E0, Aq)
            E0 = T(E0, A)
            m1.append(np.trace(E1)); m2.append(np.trace(E2))
        return np.array(m1), np.array(m2)

    def window_spectrum(self, j0, m, nev=4):
        """Leading eigenvalues of the transfer matrix of sites [j0, j0+m) (right-canonical part)."""
        Bs = self.M[j0:j0 + m]
        chi = Bs[0].shape[0]
        assert Bs[-1].shape[2] == chi
        N = chi * chi

        def mv(v):
            X = v.reshape(chi, chi)
            for B in Bs[::-1]:
                X = np.tensordot(np.tensordot(B, X, ([2], [0])), B, ([1, 2], [1, 2]))
            return X.ravel()
        if N <= 300:
            Mx = np.stack([mv(e) for e in np.eye(N)], axis=1)
            w = linalg.eigvals(Mx)
        else:
            w = eigs(LinearOperator((N, N), matvec=mv, dtype=float), k=nev, which='LM', tol=1e-10,
                     return_eigenvectors=False, maxiter=5000, ncv=max(2 * nev + 1, 24))
        return w[np.argsort(-np.abs(w))][:nev]


def run_ground(Ws, dims, chis, mps=None, sweeps_min=3, sweeps_max=40, etol=1e-11, kmax=12, seed=1, log=print,
               ortho=(), chi0=8, inter_max=None):
    dm = DMRG(Ws, dims, mps=mps, seed=seed, chi0=chi0)
    for w, P in ortho:
        dm.ortho.append((w, P))
    dm.build_envs()
    hist = []
    for ic, chi in enumerate(chis):
        Eprev = None
        nmax = sweeps_max if (ic == len(chis) - 1 or inter_max is None) else min(sweeps_max, inter_max)
        for sw in range(nmax):
            E, tr = dm.sweep(chi, kmax=kmax)
            hist.append((chi, E, tr))
            if log:
                log(f"  chi={chi} sweep={sw} E={E:.13f} trunc={tr:.1e}")
            if sw + 1 >= sweeps_min and Eprev is not None and (E - Eprev) > -etol:
                break
            Eprev = E
    return dm, E, hist
