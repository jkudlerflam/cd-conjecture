#!/usr/bin/env python3
"""Core routines for the bound attack: the gauge-optimized bond scale M, a symmetry-free
matrix-free ring ED (any Hermitian bond on (C^d)^L, periodic), Casimir fits and the
auxiliary criticality diagnostics (gap scaling, half-ring entropy).

Conventions. H = sum_i h_{i,i+1} on a ring of L sites, h a d^2 x d^2 Hermitian matrix in
the basis |a b> = a*d + b (site i first). M = min_o (lambda_max - lambda_min)/2 of
h + o (x) 1 - 1 (x) o over Hermitian o. A is the Casimir coefficient in
E_0(L) = e_0 L - pi A/(6 L) + ..., so the conjecture reads A <= pi M log2 d.

numpy and scipy only.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import numpy as np
from scipy import optimize
from scipy.sparse.linalg import LinearOperator, eigsh

LOG2 = np.log(2.0)


# ---------------------------------------------------------------- gauge-optimized M
def herm_basis(d):
    B = []
    for a in range(d):
        for b in range(a + 1, d):
            m = np.zeros((d, d), complex); m[a, b] = m[b, a] = 1 / np.sqrt(2); B.append(m)
            m = np.zeros((d, d), complex); m[a, b] = -1j / np.sqrt(2); m[b, a] = 1j / np.sqrt(2); B.append(m)
    for a in range(1, d):
        m = np.zeros((d, d), complex)
        m[np.arange(a), np.arange(a)] = 1; m[a, a] = -a
        B.append(m / np.sqrt(a * (a + 1)))
    return np.array(B)


def spread(h):
    e = np.linalg.eigvalsh((h + h.conj().T) / 2)
    return float(e[-1] - e[0])


def hnorm(h, d, return_o=False, betas=(4, 16, 64, 256, 1024, 4096), gauge_right=None):
    """M = min_o spread(h + o(x)1 - 1(x)o)/2.  The spread is convex in o; minimize a
    log-sum-exp smoothing at increasing inverse temperature (BFGS, analytic gradient),
    then report the exact spread at the optimum.  gauge_right, if given, is a d x d unitary
    R, and the gauge set becomes o(x)1 - 1(x)R o R^dag (the twisted gauge set of a bond
    (1(x)R) h (1(x)R^dag), see README, Lemma 3)."""
    h = np.asarray(h, complex); h = (h + h.conj().T) / 2
    B = herm_basis(d)
    I = np.eye(d)
    if gauge_right is None:
        KB = np.array([np.kron(b, I) - np.kron(I, b) for b in B])
    else:
        R = gauge_right
        KB = np.array([np.kron(b, I) - np.kron(I, R @ b @ R.conj().T) for b in B])
    s0 = spread(h)
    if s0 < 1e-14:
        return (0.0, np.zeros((d, d))) if return_o else 0.0

    def soft(x, beta):
        K = h + np.tensordot(x, KB, 1)
        e, U = np.linalg.eigh(K)
        bb = beta / s0
        wp = np.exp(bb * (e - e[-1])); wp /= wp.sum()
        wm = np.exp(-bb * (e - e[0])); wm /= wm.sum()
        f = (e[-1] + np.log(np.sum(np.exp(bb * (e - e[-1])))) / bb) - (e[0] - np.log(np.sum(np.exp(-bb * (e - e[0])))) / bb)
        g = np.einsum('kab,ai,bi,i->k', KB, U.conj(), U, wp - wm).real
        return f, g

    x = np.zeros(len(B))
    for beta in betas:
        r = optimize.minimize(soft, x, args=(beta,), jac=True, method='BFGS',
                              options=dict(gtol=1e-12, maxiter=2000))
        x = r.x
    best = spread(h + np.tensordot(x, KB, 1))
    if best > s0:
        best, x = s0, np.zeros(len(B))
    o = np.tensordot(x, B, 1)
    return (best / 2, o) if return_o else best / 2


# ---------------------------------------------------------------- ring ED, no symmetry
class Ring:
    """Matrix-free H = sum_i h_{i,i+1} on d^L states.  Site i is tensor axis i (site 0 most
    significant).  Works for any Hermitian bond; cost per matvec L d^(L+2)."""

    def __init__(self, h, d, L):
        self.d, self.L = d, L
        self.h = np.asarray(h, complex).reshape(d, d, d, d)  # [a', b', a, b]
        self.real = np.allclose(self.h.imag, 0)
        if self.real:
            self.h = self.h.real.copy()
        self.N = d ** L
        self.dtype = float if self.real else complex

    def matvec(self, v):
        d, L = self.d, self.L
        psi = v.reshape([d] * L)
        out = np.zeros_like(psi)
        for i in range(L - 1):
            # contract h[a'b', a b] with axes (i, i+1) of psi
            t = np.tensordot(self.h, psi, axes=([2, 3], [i, i + 1]))  # [a', b', rest...]
            out += np.moveaxis(t, [0, 1], [i, i + 1])
        # wrap bond (L-1, 0)
        t = np.tensordot(self.h, psi, axes=([2, 3], [L - 1, 0]))
        out += np.moveaxis(t, [0, 1], [L - 1, 0])
        return out.reshape(-1)

    def operator(self):
        return LinearOperator((self.N, self.N), matvec=self.matvec, dtype=self.dtype)

    def lowest(self, k=4, tol=1e-10, v0=None, return_vec=False):
        op = self.operator()
        kk = min(k, self.N - 2)
        w, v = eigsh(op, k=kk, which='SA', tol=tol, v0=v0, ncv=max(2 * kk + 10, 24), maxiter=20000)
        idx = np.argsort(w)
        w, v = w[idx], v[:, idx]
        return (w, v[:, 0]) if return_vec else w


def half_entropy(psi, d, L, cut=None):
    cut = L // 2 if cut is None else cut
    m = psi.reshape(d ** cut, -1)
    s = np.linalg.svd(m, compute_uv=False) ** 2
    s = s[s > 1e-14]
    return float(-np.sum(s * np.log(s)))


# ---------------------------------------------------------------- fits
def fit_cas(Ls, Es, e0=None, quartic=True):
    """Least squares E/L = e0 - (pi/6) A/L^2 [+ b/L^4].  Returns (A, e0, b)."""
    L = np.asarray(Ls, float); y = np.asarray(Es, float) / L
    cols = [] if e0 is not None else [np.ones_like(L)]
    cols.append(-np.pi / (6 * L ** 2))
    if quartic:
        cols.append(1 / L ** 4)
    X = np.vstack(cols).T
    yy = y - (e0 if e0 is not None else 0)
    coef = np.linalg.lstsq(X, yy, rcond=None)[0]
    if e0 is None:
        return float(coef[1]), float(coef[0]), (float(coef[2]) if quartic else 0.0)
    return float(coef[0]), e0, (float(coef[1]) if quartic else 0.0)


def ratio(A, M, d):
    return A / (np.pi * M * np.log(d) / LOG2)


def evaluate(h, d, Ls, nlev=3, with_entropy=True, M=None, tol=1e-9):
    """Ground energies, gaps and half-ring entropies on rings Ls; Casimir fit; ratio."""
    out = {'d': d, 'Ls': list(map(int, Ls)), 'E0': [], 'gap1': [], 'gap2': [], 'S': []}
    for L in Ls:
        r = Ring(h, d, L)
        if with_entropy:
            w, v = r.lowest(k=nlev, tol=tol, return_vec=True)
            out['S'].append(half_entropy(v, d, L))
        else:
            w = r.lowest(k=nlev, tol=tol)
        out['E0'].append(float(w[0]))
        out['gap1'].append(float(w[1] - w[0]) if len(w) > 1 else np.nan)
        out['gap2'].append(float(w[2] - w[0]) if len(w) > 2 else np.nan)
    A, e0, b = fit_cas(out['Ls'], out['E0'])
    out['A_fit'], out['e0_fit'], out['b_fit'] = A, e0, b
    if len(Ls) >= 3:
        A3, e03, _ = fit_cas(out['Ls'][-3:], out['E0'][-3:], quartic=True)
        out['A_fit_last3'] = A3
    A2, _, _ = fit_cas(out['Ls'][-2:], out['E0'][-2:], quartic=False)
    out['A_fit_last2'] = A2
    out['M'] = float(hnorm(h, d)) if M is None else float(M)
    out['ratio'] = float(ratio(A, out['M'], d))
    out['gapL'] = [float(g * L) for g, L in zip(out['gap1'], out['Ls'])]
    if with_entropy and len(Ls) >= 2:
        S1, S2 = out['S'][-2], out['S'][-1]
        L1, L2 = out['Ls'][-2], out['Ls'][-1]
        out['c_half'] = float(3 * (S2 - S1) / np.log(L2 / L1))
    return out
