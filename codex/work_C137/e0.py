"""Energy densities e_0 for the C137 test bonds.

Exact (literature) formulas for XXZ (Yang-Yang), TFIM (Pfeuty) and quadratic bonds;
variational upper bounds from a two-site iTEBD uniform MPS, symmetrized over the one-site shift
(a translation-invariant state on Z, so its energy density is a rigorous upper bound on e_0
up to floating point)."""
import numpy as np
import scipy.linalg as la
from scipy.integrate import quad
import scipy.sparse.linalg as sla


def e0_xxz(delta):
    """H = sum SxSx + SySy + Delta SzSz, S = sigma/2."""
    if abs(delta - 1) < 1e-14:
        return 0.25 - np.log(2)
    if -1 < delta < 1:
        g = np.arccos(delta)
        def f(x):
            if x < 1e-8:
                return (np.pi - g) / np.pi
            return 2 * np.exp(-2 * g * x) * (-np.expm1(-2 * (np.pi - g) * x)) / ((-np.expm1(-2 * np.pi * x)) * (1 + np.exp(-2 * g * x)))
        val = 2 * quad(f, 0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
        return delta / 4 - np.sin(g) / 2 * val
    if delta > 1:
        phi = np.arccosh(delta)
        s = sum(np.exp(-2 * n * phi) / (1 + np.exp(-2 * n * phi)) for n in range(1, 200000))
        return delta / 4 - np.sinh(phi) * (0.5 + 2 * s)
    if delta <= -1:
        return delta / 4
    raise ValueError


def e0_tfim(g):
    """H = -sum X X - g sum Z."""
    f = lambda k: np.sqrt(1 + g * g - 2 * g * np.cos(k))
    return -quad(f, 0, np.pi, limit=400, epsabs=1e-14, epsrel=1e-13)[0] / np.pi


# ---------------------------------------------------------------- iTEBD upper bound
def itebd_upper(h, chi=96, taus=(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001), steps=1500, seed=0):
    """Two-site iTEBD in imaginary time, then exact energy of the symmetrized uniform MPS."""
    h = np.asarray(h, complex)
    d = 2
    rng = np.random.default_rng(seed)
    G = [rng.normal(size=(1, d, 1)) + 0j, rng.normal(size=(1, d, 1)) + 0j]
    lam = [np.ones(1), np.ones(1)]
    for tau in taus:
        U = la.expm(-tau * h).reshape(d, d, d, d)
        for step in range(steps):
            for A in (0, 1):
                B = 1 - A
                la_, lb = lam[B], lam[A]
                theta = np.einsum('a,aib,b,bjc,c->aijc', la_, G[A], lam[A], G[B], la_)
                theta = np.einsum('ijkl,aklc->aijc', U, theta)
                Da, _, _, Dc = theta.shape
                X, S, Y = la.svd(theta.reshape(Da * d, d * Dc), full_matrices=False)
                k = min(chi, int(np.sum(S > 1e-13 * S[0])))
                X, S, Y = X[:, :k], S[:k], Y[:k]
                S = S / np.linalg.norm(S)
                lam[A] = S
                inv = 1 / la_
                G[A] = (inv[:, None] * X.reshape(Da, d * k)).reshape(Da, d, k)
                G[B] = (Y.reshape(k, d, Dc) * inv[None, None, :])
    # build MPS tensors A1 = lam[B] G[A], A2 = lam[A] G[B] for a two-site unit cell
    A1 = np.einsum('a,aib->aib', lam[1], G[0])
    A2 = np.einsum('a,aib->aib', lam[0], G[1])
    return mps_energy(h, A1, A2)


def mps_energy(h, A1, A2):
    """Energy density of the shift-symmetrized two-periodic MPS ...A1 A2 A1 A2..., exactly
    from transfer-matrix fixed points computed by Arnoldi (no canonical-form assumption)."""
    d = A1.shape[1]
    D = A1.shape[0]
    def right_apply(x):   # x: (D, D) env to the right of a cell; returns env to the left of the cell
        X = x.reshape(D, D)
        X = np.einsum('bjc,cd,ejd->be', A2, X, A2.conj(), optimize=True)
        X = np.einsum('aib,be,cie->ac', A1, X, A1.conj(), optimize=True)
        return X.ravel()
    def left_apply(y):
        Y = y.reshape(D, D)
        Y = np.einsum('ac,aib,cie->be', Y, A1, A1.conj(), optimize=True)
        Y = np.einsum('be,bjc,ejd->cd', Y, A2, A2.conj(), optimize=True)
        return Y.ravel()
    N = D * D
    if N <= 400:
        Tm = np.array([right_apply(e) for e in np.eye(N)]).T
        wr, vr = la.eig(Tm); wl, vl = la.eig(Tm.conj().T)
    else:
        op = lambda f: sla.LinearOperator((N, N), matvec=f, dtype=complex)
        k = min(10, N - 2)
        wr, vr = sla.eigs(op(right_apply), k=k, which='LM', tol=1e-13)
        wl, vl = sla.eigs(op(left_apply), k=k, which='LM', tol=1e-13)
    lam = np.max(wr.real[np.abs(wr.imag) < 1e-9 * np.abs(wr).max()])
    selr = np.abs(wr - lam) < 1e-7 * abs(lam); sell = np.abs(wl - np.conj(lam)) < 1e-7 * abs(lam)
    Rv = vr[:, selr]; Lv = vl[:, sell]
    if Rv.shape[1] != Lv.shape[1]:
        raise RuntimeError('fixed-point spaces of different dimension')
    # spectral projector onto the fixed-point space, P = Rv (Lv^* Rv)^{-1} Lv^*, applied to the identity:
    # P(1) and P^*(1) are positive semidefinite exact fixed points (limits of Cesaro means of T^k(1)).
    G = Lv.conj().T @ Rv
    one = np.eye(D).ravel()
    R = (Rv @ la.solve(G, Lv.conj().T @ one)).reshape(D, D)
    L = (Lv @ la.solve(G.conj().T, Rv.conj().T @ one)).reshape(D, D)
    R = (R + R.conj().T) / 2; L = (L + L.conj().T) / 2
    for X in (R, L):
        ev = la.eigvalsh(X)
        if ev[0] < -1e-8 * ev[-1]:
            raise RuntimeError('fixed point not positive: %g %g' % (ev[0], ev[-1]))
    global LAST_DIAG
    LAST_DIAG = dict(degeneracy=int(Rv.shape[1]), lam=float(lam))
    h4 = h.reshape(d, d, d, d)
    def bond(X, Y, Lm, Rm):
        theta = np.einsum('aib,bjc->aijc', X, Y)
        Ht = np.einsum('klij,aijc->aklc', h4, theta)
        num = np.einsum('ac,aijb,cijd,bd->', Lm, Ht, theta.conj(), Rm, optimize=True)
        den = np.einsum('ac,aijb,cijd,bd->', Lm, theta, theta.conj(), Rm, optimize=True)
        return (num / den).real
    e12 = bond(A1, A2, L, R)
    # bond (A2, A1 of next cell): left env = L transported through A1, right env = R transported through A2
    L2 = np.einsum('ac,aib,cie->be', L, A1, A1.conj(), optimize=True)
    R2 = np.einsum('bjc,cd,ejd->be', A2, R, A2.conj(), optimize=True)
    e21 = bond(A2, A1, L2, R2)
    return float((e12 + e21) / 2), float(e12), float(e21)
