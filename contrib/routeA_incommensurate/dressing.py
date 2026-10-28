#!/usr/bin/env python3
"""First-order dressed initial state for the full chain.

On a code-only state, prod_b (1 + s_b) with the first-order Schrieffer-Wolff amplitude
    s_b = -lam sum_a sum_{e,f} (P_e FR_a P_c) (x) (P_f FL_a P_c) / (e + f)
is an MPO of bond dimension 1 + 3 r (r = rank of the matrix 1/(e+f) over the distinct excited
cell energies), because s_{b+1} s_b = 0 on code states.  Applied to the code-embedded logical uMPS
it gives the full-chain state to first order in lam; uMPS gauge fixing by iterated QR.
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-28.
"""
import numpy as np
from scipy import linalg


def dressing_mpo(cb, lam, tol=1e-12):
    e, N, D = cb['e'], cb['N'], cb['D']
    levels = sorted(set(np.round(e[N:], 9)))
    P = [np.diag((np.abs(e - x) < 1e-6).astype(float)) for x in levels]
    Pc = np.diag((np.arange(D) < N).astype(float))
    Cm = 1.0 / (np.array(levels)[:, None] + np.array(levels)[None, :])
    U, s, Vh = linalg.svd(Cm)
    r = int(np.sum(s > tol * s[0]))
    u = U[:, :r] * s[:r]; v = Vh[:r].T
    Dw = 1 + 3 * r
    W = np.zeros((Dw, Dw, D, D), complex)
    W[0, 0] = Pc
    for a in range(3):
        for k in range(r):
            ch = 1 + a * r + k
            W[0, ch] = -lam * sum(u[i, k] * P[i] for i in range(len(levels))) @ cb['FR'][a] @ Pc
            W[ch, 0] = sum(v[i, k] * P[i] for i in range(len(levels))) @ cb['FL'][a] @ Pc
    return W


def apply_mpo(A, W):
    """A (chi, d, chi), W (w, w, d, d) -> (chi w, d, chi w)."""
    chi, d, chi2 = A.shape
    w = W.shape[0]
    T = np.einsum('xyst,atb->axsby', W, A)
    return T.reshape(chi * w, d, chi2 * w)


def _fixed_point(As, side):
    """Dominant fixed point of the cell transfer matrix: side='L': X -> sum A^dag X A (left),
    side='R': X -> sum A X A^dag (right).  Returns the PSD fixed point on bond 0, trace 1."""
    from scipy.sparse.linalg import LinearOperator, eigs
    chi = As[0].shape[0]

    def mv(v):
        X = v.reshape(chi, chi)
        if side == 'L':
            for A in As:
                X = np.tensordot(A.conj(), np.tensordot(X, A, ([1], [0])), ([0, 1], [0, 1]))
        else:
            for A in As[::-1]:
                X = np.tensordot(np.tensordot(A, X, ([2], [0])), A.conj(), ([1, 2], [1, 2]))
        return X.ravel()
    n2 = chi * chi
    if n2 <= 900:
        M = np.zeros((n2, n2), complex)
        for i in range(n2):
            e = np.zeros(n2, complex); e[i] = 1; M[:, i] = mv(e)
        w, U = linalg.eig(M); v = U[:, np.argmax(np.abs(w))]
    else:
        w, U = eigs(LinearOperator((n2, n2), matvec=mv, dtype=complex), k=1, which='LM', tol=1e-13, maxiter=20000)
        v = U[:, 0]
    X = v.reshape(chi, chi)
    X = X / np.trace(X)
    return 0.5 * (X + X.conj().T)


def canonical(As, cut=1e-12):
    """Mixed canonical form (AL, AR, C) of the uniform MPS with unit cell As (equal bond dims).
    Directions with (numerically) zero weight are removed first, so the result is injective."""
    n = len(As)
    for _ in range(3):
        # left gauge on every bond from the left fixed point
        X = _fixed_point(As, 'L')
        w, U = linalg.eigh(X)
        keep = w > cut * w.max()
        Lg = (U[:, keep] * np.sqrt(w[keep])).conj().T          # Lg^dag Lg = X restricted
        # propagate: L_k A_k = AL_k L_{k+1}
        AL, Ls = [], [Lg]
        for k in range(n):
            M = np.tensordot(Ls[-1], As[k], ([1], [0]))
            a, d, b = M.shape
            Q, R = linalg.qr(M.reshape(a * d, b), mode='economic')
            # drop null directions of R (zero-weight bond states)
            u_, s_, vh_ = linalg.svd(R)
            kk = s_ > cut * s_[0]
            Q = Q @ u_[:, kk]; R = s_[kk, None] * vh_[kk]
            AL.append(Q.reshape(a, d, -1)); Ls.append(R)
        # make the cell closed: map the bond-n space back to bond 0
        G = Ls[-1] @ linalg.pinv(Ls[0])
        AL[-1] = np.tensordot(AL[-1], G, ([2], [0]))
        # As <- AL (now left-orthonormal up to the closing map); iterate once more for safety
        As = AL
        chis = [A.shape[0] for A in As]
        if all(As[k].shape[2] == As[(k + 1) % n].shape[0] for k in range(n)):
            ok = max(np.linalg.norm(np.tensordot(A.conj(), A, ([0, 1], [0, 1])) - np.eye(A.shape[2])) for A in As)
            if ok < 1e-9:
                break
    AL = As
    # right fixed point on every bond: rho_k = E_{AL_k}(rho_{k+1})
    rho0 = _fixed_point(AL, 'R')
    rhos = [None] * n
    rhos[0] = rho0
    X = rho0
    for k in range(n - 1, 0, -1):
        X = np.tensordot(np.tensordot(AL[k], X, ([2], [0])), AL[k].conj(), ([1, 2], [1, 2]))
        X = X / np.trace(X)
        rhos[k] = 0.5 * (X + X.conj().T)
    Cs, AR, Us = [], [], []
    for k in range(n):
        w, U = linalg.eigh(rhos[k])
        w = np.clip(w, 1e-300, None)
        Cs.append((U * np.sqrt(w)) @ U.conj().T)
    for k in range(n):
        AC = np.tensordot(AL[k], Cs[(k + 1) % n], ([2], [0]))
        AR.append(np.tensordot(linalg.inv(Cs[k]), AC, ([1], [0])))
    # diagonal C gauge
    out_AL, out_AR, out_C, Uk, Vk = [], [], [], [], []
    for k in range(n):
        u_, s_, vh_ = linalg.svd(Cs[k])
        Uk.append(u_); Vk.append(vh_.conj().T); out_C.append(np.diag(s_ / np.linalg.norm(s_)) + 0j)
    for k in range(n):
        out_AL.append(np.einsum('ab,bsc,cd->asd', Uk[k].conj().T, AL[k], Uk[(k + 1) % n]))
        out_AR.append(np.einsum('ab,bsc,cd->asd', Vk[k].conj().T, AR[k], Vk[(k + 1) % n]))
    return out_AL, out_AR, out_C


def truncate(AL, AR, C, chi):
    """Keep the chi largest Schmidt values on every bond and re-canonicalize."""
    n = len(AL)
    if all(c.shape[0] <= chi for c in C):
        return AL, AR, C
    m = min(chi, min(c.shape[0] for c in C))
    As = [np.tensordot(AL[k], C[(k + 1) % n], ([2], [0]))[:m, :, :m] for k in range(n)]
    return canonical(As)
