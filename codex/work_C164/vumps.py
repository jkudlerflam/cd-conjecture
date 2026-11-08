#!/usr/bin/env python3
"""VUMPS (variational uniform MPS) with an n-site unit cell and a general
finite-range MPO, numpy/scipy only.  Written for contrib/d2_followup_220_1282
(Claude Opus 5.5, 2026-10-27).

Conventions
-----------
MPS tensors A[k] have shape (chi_k, d, chi_{k+1}); bond k sits to the left of
site k, bond n is bond 0 of the next cell.  Mixed gauge:
    A_C[k] = A_L[k] C[k+1] = C[k] A_R[k].
MPO tensors W[wl, wr, s_out, s_in] are upper triangular; channel 0 is
"nothing placed yet", channel D-1 is "term completed".
Left environments L[k] have shape (chi_k, D, chi_k) = (bra, mpo, ket);
right environments R[k] have shape (chi_k, D, chi_k) = (ket, mpo, bra), and
R[k] includes sites k, k+1, ... .

The algorithm is the parallel multi-site VUMPS of Zauner-Stauber,
Vanderstraeten, Fishman, Verstraete, Haegeman, PRB 97, 045145 (2018), with the
two-site null-space subspace expansion of the same paper to grow chi.
mode='damped' in run() replaces the eigensolves by imaginary-time steps
A_C <- exp(-tau H_AC) A_C, C <- exp(-tau H_C) C (Lanczos), with tau adapted so that
the variational energy of the A_L state never increases.  This converges for
incommensurate spirals in commensurate cells, where plain VUMPS oscillates.
Also here: local expectation values, correlators, block entropies, and the
transfer-matrix spectrum of the cell.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh, eigs, gmres

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
I2 = np.eye(2)


# ----------------------------------------------------------------------------- models
def bond(J, field):
    """Identical to d2_classification/ed.py: h = sum_ab J_ab s^a (x) s^b
    + sum_a f_a (s^a (x) 1 + 1 (x) s^a)/2 (Pauli), left site = left factor."""
    J = np.asarray(J, float); field = np.asarray(field, float)
    h = np.zeros((4, 4), complex)
    for a in range(3):
        h += field[a] * 0.5 * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a]))
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
    return h


def mpo_from_term(h, r, d=2, tol=1e-13):
    """Upper-triangular MPO for H = sum_j h_{j..j+r-1}, by successive operator SVDs.
    h is (d^r, d^r) with the leftmost site the most significant kron factor."""
    T = np.asarray(h, complex).reshape([d] * (2 * r))
    # order (s1', s1, s2', s2, ...)
    perm = [x for i in range(r) for x in (i, r + i)]
    T = T.transpose(perm)
    ops = []      # ops[i] has shape (k_{i}, d, d, k_{i+1}), k_0 = k_r = 1
    left = 1
    M = T.reshape(1, -1)
    for i in range(r - 1):
        M = M.reshape(left * d * d, -1)
        U, S, Vh = linalg.svd(M, full_matrices=False)
        keep = int(np.sum(S > tol * S[0]))
        ops.append(U[:, :keep].reshape(left, d, d, keep))
        M = S[:keep, None] * Vh[:keep]
        left = keep
    ops.append(M.reshape(left, d, d, 1))
    chans = [o.shape[3] for o in ops[:-1]]
    D = 2 + sum(chans)
    W = np.zeros((D, D, d, d), complex)
    W[0, 0] = np.eye(d)
    W[D - 1, D - 1] = np.eye(d)
    offs = [0]
    o = 1
    for c in chans:
        offs.append(o); o += c
    offs.append(D - 1)
    for i, op in enumerate(ops):
        kl, _, _, kr = op.shape
        for a in range(kl):
            for b in range(kr):
                wl = 0 if i == 0 else offs[i] + a
                wr = D - 1 if i == r - 1 else offs[i + 1] + b
                W[wl, wr] += op[a, :, :, b]
    return W


# ----------------------------------------------------------------------------- contractions
def lt(L, A, W):
    """Left environment through one site: (bra, w, ket) -> (bra', w', ket')."""
    t = np.tensordot(L, A, ([2], [0]))              # a', w, s, b
    t = np.tensordot(t, W, ([1, 2], [0, 3]))        # a', b, w', s'
    t = np.tensordot(A.conj(), t, ([0, 1], [0, 3]))  # b', b, w'
    return t.transpose(0, 2, 1)


def rt(R, B, W):
    """Right environment through one site: (ket, w, bra) -> (ket', w', bra')."""
    t = np.tensordot(B, R, ([2], [0]))              # a, s, w', b'
    t = np.tensordot(t, W, ([1, 2], [3, 1]))        # a, b', w, s'
    t = np.tensordot(t, B.conj(), ([1, 3], [2, 1]))  # a, w, a'
    return t


def tl_id(x, A):
    """x (bra, ket) -> sum_s A^s^dagger x A^s."""
    t = np.tensordot(x, A, ([1], [0]))              # a', s, b
    return np.tensordot(A.conj(), t, ([0, 1], [0, 1]))


def tr_id(x, B):
    """x (ket, bra) -> sum_s B^s x B^s^dagger."""
    t = np.tensordot(B, x, ([2], [0]))              # a, s, b'
    return np.tensordot(t, B.conj(), ([1, 2], [1, 2]))


def hac(L, W, R, x):
    t = np.tensordot(L, x, ([2], [0]))              # a', w, s, b
    t = np.tensordot(t, W, ([1, 2], [0, 3]))        # a', b, w', s'
    t = np.tensordot(t, R, ([1, 2], [0, 1]))        # a', s', b'
    return t


def hc(L, R, x):
    t = np.tensordot(L, x, ([2], [0]))              # a', w, b
    return np.tensordot(t, R, ([1, 2], [1, 0]))     # a', b'


def hac2(L, W1, W2, R, x):
    t = np.tensordot(L, x, ([2], [0]))              # a', w, s1, s2, b
    t = np.tensordot(t, W1, ([1, 2], [0, 3]))       # a', s2, b, w1, s1'
    t = np.tensordot(t, W2, ([3, 1], [0, 3]))       # a', b, s1', w2, s2'
    t = np.tensordot(t, R, ([1, 3], [0, 1]))        # a', s1', s2', b'
    return t


def polar_u(M):
    U, _, Vh = linalg.svd(M, full_matrices=False)
    return U @ Vh


def lowest(op_apply, x0, tol):
    shp = x0.shape
    n = x0.size
    if n <= 400:
        H = np.zeros((n, n), complex)
        for i in range(n):
            e = np.zeros(n, complex); e[i] = 1
            H[:, i] = op_apply(e.reshape(shp)).ravel()
        w, U = linalg.eigh((H + H.conj().T) / 2)
        return w[0], U[:, 0].reshape(shp)
    op = LinearOperator((n, n), matvec=lambda v: op_apply(v.reshape(shp)).ravel(), dtype=complex)
    v0 = x0.ravel()
    if not np.any(v0):
        v0 = np.random.default_rng(0).normal(size=n) + 0j
    w, U = eigsh(op, k=1, which='SA', v0=v0, tol=tol, ncv=min(n - 1, 30), maxiter=2000)
    return w[0], U[:, 0].reshape(shp)


def expstep(op_apply, x0, tau, m=16):
    """Normalized exp(-tau H) x0 from an m-step Lanczos basis (H Hermitian).
    tau = None returns the lowest Ritz vector (the plain VUMPS update)."""
    shp = x0.shape
    v = x0.ravel() / np.linalg.norm(x0)
    m = min(m, v.size)
    Q = np.zeros((m, v.size), complex)
    a = np.zeros(m); b = np.zeros(m)
    Q[0] = v
    k = m
    for j in range(m):
        w = op_apply(Q[j].reshape(shp)).ravel()
        a[j] = np.vdot(Q[j], w).real
        w = w - a[j] * Q[j] - (b[j - 1] * Q[j - 1] if j > 0 else 0)
        w = w - Q[:j + 1].T @ (Q[:j + 1].conj() @ w)          # full reorthogonalization
        if j < m - 1:
            b[j] = np.linalg.norm(w)
            if b[j] < 1e-12:
                k = j + 1; break
            Q[j + 1] = w / b[j]
    T = np.diag(a[:k]) + np.diag(b[:k - 1], 1) + np.diag(b[:k - 1], -1)
    ev, U = linalg.eigh(T)
    if tau is None:
        c = U[:, 0]
    else:
        c = U @ (np.exp(-tau * (ev - ev[0])) * U[0].conj())
    y = Q[:k].T @ c
    return ev[0], (y / np.linalg.norm(y)).reshape(shp)


def energy_AL(AL, W):
    """Variational energy per site of the uMPS defined by A_L alone (exact right fixed point)."""
    n = len(AL)
    chi0 = AL[0].shape[0]

    def mvr(v):
        X = v.reshape(chi0, chi0)
        for k in range(n - 1, -1, -1):
            X = tr_id(X, AL[k])
        return X.ravel()
    N = chi0 * chi0
    if N <= 400:
        M = np.zeros((N, N), complex)
        for i in range(N):
            e = np.zeros(N, complex); e[i] = 1
            M[:, i] = mvr(e)
        w, U = linalg.eig(M)
        i = np.argmax(np.abs(w)); r = U[:, i]
    else:
        w, U = eigs(LinearOperator((N, N), matvec=mvr, dtype=complex), k=1, which='LM', tol=1e-12)
        r = U[:, 0]
    rho = r.reshape(chi0, chi0)
    rho = rho / np.trace(rho)
    rho = 0.5 * (rho + rho.conj().T)
    D = W[0].shape[0]
    Lb = np.zeros((chi0, D, chi0), complex)
    Lb[:, 0, :] = np.eye(chi0)
    for _ in range(int(np.ceil((D - 2) / n)) + 1):
        x = Lb
        for k in range(n):
            x = lt(x, AL[k], W[k])
        x[:, D - 1, :] = 0
        x[:, 0, :] = np.eye(chi0)
        Lb = x
    x = Lb
    for k in range(n):
        x = lt(x, AL[k], W[k])
    return float(np.trace(x[:, D - 1, :] @ rho).real) / n


# ----------------------------------------------------------------------------- uMPS
class UMPS:
    def __init__(self, AL, AR, C):
        self.AL, self.AR, self.C = AL, AR, C
        self.n = len(AL)
        self.AC = [np.tensordot(AL[k], C[(k + 1) % self.n], ([2], [0])) for k in range(self.n)]

    @staticmethod
    def random(n, d, chi, seed=0):
        rng = np.random.default_rng(seed)
        AL, AR, C = [], [], []
        for k in range(n):
            A = rng.normal(size=(chi, d, chi)) + 1j * rng.normal(size=(chi, d, chi))
            AL.append(np.linalg.qr(A.reshape(chi * d, chi))[0].reshape(chi, d, chi))
            AR.append(np.linalg.qr(A.reshape(chi, d * chi).T)[0].T.reshape(chi, d, chi).conj())
            c = np.diag(rng.uniform(0.5, 1, chi)) + 0j
            C.append(c / np.linalg.norm(c))
        return UMPS(AL, AR, C)


def left_envs(AL, C, W, tol, x0=None):
    """L[k] for k = 0..n-1 and the energy per cell."""
    n = len(AL)
    D = W[0].shape[0]
    chi0 = AL[0].shape[0]
    Lb = np.zeros((chi0, D, chi0), complex)
    Lb[:, 0, :] = np.eye(chi0)
    npass = int(np.ceil((D - 2) / n)) + 1
    for _ in range(npass):
        x = Lb
        for k in range(n):
            x = lt(x, AL[k], W[k])
        x[:, D - 1, :] = 0
        x[:, 0, :] = np.eye(chi0)
        Lb = x
    x = Lb
    for k in range(n):
        x = lt(x, AL[k], W[k])
    Y = x[:, D - 1, :]
    rho = C[0] @ C[0].conj().T                       # (ket, bra)
    e = np.trace(Y @ rho)
    rhs = (Y - e * np.eye(chi0)).ravel()

    def mv(v):
        X = v.reshape(chi0, chi0)
        Z = X
        for k in range(n):
            Z = tl_id(Z, AL[k])
        return (X - Z + np.trace(X @ rho) * np.eye(chi0)).ravel()
    op = LinearOperator((chi0 * chi0,) * 2, matvec=mv, dtype=complex)
    sol, info = gmres(op, rhs, x0=None if x0 is None else x0.ravel(), rtol=tol, atol=0, restart=40, maxiter=200)
    Lb[:, D - 1, :] = sol.reshape(chi0, chi0)
    Ls = [Lb]
    for k in range(n - 1):
        Ls.append(lt(Ls[-1], AL[k], W[k]))
    return Ls, e.real, info


def right_envs(AR, C, W, tol, x0=None):
    """R[k] for k = 0..n-1 (R[0] = R[n]) and the energy per cell."""
    n = len(AR)
    D = W[0].shape[0]
    chi0 = AR[0].shape[0]
    Rb = np.zeros((chi0, D, chi0), complex)
    Rb[:, D - 1, :] = np.eye(chi0)
    npass = int(np.ceil((D - 2) / n)) + 1
    for _ in range(npass):
        x = Rb
        for k in range(n - 1, -1, -1):
            x = rt(x, AR[k], W[k])
        x[:, 0, :] = 0
        x[:, D - 1, :] = np.eye(chi0)
        Rb = x
    x = Rb
    for k in range(n - 1, -1, -1):
        x = rt(x, AR[k], W[k])
    Y = x[:, 0, :]
    lam = C[0].conj().T @ C[0]                       # (bra, ket)
    e = np.trace(Y @ lam)
    rhs = (Y - e * np.eye(chi0)).ravel()

    def mv(v):
        X = v.reshape(chi0, chi0)
        Z = X
        for k in range(n - 1, -1, -1):
            Z = tr_id(Z, AR[k])
        return (X - Z + np.trace(X @ lam) * np.eye(chi0)).ravel()
    op = LinearOperator((chi0 * chi0,) * 2, matvec=mv, dtype=complex)
    sol, info = gmres(op, rhs, x0=None if x0 is None else x0.ravel(), rtol=tol, atol=0, restart=40, maxiter=200)
    Rb[:, 0, :] = sol.reshape(chi0, chi0)
    Rs = [None] * n
    Rs[0] = Rb
    x = Rb
    for k in range(n - 1, 0, -1):
        x = rt(x, AR[k], W[k])
        Rs[k] = x
    return Rs, e.real, info


def vumps_step(psi, W, tol_eig, tol_env, cache, tau=None, krylov=16):
    n = psi.n
    Ls, eL, _ = left_envs(psi.AL, psi.C, W, tol_env, cache.get('L'))
    Rs, eR, _ = right_envs(psi.AR, psi.C, W, tol_env, cache.get('R'))
    D = W[0].shape[0]
    cache['L'] = Ls[0][:, D - 1, :]
    cache['R'] = Rs[0][:, 0, :]
    AC, C = [], []
    for k in range(n):
        Rn = Rs[(k + 1) % n]
        if tau is None:
            _, ac = lowest(lambda x: hac(Ls[k], W[k], Rn, x), psi.AC[k], tol_eig)
            _, c = lowest(lambda x: hc(Ls[k], Rs[k], x), psi.C[k], tol_eig)
        else:
            _, ac = expstep(lambda x: hac(Ls[k], W[k], Rn, x), psi.AC[k], tau, krylov)
            _, c = expstep(lambda x: hc(Ls[k], Rs[k], x), psi.C[k], tau, krylov)
        AC.append(ac / np.linalg.norm(ac)); C.append(c / np.linalg.norm(c))
    AL, AR = [], []
    err = 0.0
    for k in range(n):
        chil, d, chir = AC[k].shape
        Uac = polar_u(AC[k].reshape(chil * d, chir))
        Uc = polar_u(C[(k + 1) % n])
        AL.append((Uac @ Uc.conj().T).reshape(chil, d, chir))
        Uac2 = polar_u(AC[k].reshape(chil, d * chir))
        Uc2 = polar_u(C[k])
        AR.append((Uc2.conj().T @ Uac2).reshape(chil, d, chir))
        x1 = np.tensordot(AL[k], C[(k + 1) % n], ([2], [0]))
        x2 = np.tensordot(C[k], AR[k], ([1], [0]))
        for x in (x1, x2):
            ov = np.vdot(x, AC[k])
            ph = ov / abs(ov) if abs(ov) > 0 else 1.0
            err = max(err, float(np.linalg.norm(AC[k] - ph * x)))
    new = UMPS(AL, AR, C)
    new.AC = AC
    return new, 0.5 * (eL + eR) / n, err, (Ls, Rs)


def expand(psi, W, envs, chi_new):
    """Two-site null-space subspace expansion to bond dimension <= chi_new at every bond."""
    n = psi.n
    Ls, Rs = envs
    d = psi.AL[0].shape[1]
    add = []   # for bond b = k+1 (between site k and k+1): (UL (chi_k d, delta), VR (delta, d chi_{k+2}))
    for k in range(n):
        k1 = (k + 1) % n
        A2 = np.tensordot(psi.AL[k], psi.AC[k1], ([2], [0]))          # chi_k d d chi_{k+2}
        H2 = hac2(Ls[k], W[k], W[k1], Rs[(k + 2) % n], A2)
        chik, _, _, chik2 = A2.shape
        chib = psi.AL[k].shape[2]
        delta = min(chi_new - chib, chik * d - chib, d * chik2 - chib)
        if delta <= 0:
            add.append(None); continue
        AL = psi.AL[k].reshape(chik * d, chib)
        NL = linalg.null_space(AL.conj().T)                           # (chik d, chik d - chib)
        AR = psi.AR[k1].reshape(chib, d * chik2)
        NR = linalg.null_space(AR).conj().T                            # (d chik2 - chib, d chik2)
        X = NL.conj().T @ H2.reshape(chik * d, d * chik2) @ NR.conj().T
        U, S, Vh = linalg.svd(X, full_matrices=False)
        delta = min(delta, int(np.sum(S > 1e-14 * max(S[0], 1e-300))))
        if delta <= 0:
            add.append(None); continue
        add.append((NL @ U[:, :delta], Vh[:delta] @ NR))
    chis = [psi.AL[k].shape[0] for k in range(n)]
    newchi = [chis[b] + (0 if add[(b - 1) % n] is None else add[(b - 1) % n][0].shape[1]) for b in range(n)]
    AL, AR, C = [], [], []
    for k in range(n):
        k1 = (k + 1) % n
        cl, cr = chis[k], chis[k1]
        nl, nr = newchi[k], newchi[k1]
        a = np.zeros((nl, d, nr), complex)
        a[:cl, :, :cr] = psi.AL[k]
        if add[k] is not None:
            a[:cl, :, cr:] = add[k][0].reshape(cl, d, -1)
        AL.append(a)
        b = np.zeros((nl, d, nr), complex)
        b[:cl, :, :cr] = psi.AR[k]
        if add[(k - 1) % n] is not None:
            b[cl:, :, :cr] = add[(k - 1) % n][1].reshape(-1, d, cr)
        AR.append(b)
        c = np.zeros((nl, nl), complex)
        c[:cl, :cl] = psi.C[k]
        C.append(c)
    new = UMPS(AL, AR, C)
    return new


# ----------------------------------------------------------------------------- observables
def entropies(psi):
    out = []
    for c in psi.C:
        s = linalg.svdvals(c)
        p = s ** 2; p = p / p.sum(); p = p[p > 1e-30]
        out.append(float(-(p * np.log(p)).sum()))
    return out


def schmidt(psi, k=0):
    return linalg.svdvals(psi.C[k])


def transfer_spectrum(psi, nev=8, tol=1e-10):
    """Eigenvalues of the A_L cell transfer matrix, largest modulus first."""
    AL = psi.AL
    chi0 = AL[0].shape[0]

    def mv(v):
        X = v.reshape(chi0, chi0)
        for A in AL:
            X = tl_id(X, A)
        return X.ravel()
    N = chi0 * chi0
    if N <= 600:
        M = np.zeros((N, N), complex)
        for i in range(N):
            e = np.zeros(N, complex); e[i] = 1
            M[:, i] = mv(e)
        w = linalg.eigvals(M)
    else:
        op = LinearOperator((N, N), matvec=mv, dtype=complex)
        w = eigs(op, k=min(nev, N - 2), which='LM', tol=tol, return_eigenvectors=False, maxiter=5000,
                 ncv=min(N - 1, max(2 * nev + 1, 40)))
    w = w[np.argsort(-np.abs(w))]
    return w[:nev]


def run(W, n, chis, d=2, tol=1e-9, maxit=300, seed=0, log=print, minit=5, psi=None, on_chi=None,
        mode='vumps', tau0=0.5, tau_max=50.0, krylov=16, stol=1e-7, errmax_s=1e-3):
    """Grow chi through the list `chis`, converging at each value.
    mode='vumps': plain parallel VUMPS (lowest eigenvectors).
    mode='damped': A_C <- exp(-tau H_AC) A_C, C <- exp(-tau H_C) C (imaginary-time TDVP steps),
    with tau adapted so that the variational energy of the A_L state never increases."""
    if psi is None:
        psi = UMPS.random(n, d, chis[0], seed)
    cache = {}
    results = []
    first = True
    tau = tau0
    for chi in chis:
        if not first:
            _, _, _, envs = vumps_step(psi, W, 1e-6, 1e-8, cache, tau=(None if mode == 'vumps' else 1e-9), krylov=4)
            psi = expand(psi, W, envs, chi)
            cache = {}
        first = False
        err = 1.0
        e_old = np.inf
        nrej = 0
        Shist = []
        for it in range(maxit):
            te = max(min(err * 1e-2, 1e-4), 1e-13)
            tenv = max(min(err * 1e-3, 1e-6), 1e-13)
            if mode == 'vumps':
                psi, e, err, envs = vumps_step(psi, W, te, tenv, cache)
            else:
                new, e, err_new, envs = vumps_step(psi, W, te, tenv, cache, tau=tau, krylov=krylov)
                e_new = energy_AL(new.AL, W)
                if e_new <= e_old + 1e-13 or it == 0:
                    psi, err, e_old = new, err_new, e_new
                    tau = min(tau * 1.25, tau_max)
                else:
                    nrej += 1
                    tau = max(tau * 0.5, 1e-3)
                    cache = {}
                e = e_old
            Shist.append(float(np.mean(entropies(psi))))
            if it >= minit and err < tol:
                break
            # stationarity of the entropy: the sliding mode of a spiral relaxes very slowly in the
            # gradient norm while S (and xi) have long converged
            if it >= max(minit, 20) and err < errmax_s and abs(Shist[-1] - Shist[-21]) < stol:
                break
        e_var = energy_AL(psi.AL, W)
        S = entropies(psi)
        lam = transfer_spectrum(psi)
        rec = dict(chi=chi, chis_actual=[a.shape[0] for a in psi.AL], e=float(e_var), err=float(err), iters=it + 1,
                   S=S, lam=[[float(z.real), float(z.imag)] for z in lam], mode=mode, tau=float(tau),
                   xi=float(-n / np.log(abs(lam[1]) / abs(lam[0]))))
        results.append(rec)
        log(f"n={n} chi={chi} e={e_var:.12f} err={err:.1e} it={it+1} tau={tau:.2g} S={np.round(S,5).tolist()} "
            f"xi={rec['xi']:.3f} |l1|={abs(lam[1]):.6f} arg/pi={np.angle(lam[1])/np.pi:.4f}", flush=True) if log is print else log(rec)
        if on_chi is not None:
            on_chi(psi, rec)
    return psi, results


def local_expect(psi, ops):
    """<O> on every site of the cell, from A_C[k]: array (len(ops), n)."""
    out = np.zeros((len(ops), psi.n))
    for k in range(psi.n):
        A = psi.AL[k]
        AC = np.tensordot(A, psi.C[(k + 1) % psi.n], ([2], [0]))
        nrm = np.vdot(AC, AC).real
        for i, O in enumerate(ops):
            out[i, k] = (np.vdot(AC, np.tensordot(O, AC, ([1], [1])).transpose(1, 0, 2)) / nrm).real
    return out


def correlator(psi, O1, O2, k0, rmax):
    """<O1_{k0} O2_{k0+r}> - <O1><O2> for r = 1..rmax (sites counted along the infinite chain)."""
    n = psi.n
    AL = psi.AL
    AC0 = np.tensordot(AL[k0], psi.C[(k0 + 1) % n], ([2], [0]))
    nrm = np.vdot(AC0, AC0).real
    loc = local_expect(psi, [O1, O2])
    # left part: (b', b) after applying O1 on site k0 in the A_C gauge, then propagate with A_R
    t = np.tensordot(O1, AC0, ([1], [1])).transpose(1, 0, 2)           # a, s', b
    Y = np.tensordot(AC0.conj(), t, ([0, 1], [0, 1])) / nrm            # b', b
    # A_R from A_C = C A_R  ->  use A_R = polar construction consistent with A_L, C
    out = []
    for r in range(1, rmax + 1):
        k = (k0 + r) % n
        B = psi.AR[k]
        tb = np.tensordot(Y, B, ([1], [0]))                              # b', s, c
        t2 = np.tensordot(tb, O2, ([1], [1]))                            # b', c, s'
        val = np.tensordot(B.conj(), t2, ([0, 1, 2], [0, 2, 1]))
        out.append(val - loc[0, k0] * loc[1, k])
        Y = np.tensordot(B.conj(), tb, ([0, 1], [0, 1]))
    return np.array(out)


def block_entropies(psi, k0, lmax, every=1):
    """Von Neumann entropy of a block of l consecutive sites (l = 1..lmax) of the infinite uMPS,
    starting at site k0 of the cell.  The block state is |Phi_ab> = (A_L ... A_L C)_ab, whose Gram
    matrix (chi^2 x chi^2) has the spectrum of the block density matrix."""
    n = psi.n
    chi0 = psi.AL[k0].shape[0]
    T = np.einsum('ac,bd->abcd', np.eye(chi0), np.eye(chi0)).astype(complex)   # (a, a', c, c')
    out = []
    for l in range(1, lmax + 1):
        A = psi.AL[(k0 + l - 1) % n]
        # T[a,a',c,c'] -> sum_s conj(A[c,s,e]) A[c',s,e'] T[a,a',c,c']
        t = np.tensordot(T, A.conj(), ([2], [0]))          # a a' c' s e
        t = np.tensordot(t, A, ([2, 3], [0, 1]))           # a a' e e'
        T = t
        if l % every:
            continue
        C = psi.C[(k0 + l) % n]
        G = np.tensordot(T, C.conj(), ([2], [0]))          # a a' c' b
        G = np.tensordot(G, C, ([2], [0]))                 # a a' b b'
        chi = G.shape[0]; chib = G.shape[2]
        G = G.transpose(0, 2, 1, 3).reshape(chi * chib, chi * chib)
        w = linalg.eigvalsh(0.5 * (G + G.conj().T))
        w = w[w > 1e-14]; w = w / w.sum()
        out.append((l, float(-(w * np.log(w)).sum())))
    return out
