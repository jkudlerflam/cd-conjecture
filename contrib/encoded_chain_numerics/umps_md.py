#!/usr/bin/env python3
"""Multi-site VUMPS with site-dependent physical dimensions and MPO bond dimensions.

Thin layer over ../d2_followup_220_1282/vumps.py (whose contractions, environments and
eigensolvers are shape-agnostic): this file replaces the random initial state, the two-site
null-space expansion and the driver so that each site k may have its own d_k and each bond
its own chi target.  Conventions as in vumps.py.

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import os, sys
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from scipy import linalg
import vumps as V                                     # noqa: E402
from vumps import UMPS, vumps_step, energy_AL, entropies, transfer_spectrum, hac2   # noqa: E402


def random_md(dims, chis, seed=0):
    """chis[k] = bond dimension of bond k (left of site k)."""
    n = len(dims)
    rng = np.random.default_rng(seed)
    AL, AR, C = [], [], []
    for k in range(n):
        cl, d, cr = chis[k], dims[k], chis[(k + 1) % n]
        A = rng.normal(size=(cl, d, cr)) + 1j * rng.normal(size=(cl, d, cr))
        AL.append(np.linalg.qr(A.reshape(cl * d, cr))[0].reshape(cl, d, cr))
        AR.append(np.linalg.qr(A.reshape(cl, d * cr).T)[0].T.reshape(cl, d, cr).conj())
        c = np.diag(rng.uniform(0.5, 1, cl)) + 0j
        C.append(c / np.linalg.norm(c))
    return UMPS(AL, AR, C)


def feasible_chis(dims, target):
    """Clip bond targets so that chi_{k+1} <= d_k chi_k and chi_k <= d_k chi_{k+1} (cyclically)."""
    n = len(dims)
    ch = list(target)
    for _ in range(4 * n):
        for k in range(n):
            k1 = (k + 1) % n
            ch[k1] = min(ch[k1], dims[k] * ch[k])
            ch[k] = min(ch[k], dims[k] * ch[k1])
    return ch


def expand_md(psi, W, envs, chi_new):
    """Two-site null-space expansion; chi_new[b] is the target for bond b (left of site b)."""
    n = psi.n
    Ls, Rs = envs
    add = []
    for k in range(n):
        k1 = (k + 1) % n
        A2 = np.tensordot(psi.AL[k], psi.AC[k1], ([2], [0]))
        H2 = hac2(Ls[k], W[k], W[k1], Rs[(k + 2) % n], A2)
        chik, dk, dk1, chik2 = A2.shape
        chib = psi.AL[k].shape[2]
        delta = min(chi_new[k1] - chib, chik * dk - chib, dk1 * chik2 - chib)
        if delta <= 0:
            add.append(None); continue
        ALm = psi.AL[k].reshape(chik * dk, chib)
        NL = linalg.null_space(ALm.conj().T)
        ARm = psi.AR[k1].reshape(chib, dk1 * chik2)
        NR = linalg.null_space(ARm).conj().T
        X = NL.conj().T @ H2.reshape(chik * dk, dk1 * chik2) @ NR.conj().T
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
        d = psi.AL[k].shape[1]
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
    return UMPS(AL, AR, C)


def two_site_expect(psi, k, O):
    """<O> for an operator O on sites (k, k+1) of the cell, O of shape (d_k d_k1, d_k d_k1)."""
    n = psi.n
    k1 = (k + 1) % n
    th = np.tensordot(psi.AC[k], psi.AR[k1], ([2], [0]))     # a s1 s2 b
    a, d1, d2, b = th.shape
    t = th.transpose(1, 2, 0, 3).reshape(d1 * d2, a * b)
    return float(np.real(np.vdot(t, O @ t)) / np.real(np.vdot(t, t)))


def one_site_expect(psi, k, O):
    AC = psi.AC[k]
    t = np.tensordot(O, AC, ([1], [1])).transpose(1, 0, 2)
    return float(np.real(np.vdot(AC, t)) / np.real(np.vdot(AC, AC)))


def run_md(W, dims, chis, mults=None, tol=1e-9, maxit=300, seed=0, log=print, minit=5, psi=None,
           on_chi=None, mode='vumps', tau0=0.5, tau_max=50.0, krylov=16, stol=1e-7, errmax_s=1e-3,
           expand=None, eig_factor=1e-2, verbose=0):
    """Grow the boundary bond dimension through `chis`; bond b gets target mults[b]*chi.
    Logic copied from vumps.run (damped mode included)."""
    n = len(dims)
    mults = mults or [1] * n
    tgt = lambda chi: feasible_chis(dims, [int(m * chi) for m in mults])
    if psi is None:
        psi = random_md(dims, tgt(chis[0]), seed)
    cache = {}
    results = []
    first = True
    tau = tau0
    for chi in chis:
        if not first:
            _, _, _, envs = vumps_step(psi, W, 1e-6, 1e-8, cache, tau=(None if mode == 'vumps' else 1e-9), krylov=4)
            psi = (expand or expand_md)(psi, W, envs, tgt(chi))
            cache = {}
        first = False
        err = 1.0
        e_old = np.inf
        Shist = []
        for it in range(maxit):
            te = max(min(err * eig_factor, 1e-4), 1e-13)
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
                    tau = max(tau * 0.5, 1e-3)
                    cache = {}
                e = e_old
            Shist.append(float(np.mean(entropies(psi))))
            if verbose and it % verbose == 0:
                print(f'    it={it} e={e:.12f} err={err:.2e} S={Shist[-1]:.6f}', file=sys.stderr, flush=True)
            if it >= minit and err < tol:
                break
            if it >= max(minit, 20) and err < errmax_s and abs(Shist[-1] - Shist[-21]) < stol:
                break
        e_var = energy_AL(psi.AL, W)
        S = entropies(psi)
        lam = transfer_spectrum(psi)
        rec = dict(chi=chi, chis_actual=[a.shape[0] for a in psi.AL], e=float(e_var), err=float(err), iters=it + 1,
                   S=S, lam=[[float(z.real), float(z.imag)] for z in lam], mode=mode, tau=float(tau),
                   xi=float(-1.0 / np.log(abs(lam[1]) / abs(lam[0]))))   # in units of the unit cell
        results.append(rec)
        if on_chi is not None:
            on_chi(psi, rec)
        log(rec)
    return psi, results
