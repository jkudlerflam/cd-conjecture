#!/usr/bin/env python3
"""Sparse-MPO back end for vumps.py, so that the literal big site (D = q^4 2^n, 324 for q = 3)
can be used directly.

An MPO tensor is an SpW object: shape = (D_left, D_right) and a list of terms (wl, wr, op) with op
a scipy.sparse (d x d) matrix acting as W[wl, wr, s', s] = op[s', s].  install() redirects the
contractions lt, rt, hac of vumps.py to sparse versions whenever they receive an SpW (so the
environment solvers and the VUMPS step of vumps.py are reused unchanged), and expand_sp() replaces
the two-site null-space expansion by its low-rank form: only MPO paths that are open across the
bond survive the null-space projections, so P_L H2 P_R = sum_m (P_L U_m)(V_m P_R) over the open
middle channels m, and its SVD needs QRs of (chi d) x (M chi) matrices instead of the
(chi d) x (d chi) two-site tensor.

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import numpy as np
import scipy.sparse as sp
from scipy import linalg
import umps_md as U
import vumps as V

_lt, _rt, _hac = V.lt, V.rt, V.hac


class SpW:
    def __init__(self, Dl, Dr, d, terms):
        self.shape = (Dl, Dr, d, d)
        self.terms = [(wl, wr, (None if op is None else sp.csr_matrix(op))) for (wl, wr, op) in terms]


def apply_phys(op, x):
    if op is None:
        return x
    a, d, b = x.shape
    y = op @ x.transpose(1, 0, 2).reshape(d, a * b)
    return np.asarray(y).reshape(d, a, b).transpose(1, 0, 2)


def lt(L, A, W):
    if not isinstance(W, SpW):
        return _lt(L, A, W)
    chi = A.shape[2]
    out = np.zeros((chi, W.shape[1], chi), complex)
    cache = {}
    for wl, wr, op in W.terms:
        if wl not in cache:
            cache[wl] = np.tensordot(L[:, wl, :], A, ([1], [0]))          # a', s, b
        t = apply_phys(op, cache[wl])
        out[:, wr, :] += np.tensordot(A.conj(), t, ([0, 1], [0, 1]))
    return out


def rt(R, B, W):
    if not isinstance(W, SpW):
        return _rt(R, B, W)
    chi = B.shape[0]
    out = np.zeros((chi, W.shape[0], chi), complex)
    cache = {}
    for wl, wr, op in W.terms:
        if wr not in cache:
            cache[wr] = np.tensordot(B, R[:, wr, :], ([2], [0]))           # a, s, b'
        t = apply_phys(op, cache[wr])
        out[:, wl, :] += np.tensordot(t, B.conj(), ([1, 2], [1, 2]))
    return out


def hac(L, W, R, x):
    if not isinstance(W, SpW):
        return _hac(L, W, R, x)
    y = np.zeros(x.shape, complex)
    Dr = W.shape[1]
    for wl, wr, op in W.terms:
        t = apply_phys(op, x)
        if wl != 0:                       # channel 0 of a left environment is the identity
            t = np.tensordot(L[:, wl, :], t, ([1], [0]))
        if wr != Dr - 1:                  # last channel of a right environment is the identity
            t = np.tensordot(t, R[:, wr, :], ([2], [0]))
        y += t
    return y


def install():
    V.lt, V.rt, V.hac = lt, rt, hac


def expand_sp(psi, W, envs, chi_new):
    """Low-rank null-space expansion; chi_new[b] = target for bond b."""
    n = psi.n
    Ls, Rs = envs
    add = []
    for k in range(n):
        k1 = (k + 1) % n
        AL, AC = psi.AL[k], psi.AC[k1]
        chik, dk, chib = AL.shape
        _, dk1, chik2 = AC.shape
        delta = min(chi_new[k1] - chib, chik * dk - chib, dk1 * chik2 - chib)
        if delta <= 0:
            add.append(None); continue
        Dmid = W[k].shape[1]
        mids = range(1, Dmid - 1)
        Lk, Rk2 = Ls[k], Rs[(k + 2) % n]
        Us, Vs = [], []
        for m in mids:
            u = np.zeros((chik, dk, chib), complex)
            for wl, wr, op in W[k].terms:
                if wr == m:
                    t = apply_phys(op, AL)
                    u += t if wl == 0 else np.tensordot(Lk[:, wl, :], t, ([1], [0]))
            v = np.zeros((chib, dk1, chik2), complex)
            for wl, wr, op in W[k1].terms:
                if wl == m:
                    t = apply_phys(op, AC)
                    v += t if wr == W[k1].shape[1] - 1 else np.tensordot(t, Rk2[:, wr, :], ([2], [0]))
            if np.any(u) and np.any(v):
                Us.append(u.reshape(chik * dk, chib)); Vs.append(v.reshape(chib, dk1 * chik2))
        if not Us:
            add.append(None); continue
        ALm = AL.reshape(chik * dk, chib)
        ARm = psi.AR[k1].reshape(chib, dk1 * chik2)
        Ut = np.hstack(Us); Ut = Ut - ALm @ (ALm.conj().T @ Ut)
        Vt = np.vstack(Vs); Vt = Vt - (Vt @ ARm.conj().T) @ ARm
        Q1, R1 = linalg.qr(Ut, mode='economic')
        Q2, R2 = linalg.qr(Vt.conj().T, mode='economic')
        u_, s_, vh_ = linalg.svd(R1 @ R2.conj().T)
        delta = min(delta, int(np.sum(s_ > 1e-14 * max(s_[0], 1e-300))))
        if delta <= 0:
            add.append(None); continue
        left = Q1 @ u_[:, :delta]
        right = (Q2 @ vh_[:delta].conj().T).conj().T
        # re-orthogonalize against the existing isometries (round-off)
        left = left - ALm @ (ALm.conj().T @ left); left = linalg.qr(left, mode='economic')[0]
        right = right - (right @ ARm.conj().T) @ ARm; right = linalg.qr(right.conj().T, mode='economic')[0].conj().T
        add.append((left, right))
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
    return V.UMPS(AL, AR, C)


def bigsite_spmpo(c, lam, ncells, bigsite_H, bigsite_ports):
    """Sparse MPO of the literal big-site chain: channels 0 (start), 1..r (F_a(pR) placed), r+1 (done)."""
    H = bigsite_H(c)
    FL, FR = bigsite_ports(c)
    r = len(FL)
    D = H.shape[0]
    Dw = r + 2
    terms = [(0, 0, None), (Dw - 1, Dw - 1, None), (0, Dw - 1, H)]
    for a in range(r):
        terms.append((0, 1 + a, lam * FR[a]))
        terms.append((1 + a, Dw - 1, FL[a]))
    return [D] * ncells, [SpW(Dw, Dw, D, terms) for _ in range(ncells)]


# ----------------------------------------------------------------------------- two-site (finite DMRG)
def from_dense(W, tol=1e-14):
    Dl, Dr, d, _ = W.shape
    terms = []
    for wl in range(Dl):
        for wr in range(Dr):
            blk = W[wl, wr]
            if np.max(np.abs(blk)) > tol:
                terms.append((wl, wr, blk))
    return SpW(Dl, Dr, d, terms)


def _apply_axis(op, x, axis):
    if op is None:
        return x
    xm = np.moveaxis(x, axis, 0)
    shp = xm.shape
    y = op @ xm.reshape(shp[0], -1)
    return np.moveaxis(np.asarray(y).reshape(shp), 0, axis)


def hac2(L, W1, W2, R, x):
    """Two-site effective Hamiltonian with sparse MPO tensors; x has shape (a, s1, s2, b).
    No identity shortcuts here (finite-DMRG environments are not normalized that way)."""
    y = np.zeros(x.shape, complex)
    Tm = {}
    Lx = {}
    for wl, m, op1 in W1.terms:
        if wl not in Lx:
            Lx[wl] = np.tensordot(L[:, wl, :], x, ([1], [0]))        # a', s1, s2, b
        t = _apply_axis(op1, Lx[wl], 1)
        Tm[m] = Tm[m] + t if m in Tm else t
    for m, wr, op2 in W2.terms:
        if m not in Tm:
            continue
        t = _apply_axis(op2, Tm[m], 2)
        y += np.tensordot(t, R[:, wr, :], ([3], [0]))
    return y


def lt_full(L, A, W):
    return lt(L, A, W)


def rt_full(R, B, W):
    return rt(R, B, W)


# ----------------------------------------------------------------------------- eigensolver
_lowest_orig = V.lowest


def lanczos_lowest(op_apply, x0, tol, maxm=30, restarts=6):
    """Lowest eigenpair by restarted Lanczos with full reorthogonalization (cheaper than ARPACK
    for the moderate accuracy VUMPS needs).  Stops when the residual norm is below tol."""
    shp = x0.shape
    n = x0.size
    if n <= 400:
        return _lowest_orig(op_apply, x0, tol)
    v = x0.ravel().astype(complex)
    if not np.any(v):
        v = np.random.default_rng(0).normal(size=n) + 0j
    v = v / np.linalg.norm(v)
    e = None
    for _ in range(restarts):
        m = min(maxm, n)
        Q = np.zeros((m, n), complex)
        a = np.zeros(m); b = np.zeros(m)
        Q[0] = v
        k = m
        for j in range(m):
            w = op_apply(Q[j].reshape(shp)).ravel()
            a[j] = np.vdot(Q[j], w).real
            w = w - Q[:j + 1].T @ (Q[:j + 1].conj() @ w)
            w = w - Q[:j + 1].T @ (Q[:j + 1].conj() @ w)
            if j < m - 1:
                b[j] = np.linalg.norm(w)
                if b[j] < 1e-14:
                    k = j + 1; break
                Q[j + 1] = w / b[j]
        T = np.diag(a[:k]) + np.diag(b[:k - 1], 1) + np.diag(b[:k - 1], -1)
        ev, Uv = linalg.eigh(T)
        y = Q[:k].T @ Uv[:, 0]
        y = y / np.linalg.norm(y)
        e = ev[0]
        r = op_apply(y.reshape(shp)).ravel() - e * y
        v = y
        if np.linalg.norm(r) < max(tol, 1e-13) or k < m:
            break
    return e, v.reshape(shp)


def install_lanczos():
    V.lowest = lanczos_lowest
