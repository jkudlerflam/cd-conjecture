#!/usr/bin/env python3
"""The full route-A chain: every cell is ONE site of dimension D = 64 carrying the exact cell
Hamiltonian H0 of cell.py, neighbours coupled by lam sum_a FR_a(j) FL_a(j+1).  No projection and no
perturbation theory:
    H = sum_j H0_j + lam sum_j sum_{a = X,Y,Z} FR_a(j) FL_a(j+1),      Delta = 1.

The logical hopping J0 lam^2 ~ 0.022 lam^2 is 1e-4 to 1e-3 of the cell gap, so every local
eigenproblem is solved by Davidson with the preconditioner (H0 on the physical leg + delta)^{-1}
(H0 is diagonal in the cell eigenbasis used here), as in contrib/routeA_numerics/bigsite.py, and the
MPO is applied term by term (sparse contractions and the low-rank null-space expansion follow
contrib/encoded_chain_numerics/sparse_mpo.py; both copied here, not imported).

Modes
  fes   VUMPS on the infinite chain with an n-cell unit cell over a chi list.
  ring  Lanczos on rings of L = 3, 4 cells (full Hilbert space D^L), lowest levels per parity sector,
        for comparison with the effective logical chain on the same rings.
usage: bigchain.py fes --target INC --lam 0.1 --n 2 --chis 16 24 32 48 64 96 --out data/big_INC_0.1.json
       bigchain.py ring --target INC --lam 0.1 --Ls 3 4 --nev 8 --out data/bigring_INC_0.1.json
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-27.
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh
import vumps as V
import cell as CE
from ed_check import eigbasis

_lt, _rt, _hac = V.lt, V.rt, V.hac
DELTA = 1e-3
EVALS = None
PREC = 'env'


# ----------------------------------------------------------------------------- sparse MPO
class SpW:
    def __init__(self, Dl, Dr, d, terms):
        self.shape = (Dl, Dr, d, d)
        self.terms = terms              # (wl, wr, op) with op a dense d x d matrix or None (identity)


def apply_phys(op, x):
    if op is None:
        return x
    a, d, b = x.shape
    y = op @ x.transpose(1, 0, 2).reshape(d, a * b)
    return y.reshape(d, a, b).transpose(1, 0, 2)


def lt(L, A, W):
    if not isinstance(W, SpW):
        return _lt(L, A, W)
    chi = A.shape[2]
    out = np.zeros((chi, W.shape[1], chi), complex)
    cache = {}
    for wl, wr, op in W.terms:
        if wl not in cache:
            cache[wl] = np.tensordot(L[:, wl, :], A, ([1], [0]))
        out[:, wr, :] += np.tensordot(A.conj(), apply_phys(op, cache[wl]), ([0, 1], [0, 1]))
    return out


def rt(R, B, W):
    if not isinstance(W, SpW):
        return _rt(R, B, W)
    chi = B.shape[0]
    out = np.zeros((chi, W.shape[0], chi), complex)
    cache = {}
    for wl, wr, op in W.terms:
        if wr not in cache:
            cache[wr] = np.tensordot(B, R[:, wr, :], ([2], [0]))
        out[:, wl, :] += np.tensordot(apply_phys(op, cache[wr]), B.conj(), ([1, 2], [1, 2]))
    return out


LAST = {}


def hac(L, W, R, x):
    if not isinstance(W, SpW):
        return _hac(L, W, R, x)
    LAST['L'], LAST['R'] = L, R
    y = np.zeros(x.shape, complex)
    Dr = W.shape[1]
    for wl, wr, op in W.terms:
        t = apply_phys(op, x)
        if wl != 0:
            t = np.tensordot(L[:, wl, :], t, ([1], [0]))
        if wr != Dr - 1:
            t = np.tensordot(t, R[:, wr, :], ([2], [0]))
        y += t
    return y


V.lt, V.rt, V.hac = lt, rt, hac


def expand_sp(psi, W, envs, chi_new):
    """Low-rank two-site null-space expansion (only MPO channels open across the bond contribute)."""
    n = psi.n
    Ls, Rs = envs
    add = []
    for k in range(n):
        k1 = (k + 1) % n
        AL, AC = psi.AL[k], psi.AC[k1]
        chik, dk, chib = AL.shape
        _, dk1, chik2 = AC.shape
        delta = min(chi_new - chib, chik * dk - chib, dk1 * chik2 - chib)
        if delta <= 0:
            add.append(None); continue
        Dmid = W[k].shape[1]
        Lk, Rk2 = Ls[k], Rs[(k + 2) % n]
        Us, Vs = [], []
        for m in range(1, Dmid - 1):
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
        a = np.zeros((nl, d, nr), complex); a[:cl, :, :cr] = psi.AL[k]
        if add[k] is not None:
            a[:cl, :, cr:] = add[k][0].reshape(cl, d, -1)
        AL.append(a)
        b = np.zeros((nl, d, nr), complex); b[:cl, :, :cr] = psi.AR[k]
        if add[(k - 1) % n] is not None:
            b[cl:, :, :cr] = add[(k - 1) % n][1].reshape(-1, d, cr)
        AR.append(b)
        c = np.zeros((nl, nl), complex); c[:cl, :cl] = psi.C[k]
        C.append(c)
    return V.UMPS(AL, AR, C)


# ----------------------------------------------------------------------------- eigensolver
def davidson(apply, x0, prec, tol=1e-10, maxit=400, mmax=40, keep=4):
    """Lowest eigenpair, preconditioned Davidson with thick restart and Olsen correction."""
    shp = x0.shape
    n = x0.size
    x = x0.ravel().astype(complex)
    if np.linalg.norm(x) < 1e-300:
        x = np.random.default_rng(0).normal(size=n) + 0j
    Vm = np.zeros((mmax + 1, n), complex)
    Am = np.zeros((mmax + 1, n), complex)
    Hs = np.zeros((mmax + 1, mmax + 1), complex)
    Vm[0] = x / np.linalg.norm(x)
    Am[0] = apply(Vm[0].reshape(shp)).ravel()
    Hs[0, 0] = np.vdot(Vm[0], Am[0])
    m = 1
    for it in range(maxit):
        w, s = linalg.eigh(0.5 * (Hs[:m, :m] + Hs[:m, :m].conj().T))
        theta = w[0]
        u = s[:, 0] @ Vm[:m]
        r = s[:, 0] @ Am[:m] - theta * u
        rn = np.linalg.norm(r)
        if rn < tol:
            break
        if m >= mmax:
            S = s[:, :keep]
            Vm[:keep] = S.T @ Vm[:m]; Am[:keep] = S.T @ Am[:m]
            Hs[:keep, :keep] = np.diag(w[:keep]); Hs[:keep, keep:] = 0; Hs[keep:, :keep] = 0
            m = keep
        Mr = prec(r.reshape(shp)).ravel()
        Mu = prec(u.reshape(shp)).ravel()
        t = Mr - (np.vdot(u, Mr) / np.vdot(u, Mu)) * Mu
        for _ in range(2):
            t -= Vm[:m].T @ (Vm[:m].conj() @ t)
        nt = np.linalg.norm(t)
        if nt < 1e-14:
            break
        Vm[m] = t / nt
        Am[m] = apply(Vm[m].reshape(shp)).ravel()
        Hs[:m + 1, m] = Vm[:m + 1].conj() @ Am[m]
        Hs[m, :m] = Hs[:m, m].conj()
        m += 1
    return theta, (u / np.linalg.norm(u)).reshape(shp), it, rn


STATS = dict(calls=0, iters=0, last_res=0.0)


def prec_phys(r):
    return r / (EVALS[None, :, None] + DELTA)


def env_prec(L, R, theta, phys=True):
    """(H_L + [H0] + H_R - theta + delta)^{-1}, positive: H_L = L[:, -1, :], H_R = R[:, 0, :]."""
    Dw = L.shape[1]
    AL_ = L[:, Dw - 1, :]; AL_ = 0.5 * (AL_ + AL_.conj().T)
    AR_ = R[:, 0, :]; AR_ = 0.5 * (AR_ + AR_.conj().T)
    l, UL = linalg.eigh(AL_)
    rr, VR = linalg.eigh(AR_)
    if phys:
        den = l[:, None, None] + EVALS[None, :, None] + rr[None, None, :] - theta
        ax = 2
    else:
        den = l[:, None] + rr[None, :] - theta
        ax = 1
    den = np.abs(den) + DELTA

    def p(x):
        t = np.tensordot(UL.conj().T, x, ([1], [0]))
        t = np.tensordot(t, VR, ([ax], [0]))
        t = t / den
        t = np.tensordot(UL, t, ([1], [0]))
        return np.tensordot(t, VR.conj().T, ([ax], [0]))
    return p


def pgd_step(psi, W, eta, cache, tol_env):
    """Preconditioned gradient (imaginary-time-like) update of every A_C and C, then the usual
    polar decompositions.  Replaces the eigensolves of VUMPS, which oscillate here."""
    n = psi.n
    Ls, eL, _ = V.left_envs(psi.AL, psi.C, W, tol_env, cache.get('L'))
    Rs, eR, _ = V.right_envs(psi.AR, psi.C, W, tol_env, cache.get('R'))
    Dw = W[0].shape[0]
    cache['L'] = Ls[0][:, Dw - 1, :]
    cache['R'] = Rs[0][:, 0, :]
    AC, C = [], []
    gn = 0.0
    for k in range(n):
        Rn = Rs[(k + 1) % n]
        x = psi.AC[k] / np.linalg.norm(psi.AC[k])
        y = V.hac(Ls[k], W[k], Rn, x)
        th = np.real(np.vdot(x, y))
        g = y - th * x
        gn = max(gn, float(np.linalg.norm(g)))
        ac = x - eta * env_prec(Ls[k], Rn, th, True)(g)
        AC.append(ac / np.linalg.norm(ac))
        c = psi.C[k] / np.linalg.norm(psi.C[k])
        yc = V.hc(Ls[k], Rs[k], c)
        thc = np.real(np.vdot(c, yc))
        gc = yc - thc * c
        cn = c - eta * env_prec(Ls[k], Rs[k], thc, False)(gc)
        C.append(cn / np.linalg.norm(cn))
    AL, AR = [], []
    err = 0.0
    for k in range(n):
        chil, d, chir = AC[k].shape
        Uac = V.polar_u(AC[k].reshape(chil * d, chir))
        Uc = V.polar_u(C[(k + 1) % n])
        AL.append((Uac @ Uc.conj().T).reshape(chil, d, chir))
        Uac2 = V.polar_u(AC[k].reshape(chil, d * chir))
        Uc2 = V.polar_u(C[k])
        AR.append((Uc2.conj().T @ Uac2).reshape(chil, d, chir))
        x1 = np.tensordot(AL[k], C[(k + 1) % n], ([2], [0]))
        x2 = np.tensordot(C[k], AR[k], ([1], [0]))
        for xx in (x1, x2):
            ov = np.vdot(xx, AC[k])
            ph = ov / abs(ov) if abs(ov) > 0 else 1.0
            err = max(err, float(np.linalg.norm(AC[k] - ph * xx)))
    new = V.UMPS(AL, AR, C)
    new.AC = AC
    return new, 0.5 * (eL + eR) / n, err, gn, (Ls, Rs)


def mix_step(psi, W, beta, cache, tol_eig, tol_env):
    """Under-relaxed VUMPS: A_C <- A_C + beta (A_C^eig - A_C), C likewise (phases aligned)."""
    n = psi.n
    Ls, eL, _ = V.left_envs(psi.AL, psi.C, W, tol_env, cache.get('L'))
    Rs, eR, _ = V.right_envs(psi.AR, psi.C, W, tol_env, cache.get('R'))
    Dw = W[0].shape[0]
    cache['L'] = Ls[0][:, Dw - 1, :]
    cache['R'] = Rs[0][:, 0, :]
    AC, C = [], []
    for k in range(n):
        Rn = Rs[(k + 1) % n]
        x = psi.AC[k] / np.linalg.norm(psi.AC[k])
        _, ac = V.lowest(lambda v: V.hac(Ls[k], W[k], Rn, v), x, tol_eig)
        ov = np.vdot(ac, x); ac = ac * (ov / abs(ov) if abs(ov) > 0 else 1)
        ac = x + beta * (ac - x)
        AC.append(ac / np.linalg.norm(ac))
        c0 = psi.C[k] / np.linalg.norm(psi.C[k])
        _, c = V.lowest(lambda v: V.hc(Ls[k], Rs[k], v), c0, tol_eig)
        ov = np.vdot(c, c0); c = c * (ov / abs(ov) if abs(ov) > 0 else 1)
        c = c0 + beta * (c - c0)
        C.append(c / np.linalg.norm(c))
    AL, AR = [], []
    err = 0.0
    for k in range(n):
        chil, d, chir = AC[k].shape
        Uac = V.polar_u(AC[k].reshape(chil * d, chir))
        Uc = V.polar_u(C[(k + 1) % n])
        AL.append((Uac @ Uc.conj().T).reshape(chil, d, chir))
        Uac2 = V.polar_u(AC[k].reshape(chil, d * chir))
        Uc2 = V.polar_u(C[k])
        AR.append((Uc2.conj().T @ Uac2).reshape(chil, d, chir))
        x1 = np.tensordot(AL[k], C[(k + 1) % n], ([2], [0]))
        x2 = np.tensordot(C[k], AR[k], ([1], [0]))
        for xx in (x1, x2):
            ov = np.vdot(xx, AC[k])
            ph = ov / abs(ov) if abs(ov) > 0 else 1.0
            err = max(err, float(np.linalg.norm(AC[k] - ph * xx)))
    new = V.UMPS(AL, AR, C)
    new.AC = AC
    return new, 0.5 * (eL + eR) / n, err, (Ls, Rs)


def make_prec_env(theta_ref):
    """(H_L (x) 1 (x) 1 + 1 (x) H0 (x) 1 + 1 (x) 1 (x) H_R - theta + delta)^{-1}, with H_L, H_R the
    energy blocks of the current environments (they carry the ~Delta energies of leakage directions
    in the bond spaces, which the physical-leg preconditioner misses).  Bond terms are ignored."""
    L, R = LAST['L'], LAST['R']
    Dw = L.shape[1]
    AL_ = L[:, Dw - 1, :]; AL_ = 0.5 * (AL_ + AL_.conj().T)
    AR_ = R[:, 0, :]; AR_ = 0.5 * (AR_ + AR_.conj().T)
    l, UL = linalg.eigh(AL_)
    rr, VR = linalg.eigh(AR_)
    den = l[:, None, None] + EVALS[None, :, None] + rr[None, None, :] - theta_ref
    den = np.where(np.abs(den) < DELTA, DELTA, np.abs(den) + DELTA)

    def p(x):
        t = np.tensordot(UL.conj().T, x, ([1], [0]))
        t = np.tensordot(t, VR, ([2], [0]))             # the right block acts as t -> t A_R, A_R = VR r VR^dag
        t = t / den
        t = np.tensordot(UL, t, ([1], [0]))
        return np.tensordot(t, VR.conj().T, ([2], [0]))
    return p


def lowest_patched(op_apply, x0, tol):
    n = x0.size
    if x0.ndim == 3 and x0.shape[1] == len(EVALS):
        y0 = op_apply(x0)
        th0 = np.real(np.vdot(x0, y0) / np.vdot(x0, x0))
        prec = make_prec_env(th0) if PREC == 'env' else prec_phys
        e, v, it, rn = davidson(op_apply, x0, prec, tol=max(tol, 1e-14))
        STATS['calls'] += 1; STATS['iters'] += it; STATS['last_res'] = float(rn)
        return e, v
    if n <= 4096:
        H = np.zeros((n, n), complex)
        for i in range(n):
            ee = np.zeros(n, complex); ee[i] = 1
            H[:, i] = op_apply(ee.reshape(x0.shape)).ravel()
        w, U = linalg.eigh((H + H.conj().T) / 2)
        return w[0], U[:, 0].reshape(x0.shape)
    e, v, it, rn = davidson(op_apply, x0, lambda r: r, tol=max(tol, 1e-12), maxit=2000)
    return e, v


V.lowest = lowest_patched


# ----------------------------------------------------------------------------- model
def mpo(cb, lam):
    D = cb['D']
    terms = [(0, 0, None), (4, 4, None), (0, 4, np.diag(cb['e']).astype(complex))]
    for a in range(3):
        terms.append((0, 1 + a, lam * cb['FR'][a]))
        terms.append((1 + a, 4, cb['FL'][a]))
    return SpW(5, 5, D, terms)


def bond_expect(psi, j, FR, FL):
    n = psi.n
    AC, AR = psi.AC[j], psi.AR[(j + 1) % n]
    nrm = np.real(np.vdot(AC, AC))
    tot = 0.0
    for fr, fl in zip(FR, FL):
        t = apply_phys(fr, AC)
        Xm = np.tensordot(AC.conj(), t, ([0, 1], [0, 1]))
        t2 = np.tensordot(Xm, apply_phys(fl, AR), ([1], [0]))
        tot += np.real(np.tensordot(AR.conj(), t2, ([0, 1, 2], [0, 1, 2])))
    return float(tot / nrm)


def run_fes(a, cb, cell):
    global DELTA
    D, N = cb['D'], cb['N']
    J = cell['pred']['J0'] * a.lam ** 2
    DELTA = a.delta_mult * J
    W = [mpo(cb, a.lam)] * a.n
    nops = []
    for x in cell['a']:
        O = np.zeros((D, D), complex); O[:N, :N] = x.conj().T @ x; nops.append(O)
    Pc = np.zeros((D, D), complex); Pc[:N, :N] = np.eye(N)
    data = dict(args=vars(a), J=J, pred=cell['pred'], results=[])
    t0 = time.time()
    # initial state: ground state of the second-order logical chain at chi_init (VUMPS on C^N),
    # embedded in the code space, then a two-site null-space expansion to chis[0] that adds the
    # bond-leakage directions lam FR (x) FL acting on the code (a code-only uMPS is a degenerate
    # stationary point of the single-site update)
    heff = CE.expected_h2(cell)
    heff = (heff - np.trace(heff) / N ** 2 * np.eye(N ** 2)) / cell['pred']['J0']
    chi_i = a.chi_init or max(2, a.chis[0] // (10 if a.init == 'dress' else 4))
    p0, _ = V.run([V.mpo_from_term(heff, 2, d=N)] * a.n, a.n, [chi_i], d=N, tol=1e-10, maxit=400,
                  mode='damped', log=lambda r_: None, seed=a.seed)
    if a.init == 'dress':
        import dressing as DR
        Wd = DR.dressing_mpo(cb, a.lam)
        As = [DR.apply_mpo(np.concatenate([x, np.zeros((x.shape[0], D - N, x.shape[2]))], 1), Wd) for x in p0.AL]
        ALd, ARd, Cd = DR.canonical(As)
        ALd, ARd, Cd = DR.truncate(ALd, ARd, Cd, a.chis[0])
        psi = V.UMPS(ALd, ARd, Cd)
        print(f"dressed start: chi={ALd[0].shape[0]} E={V.energy_AL(psi.AL, W):.12f} S={np.round(V.entropies(psi), 4).tolist()}", flush=True)
    else:
        AL = [np.concatenate([x, np.zeros((x.shape[0], D - N, x.shape[2]))], 1) for x in p0.AL]
        AR = [np.concatenate([x, np.zeros((x.shape[0], D - N, x.shape[2]))], 1) for x in p0.AR]
        psi = V.UMPS(AL, AR, p0.C)
        Ls, _, _ = V.left_envs(psi.AL, psi.C, W, 1e-12)
        Rs, _, _ = V.right_envs(psi.AR, psi.C, W, 1e-12)
        psi = expand_sp(psi, W, (Ls, Rs), a.chis[0])
    cache = {}
    first = True
    for chi in a.chis:
        if not first:
            Ls_, _, _ = V.left_envs(psi.AL, psi.C, W, 1e-10)
            Rs_, _, _ = V.right_envs(psi.AR, psi.C, W, 1e-10)
            psi = expand_sp(psi, W, (Ls_, Rs_), chi)
            cache = {}
        first = False
        err = 1.0
        Sh = []
        eta = a.eta0
        e_best = np.inf
        for it in range(a.maxit):
            te = max(min(err * 1e-2, 1e-3) * a.tol_scale * J, 1e-14)   # residual tolerance relative to J
            tenv = max(min(err * 1e-3, 1e-7), 1e-13)
            t1 = time.time()
            if a.opt == 'pgd':
                new, e, err_new, gnorm, envs = pgd_step(psi, W, eta, cache, tenv)
                e_new = V.energy_AL(new.AL, W)
                if e_new <= e_best + 1e-15 or it == 0:
                    psi, err, e_best = new, err_new, e_new
                    eta = min(eta * 1.2, a.eta_max)
                else:
                    eta = max(eta * 0.5, 1e-3)
                    cache = {}
                e = e_best
            elif a.opt == 'mix':
                new, e, err_new, envs = mix_step(psi, W, eta, cache, te, tenv)
                e_new = V.energy_AL(new.AL, W)
                if e_new <= e_best + 1e-15 or it == 0:
                    psi, err, e_best = new, err_new, e_new
                    eta = min(eta * 1.2, 1.0)
                else:
                    eta = max(eta * 0.5, 0.02)
                    cache = {}
                e = e_best
            else:
                psi, e, err, envs = V.vumps_step(psi, W, te, tenv, cache)
            Sh.append(float(np.mean(V.entropies(psi))))
            if a.verbose:
                print(f"   it={it} e={e:.14f} err={err:.2e} S={Sh[-1]:.5f} eta={eta:.3g} dav_it={STATS['iters']} dt={time.time()-t1:.1f}s",
                      flush=True)
            if it >= 5 and err < a.tol:
                break
            if it >= 30 and err < 1e-4 and abs(Sh[-1] - Sh[-21]) < a.stol:
                break
        e_var = V.energy_AL(psi.AL, W)
        lam_t = V.transfer_spectrum(psi)
        rec = dict(chi=chi, chis_actual=[x.shape[0] for x in psi.AL], e=float(e_var), err=float(err), iters=it + 1,
                   S=V.entropies(psi), lam=[[float(z.real), float(z.imag)] for z in lam_t],
                   xi=float(-a.n / np.log(abs(lam_t[1]) / abs(lam_t[0]))),
                   ebond=[a.lam * bond_expect(psi, j, cb['FR'], cb['FL']) for j in range(a.n)],
                   ecell=[float(np.real(np.vdot(psi.AC[j], apply_phys(np.diag(cb['e']).astype(complex), psi.AC[j])))
                                / np.real(np.vdot(psi.AC[j], psi.AC[j]))) for j in range(a.n)],
                   code=V.local_expect(psi, [Pc])[0].tolist(),
                   dens=V.local_expect(psi, nops).tolist(),
                   dav=dict(STATS), time_s=time.time() - t0)
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
        print(f"lam={a.lam} chi={chi} e={e_var:.14f} err={err:.1e} it={it+1} S={np.round(rec['S'], 5).tolist()} "
              f"xi={rec['xi']:.3f} eb={np.round(rec['ebond'], 10).tolist()} code={np.round(rec['code'], 5).tolist()} "
              f"dens={np.round(rec['dens'], 4).tolist()} dav={STATS} t={rec['time_s']:.0f}s", flush=True)


def run_ring(a, cb, cell):
    """Rings of L cells in the full Hilbert space D^L, each parity sector of U^{(x)L} separately.
    Block Davidson with the preconditioner (sum_j H0_j - w + delta)^{-1}, started from the lowest
    eigenvectors of the second-order logical ring embedded in the code space.  Momenta from the
    cyclic shift of cells.  The same rings are diagonalized for the second-order logical chain."""
    D, N = cb['D'], cb['N']
    udiag = np.rint(np.real(np.diag(cb['U'])))          # exact +-1 (the eigenbasis carries round-off)
    assert np.allclose(cb['U'], np.diag(udiag))
    h2p = CE.expected_h2(cell) * a.lam ** 2
    data = dict(args=vars(a), J=cell['pred']['J0'] * a.lam ** 2, results=[])
    k = a.nev
    for L in a.Ls:
        dim = D ** L
        e = cb['e']
        diag = np.zeros([D] * L)
        for j in range(L):
            sh = [1] * L; sh[j] = D
            diag = diag + e.reshape(sh)
        diag = diag.ravel()

        def mv(x):
            x = x.reshape([D] * L + [-1])
            y = diag.reshape([D] * L + [1]) * x
            for j in range(L):
                j1 = (j + 1) % L
                for c in range(3):
                    t = np.moveaxis(np.tensordot(cb['FR'][c], x, ([1], [j])), 0, j)
                    t = np.moveaxis(np.tensordot(cb['FL'][c], t, ([1], [j1])), 0, j1)
                    y += a.lam * t
            return y.reshape(dim, -1)

        def shift(x):
            x = x.reshape([D] * L + [-1])
            return np.moveaxis(x, 0, L - 1).reshape(dim, -1)
        par = np.ones([D] * L)
        for j in range(L):
            sh = [1] * L; sh[j] = D
            par = par * udiag.reshape(sh)
        par = par.ravel()
        # second-order logical ring (dense)
        hb = h2p.reshape(N, N, N, N)
        Hl = np.zeros((N ** L, N ** L), complex)
        Il = np.eye(N ** L).reshape([N] * L + [N ** L])
        for j in range(L):
            j1 = (j + 1) % L
            t = np.tensordot(hb, Il, ([2, 3], [j, j1]))
            Hl += np.moveaxis(t, [0, 1], [j, j1]).reshape(N ** L, N ** L)
        Hl = 0.5 * (Hl + Hl.conj().T)
        pl = np.ones([N] * L)
        for j in range(L):
            sh = [1] * L; sh[j] = N
            pl = pl * udiag[:N].reshape(sh)
        pl = pl.ravel()
        out = []
        for Q in (1, -1):
            sl = np.where(pl == Q)[0]
            wl, Ul = linalg.eigh(Hl[np.ix_(sl, sl)])
            sec = np.where(par == Q)[0]
            ds = len(sec)
            # embed the logical eigenvectors (code = first N cell states)
            Qm = np.zeros((dim, k), complex)
            idx = np.array(np.unravel_index(sl, [N] * L))
            flat = np.ravel_multi_index(tuple(idx), [D] * L)
            Qm[flat, :] = Ul[:, :k]
            Qm = Qm[sec]
            dsec = diag[sec]

            def mvs(v):
                full = np.zeros((dim, v.shape[1]), complex)
                full[sec] = v
                return mv(full)[sec]
            res = 1.0
            for it in range(200):
                HQ = mvs(Qm)
                Hs = Qm.conj().T @ HQ
                w, s = linalg.eigh(0.5 * (Hs + Hs.conj().T))
                Qm = Qm @ s; HQ = HQ @ s
                R = HQ - Qm * w[None, :]
                res = float(np.linalg.norm(R, axis=0).max())
                if res < a.ring_tol:
                    break
                T = R / np.where(dsec[:, None] > 0.5, dsec[:, None] - w[None, :], 1.0 + a.lam)
                for _ in range(2):
                    T = T - Qm @ (Qm.conj().T @ T)
                Qb, _ = np.linalg.qr(np.hstack([Qm, T]))
                HQb = mvs(Qb)
                Hs = Qb.conj().T @ HQb
                w2, s2 = linalg.eigh(0.5 * (Hs + Hs.conj().T))
                Qm = Qb @ s2[:, :k]
            full = np.zeros((dim, k), complex); full[sec] = Qm
            TQ = shift(full)[sec]
            Mt = Qm.conj().T @ TQ
            # momenta: diagonalize T inside clusters of (near) degenerate levels
            ks = np.zeros(k); i = 0
            while i < k:
                j = i
                while j + 1 < k and abs(w[j + 1] - w[i]) < 1e-7 * data['J']:
                    j += 1
                ks[i:j + 1] = np.sort(np.angle(linalg.eigvals(Mt[i:j + 1, i:j + 1])))
                i = j + 1
            code = np.real(np.einsum('ij,ij->j', Qm.conj()[np.isin(sec, flat)], Qm[np.isin(sec, flat)]))
            out.append(dict(Q=Q, E=w.tolist(), k=ks.tolist(), res=res, iters=it + 1, code_weight=code.tolist(),
                            E_log2=wl[:k].tolist()))
            print(f"ring L={L} Q={Q} res={res:.1e} it={it+1} E={np.round(w, 12).tolist()} k/pi={np.round(ks/np.pi, 3).tolist()} "
                  f"E2nd={np.round(wl[:k], 12).tolist()}", flush=True)
        data['results'].append(dict(L=L, sectors=out))
        json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)


def main():
    global EVALS
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['fes', 'ring'])
    ap.add_argument('--target', default='INC')
    ap.add_argument('--lam', type=float, default=0.1)
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64])
    ap.add_argument('--Ls', type=int, nargs='+', default=[3, 4])
    ap.add_argument('--nev', type=int, default=8)
    ap.add_argument('--maxit', type=int, default=300)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--stol', type=float, default=1e-6)
    ap.add_argument('--delta_mult', type=float, default=10.0)
    ap.add_argument('--tol_scale', type=float, default=1.0, help='eigensolver residual tol = err 1e-2 J tol_scale')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--ring_tol', type=float, default=1e-11)
    ap.add_argument('--prec', default='env', choices=['env', 'phys'])
    ap.add_argument('--init', default='dress', choices=['dress', 'expand'])
    ap.add_argument('--verbose', type=int, default=1)
    ap.add_argument('--opt', default='mix', choices=['pgd', 'vumps', 'mix'])
    ap.add_argument('--eta0', type=float, default=0.3)
    ap.add_argument('--eta_max', type=float, default=4.0)
    ap.add_argument('--chi_init', type=int, default=0, help='logical chi of the start (default chis[0]//4)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    par = getattr(CE, a.target)
    cell = CE.build(par, extra=6 if a.target == 'INC3' else 0)
    cb = eigbasis(cell)
    EVALS = cb['e'].copy()
    global PREC
    PREC = a.prec
    print(json.dumps(dict(target=a.target, lam=a.lam, J=cell['pred']['J0'] * a.lam ** 2, D=cb['D'])), flush=True)
    if a.mode == 'fes':
        run_fes(a, cb, cell)
    else:
        run_ring(a, cb, cell)


if __name__ == '__main__':
    main()
