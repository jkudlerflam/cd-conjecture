"""C196 (Claude Opus 5.5 agent, 2026-11-10): one-site VUMPS (Zauner-Stauber et al. 2018) for the nearest-neighbour
qubit chain with exact Z2 parity labels, numpy or cupy.

Why not iDMRG (PROOF 2): the growing-window iDMRG has no fixed point for a (nearly) number-conserving chain at
incommensurate filling, since each step must add a non-integer number of particles.  A uniform MPS has no such
constraint.

Labels.  One bond basis, even states first (n0 of them), D in total.  A_L[a, s, b] = 0 unless q_b = q_a ^ s; C is
block diagonal.  Environments use the C190 convention Le[a_bra, w, a_ket], Re[c_bra, u, c_ket], MPO channels
0 = start, 2 = end (z2dmrg.MPO).  The half-infinite Hamiltonians sit in Le[:, 2, :] and Re[:, 0, :], with the
energy per site subtracted (geometric sum solved by GMRES with the fixed-point projector).

Eigenproblem.  A_C is the lowest eigenvector of H_AC, found with C190's Davidson (z2dmrg_xp.davidson) and the
single-site analogue of C190's (s1, s2)-block Kronecker preconditioner (PROOF 3):
    H0 = H_L^(s) (x) 1 + 1 (x) H_R^(s) - h1(s),  H_L^(s) = sum_{w even} Le[:, w, :] W[w, end, s, s],
    H_R^(s) = sum_{u even} W[start, u, s, s] Re[:, u, :],
which keeps every term diagonal in s (fields, both ZZ couplings of the site, the two half-chain Hamiltonians) and
drops only the terms that flip s.  C is the lowest eigenvector of H_C (Davidson, preconditioned with the bond
analogue Le[:, 2, :] (x) 1 + 1 (x) Re[:, 0, :] + c_zz Z_L (x) Z_R in the diagonal of the ZZ-channel, see PROOF 3).
Update A_L = polar(A_C) polar(C)^dag, A_R = polar(C)^dag polar(A_C) (VUMPS).
Errors: eps_L = |A_C - A_L C|, eps_R = |A_C - C A_R| (the gradient norm of the registration is max of the two).
"""
import os
import sys
import time
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg as Z0          # noqa: E402
import z2dmrg_xp as X        # noqa: E402

bsl, bdim = X.bsl, X.bdim


def label_vec(n0, D):
    q = np.zeros(D, dtype=int)
    q[n0:] = 1
    return q


# ----------------------------------------------------------------------------------------------------------
# Linear algebra helpers
# ----------------------------------------------------------------------------------------------------------
def polar_tall(M):
    xm = X.amod(M)
    U, s, Vh = xm.linalg.svd(M, full_matrices=False)
    return U @ Vh


def gmres(op, b, x0, tol=1e-12, m=40, maxrestart=50):
    """Restarted GMRES on vectors of module of b.  Returns x, relative residual, matvecs."""
    xm = X.amod(b)
    nb = float(xm.linalg.norm(b))
    if nb == 0.0:
        return xm.zeros_like(b), 0.0, 0
    x = x0.copy() if x0 is not None else xm.zeros_like(b)
    nmv = 0
    rel = np.inf
    for rs in range(maxrestart):
        r = b - op(x)
        nmv += 1
        beta = float(xm.linalg.norm(r))
        rel = beta / nb
        if rel < tol:
            break
        V = xm.zeros((m + 1, b.size), dtype=b.dtype)
        H = np.zeros((m + 1, m), dtype=np.complex128)
        V[0] = r / beta
        k_used = m
        for j in range(m):
            w = op(V[j])
            nmv += 1
            h = V[:j + 1].conj() @ w
            w = w - V[:j + 1].T @ h
            h2 = V[:j + 1].conj() @ w
            w = w - V[:j + 1].T @ h2
            H[:j + 1, j] = X.asnp(h + h2)
            hn = float(xm.linalg.norm(w))
            H[j + 1, j] = hn
            # residual estimate from the small least-squares problem
            e1 = np.zeros(j + 2, dtype=np.complex128); e1[0] = beta
            y, *_ = np.linalg.lstsq(H[:j + 2, :j + 1], e1, rcond=None)
            res = np.linalg.norm(H[:j + 2, :j + 1] @ y - e1) / nb
            if hn < 1e-300 or res < tol:
                k_used = j + 1
                break
            V[j + 1] = w / hn
        e1 = np.zeros(k_used + 1, dtype=np.complex128); e1[0] = beta
        y, *_ = np.linalg.lstsq(H[:k_used + 1, :k_used], e1, rcond=None)
        x = x + V[:k_used].T @ X.to_mod(y, xm)
    r = b - op(x)
    rel = float(xm.linalg.norm(r)) / nb
    return x, rel, nmv


# ----------------------------------------------------------------------------------------------------------
# Transfer maps with an operator on the site
# ----------------------------------------------------------------------------------------------------------
def tl(A, Xm, O):
    """sum_{s', s} O[s', s] A_{s'}^dag X A_s  (X[a_bra, a_ket])."""
    xm = X.amod(A)
    XA = [Xm @ A[:, s, :] for s in (0, 1)]
    out = None
    for sp in (0, 1):
        acc = None
        for s in (0, 1):
            o = complex(O[sp, s])
            if o != 0:
                acc = o * XA[s] if acc is None else acc + o * XA[s]
        if acc is not None:
            t = A[:, sp, :].conj().T @ acc
            out = t if out is None else out + t
    return out if out is not None else xm.zeros_like(Xm)


def trr(B, Xm, O):
    """sum_{s', s} O[s', s] conj(B_{s'}) X B_s^T  (X[c_bra, c_ket])."""
    xm = X.amod(B)
    XB = [Xm @ B[:, s, :].T for s in (0, 1)]
    out = None
    for sp in (0, 1):
        acc = None
        for s in (0, 1):
            o = complex(O[sp, s])
            if o != 0:
                acc = o * XB[s] if acc is None else acc + o * XB[s]
        if acc is not None:
            t = B[:, sp, :].conj() @ acc
            out = t if out is None else out + t
    return out if out is not None else xm.zeros_like(Xm)


def left_env(AL, rho, Wh, L2_guess=None, tol=1e-12):
    """Infinite left environment of the uniform left-orthonormal AL.  rho = C C^dag (right fixed point of T_L).
    Returns Le (D, w, D), energy per site e, gmres residual and matvecs."""
    xm = X.amod(AL)
    D = AL.shape[0]
    w = Wh.shape[0]
    nz = lambda a, b: np.max(np.abs(Wh[a, b])) > 0
    I_ = xm.eye(D, dtype=np.complex128)
    L = {0: I_}
    order = [c for c in (0, 1, 3, 4, 5, 6, 2) if c < w]
    info = {}
    for b in order:
        if b == 0:
            continue
        Y = xm.zeros((D, D), dtype=np.complex128)
        for a in L:
            if a != b and nz(a, b):
                Y = Y + tl(AL, L[a], Wh[a, b])
        if b == 2:
            e = complex(xm.sum(Y * rho.T)).real
            Yp = Y - e * I_
            op = lambda v: (v.reshape(D, D) - tl(AL, v.reshape(D, D), np.eye(2))
                            + xm.sum(v.reshape(D, D) * rho.T) * I_).reshape(-1)
            x0 = L2_guess.reshape(-1) if L2_guess is not None else None
            x, rel, nmv = gmres(op, Yp.reshape(-1), x0, tol=tol)
            L[2] = x.reshape(D, D)
            info = dict(e=e, gm_res=rel, gm_mv=nmv)
        else:
            L[b] = Y
    Le = xm.zeros((D, w, D), dtype=np.complex128)
    for b in L:
        Le[:, b, :] = L[b]
    return Le, info


def right_env(AR, rhoR, Wh, R0_guess=None, tol=1e-12):
    """Infinite right environment of the uniform right-orthonormal AR.  rhoR = C^dag C."""
    xm = X.amod(AR)
    D = AR.shape[0]
    w = Wh.shape[0]
    nz = lambda a, b: np.max(np.abs(Wh[a, b])) > 0
    I_ = xm.eye(D, dtype=np.complex128)
    R = {2: I_}
    order = [c for c in (2, 1, 5, 6, 3, 4, 0) if c < w]
    info = {}
    for b in order:
        if b == 2:
            continue
        Y = xm.zeros((D, D), dtype=np.complex128)
        for v in list(R.keys()):
            if v != b and nz(b, v):
                Y = Y + trr(AR, R[v], Wh[b, v])
        if b == 0:
            e = complex(xm.sum(Y * rhoR.T)).real
            Yp = Y - e * I_
            op = lambda v: (v.reshape(D, D) - trr(AR, v.reshape(D, D), np.eye(2))
                            + xm.sum(v.reshape(D, D) * rhoR.T) * I_).reshape(-1)
            x0 = R0_guess.reshape(-1) if R0_guess is not None else None
            x, rel, nmv = gmres(op, Yp.reshape(-1), x0, tol=tol)
            R[0] = x.reshape(D, D)
            info = dict(e=e, gm_res=rel, gm_mv=nmv)
        else:
            R[b] = Y
    Re = xm.zeros((D, w, D), dtype=np.complex128)
    for b in R:
        Re[:, b, :] = R[b]
    return Re, info


# ----------------------------------------------------------------------------------------------------------
# One-site and bond effective Hamiltonians on packed Z2 blocks
# ----------------------------------------------------------------------------------------------------------
class OneSite:
    def __init__(self, Le, W, Re, n0, D, Wh):
        self.xm = X.amod(Le)
        self.Le, self.W, self.Re, self.n0, self.D, self.Wh = Le, W, Re, n0, D, Wh
        self.blocks = []
        off = 0
        for qa in (0, 1):
            for s in (0, 1):
                qc = qa ^ s
                shp = (bdim(n0, D, qa), bdim(n0, D, qc))
                n = shp[0] * shp[1]
                self.blocks.append((qa, s, qc, shp, off, n))
                off += n
        self.size = off

    def pack(self, A):
        out = self.xm.empty(self.size, dtype=np.complex128)
        for qa, s, qc, shp, off, n in self.blocks:
            if n:
                out[off:off + n] = A[bsl(self.n0, self.D, qa), s, bsl(self.n0, self.D, qc)].reshape(-1)
        return out

    def unpack(self, x):
        A = self.xm.zeros((self.D, 2, self.D), dtype=np.complex128)
        for qa, s, qc, shp, off, n in self.blocks:
            if n:
                A[bsl(self.n0, self.D, qa), s, bsl(self.n0, self.D, qc)] = x[off:off + n].reshape(shp)
        return A

    def matvec(self, x):
        xm = self.xm
        A = self.unpack(x)
        T = xm.tensordot(self.Le, A, axes=([2], [0]))             # (a', w, s, c)
        T = xm.tensordot(T, self.W, axes=([1, 2], [0, 3]))         # (a', c, u, s')
        T = xm.tensordot(T, self.Re, axes=([1, 2], [2, 1]))       # (a', s', c')
        return self.pack(T)

    def build_prec(self):
        xm = self.xm
        Le, Re, Wh, n0, D = self.Le, self.Re, self.Wh, self.n0, self.D
        ev = [0, 1, 2]
        self.kL, self.kR = {}, {}
        for s in (0, 1):
            HL = sum(Le[:, w, :] * complex(Wh[w, 2, s, s]) for w in ev)
            HR = sum(Re[:, u, :] * complex(Wh[0, u, s, s]) for u in ev)
            for q in (0, 1):
                sl = bsl(n0, D, q)
                if bdim(n0, D, q) == 0:
                    continue
                H = HL[sl, sl]
                lam, U = xm.linalg.eigh((H + H.conj().T) / 2)
                self.kL[(s, q)] = (lam, U, float(lam[0]))
                H = HR[sl, sl]
                mu, V = xm.linalg.eigh((H + H.conj().T) / 2)
                self.kR[(s, q)] = (mu, V, float(mu[0]))
        self.kc = {s: -float(np.real(Wh[0, 2, s, s])) for s in (0, 1)}
        e0 = np.inf
        for qa, s, qc, shp, off, n in self.blocks:
            if n:
                e0 = min(e0, self.kL[(s, qa)][2] + self.kR[(s, qc)][2] + self.kc[s])
        self.e0 = e0

    def prec(self, r, theta, sigma):
        xm = self.xm
        out = xm.empty_like(r)
        sh = min(theta, self.e0) - sigma
        for qa, s, qc, shp, off, n in self.blocks:
            if n == 0:
                continue
            lam, U, _ = self.kL[(s, qa)]
            mu, V, _ = self.kR[(s, qc)]
            Xb = U.conj().T @ r[off:off + n].reshape(shp) @ V.conj()
            Xb = Xb / (lam[:, None] + mu[None, :] + (self.kc[s] - sh))
            out[off:off + n] = (U @ Xb @ V.T).reshape(-1)
        return out


class Bond:
    """H_C on block-diagonal C (q_a = q_c)."""

    def __init__(self, Le, Re, n0, D):
        self.xm = X.amod(Le)
        self.Le, self.Re, self.n0, self.D = Le, Re, n0, D
        self.blocks = []
        off = 0
        for q in (0, 1):
            n = bdim(n0, D, q) ** 2
            self.blocks.append((q, (bdim(n0, D, q),) * 2, off, n))
            off += n
        self.size = off

    def pack(self, C):
        out = self.xm.empty(self.size, dtype=np.complex128)
        for q, shp, off, n in self.blocks:
            if n:
                sl = bsl(self.n0, self.D, q)
                out[off:off + n] = C[sl, sl].reshape(-1)
        return out

    def unpack(self, x):
        C = self.xm.zeros((self.D, self.D), dtype=np.complex128)
        for q, shp, off, n in self.blocks:
            if n:
                sl = bsl(self.n0, self.D, q)
                C[sl, sl] = x[off:off + n].reshape(shp)
        return C

    def matvec(self, x):
        xm = self.xm
        C = self.unpack(x)
        T = xm.tensordot(self.Le, C, axes=([2], [0]))              # (a', w, c)
        T = xm.tensordot(T, self.Re, axes=([1, 2], [1, 2]))        # (a', c')
        return self.pack(T)

    def build_prec(self):
        xm = self.xm
        Le, Re, n0, D = self.Le, self.Re, self.n0, self.D
        self.kL, self.kR = {}, {}
        for q in (0, 1):
            sl = bsl(n0, D, q)
            if bdim(n0, D, q) == 0:
                continue
            H = Le[sl, 2, sl]
            lam, U = xm.linalg.eigh((H + H.conj().T) / 2)
            self.kL[q] = (lam, U)
            H = Re[sl, 0, sl]
            mu, V = xm.linalg.eigh((H + H.conj().T) / 2)
            self.kR[q] = (mu, V)
        self.e0 = min(float(self.kL[q][0][0] + self.kR[q][0][0]) for q in self.kL)

    def prec(self, r, theta, sigma):
        xm = self.xm
        out = xm.empty_like(r)
        sh = min(theta, self.e0) - sigma
        for q, shp, off, n in self.blocks:
            if n == 0:
                continue
            lam, U = self.kL[q]
            mu, V = self.kR[q]
            Xb = U.conj().T @ r[off:off + n].reshape(shp) @ V.conj()
            Xb = Xb / (lam[:, None] + mu[None, :] - sh)
            out[off:off + n] = (U @ Xb @ V.T).reshape(-1)
        return out


def AL_old_proj(AL, AC):
    """C = sum_s A_L(s)^dag A_C(s)."""
    return AL[:, 0, :].conj().T @ AC[:, 0, :] + AL[:, 1, :].conj().T @ AC[:, 1, :]


# ----------------------------------------------------------------------------------------------------------
# VUMPS state
# ----------------------------------------------------------------------------------------------------------
class VUMPS:
    def __init__(self, Wh, mpo, D=16, seed=0, cmode='eig'):
        self.cmode = cmode
        self.damp = float(os.environ.get('C196_DAMP', '1.0'))
        self.Wh = np.ascontiguousarray(Wh.astype(np.complex128))
        self.mpo = mpo
        self.rng = np.random.default_rng(seed)
        self.xm = np
        self.n0, self.D = D // 2, D
        AC = self._random_block_tensor(self.n0, D)
        C = np.diag(np.linspace(1.0, 0.1, D)).astype(np.complex128)
        self.AC, self.C = AC / np.linalg.norm(AC), C / np.linalg.norm(C)
        self.AL, self.AR = self._gauge(self.AC, self.C)
        self.Le = self.Re = None
        self.it = 0

    def _random_block_tensor(self, n0, D, scale=1.0):
        A = np.zeros((D, 2, D), dtype=np.complex128)
        q = label_vec(n0, D)
        for s in (0, 1):
            mask = (q[:, None] ^ s) == q[None, :]
            R = self.rng.standard_normal((D, D)) + 1j * self.rng.standard_normal((D, D))
            A[:, s, :] = R * mask * scale
        return A

    def to(self, xm):
        if xm is self.xm:
            return
        f = lambda a: X.to_mod(a, xm) if a is not None else None
        self.AC, self.C, self.AL, self.AR = f(self.AC), f(self.C), f(self.AL), f(self.AR)
        self.Le, self.Re = f(self.Le), f(self.Re)
        self.xm = xm

    def _gauge(self, AC, C):
        D = C.shape[0]
        UC = polar_tall(C)
        AL = (polar_tall(AC.reshape(D * 2, D)) @ UC.conj().T).reshape(D, 2, D)
        URt = polar_tall(AC.reshape(D, 2 * D).T)                     # tall (2D x D), transpose of the wide polar
        AR = (UC.conj().T @ URt.T).reshape(D, 2, D)
        return AL, AR

    def errors(self):
        xm = self.xm
        epsL = float(xm.linalg.norm(self.AC - xm.tensordot(self.AL, self.C, axes=([2], [0]))))
        epsR = float(xm.linalg.norm(self.AC - xm.tensordot(self.C, self.AR, axes=([1], [0]))))
        return epsL, epsR

    def step(self, tol_eig=None, sigma=1e-2, maxit=300, env_tol=1e-12):
        """One VUMPS iteration.  Returns a record."""
        xm = self.xm
        D, n0 = self.D, self.n0
        Wd = X.to_mod(self.Wh, xm)
        rhoL = self.C @ self.C.conj().T
        rhoR = (self.C.conj().T @ self.C).T          # right_env uses sum(Y * rhoR.T) = tr(Y C^dag C)-form, PROOF 3
        t0 = time.time()
        self.Le, il = left_env(self.AL, rhoL, self.Wh,
                               None if self.Le is None else self.Le[:, 2, :], tol=env_tol)
        self.Re, ir = right_env(self.AR, rhoR, self.Wh,
                                None if self.Re is None else self.Re[:, 0, :], tol=env_tol)
        t_env = time.time() - t0
        epsL0, epsR0 = self.errors()
        eps0 = max(epsL0, epsR0)
        tol = tol_eig if tol_eig is not None else max(1e-11, min(1e-6, eps0 * 1e-2))
        if self.cmode == 'two':
            return self._step_two(Wd, tol, sigma, maxit, eps0, il, ir, t_env)
        # A_C
        os_ = OneSite(self.Le, Wd, self.Re, n0, D, self.Wh)
        x0 = os_.pack(self.AC)
        x0 = x0 / xm.linalg.norm(x0)
        os_.build_prec()
        eAC, vAC, nmvA, rA = X.davidson(os_.matvec, x0, None, res_tol=tol, maxit=maxit, mmax=32, keep=6,
                                         prec=lambda r, th: os_.prec(r, th, sigma))
        AC = os_.unpack(vAC)
        # C
        if self.cmode == 'eig':
            bd = Bond(self.Le, self.Re, n0, D)
            y0 = bd.pack(self.C)
            y0 = y0 / xm.linalg.norm(y0)
            AL = self.AL

            def cprec(r, th):
                # lift with A_L, apply the one-site preconditioner, project back (PROOF 3)
                Cr = bd.unpack(r)
                R = xm.tensordot(AL, Cr, axes=([2], [0]))
                Pr = os_.unpack(os_.prec(os_.pack(R), th, sigma))
                return bd.pack(xm.tensordot(AL.conj(), Pr, axes=([0, 1], [0, 1])))
            eC, vC, nmvC, rC = X.davidson(bd.matvec, y0, None, res_tol=tol, maxit=maxit, mmax=32, keep=6,
                                           prec=cprec)
            C = bd.unpack(vC)
        else:
            # projection of the new A_C on the old A_L (equal to the H_C eigenvector at the fixed point, PROOF 3)
            C = AL_old_proj(self.AL, AC)
            C = C / xm.linalg.norm(C)
            eC, nmvC, rC = 0.0, 0, 0.0
        # fix the phase of C to that of the previous C (A_L, A_R are insensitive, eps is not)
        ph = complex(xm.vdot(self.C, C))
        if abs(ph) > 0:
            C = C * (abs(ph) / ph)
        ph = complex(xm.vdot(self.AC, AC))
        if abs(ph) > 0:
            AC = AC * (abs(ph) / ph)
        if self.damp < 1.0:
            AC = self.AC + self.damp * (AC - self.AC)
            C = self.C + self.damp * (C - self.C)
            AC, C = AC / xm.linalg.norm(AC), C / xm.linalg.norm(C)
        self.AC, self.C = AC, C
        self.AL, self.AR = self._gauge(AC, C)
        epsL, epsR = self.errors()
        self.it += 1
        dens = float(X.asnp(xm.sum(xm.abs(self.AC[:, 1, :]) ** 2)).real)
        return dict(it=self.it, e=il['e'], eR=ir['e'], eps_prev=eps0, epsL=epsL, epsR=epsR, eps=max(epsL, epsR),
                    tol=tol, nmvA=int(nmvA), rA=float(rA), nmvC=int(nmvC), rC=float(rC), eAC=float(eAC),
                    eC=float(eC), gmL=il['gm_res'], gmR=ir['gm_res'], gmvL=il['gm_mv'], gmvR=ir['gm_mv'],
                    n=dens, S=self.entropy(), t_env=t_env)

    def _step_two(self, Wd, tol, sigma, maxit, eps0, il, ir, t_env):
        """Two-site centre: theta = lowest eigenvector of H_AC2 (C190 TwoSite, kron4 preconditioner), started from
        A_L C A_R; A_C = theta . A_R^dag, C = A_L^dag . theta . A_R^dag."""
        xm = self.xm
        D, n0 = self.D, self.n0
        ts = X.TwoSite(self.Le, Wd, Wd, self.Re, n0, D, n0, D, self.mpo)
        th0 = xm.tensordot(xm.tensordot(self.AL, self.C, axes=([2], [0])), self.AR, axes=([2], [0]))
        x0 = ts.pack_dense(th0)
        x0 = x0 / xm.linalg.norm(x0)
        ts.build_kron4()
        e2, v, nmv, r = X.davidson(ts.matvec, x0, None, res_tol=tol, maxit=maxit, mmax=32, keep=6,
                                   prec=lambda rr, thv: ts.precond4(rr, thv, sigma))
        th = ts.unpack_dense(v)
        AC = xm.tensordot(th, self.AR.conj(), axes=([2, 3], [1, 2]))          # (a, s1, b)
        C = xm.tensordot(self.AL.conj(), AC, axes=([0, 1], [0, 1]))           # (b', b)
        nC = xm.linalg.norm(C)
        AC, C = AC / nC, C / nC
        ph = complex(xm.vdot(self.C, C))
        if abs(ph) > 0:
            C, AC = C * (abs(ph) / ph), AC * (abs(ph) / ph)
        if self.damp < 1.0:
            AC = self.AC + self.damp * (AC - self.AC)
            C = self.C + self.damp * (C - self.C)
            AC, C = AC / xm.linalg.norm(AC), C / xm.linalg.norm(C)
        self.AC, self.C = AC, C
        self.AL, self.AR = self._gauge(AC, C)
        epsL, epsR = self.errors()
        self.it += 1
        dens = float(X.asnp(xm.sum(xm.abs(self.AC[:, 1, :]) ** 2)).real) / float(xm.linalg.norm(self.AC)) ** 2
        return dict(it=self.it, e=il['e'], eR=ir['e'], eps_prev=eps0, epsL=epsL, epsR=epsR, eps=max(epsL, epsR),
                    tol=tol, nmvA=int(nmv), rA=float(r), nmvC=0, rC=0.0, eAC=float(e2), eC=0.0,
                    gmL=il['gm_res'], gmR=ir['gm_res'], gmvL=il['gm_mv'], gmvR=ir['gm_mv'], n=dens, S=self.entropy(),
                    t_env=t_env)

    def entropy(self):
        s = np.linalg.svd(X.asnp(self.C), compute_uv=False)
        return Z0.entropy(s ** 2)

    def schmidt(self):
        return np.linalg.svd(X.asnp(self.C), compute_uv=False)

    # ------------------------------------------------------------------------------------------ growth
    def expand(self, Dnew, scale=1e-3):
        """Pad to Dnew (new states split equally between the parities) with small random entries, then regauge."""
        xm = self.xm
        D, n0 = self.D, self.n0
        add = Dnew - D
        if add <= 0:
            return
        a0 = add // 2
        n0n = n0 + a0
        qn = label_vec(n0n, Dnew)
        idx_old = np.r_[np.arange(n0), n0n + np.arange(D - n0)]   # positions of the old states
        AC = X.asnp(self.AC)
        C = X.asnp(self.C)
        ACn = self._random_block_tensor(n0n, Dnew, scale=scale * np.abs(AC).max())
        ACn[np.ix_(idx_old, [0, 1], idx_old)] = AC
        smin = np.linalg.svd(C, compute_uv=False).min()
        Cn = np.zeros((Dnew, Dnew), dtype=np.complex128)
        new = np.setdiff1d(np.arange(Dnew), idx_old)
        Cn[np.ix_(idx_old, idx_old)] = C
        Cn[new, new] = smin * 0.1 * (1 + 0.1 * self.rng.random(new.size))
        ACn /= np.linalg.norm(ACn)
        Cn /= np.linalg.norm(Cn)
        self.D, self.n0 = Dnew, n0n
        self.AC, self.C = X.to_mod(ACn, xm), X.to_mod(Cn, xm)
        self.AL, self.AR = self._gauge(self.AC, self.C)
        self.Le = self.Re = None

    # ------------------------------------------------------------------------------------------ transfer matrix
    def transfer_spectrum(self, k=6, m=60, tol=1e-10, seed=1):
        """Leading |eigenvalues| of T(X) = sum_s A_L^dag X A_L in the even and odd sectors (one-site cell)."""
        import idmrg196 as I
        xm = self.xm
        AL = self.AL
        n0, D = self.n0, self.D
        out = {}
        for r in (0, 1):
            blocks = [(q, q ^ r) for q in (0, 1) if bdim(n0, D, q) and bdim(n0, D, q ^ r)]
            shapes = [(bdim(n0, D, a), bdim(n0, D, b)) for a, b in blocks]
            sizes = [a * b for a, b in shapes]
            ntot = int(sum(sizes))

            def embed(v):
                Xf = xm.zeros((D, D), dtype=np.complex128)
                off = 0
                for (a, b), shp, sz in zip(blocks, shapes, sizes):
                    Xf[bsl(n0, D, a), bsl(n0, D, b)] = v[off:off + sz].reshape(shp)
                    off += sz
                return Xf

            def extract(Xf):
                return xm.concatenate([Xf[bsl(n0, D, a), bsl(n0, D, b)].reshape(-1) for a, b in blocks])

            def mv(v):
                Xf = embed(v)
                Y = AL[:, 0, :].conj().T @ Xf @ AL[:, 0, :] + AL[:, 1, :].conj().T @ Xf @ AL[:, 1, :]
                return extract(Y)

            rng = np.random.default_rng(seed + r)
            kk = min(k, max(1, ntot - 1))
            if ntot <= 400:
                Tm = xm.stack([mv(e) for e in xm.eye(ntot, dtype=np.complex128)]).T
                lam = np.linalg.eigvals(X.asnp(Tm))
                lam = lam[np.argsort(-np.abs(lam))][:kk]
                out[r] = dict(abs=np.abs(lam).tolist(), arg=np.angle(lam).tolist(), nmv=ntot, conv=True)
                continue
            v0 = X.to_mod(rng.standard_normal(ntot) + 1j * rng.standard_normal(ntot), xm)
            ev, res, nmv, conv = I.arnoldi_top(mv, v0, k=kk, m=min(m, ntot - 1), tol=tol)
            out[r] = dict(abs=np.abs(ev[:kk]).tolist(), arg=np.angle(ev[:kk]).tolist(), nmv=int(nmv),
                          conv=bool(conv), res=[float(x) for x in res[:kk]])
        lam0 = out[0]['abs'][0]
        ev_even = out[0]['abs'][1] if len(out[0]['abs']) > 1 else 0.0
        ev_odd = out[1]['abs'][0]
        xi_e = float(-1.0 / np.log(ev_even / lam0)) if ev_even > 0 else 0.0
        xi_o = float(-1.0 / np.log(ev_odd / lam0)) if ev_odd > 0 else 0.0
        return dict(lam0=lam0, xi_even=xi_e, xi_odd=xi_o, xi=max(xi_e, xi_o), spec=out)

    # ------------------------------------------------------------------------------------------ checkpoint
    def save(self, fn, meta):
        d = dict(AC=X.asnp(self.AC), C=X.asnp(self.C), ints=np.array([self.n0, self.D, self.it]),
                 meta=np.array(json.dumps(meta)))
        if self.Le is not None:
            d['L2'] = X.asnp(self.Le[:, 2, :])
            d['R0'] = X.asnp(self.Re[:, 0, :])
        tmp = fn + '.tmp.npz'
        np.savez(tmp, **d)
        os.replace(tmp, fn)

    def load(self, fn):
        z = np.load(fn, allow_pickle=False)
        self.xm = np
        self.AC, self.C = z['AC'], z['C']
        self.n0, self.D, self.it = [int(x) for x in z['ints']]
        self.AL, self.AR = self._gauge(self.AC, self.C)
        self.Le = self.Re = None
        return json.loads(str(z['meta']))
