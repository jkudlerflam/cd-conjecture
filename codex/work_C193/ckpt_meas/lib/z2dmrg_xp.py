"""C190g (Claude Opus 5.5 agent, 2026-11-09): array-module port of z2dmrg.py (numpy or cupy).

Same algorithm, same operation order and same tolerances as z2dmrg.py.  Every sweep-time routine works on
whichever array module its inputs live in (numpy on the host, cupy on a GPU), so the numpy path of this file is
the CPU code line by line, and the cupy path keeps the whole local solve on the device: the two-site tensor,
both environments, the MPO tensors, the preconditioner eigenbases and the Davidson basis.  Inside the Davidson
loop the only host traffic is the m x m Rayleigh-Ritz matrix (m <= 32), its eigenvector and the residual norm.

Options (module globals, set by the driver):
  USE_GPU      run sweeps on the GPU (needs cupy).
  GPU_MIN_CHI  sweeps with chi below this run on the host (small blocks are launch-bound on the GPU).
  SVD_HOST     do the truncation SVD / mixer eigh and the entropy SVD on the host (numpy/LAPACK) instead of cusolver.
Checkpoints are numpy arrays (MPS.save converts), so a run can stop on one backend and resume on the other.
The Hamiltonian classes and the final measurements (measure, entropies, variance) are the unchanged host code of
z2dmrg.py, imported from it.
"""
import os
import sys
import time
import json
import numpy as np

try:
    import z2dmrg as _Z0
except ImportError:                                   # local layout: z2dmrg.py lives one folder up
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import z2dmrg as _Z0

try:
    import cupy as _cp
except Exception:                                     # no cupy (local machine): numpy only
    _cp = None

# unchanged host pieces
I2, ZZ_, SP, SM, PX, PY, OPS, OPAR = _Z0.I2, _Z0.ZZ_, _Z0.SP, _Z0.SM, _Z0.PX, _Z0.PY, _Z0.OPS, _Z0.OPAR
pauli_bond, decompose, MPO, gauge_phases, rotate_bond = (_Z0.pauli_bond, _Z0.decompose, _Z0.MPO,
                                                         _Z0.gauge_phases, _Z0.rotate_bond)
entropy, sbar = _Z0.entropy, _Z0.sbar
measure, entropies, parity_expectation, energy, variance = (_Z0.measure, _Z0.entropies, _Z0.parity_expectation,
                                                            _Z0.energy, _Z0.variance)

USE_GPU = False
GPU_MIN_CHI = 0
SVD_HOST = True


def amod(a):
    """Array module of a (cupy for device arrays, numpy otherwise)."""
    if _cp is not None and isinstance(a, _cp.ndarray):
        return _cp
    return np


def asnp(a):
    if _cp is not None and isinstance(a, _cp.ndarray):
        return _cp.asnumpy(a)
    return a


def to_mod(a, m):
    if a is None:
        return None
    if m is np:
        return asnp(a)
    return _cp.asarray(a)


def sync():
    if _cp is not None and USE_GPU:
        _cp.cuda.Device().synchronize()


# ----------------------------------------------------------------------------------------------------------
# Block helpers
# ----------------------------------------------------------------------------------------------------------
def bsl(n0, D, q):
    return slice(0, n0) if q == 0 else slice(n0, D)


def bdim(n0, D, q):
    return n0 if q == 0 else D - n0


class MPS:
    def __init__(self, M, n0, P):
        self.M, self.n0, self.P = M, n0, P
        self.L = len(M)

    @property
    def dims(self):
        return [1] + [m.shape[2] for m in self.M]

    def check(self, tol=1e-13):
        L = self.L
        D = self.dims
        for j in range(L):
            M = self.M[j]
            xm = amod(M)
            for qa in (0, 1):
                for s in (0, 1):
                    qb = 1 - (qa ^ s)
                    blk = M[bsl(self.n0[j], D[j], qa), s, bsl(self.n0[j + 1], D[j + 1], qb)]
                    if blk.size and float(xm.max(xm.abs(blk))) > tol:
                        return False
        return True

    def copy(self):
        return MPS([m.copy() for m in self.M], list(self.n0), self.P)

    def to(self, m):
        self.M = [to_mod(x, m) for x in self.M]
        return self

    def host(self):
        """A numpy copy as a z2dmrg.MPS (for the unchanged host measurements)."""
        return _Z0.MPS([np.array(asnp(x)) for x in self.M], list(self.n0), self.P)

    def save(self, fn, meta=None):
        d = {f'M{j}': asnp(m) for j, m in enumerate(self.M)}
        d['n0'] = np.array(self.n0)
        d['P'] = np.array(self.P)
        if meta is not None:
            d['meta'] = np.array(json.dumps(meta))
        np.savez(fn, **d)

    @staticmethod
    def load(fn):
        z = np.load(fn, allow_pickle=False)
        L = len([k for k in z.files if k.startswith('M')])
        return MPS([z[f'M{j}'] for j in range(L)], [int(x) for x in z['n0']], int(z['P']))


def random_mps(L, P, D0=8, seed=0, dtype=float):
    """Identical to z2dmrg.random_mps (host, same generator calls)."""
    p = _Z0.random_mps(L, P, D0=D0, seed=seed, dtype=dtype)
    return MPS(p.M, p.n0, p.P)


def product_mps(L, occ, dtype=float):
    p = _Z0.product_mps(L, occ, dtype=dtype)
    return MPS(p.M, p.n0, p.P)


def right_canonicalize(psi):
    L = psi.L
    for j in range(L - 1, 0, -1):
        _rq_move_left(psi, j)
    xm = amod(psi.M[0])
    psi.M[0] /= xm.linalg.norm(psi.M[0])


def _rq_move_left(psi, j):
    M = psi.M[j]
    xm = amod(M)
    Dl, _, Dr = M.shape
    n0l, n0r = psi.n0[j], psi.n0[j + 1]
    newQ, newR, newdims = {}, {}, {}
    for ql in (0, 1):
        rows = bsl(n0l, Dl, ql)
        nr = bdim(n0l, Dl, ql)
        c0, c1 = bsl(n0r, Dr, ql), bsl(n0r, Dr, 1 - ql)
        A = xm.hstack([M[rows, 0, c0], M[rows, 1, c1]])
        if nr == 0 or A.shape[1] == 0:
            newQ[ql] = xm.zeros((0, A.shape[1]), dtype=M.dtype)
            newR[ql] = xm.zeros((nr, 0), dtype=M.dtype)
            newdims[ql] = 0
            continue
        Q, R = xm.linalg.qr(A.conj().T)
        newQ[ql] = Q.conj().T
        newR[ql] = R.conj().T
        newdims[ql] = Q.shape[1]
    k0, k1 = newdims[0], newdims[1]
    Dn = k0 + k1
    Mn = xm.zeros((Dn, 2, Dr), dtype=M.dtype)
    for ql in (0, 1):
        rr = slice(0, k0) if ql == 0 else slice(k0, Dn)
        c0, c1 = bsl(n0r, Dr, ql), bsl(n0r, Dr, 1 - ql)
        n_c0 = bdim(n0r, Dr, ql)
        Q = newQ[ql]
        Mn[rr, 0, c0] = Q[:, :n_c0]
        Mn[rr, 1, c1] = Q[:, n_c0:]
    Rfull = xm.zeros((Dl, Dn), dtype=M.dtype)
    Rfull[bsl(n0l, Dl, 0), 0:k0] = newR[0]
    Rfull[bsl(n0l, Dl, 1), k0:Dn] = newR[1]
    psi.M[j] = Mn
    psi.M[j - 1] = xm.tensordot(psi.M[j - 1], Rfull, axes=([2], [0]))
    psi.n0[j] = k0


def _qr_move_right(psi, j):
    M = psi.M[j]
    xm = amod(M)
    Dl, _, Dr = M.shape
    n0l, n0r = psi.n0[j], psi.n0[j + 1]
    newQ, newR, newdims = {}, {}, {}
    for qr in (0, 1):
        r0, r1 = bsl(n0l, Dl, qr), bsl(n0l, Dl, 1 - qr)
        cols = bsl(n0r, Dr, qr)
        A = xm.vstack([M[r0, 0, cols], M[r1, 1, cols]])
        nc = bdim(n0r, Dr, qr)
        if nc == 0 or A.shape[0] == 0:
            newQ[qr] = xm.zeros((A.shape[0], 0), dtype=M.dtype)
            newR[qr] = xm.zeros((0, nc), dtype=M.dtype)
            newdims[qr] = 0
            continue
        Q, R = xm.linalg.qr(A)
        newQ[qr], newR[qr], newdims[qr] = Q, R, Q.shape[1]
    k0, k1 = newdims[0], newdims[1]
    Dn = k0 + k1
    Mn = xm.zeros((Dl, 2, Dn), dtype=M.dtype)
    for qr in (0, 1):
        cc = slice(0, k0) if qr == 0 else slice(k0, Dn)
        r0, r1 = bsl(n0l, Dl, qr), bsl(n0l, Dl, 1 - qr)
        n_r0 = bdim(n0l, Dl, qr)
        Mn[r0, 0, cc] = newQ[qr][:n_r0]
        Mn[r1, 1, cc] = newQ[qr][n_r0:]
    Rfull = xm.zeros((Dn, Dr), dtype=M.dtype)
    Rfull[0:k0, bsl(n0r, Dr, 0)] = newR[0]
    Rfull[k0:Dn, bsl(n0r, Dr, 1)] = newR[1]
    psi.M[j] = Mn
    psi.M[j + 1] = xm.tensordot(Rfull, psi.M[j + 1], axes=([1], [0]))
    psi.n0[j + 1] = k0


def flip_site(psi, j):
    p = _Z0.flip_site(psi.host(), j)
    return MPS(p.M, p.n0, p.P)


# ----------------------------------------------------------------------------------------------------------
# Environments
# ----------------------------------------------------------------------------------------------------------
def env_left(Le, A, W, n0a, n0b, mpo):
    xm = amod(A)
    Da, _, Db = A.shape
    w = mpo.w
    dt = np.result_type(Le.dtype, A.dtype, W.dtype)
    T1 = xm.zeros((Da, w, 2, Db), dtype=dt)
    for qa in (0, 1):
        a = bsl(n0a, Da, qa)
        na = bdim(n0a, Da, qa)
        if na == 0:
            continue
        cols = []
        for s in (0, 1):
            cols.append(A[a, s, bsl(n0b, Db, qa ^ s)])
        C = xm.hstack(cols)
        for m in (0, 1):
            ap = bsl(n0a, Da, qa ^ m)
            nap = bdim(n0a, Da, qa ^ m)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if nap == 0 or nm == 0 or C.shape[1] == 0:
                continue
            Lb = Le[ap, ch, a].reshape(nap * nm, na)
            R = (Lb @ C).reshape(nap, nm, -1)
            off = 0
            for s in (0, 1):
                nb = bdim(n0b, Db, qa ^ s)
                T1[ap, ch, s, bsl(n0b, Db, qa ^ s)] = R[:, :, off:off + nb]
                off += nb
    T2 = xm.tensordot(T1, W, axes=([1, 2], [0, 3]))
    Ln = xm.zeros((Db, w, Db), dtype=dt)
    for qbp in (0, 1):
        bp = bsl(n0b, Db, qbp)
        nbp = bdim(n0b, Db, qbp)
        if nbp == 0:
            continue
        for m in (0, 1):
            qb = qbp ^ m
            bk = bsl(n0b, Db, qb)
            nb = bdim(n0b, Db, qb)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if nb == 0 or nm == 0:
                continue
            acc = xm.zeros((nbp, nb * nm), dtype=dt)
            for sp in (0, 1):
                qap = qbp ^ sp
                ap = bsl(n0a, Da, qap)
                if bdim(n0a, Da, qap) == 0:
                    continue
                Ab = A[ap, sp, bp]
                Tb = T2[ap, bk, ch, sp]
                acc += Ab.conj().T @ Tb.reshape(Tb.shape[0], -1)
            Ln[bp, ch, bk] = acc.reshape(nbp, nb, nm).transpose(0, 2, 1)
    return Ln


def env_right(Re, B, W, n0a, n0b, mpo):
    xm = amod(B)
    Da, _, Db = B.shape
    w = mpo.w
    dt = np.result_type(Re.dtype, B.dtype, W.dtype)
    T1 = xm.zeros((Da, 2, Db, w), dtype=dt)
    for qb in (0, 1):
        b = bsl(n0b, Db, qb)
        nb = bdim(n0b, Db, qb)
        if nb == 0:
            continue
        rows = []
        for s in (0, 1):
            rows.append(B[bsl(n0a, Da, qb ^ s), s, b])
        Rw = xm.vstack(rows)
        for m in (0, 1):
            qbp = qb ^ m
            bp = bsl(n0b, Db, qbp)
            nbp = bdim(n0b, Db, qbp)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if nbp == 0 or nm == 0 or Rw.shape[0] == 0:
                continue
            Rb = Re[bp, ch, b].reshape(nbp * nm, nb)
            X = (Rw @ Rb.T).reshape(-1, nbp, nm)
            off = 0
            for s in (0, 1):
                na = bdim(n0a, Da, qb ^ s)
                T1[bsl(n0a, Da, qb ^ s), s, bp, ch] = X[off:off + na]
                off += na
    T2 = xm.tensordot(T1, W, axes=([1, 3], [3, 1]))
    Rn = xm.zeros((Da, w, Da), dtype=dt)
    for qap in (0, 1):
        ap = bsl(n0a, Da, qap)
        nap = bdim(n0a, Da, qap)
        if nap == 0:
            continue
        for m in (0, 1):
            qa = qap ^ m
            ak = bsl(n0a, Da, qa)
            na = bdim(n0a, Da, qa)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if na == 0 or nm == 0:
                continue
            acc = xm.zeros((nap, na * nm), dtype=dt)
            for sp in (0, 1):
                qbp = qap ^ sp
                bp = bsl(n0b, Db, qbp)
                if bdim(n0b, Db, qbp) == 0:
                    continue
                Bb = B[ap, sp, bp]
                Tb = T2[ak, bp, ch, sp]
                acc += Bb.conj() @ Tb.transpose(1, 2, 0).reshape(Tb.shape[1], -1)
            Rn[ap, ch, ak] = acc.reshape(nap, nm, na)
    return Rn


# ----------------------------------------------------------------------------------------------------------
# Two-site effective Hamiltonian
# ----------------------------------------------------------------------------------------------------------
class TwoSite:
    def __init__(self, Le, W1, W2, Re, n0a, Da, n0c, Dc, mpo, P_left_check=None):
        xm = amod(Le)
        self.xm = xm
        self.n0a, self.Da, self.n0c, self.Dc, self.mpo = n0a, Da, n0c, Dc, mpo
        self.dt = np.result_type(Le.dtype, W1.dtype, Re.dtype)
        self.blocks = []
        off = 0
        for qa in (0, 1):
            for s1 in (0, 1):
                for s2 in (0, 1):
                    qc = qa ^ s1 ^ s2
                    shp = (bdim(n0a, Da, qa), bdim(n0c, Dc, qc))
                    n = shp[0] * shp[1]
                    self.blocks.append((qa, s1, s2, qc, shp, off, n))
                    off += n
        self.size = off
        self.Lb = {}
        for qa in (0, 1):
            na = bdim(n0a, Da, qa)
            for m in (0, 1):
                ch = mpo.chs(m)
                nm = ch.stop - ch.start
                nap = bdim(n0a, Da, qa ^ m)
                if na and nap and nm:
                    self.Lb[(qa, m)] = xm.ascontiguousarray(
                        Le[bsl(n0a, Da, qa ^ m), ch, bsl(n0a, Da, qa)].reshape(nap * nm, na))
        self.RbT = {}
        for qcp in (0, 1):
            ncp = bdim(n0c, Dc, qcp)
            for m in (0, 1):
                ch = mpo.chs(m)
                nm = ch.stop - ch.start
                nc = bdim(n0c, Dc, qcp ^ m)
                if ncp and nc and nm:
                    Rb = Re[bsl(n0c, Dc, qcp), ch, bsl(n0c, Dc, qcp ^ m)]
                    self.RbT[(qcp, m)] = xm.ascontiguousarray(Rb.reshape(ncp, nm * nc).T)
        self.W12 = xm.einsum('wvas,vubt->wuabst', W1, W2)
        wd = mpo.w
        self.W12m = xm.ascontiguousarray(self.W12.transpose(1, 2, 3, 0, 4, 5).reshape(wd * 4, wd * 4))
        self.W1 = W1
        self.W2 = W2
        self.W1h, self.W2h = asnp(W1), asnp(W2)
        self.Le, self.Re = Le, Re
        self.diag = None            # only used without a preconditioner; built on demand (diag_vec)

    def diag_vec(self):
        """Diagonal of the two-site problem (z2dmrg.TwoSite.diag), for the diagonal-preconditioned path."""
        if self.diag is None:
            xm = self.xm
            l = xm.einsum('awa->aw', self.Le)
            r = xm.einsum('cuc->cu', self.Re)
            w1 = xm.einsum('wvss->wvs', self.W1)
            w2 = xm.einsum('vuss->vus', self.W2)
            t = xm.einsum('aw,wvs->avs', l, w1)
            t = xm.einsum('avs,vut->astu', t, w2)
            dg = xm.einsum('astu,cu->astc', t, r)
            self.diag = self.pack_dense(dg).real.copy()
        return self.diag

    def build_kron4(self):
        """(s1, s2)-block Kronecker preconditioner, PROOF eq. (5.1); see z2dmrg.TwoSite.build_kron4."""
        xm = self.xm
        n0a, Da, n0c, Dc, mpo = self.n0a, self.Da, self.n0c, self.Dc, self.mpo
        W1, W2, Le, Re = self.W1, self.W2, self.Le, self.Re
        ev = mpo.chs(0)
        self.k4L, self.k4R = {}, {}
        for s1 in (0, 1):
            for qa in (0, 1):
                a = bsl(n0a, Da, qa)
                if bdim(n0a, Da, qa) == 0:
                    continue
                H = xm.tensordot(Le[a, ev, a], W1[ev, 2, s1, s1], axes=([1], [0]))
                lam, U = xm.linalg.eigh((H + H.conj().T) / 2)
                self.k4L[(s1, qa)] = (lam, U, float(lam[0]))
        for s2 in (0, 1):
            for qc in (0, 1):
                c = bsl(n0c, Dc, qc)
                if bdim(n0c, Dc, qc) == 0:
                    continue
                H = xm.tensordot(Re[c, ev, c], W2[0, ev, s2, s2], axes=([1], [0]))
                mu, V = xm.linalg.eigh((H + H.conj().T) / 2)
                self.k4R[(s2, qc)] = (mu, V, float(mu[0]))
        W1h, W2h = self.W1h, self.W2h
        self.k4c = {(s1, s2): float(np.real(W1h[0, 1, s1, s1] * W2h[1, 2, s2, s2]))
                    for s1 in (0, 1) for s2 in (0, 1)}
        e0 = np.inf
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n:
                e0 = min(e0, self.k4L[(s1, qa)][2] + self.k4R[(s2, qc)][2] + self.k4c[(s1, s2)])
        self.k4e0 = e0

    def precond4(self, r_vec, theta, sigma):
        xm = self.xm
        out = xm.empty_like(r_vec)
        sh = min(theta, self.k4e0) - sigma
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n == 0:
                continue
            lam, U, _ = self.k4L[(s1, qa)]
            mu, V, _ = self.k4R[(s2, qc)]
            X = U.conj().T @ r_vec[off:off + n].reshape(shp) @ V.conj()
            X = X / (lam[:, None] + mu[None, :] + (self.k4c[(s1, s2)] - sh))
            out[off:off + n] = (U @ X @ V.T).reshape(-1)
        return out

    def pack_dense(self, th):
        out = self.xm.empty(self.size, dtype=th.dtype)
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n:
                out[off:off + n] = th[bsl(self.n0a, self.Da, qa), s1, s2, bsl(self.n0c, self.Dc, qc)].reshape(-1)
        return out

    def unpack_dense(self, x):
        th = self.xm.zeros((self.Da, 2, 2, self.Dc), dtype=x.dtype)
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n:
                th[bsl(self.n0a, self.Da, qa), s1, s2, bsl(self.n0c, self.Dc, qc)] = x[off:off + n].reshape(shp)
        return th

    def step1(self, x):
        xm = self.xm
        n0a, Da, n0c, Dc, mpo = self.n0a, self.Da, self.n0c, self.Dc, self.mpo
        dt = np.result_type(self.dt, x.dtype)
        Y = xm.zeros((Da, mpo.w, 2, 2, Dc), dtype=dt)
        for qa in (0, 1):
            bl = [b for b in self.blocks if b[0] == qa and b[6] > 0]
            if not bl:
                continue
            C = xm.hstack([x[b[5]:b[5] + b[6]].reshape(b[4]) for b in bl])
            for m in (0, 1):
                if (qa, m) not in self.Lb:
                    continue
                ch = mpo.chs(m)
                nm = ch.stop - ch.start
                nap = bdim(n0a, Da, qa ^ m)
                R = (self.Lb[(qa, m)] @ C).reshape(nap, nm, -1)
                off = 0
                for b in bl:
                    nc = b[4][1]
                    Y[bsl(n0a, Da, qa ^ m), ch, b[1], b[2], bsl(n0c, Dc, b[3])] = R[:, :, off:off + nc]
                    off += nc
        return Y

    def matvec(self, x):
        xm = self.xm
        n0a, Da, n0c, Dc, mpo = self.n0a, self.Da, self.n0c, self.Dc, self.mpo
        Y = self.step1(x)
        wd = mpo.w
        Zt = xm.matmul(self.W12m.astype(np.result_type(self.W12m.dtype, Y.dtype), copy=False),
                       Y.reshape(Da, wd * 4, Dc))
        Zt = Zt.reshape(Da, wd, 2, 2, Dc)
        out = xm.empty(self.size, dtype=Zt.dtype)
        for qap, s1, s2, qcp, shp, off, n in self.blocks:
            if n == 0:
                continue
            acc = xm.zeros(shp, dtype=Zt.dtype)
            for m in (0, 1):
                if (qcp, m) not in self.RbT:
                    continue
                ch = mpo.chs(m)
                Zb = Zt[bsl(n0a, Da, qap), ch, s1, s2, bsl(n0c, Dc, qcp ^ m)]
                Zb = Zb.reshape(shp[0], Zb.shape[1] * Zb.shape[2])
                acc += Zb @ self.RbT[(qcp, m)]
            out[off:off + n] = acc.reshape(-1)
        return out

    def mixer_rho_left(self, x, r):
        xm = self.xm
        n0a, Da = self.n0a, self.Da
        Y = self.step1(x)
        Xt = xm.tensordot(Y, self.W1, axes=([1, 2], [0, 3]))
        keep = [v for v in range(self.mpo.w) if v not in (0, 2)]
        Xt = Xt[:, :, :, keep, :]
        ncol = int(np.prod(Xt.shape[1:4]))
        rows = [Xt[bsl(n0a, Da, r), :, :, :, 0].reshape(bdim(n0a, Da, r), ncol),
                Xt[bsl(n0a, Da, 1 - r), :, :, :, 1].reshape(bdim(n0a, Da, 1 - r), ncol)]
        X = xm.vstack(rows)
        return X @ X.conj().T

    def mixer_rho_right(self, x, r):
        xm = self.xm
        n0c, Dc = self.n0c, self.Dc
        th = self.unpack_dense(x)
        T = xm.tensordot(th, self.Re, axes=([3], [2]))
        T = xm.tensordot(T, self.W2, axes=([2, 4], [3, 1]))
        keep = [v for v in range(self.mpo.w) if v not in (0, 2)]
        T = T[:, :, :, keep, :]
        cols = [T[:, :, bsl(n0c, Dc, r), :, 0], T[:, :, bsl(n0c, Dc, 1 - r), :, 1]]
        cols = [c.transpose(2, 0, 1, 3).reshape(c.shape[2], c.shape[0] * c.shape[1] * c.shape[3]) for c in cols]
        X = xm.vstack(cols)
        return X.conj() @ X.T


def davidson(mv, x0, diag, res_tol=1e-8, maxit=400, mmax=32, keep=4, sigma=1e-4, prec=None):
    """z2dmrg.davidson with the basis on the device of x0; the m x m Ritz problem is solved on the host."""
    xm = amod(x0)
    n = x0.size
    dt = np.result_type(x0.dtype, float)
    x = x0 / xm.linalg.norm(x0)
    mmax = min(mmax, n)
    keep = min(keep, mmax - 1)
    V = xm.zeros((mmax, n), dtype=dt)
    AV = xm.zeros((mmax, n), dtype=dt)
    V[0] = x
    AV[0] = mv(x)
    m = 1
    nmv = 1
    theta, u, rn = None, x, np.inf
    for it in range(maxit):
        Hs = asnp(V[:m].conj() @ AV[:m].T)
        Hs = (Hs + Hs.conj().T) / 2
        th, Y = np.linalg.eigh(Hs)
        Yd = to_mod(Y, xm)
        y = Yd[:, 0]
        theta = th[0]
        u = y @ V[:m]
        Au = y @ AV[:m]
        r = Au - theta * u
        rn = float(xm.linalg.norm(r))
        if rn < res_tol or m >= n:
            break
        if prec is not None:
            t = prec(r, theta)
        else:
            den = xm.maximum(diag - theta, sigma)
            t = r / den
        if m == mmax:
            Vn = Yd[:, :keep].T @ V[:m]
            AVn = Yd[:, :keep].T @ AV[:m]
            V[:keep] = Vn
            AV[:keep] = AVn
            m = keep
        for _ in range(2):
            t = t - V[:m].T @ (V[:m].conj() @ t)
        tn = float(xm.linalg.norm(t))
        if tn < 1e-14:
            t = r - V[:m].T @ (V[:m].conj() @ r)
            tn = float(xm.linalg.norm(t))
            if tn < 1e-14:
                break
        V[m] = t / tn
        AV[m] = mv(V[m])
        m += 1
        nmv += 1
    return theta, u / xm.linalg.norm(u), nmv, rn


# ----------------------------------------------------------------------------------------------------------
# Truncation
# ----------------------------------------------------------------------------------------------------------
def _theta_block(th, n0a, Da, n0c, Dc, r):
    xm = amod(th)
    a0, a1 = bsl(n0a, Da, r), bsl(n0a, Da, 1 - r)
    c0, c1 = bsl(n0c, Dc, r), bsl(n0c, Dc, 1 - r)
    top = xm.hstack([th[a0, 0, 0, c0], th[a0, 0, 1, c1]])
    bot = xm.hstack([th[a1, 1, 0, c0], th[a1, 1, 1, c1]])
    return xm.vstack([top, bot])


def _svd(Bm):
    """Thin SVD on the host (SVD_HOST) or on the device of Bm; returns (U, s_host, Vh) with U, Vh on that device."""
    xm = amod(Bm)
    if xm is np or SVD_HOST:
        U, s, Vh = np.linalg.svd(asnp(Bm), full_matrices=False)
        return to_mod(U, xm), s, to_mod(Vh, xm)
    U, s, Vh = xm.linalg.svd(Bm, full_matrices=False)
    return U, asnp(s), Vh


def _eigh(H):
    xm = amod(H)
    if xm is np or SVD_HOST:
        ev, U = np.linalg.eigh(asnp(H))
        return ev, to_mod(U, xm)
    ev, U = xm.linalg.eigh(H)
    return asnp(ev), U


def _svdvals(M):
    xm = amod(M)
    if xm is np or SVD_HOST:
        return np.linalg.svd(asnp(M), compute_uv=False)
    return asnp(xm.linalg.svd(M, compute_uv=False))


def split(th, n0a, Da, n0c, Dc, chi, cutoff, direction, alpha=0.0, ts=None, x=None):
    """z2dmrg.split; the decomposition backend follows SVD_HOST, the assembly stays on the device of th."""
    xm = amod(th)
    blk = {r: _theta_block(th, n0a, Da, n0c, Dc, r) for r in (0, 1)}
    cand = []
    dec = {}
    for r in (0, 1):
        Bm = blk[r]
        if Bm.size == 0:
            dec[r] = None
            continue
        if alpha > 0:
            if direction == 'right':
                rho = Bm @ Bm.conj().T
                pert = ts.mixer_rho_left(x, r)
            else:
                rho = Bm.conj().T @ Bm
                pert = ts.mixer_rho_right(x, r)
            tp = float(xm.trace(pert).real)
            rho = rho + (alpha * pert / tp if tp > 1e-300 else 0)
            ev, U = _eigh((rho + rho.conj().T) / 2)
            ev, U = ev[::-1], U[:, ::-1]
            dec[r] = ('eig', ev, U)
            for i, e in enumerate(ev):
                cand.append((max(e, 0.0), r, i))
        else:
            U, s, Vh = _svd(Bm)
            dec[r] = ('svd', s, U, Vh)
            for i, e in enumerate(s):
                cand.append((e * e, r, i))
    cand.sort(key=lambda t: -t[0])
    nk = min(chi, len(cand))
    while nk > 1 and cand[nk - 1][0] <= cutoff * cand[0][0]:
        nk -= 1
    kept = {0: [], 1: []}
    for c in cand[:nk]:
        kept[c[1]].append(c[2])
    for r in (0, 1):
        kept[r].sort()
    k0, k1 = len(kept[0]), len(kept[1])
    Dn = k0 + k1
    dt = th.dtype
    A = xm.zeros((Da, 2, Dn), dtype=dt)
    B = xm.zeros((Dn, 2, Dc), dtype=dt)
    for r in (0, 1):
        if not kept[r]:
            continue
        idx = kept[r]
        rr = slice(0, k0) if r == 0 else slice(k0, Dn)
        na0 = bdim(n0a, Da, r)
        nc0 = bdim(n0c, Dc, r)
        Bm = blk[r]
        if dec[r][0] == 'svd':
            _, s, U, Vh = dec[r]
            idd = idx if xm is np else xm.asarray(idx)
            Uk, sk, Vk = U[:, idd], to_mod(s[idx], xm), Vh[idd]
            if direction == 'right':
                L_, R_ = Uk, sk[:, None] * Vk
            else:
                L_, R_ = Uk * sk[None, :], Vk
        else:
            _, ev, U = dec[r]
            idd = idx if xm is np else xm.asarray(idx)
            if direction == 'right':
                L_ = U[:, idd]
                R_ = L_.conj().T @ Bm
            else:
                R_ = U[:, idd].conj().T
                L_ = Bm @ R_.conj().T
        A[bsl(n0a, Da, r), 0, rr] = L_[:na0]
        A[bsl(n0a, Da, 1 - r), 1, rr] = L_[na0:]
        B[rr, 0, bsl(n0c, Dc, r)] = R_[:, :nc0]
        B[rr, 1, bsl(n0c, Dc, 1 - r)] = R_[:, nc0:]
    nrm2 = float(xm.linalg.norm(A if direction == 'left' else B)) ** 2
    tot_th = float(xm.linalg.norm(th)) ** 2
    disc = max(0.0, 1.0 - nrm2 / tot_th)
    if direction == 'right':
        B /= np.sqrt(nrm2)
    else:
        A /= np.sqrt(nrm2)
    return A, B, k0, disc


def schmidt_probs_center(C, n0a, Da, n0b, Db):
    xm = amod(C)
    ps = []
    for qb in (0, 1):
        cols = bsl(n0b, Db, qb)
        if bdim(n0b, Db, qb) == 0:
            continue
        Mrow = xm.vstack([C[bsl(n0a, Da, qb), 0, cols], C[bsl(n0a, Da, 1 - qb), 1, cols]])
        if Mrow.size:
            ps.append(_svdvals(Mrow) ** 2)
    p = np.concatenate(ps)
    return p / p.sum()


# ----------------------------------------------------------------------------------------------------------
# DMRG driver
# ----------------------------------------------------------------------------------------------------------
TIMERS = {}


def _tick(key, t0):
    sync()
    TIMERS[key] = TIMERS.get(key, 0.0) + time.time() - t0
    return time.time()


def dmrg(psi, mpo, sweeps, cutoff=1e-14, log=None, ckpt=None, stop=None, energy_shift_per_site=0.0,
         meas_every=None, profile=False):
    """z2dmrg.dmrg.  Each sweep runs on the GPU if USE_GPU and chi >= GPU_MIN_CHI, otherwise on the host;
    tensors and environments move when the backend changes between sweeps.  With profile=True the time is
    split into TIMERS (synchronized, so profiling slows the GPU path slightly)."""
    L = psi.L
    dt = np.result_type(mpo.dtype, psi.M[0].dtype)

    def mod_for(sw):
        return _cp if (USE_GPU and sw['chi'] >= GPU_MIN_CHI) else np

    cur = mod_for(sweeps[0])
    psi.M = [to_mod(m, cur).astype(dt) for m in psi.M]
    right_canonicalize(psi)
    Wh = mpo.W
    W = [to_mod(x, cur) for x in Wh]
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    w = mpo.w
    Lenv[0] = cur.zeros((1, w, 1), dtype=dt); Lenv[0][0, 0, 0] = 1.0
    Renv[L] = cur.zeros((1, w, 1), dtype=dt); Renv[L][0, 2, 0] = 1.0
    for j in range(L - 1, 0, -1):
        Renv[j] = env_right(Renv[j + 1], psi.M[j], W[j], psi.n0[j], psi.n0[j + 1], mpo)
    hist = []
    Scut = np.zeros(L - 1)
    Zprof = np.zeros(L)
    t_start = time.time()
    tick = _tick if profile else (lambda k, t: t)

    def local(j, sw):
        xm = cur
        D = psi.dims
        t0 = time.time()
        ts = TwoSite(Lenv[j], W[j], W[j + 1], Renv[j + 2], psi.n0[j], D[j], psi.n0[j + 2], D[j + 2], mpo)
        th = xm.tensordot(psi.M[j], psi.M[j + 1], axes=([2], [0]))
        x0 = ts.pack_dense(th)
        t0 = tick('setup', t0)
        if ts.size <= 128:
            H = xm.stack([ts.matvec(e) for e in xm.eye(ts.size, dtype=dt)]).T
            ev, vv = xm.linalg.eigh((H + H.conj().T) / 2)
            e, v = ev[0], vv[:, 0]
            r = float(xm.linalg.norm(ts.matvec(v) - e * v))
            e = float(e)
            nmv = ts.size
            tick('small', t0)
        else:
            if not bool(xm.any(x0)):
                x0 = to_mod(np.random.default_rng(j).standard_normal(ts.size).astype(dt), xm)
            prec = None
            if sw.get('prec', 'kron4') == 'kron4':
                ts.build_kron4()
                t0 = tick('kron_eigh', t0)
                sg = sw.get('sigma', 1e-2)
                if profile:
                    def prec(rr, thv):
                        t1 = time.time()
                        out = ts.precond4(rr, thv, sg)
                        _tick('precond', t1)
                        return out
                else:
                    prec = lambda rr, thv: ts.precond4(rr, thv, sg)
            if profile:
                def mvf(xx):
                    t1 = time.time()
                    out = ts.matvec(xx)
                    _tick('matvec', t1)
                    return out
            else:
                mvf = ts.matvec
            dg = ts.diag_vec() if prec is None else None
            e, v, nmv, r = davidson(mvf, x0, dg, res_tol=sw['res_tol'], maxit=sw.get('maxit', 400),
                                    mmax=sw.get('mmax', 32), keep=sw.get('keep', 6), sigma=sw.get('sigma', 1e-3),
                                    prec=prec)
            tick('davidson_total', t0)
        return ts, e, v, nmv, r

    for isw, sw in enumerate(sweeps):
        new = mod_for(sw)
        if new is not cur:
            cur = new
            psi.M = [to_mod(m, cur) for m in psi.M]
            W = [to_mod(x, cur) for x in Wh]
            Lenv = [to_mod(x, cur) for x in Lenv]
            Renv = [to_mod(x, cur) for x in Renv]
        t0 = time.time()
        chi, alpha = sw['chi'], sw.get('alpha', 0.0)
        trunc, rmax, nmvs, E = 0.0, 0.0, [], None
        for j in range(L - 1):
            ts, E, v, nmv, r = local(j, sw)
            nmvs.append(nmv); rmax = max(rmax, r)
            tt = time.time()
            th = ts.unpack_dense(v)
            A, B, k0, disc = split(th, psi.n0[j], psi.dims[j], psi.n0[j + 2], psi.dims[j + 2], chi, cutoff,
                                   'right', alpha=alpha, ts=ts, x=v)
            tt = tick('split', tt)
            trunc = max(trunc, disc)
            psi.M[j], psi.M[j + 1] = A, B
            psi.n0[j + 1] = k0
            Lenv[j + 1] = env_left(Lenv[j], A, W[j], psi.n0[j], psi.n0[j + 1], mpo)
            Renv[j + 2] = None if j + 2 < L else Renv[j + 2]
            tick('env', tt)
        for j in range(L - 2, -1, -1):
            ts, E, v, nmv, r = local(j, sw)
            nmvs.append(nmv); rmax = max(rmax, r)
            tt = time.time()
            th = ts.unpack_dense(v)
            ab = cur.abs(th) ** 2
            p1 = asnp(ab.sum(axis=(0, 2, 3)))
            Zprof[j] = p1[0] - p1[1]
            if j == L - 2:
                p2 = asnp(ab.sum(axis=(0, 1, 3)))
                Zprof[L - 1] = p2[0] - p2[1]
            A, B, k0, disc = split(th, psi.n0[j], psi.dims[j], psi.n0[j + 2], psi.dims[j + 2], chi, cutoff,
                                   'left', alpha=alpha, ts=ts, x=v)
            tt = tick('split', tt)
            trunc = max(trunc, disc)
            psi.M[j], psi.M[j + 1] = A, B
            psi.n0[j + 1] = k0
            Renv[j + 1] = env_right(Renv[j + 2], B, W[j + 1], psi.n0[j + 1], psi.n0[j + 2], mpo)
            if j + 1 < L:
                Lenv[j + 1] = None
            tt = tick('env', tt)
            if alpha == 0.0:
                Scut[j] = entropy(schmidt_probs_center(psi.M[j], psi.n0[j], psi.dims[j], psi.n0[j + 1],
                                                       psi.dims[j + 1]))
                tick('entropy_svd', tt)
        rec = dict(sweep=isw, chi=chi, alpha=alpha, E=float(np.real(E)), trunc=trunc, res_max=float(rmax),
                   res_tol=sw['res_tol'], mv_mean=float(np.mean(nmvs)), mv_max=int(np.max(nmvs)),
                   Dmax=int(max(psi.dims)), t=time.time() - t0, t_tot=time.time() - t_start)
        dens = (1 - Zprof) / 2
        rec['n_mid'] = float(np.mean(dens[L // 4: 3 * L // 4]))
        rec['N_tot'] = float(np.sum(dens))
        rec['Zprof'] = Zprof.tolist()
        if alpha == 0.0:
            rec['Sbar1'] = sbar(Scut, L, 1); rec['Sbar6'] = sbar(Scut, L, 6); rec['Sbar12'] = sbar(Scut, L, 12)
            rec['Sprof'] = Scut.tolist()
        if meas_every and alpha == 0.0:
            rec.update(meas_every(psi))
        hist.append(rec)
        if log:
            log(rec)
        if ckpt:
            ckpt(psi, hist)
        if stop is not None and alpha == 0.0 and stop(hist):
            break
    return psi, hist
