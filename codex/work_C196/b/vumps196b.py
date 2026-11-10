"""C196b route 2 (Claude Opus 5.5 agent, 2026-11-10): VUMPS with an Lc-site unit cell (Zauner-Stauber et al. 2018,
multi-site version, parallel update) for the nearest-neighbour qubit chain with exact Z2 labels, numpy or cupy.

Start. The uniform MPS of period Lc = 2q closed from the last iDMRG iteration (idmrg196b),
    M = A[0] ... A[q-1] (Lam_t B[q]) B[q+1] ... (B[2q-1] Lam_{t-1}^{-1}),
brought to the mixed canonical form (A_L^j, A_R^j, C[j]) by fixed-point QR / LQ sweeps around the cell. Bond j is
the bond left of site j; bond Lc is bond 0. A_L^j C[j+1] = C[j] A_R^j = A_C^j.

One step. Infinite environments at bond 0 (C190 MPO convention: channel 0 = start, 2 = end; the completed
half-chain Hamiltonians in Le[:, 2, :] and Re[:, 0, :], each solved by GMRES with the fixed-point projector, the
energy per cell subtracted), propagated to every bond with z2dmrg_xp.env_left / env_right. Then for every site j,
A_C^j = lowest eigenvector of H_AC (C190 Davidson, single-site Kronecker preconditioner of the first attempt,
generalized to different left and right bonds) and for every bond C[j] = lowest eigenvector of H_C (Davidson,
preconditioned through the lift with A_L as in the first attempt). A_L^j = polar(A_C^j) polar(C[j+1])^dag,
A_R^j = polar(C[j])^dag polar(A_C^j). Gradient norm eps = max_j max(|A_C^j - A_L^j C[j+1]|, |A_C^j - C[j] A_R^j|).
Davidson residual target 1e-9 and at most 300 iterations (the iDMRG settings), GMRES tolerance 1e-12.
"""
import os
import sys
import time
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg as Z0          # noqa: E402
import z2dmrg_xp as X        # noqa: E402
import vumps196 as V1        # noqa: E402  (first attempt: gmres, tl, trr)

bsl, bdim = X.bsl, X.bdim


def lab(n0, D):
    q = np.zeros(D, dtype=int)
    q[n0:] = 1
    return q


def site_mask(n0a, Da, n0c, Dc):
    qa, qc = lab(n0a, Da), lab(n0c, Dc)
    m = np.zeros((Da, 2, Dc))
    for s in (0, 1):
        m[:, s, :] = ((qa[:, None] ^ s) == qc[None, :])
    return m


def bond_mask(n0, D):
    q = lab(n0, D)
    return (q[:, None] == q[None, :]).astype(float)


def polar(M):
    xm = X.amod(M)
    U, s, Vh = xm.linalg.svd(M, full_matrices=False)
    return U @ Vh


class CellVUMPS:
    def __init__(self, mpo, Wh, e_off=0.0):
        self.mpo, self.Wh, self.e_off = mpo, np.ascontiguousarray(Wh), e_off
        self.xm = np
        self.Le = self.Re = None
        self.it = 0

    # ------------------------------------------------------------------------------------------ start
    def from_cell(self, Ms, n0, maxit=500, tol=1e-13):
        """Ms: Lc tensors of a uniform MPS of period Lc (bond Lc identified with bond 0), n0: Lc + 1 labels."""
        xm = np
        Lc = len(Ms)
        assert n0[Lc] == n0[0] and Ms[0].shape[0] == Ms[-1].shape[2]
        self.Lc = Lc
        self.n0 = list(n0[:Lc])
        self.D = [Ms[j].shape[0] for j in range(Lc)]
        Ms = [X.asnp(m) for m in Ms]
        # left gauge: Y_j M_j = A_L^j Y_{j+1}, iterated around the cell to the fixed point
        Y = np.eye(self.D[0], dtype=np.complex128)
        for it in range(maxit):
            Y0 = Y
            AL = []
            for j in range(Lc):
                T = np.tensordot(Y, Ms[j], axes=([1], [0]))
                Dl, _, Dr = T.shape
                Q, R = np.linalg.qr(T.reshape(Dl * 2, Dr))
                ph = np.sign(np.real(np.diag(R))); ph[ph == 0] = 1
                Q, R = Q * ph[None, :], ph[:, None] * R
                AL.append(Q.reshape(Dl, 2, Dr))
                Y = R / np.linalg.norm(R)
            if np.linalg.norm(Y - Y0) < tol:
                break
        self.YL = Y
        # right gauge: M_j Z_{j+1} = Z_j A_R^j
        Z = np.eye(self.D[0], dtype=np.complex128)
        for it in range(maxit):
            Z0_ = Z
            AR = [None] * Lc
            Zs = [None] * (Lc + 1)
            Zs[Lc] = Z
            for j in range(Lc - 1, -1, -1):
                T = np.tensordot(Ms[j], Z, axes=([2], [0]))
                Dl, _, Dr = T.shape
                Q, R = np.linalg.qr(T.reshape(Dl, 2 * Dr).T)        # T^T = Q R  ->  T = R^T Q^T
                ph = np.sign(np.real(np.diag(R))); ph[ph == 0] = 1
                Q, R = Q * ph[None, :], ph[:, None] * R
                AR[j] = Q.T.reshape(Dl, 2, Dr)
                Z = R.T / np.linalg.norm(R)
                Zs[j] = Z
            if np.linalg.norm(Z - Z0_) < tol:
                break
        # left gauges at every bond with the converged Y
        Ys = [None] * (Lc + 1)
        Ys[0] = Y
        for j in range(Lc):
            T = np.tensordot(Ys[j], Ms[j], axes=([1], [0]))
            Dl, _, Dr = T.shape
            Q, R = np.linalg.qr(T.reshape(Dl * 2, Dr))
            ph = np.sign(np.real(np.diag(R))); ph[ph == 0] = 1
            AL[j] = (Q * ph[None, :]).reshape(Dl, 2, Dr)
            Ys[j + 1] = ph[:, None] * R
        C = []
        for j in range(Lc):
            Cj = Ys[j] @ Zs[j]
            C.append(Cj / np.linalg.norm(Cj))
        # masks (exact Z2 structure) and the gauge-consistent A_C, then regauge with the polar rule
        self.mask = [site_mask(self.n0[j], self.D[j], self.n0[(j + 1) % Lc], self.D[(j + 1) % Lc]) for j in range(Lc)]
        self.bmask = [bond_mask(self.n0[j], self.D[j]) for j in range(Lc)]
        self.C = [C[j] * self.bmask[j] for j in range(Lc)]
        self.AC = [np.tensordot(AL[j], self.C[(j + 1) % Lc], axes=([2], [0])) * self.mask[j] for j in range(Lc)]
        self.AL, self.AR = AL, AR
        self.regauge()
        return self.errors()

    # ------------------------------------------------------------------------------------------ backend
    def to(self, xm):
        if xm is self.xm:
            return
        f = lambda a: X.to_mod(a, xm) if a is not None else None
        for nm in ('AL', 'AR', 'AC', 'C', 'mask', 'bmask'):
            setattr(self, nm, [f(a) for a in getattr(self, nm)])
        self.Le = [f(a) for a in self.Le] if self.Le is not None else None
        self.Re = [f(a) for a in self.Re] if self.Re is not None else None
        self.xm = xm

    # ------------------------------------------------------------------------------------------ gauge
    def regauge(self):
        xm = self.xm
        Lc = self.Lc
        for j in range(Lc):
            Dl, _, Dr = self.AC[j].shape
            UC_r = polar(self.C[(j + 1) % Lc])
            UC_l = polar(self.C[j])
            AL = (polar(self.AC[j].reshape(Dl * 2, Dr)) @ UC_r.conj().T).reshape(Dl, 2, Dr)
            AR = (UC_l.conj().T @ polar(self.AC[j].reshape(Dl, 2 * Dr).T).T).reshape(Dl, 2, Dr)
            self.AL[j] = AL * self.mask[j]
            self.AR[j] = AR * self.mask[j]

    def errors(self):
        xm = self.xm
        Lc = self.Lc
        e = 0.0
        for j in range(Lc):
            eL = float(xm.linalg.norm(self.AC[j] - xm.tensordot(self.AL[j], self.C[(j + 1) % Lc], axes=([2], [0]))))
            eR = float(xm.linalg.norm(self.AC[j] - xm.tensordot(self.C[j], self.AR[j], axes=([1], [0]))))
            e = max(e, eL, eR)
        return e

    # ------------------------------------------------------------------------------------------ environments
    def _cell_left(self, L):
        W = X.to_mod(self.Wh, self.xm)
        out = [L]
        for j in range(self.Lc):
            L = X.env_left(L, self.AL[j], W, self.n0[j], self.n0[(j + 1) % self.Lc], self.mpo)
            out.append(L)
        return out

    def _cell_right(self, R):
        W = X.to_mod(self.Wh, self.xm)
        out = [None] * (self.Lc + 1)
        out[self.Lc] = R
        for j in range(self.Lc - 1, -1, -1):
            R = X.env_right(R, self.AR[j], W, self.n0[j], self.n0[(j + 1) % self.Lc], self.mpo)
            out[j] = R
        return out

    def environments(self, tol=1e-12):
        xm = self.xm
        Lc, w = self.Lc, self.mpo.w
        D0 = self.D[0]
        I_ = xm.eye(D0, dtype=np.complex128)
        # left
        L = xm.zeros((D0, w, D0), dtype=np.complex128); L[:, 0, :] = I_
        for _ in range(2):
            L = self._cell_left(L)[-1]
            L[:, 0, :] = I_
            L[:, 2, :] = 0
        Lp = self._cell_left(L)[-1]
        Y = Lp[:, 2, :]
        rho = self.C[0] @ self.C[0].conj().T
        e = complex(xm.sum(Y * rho.T)).real
        Ip = np.eye(2)

        def opL(v):
            Xm = v.reshape(D0, D0)
            Tx = Xm
            for j in range(Lc):
                Tx = V1.tl(self.AL[j], Tx, Ip)
            return (Xm - Tx + xm.sum(Xm * rho.T) * I_).reshape(-1)
        x0 = self.Le[0][:, 2, :].reshape(-1) if self.Le is not None else None
        x, relL, mvL = V1.gmres(opL, (Y - e * I_).reshape(-1), x0, tol=tol)
        L[:, 2, :] = x.reshape(D0, D0)
        self.Le = self._cell_left(L)[:Lc]
        # right
        R = xm.zeros((D0, w, D0), dtype=np.complex128); R[:, 2, :] = I_
        for _ in range(2):
            R = self._cell_right(R)[0]
            R[:, 2, :] = I_
            R[:, 0, :] = 0
        Rp = self._cell_right(R)[0]
        Yr = Rp[:, 0, :]
        rhoR = (self.C[0].conj().T @ self.C[0]).T
        eR = complex(xm.sum(Yr * rhoR.T)).real

        def opR(v):
            Xm = v.reshape(D0, D0)
            Tx = Xm
            for j in range(Lc - 1, -1, -1):
                Tx = V1.trr(self.AR[j], Tx, Ip)
            return (Xm - Tx + xm.sum(Xm * rhoR.T) * I_).reshape(-1)
        x0 = self.Re[0][:, 0, :].reshape(-1) if self.Re is not None else None
        x, relR, mvR = V1.gmres(opR, (Yr - eR * I_).reshape(-1), x0, tol=tol)
        R[:, 0, :] = x.reshape(D0, D0)
        self.Re = self._cell_right(R)[:Lc]
        return dict(e_cell=e, e_cellR=eR, gmL=relL, gmR=relR, gmvL=mvL, gmvR=mvR)

    # ------------------------------------------------------------------------------------------ local problems
    def _onesite(self, j):
        Lc = self.Lc
        return OneSite2(self.Le[j], X.to_mod(self.Wh, self.xm), self.Re[(j + 1) % Lc], self.n0[j], self.D[j],
                        self.n0[(j + 1) % Lc], self.D[(j + 1) % Lc], self.Wh)

    def step(self, res_tol=1e-9, maxit=300, sigma=1e-2, env_tol=1e-12):
        xm = self.xm
        Lc = self.Lc
        t0 = time.time()
        eps0 = self.errors()
        info = self.environments(tol=env_tol)
        t_env = time.time() - t0
        newAC, newC = [None] * Lc, [None] * Lc
        rmax, nmv = 0.0, []
        for j in range(Lc):
            os_ = self._onesite(j)
            os_.build_prec()
            x0 = os_.pack(self.AC[j]); x0 = x0 / xm.linalg.norm(x0)
            e1, v, nm, r = X.davidson(os_.matvec, x0, None, res_tol=res_tol, maxit=maxit, mmax=32, keep=6,
                                      prec=lambda rr, th: os_.prec(rr, th, sigma))
            AC = os_.unpack(v)
            ph = complex(xm.vdot(self.AC[j], AC))
            newAC[j] = AC * (abs(ph) / ph) if abs(ph) > 0 else AC
            rmax = max(rmax, float(r)); nmv.append(nm)
            # C on bond j (left of site j): environments Le[j], Re[j]; preconditioned through the lift with A_L^{j-1}
            jm = (j - 1) % Lc
            bd = V1.Bond(self.Le[j], self.Re[j], self.n0[j], self.D[j])
            osm = self._onesite(jm)
            osm.build_prec()
            ALm = self.AL[jm]

            def cprec(rr, th, bd=bd, osm=osm, ALm=ALm):
                Cr = bd.unpack(rr)
                Rr = xm.tensordot(ALm, Cr, axes=([2], [0]))
                Pr = osm.unpack(osm.prec(osm.pack(Rr), th, sigma))
                return bd.pack(xm.tensordot(ALm.conj(), Pr, axes=([0, 1], [0, 1])))
            y0 = bd.pack(self.C[j]); y0 = y0 / xm.linalg.norm(y0)
            eC, vC, nmC, rC = X.davidson(bd.matvec, y0, None, res_tol=res_tol, maxit=maxit, mmax=32, keep=6,
                                         prec=cprec)
            C = bd.unpack(vC)
            ph = complex(xm.vdot(self.C[j], C))
            newC[j] = C * (abs(ph) / ph) if abs(ph) > 0 else C
            rmax = max(rmax, float(rC)); nmv.append(nmC)
        self.AC, self.C = newAC, newC
        self.regauge()
        eps = self.errors()
        self.it += 1
        return dict(it=self.it, e_site=info['e_cell'] / Lc + self.e_off, eR=info['e_cellR'] / Lc + self.e_off,
                    eps=eps, eps_prev=eps0, res_max=rmax, mv_mean=float(np.mean(nmv)), gmL=info['gmL'],
                    gmR=info['gmR'], S=[self.entropy(j) for j in range(Lc)], t_env=t_env, secs=time.time() - t0)

    # ------------------------------------------------------------------------------------------ measurements
    def entropy(self, j):
        s = np.linalg.svd(X.asnp(self.C[j]), compute_uv=False)
        return Z0.entropy(s ** 2)

    def density(self):
        """<n_j> = sum_s=1 |A_C^j|^2 (normalized)."""
        xm = self.xm
        out = []
        for j in range(self.Lc):
            a = self.AC[j]
            out.append(float(xm.sum(xm.abs(a[:, 1, :]) ** 2).real / xm.sum(xm.abs(a) ** 2).real))
        return out

    def transfer_spectrum(self, k=6, m=60, tol=1e-10, seed=1):
        """Leading |eigenvalues| of the cell transfer map of A_L per sector r; xi_r = -Lc / ln|lambda_r|."""
        xm = self.xm
        Lc = self.Lc
        n0, Dd = self.n0[0], self.D[0]
        Ip = np.eye(2)

        def apply(Xf):
            for j in range(Lc):
                Xf = V1.tl(self.AL[j], Xf, Ip)
            return Xf
        out = {}
        import idmrg196 as I1
        for r in (0, 1):
            blocks = [(qq, qq ^ r) for qq in (0, 1) if bdim(n0, Dd, qq) and bdim(n0, Dd, qq ^ r)]
            shapes = [(bdim(n0, Dd, a), bdim(n0, Dd, b)) for a, b in blocks]
            sizes = [a * b for a, b in shapes]
            ntot = int(sum(sizes))

            def embed(v):
                Xf = xm.zeros((Dd, Dd), dtype=np.complex128)
                off = 0
                for (a, b), shp, sz in zip(blocks, shapes, sizes):
                    Xf[bsl(n0, Dd, a), bsl(n0, Dd, b)] = v[off:off + sz].reshape(shp)
                    off += sz
                return Xf

            def extract(Xf):
                return xm.concatenate([Xf[bsl(n0, Dd, a), bsl(n0, Dd, b)].reshape(-1) for a, b in blocks])
            mv = lambda v: extract(apply(embed(v)))
            rng = np.random.default_rng(seed + r)
            v0 = X.to_mod(rng.standard_normal(ntot) + 1j * rng.standard_normal(ntot), xm)
            kk = min(k, max(1, ntot - 1))
            ev, res, nmv, conv = I1.arnoldi_top(mv, v0, k=kk, m=min(m, ntot - 1), tol=tol)
            out[r] = dict(abs=np.abs(ev[:kk]).tolist(), arg=np.angle(ev[:kk]).tolist(), nmv=int(nmv), conv=bool(conv))
        lam0 = out[0]['abs'][0]
        ee = out[0]['abs'][1] if len(out[0]['abs']) > 1 else 0.0
        eo = out[1]['abs'][0] if out[1]['abs'] else 0.0
        xe = float(-Lc / np.log(ee / lam0)) if 0 < ee < lam0 else 0.0
        xo = float(-Lc / np.log(eo / lam0)) if 0 < eo < lam0 else 0.0
        return dict(lam0=lam0, xi_even=xe, xi_odd=xo, xi=max(xe, xo), spec=out)


class OneSite2(V1.OneSite):
    """The first attempt's OneSite with different left (n0a, Da) and right (n0c, Dc) bonds."""

    def __init__(self, Le, W, Re, n0a, Da, n0c, Dc, Wh):
        self.xm = X.amod(Le)
        self.Le, self.W, self.Re, self.Wh = Le, W, Re, Wh
        self.n0a, self.Da, self.n0c, self.Dc = n0a, Da, n0c, Dc
        self.blocks = []
        off = 0
        for qa in (0, 1):
            for s in (0, 1):
                qc = qa ^ s
                shp = (bdim(n0a, Da, qa), bdim(n0c, Dc, qc))
                n = shp[0] * shp[1]
                self.blocks.append((qa, s, qc, shp, off, n))
                off += n
        self.size = off

    def pack(self, A):
        out = self.xm.empty(self.size, dtype=np.complex128)
        for qa, s, qc, shp, off, n in self.blocks:
            if n:
                out[off:off + n] = A[bsl(self.n0a, self.Da, qa), s, bsl(self.n0c, self.Dc, qc)].reshape(-1)
        return out

    def unpack(self, x):
        A = self.xm.zeros((self.Da, 2, self.Dc), dtype=np.complex128)
        for qa, s, qc, shp, off, n in self.blocks:
            if n:
                A[bsl(self.n0a, self.Da, qa), s, bsl(self.n0c, self.Dc, qc)] = x[off:off + n].reshape(shp)
        return A

    def build_prec(self):
        xm = self.xm
        Le, Re, Wh = self.Le, self.Re, self.Wh
        ev = [0, 1, 2]
        self.kL, self.kR = {}, {}
        for s in (0, 1):
            HL = sum(Le[:, w, :] * complex(Wh[w, 2, s, s]) for w in ev)
            HR = sum(Re[:, u, :] * complex(Wh[0, u, s, s]) for u in ev)
            for q in (0, 1):
                if bdim(self.n0a, self.Da, q):
                    sl = bsl(self.n0a, self.Da, q)
                    H = HL[sl, sl]
                    lam, U = xm.linalg.eigh((H + H.conj().T) / 2)
                    self.kL[(s, q)] = (lam, U, float(lam[0]))
                if bdim(self.n0c, self.Dc, q):
                    sl = bsl(self.n0c, self.Dc, q)
                    H = HR[sl, sl]
                    mu, Vv = xm.linalg.eigh((H + H.conj().T) / 2)
                    self.kR[(s, q)] = (mu, Vv, float(mu[0]))
        self.kc = {s: -float(np.real(Wh[0, 2, s, s])) for s in (0, 1)}
        e0 = np.inf
        for qa, s, qc, shp, off, n in self.blocks:
            if n:
                e0 = min(e0, self.kL[(s, qa)][2] + self.kR[(s, qc)][2] + self.kc[s])
        self.e0 = e0


def cell_from_idmrg(st):
    """The 2q-site closing cell of the last iDMRG iteration as (Ms, n0) on the host."""
    L = st.last
    A, B = [X.asnp(a) for a in L['A']], [X.asnp(b) for b in L['B']]
    lam, lp = np.asarray(L['lam']), np.asarray(L['lam_prev'])
    linv = np.where(lp > 1e-10 * lp.max(), 1.0 / np.maximum(lp, 1e-300), 0.0)
    B = list(B)
    B[0] = B[0] * lam[:, None, None]
    B[-1] = B[-1] * linv[None, None, :]
    n0 = list(L['n0'])
    return A + B, n0
