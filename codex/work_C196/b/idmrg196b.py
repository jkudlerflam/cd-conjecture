"""C196b (Claude Opus 5.5 agent, 2026-11-10): infinite DMRG with a q-site unit cell (McCulloch 2008, window of two
cells, the growth scheme of C189's idmrg.py) for the nearest-neighbour qubit chain, on the C190 solver
(z2dmrg_xp.py, copied unchanged): exact Z2 parity labels, the (s1, s2)-block Kronecker preconditioned Davidson
of C190 PROOF section 5, residual target 1e-9, at most 300 Davidson iterations, sigma = 1e-2 (the settings of the
first attempt's idmrg196.py).

State at the start of iteration t (window of n = 2q sites between LP and RP):
    LP  M[0] M[1] ... M[n-1]  RP,   M[0] the centre, M[1..n-1] right-orthonormal,
bond 0 and bond n carry the basis of Lam_{t-1} (left Schmidt basis at bond 0, right Schmidt basis at bond n, the
same index and the same Z2 label, because total parity is even on every window).
One iteration: two-site updates to the right over the whole window (2q - 1), back to the left down to the pair
(q - 1, q) (q updates), whose split gives the diagonal Lam_t at the middle bond q. The left cell A[0..q-1] is
absorbed into LP, the right cell B[q..n-1] into RP, and the next window is McCulloch's prediction
    Lam_t B[q] ... B[n-1] Lam_{t-1}^{-1} A[0] ... A[q-1] Lam_t,
brought to right-canonical form by RQ steps (no truncation). LP's completed channel is shifted by the window
eigenvalue, so the next eigenvalue is the energy of the 2q added sites (PROOF 2).

Bond labels are C190's: the parity of the number of flipped spins (s = 1) to the left of the bond, with the
blocks of label 0 first. The unit cell is meant to hold an even number p of particles (PROOF 1).
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg as Z0          # noqa: E402
import z2dmrg_xp as X        # noqa: E402
import idmrg196 as I1        # noqa: E402  (first attempt: build_model, svd_split, arnoldi_top, exact_energy)

bsl, bdim = X.bsl, X.bdim
build_model, svd_split, arnoldi_top, exact_energy = I1.build_model, I1.svd_split, I1.arnoldi_top, I1.exact_energy


class _Seg:
    """Minimal MPS-like holder so that z2dmrg_xp._rq_move_left can be reused on a window (M, n0 lists)."""

    def __init__(self, M, n0):
        self.M, self.n0 = M, n0


def entropy_lam(lam):
    p = np.asarray(lam, dtype=float) ** 2
    p = p[p > 1e-30]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


class CellIDMRG:
    def __init__(self, mpo, W, q, cutoff=1e-14, seed=0, noise=0.0, e_off=0.0):
        self.mpo, self.Wh, self.q, self.n = mpo, W, q, 2 * q
        self.cutoff, self.noise, self.e_off = cutoff, noise, e_off
        self.rng = np.random.default_rng(seed)
        self.xm = np
        self.t = 0
        self.Ehist = []
        self.last = None

    # ------------------------------------------------------------------------------------------ backend
    def W(self):
        return X.to_mod(self.Wh, self.xm)

    def to(self, xm):
        if xm is self.xm:
            return
        f = lambda a: X.to_mod(a, xm) if a is not None else None
        self.M = [f(m) for m in self.M]
        self.LP, self.RP = f(self.LP), f(self.RP)
        self.xm = xm

    def dims(self):
        return [self.M[0].shape[0]] + [m.shape[2] for m in self.M]

    # ------------------------------------------------------------------------------------------ initial window
    def set_window(self, Ms, n0):
        """Open window (bond 0 and bond n of dimension 1, label 0). Ms right-canonical with the centre at site 0
        is not required: the window is right-canonicalized here."""
        assert len(Ms) == self.n and Ms[0].shape[0] == 1 and Ms[-1].shape[2] == 1
        self.M = [np.ascontiguousarray(np.asarray(m, dtype=np.complex128)) for m in Ms]
        self.n0 = list(n0)
        seg = _Seg(self.M, self.n0)
        for j in range(self.n - 1, 0, -1):
            X._rq_move_left(seg, j)
        self.M[0] = self.M[0] / np.linalg.norm(self.M[0])
        w = self.mpo.w
        self.LP = np.zeros((1, w, 1), dtype=np.complex128); self.LP[0, 0, 0] = 1.0
        self.RP = np.zeros((1, w, 1), dtype=np.complex128); self.RP[0, 2, 0] = 1.0
        self.lam_prev = np.ones(1)
        self.xm = np

    # ------------------------------------------------------------------------------------------ local solve
    def _solve(self, Le, Re, j, th0, res_tol, maxit, sigma, mmax=32, keep=6):
        xm = self.xm
        W = self._Wd
        D = self.dims()
        n0 = self.n0
        ts = X.TwoSite(Le, W, W, Re, n0[j], D[j], n0[j + 2], D[j + 2], self.mpo)
        x0 = ts.pack_dense(th0)
        nx = float(xm.linalg.norm(x0))
        if nx == 0.0:
            x0 = X.to_mod(self.rng.standard_normal(ts.size) + 0j, xm)
            nx = float(xm.linalg.norm(x0))
        x0 = x0 / nx
        if self.noise > 0:
            z = self.rng.standard_normal(ts.size) + 1j * self.rng.standard_normal(ts.size)
            x0 = x0 + self.noise * X.to_mod(z / np.linalg.norm(z), xm)
            x0 = x0 / xm.linalg.norm(x0)
        if ts.size <= 128:
            Hd = X.asnp(xm.stack([ts.matvec(e) for e in xm.eye(ts.size, dtype=np.complex128)]).T)
            evs, vv = np.linalg.eigh((Hd + Hd.conj().T) / 2)
            E, v = float(evs[0]), X.to_mod(vv[:, 0], xm)
            r = float(xm.linalg.norm(ts.matvec(v) - E * v))
            nmv = ts.size
        else:
            ts.build_kron4()
            prec = lambda rr, thv: ts.precond4(rr, thv, sigma)
            E, v, nmv, r = X.davidson(ts.matvec, x0, None, res_tol=res_tol, maxit=maxit, mmax=mmax, keep=keep,
                                      sigma=sigma, prec=prec)
            E = float(E)
        fid = abs(complex(xm.vdot(x0, v)))
        return ts, E, ts.unpack_dense(v), int(nmv), float(r), fid, x0

    def _pred_residual(self, Le, Re):
        """Residual of the predicted window at the first update, ||H x - <H> x|| (x normalized)."""
        xm = self.xm
        W = self._Wd
        D = self.dims()
        ts = X.TwoSite(Le, W, W, Re, self.n0[0], D[0], self.n0[2], D[2], self.mpo)
        x = ts.pack_dense(xm.tensordot(self.M[0], self.M[1], axes=([2], [0])))
        x = x / xm.linalg.norm(x)
        hx = ts.matvec(x)
        rq = complex(xm.vdot(x, hx)).real
        return float(xm.linalg.norm(hx - rq * x)), rq

    # ------------------------------------------------------------------------------------------ sweeps
    def _build_renv(self, upto):
        """Renv[j] for j = n down to upto (Renv[n] = RP)."""
        Wd = self._Wd
        Renv = [None] * (self.n + 1)
        Renv[self.n] = self.RP
        for j in range(self.n - 1, upto - 1, -1):
            Renv[j] = X.env_right(Renv[j + 1], self.M[j], Wd, self.n0[j], self.n0[j + 1], self.mpo)
        return Renv

    def sweep_window(self, chi, nsw=1, res_tol=1e-9, maxit=300, sigma=1e-2, alpha=0.0, log=None):
        """Full sweeps of the window between fixed LP and RP (used on the first, open window). The centre returns
        to site 0 with M[1..n-1] right-orthonormal."""
        n = self.n
        self._Wd = self.W()
        out = []
        for isw in range(nsw):
            Renv = self._build_renv(2)
            Lenv = [None] * (n + 1)
            Lenv[0] = self.LP
            rmax, E, disc = 0.0, None, 0.0
            for j in range(n - 1):
                th0 = self.xm.tensordot(self.M[j], self.M[j + 1], axes=([2], [0]))
                ts, E, th, nmv, r, fid, x = self._solve(Lenv[j], Renv[j + 2], j, th0, res_tol, maxit, sigma)
                rmax = max(rmax, r)
                D = self.dims()
                A, B, k0, dd = X.split(th, self.n0[j], D[j], self.n0[j + 2], D[j + 2], chi, self.cutoff, 'right',
                                       alpha=alpha, ts=ts, x=ts.pack_dense(th))
                disc = max(disc, dd)
                self.M[j], self.M[j + 1] = A, B
                self.n0[j + 1] = k0
                Lenv[j + 1] = X.env_left(Lenv[j], A, self._Wd, self.n0[j], k0, self.mpo)
            for j in range(n - 2, -1, -1):
                th0 = self.xm.tensordot(self.M[j], self.M[j + 1], axes=([2], [0]))
                ts, E, th, nmv, r, fid, x = self._solve(Lenv[j], Renv[j + 2], j, th0, res_tol, maxit, sigma)
                rmax = max(rmax, r)
                D = self.dims()
                A, B, k0, dd = X.split(th, self.n0[j], D[j], self.n0[j + 2], D[j + 2], chi, self.cutoff, 'left',
                                       alpha=alpha, ts=ts, x=ts.pack_dense(th))
                disc = max(disc, dd)
                self.M[j], self.M[j + 1] = A, B
                self.n0[j + 1] = k0
                Renv[j + 1] = X.env_right(Renv[j + 2], B, self._Wd, k0, self.n0[j + 2], self.mpo)
            rec = dict(sweep=isw, E=E, res_max=rmax, disc=disc, D=max(self.dims()), dens=self.window_density())
            out.append(rec)
            if log:
                log(rec)
        return out

    def iterate(self, chi, res_tol=1e-9, maxit=300, sigma=1e-2, nfull=0, alpha=0.0):
        """One iDMRG iteration (see the module docstring). nfull > 0 adds that many complete right-left sweeps of
        the window before the standard pass (used in G1, where the first windows grow from a product across the
        middle bond and one pass cannot reach the exact window ground state). alpha > 0 uses White's density-matrix
        mixer (C190's perturbation, through idmrg196.mixer_split) in every split of this iteration; it is used only
        in warm-up iterations, never in a converged point. Returns a record."""
        xm = self.xm
        n, q = self.n, self.q
        self._Wd = self.W()
        Wd = self._Wd
        t0 = time.time()
        Renv = self._build_renv(2)
        Lenv = [None] * (n + 1)
        Lenv[0] = self.LP
        r_pred, e_pred = self._pred_residual(Lenv[0], Renv[2])
        st = dict(rmax=0.0, disc=0.0, nmvs=[], fid0=None, E=None)

        def upd(j):
            th0 = xm.tensordot(self.M[j], self.M[j + 1], axes=([2], [0]))
            ts, E, th, nmv, r, fid, x = self._solve(Lenv[j], Renv[j + 2], j, th0, res_tol, maxit, sigma)
            if st['fid0'] is None:
                st['fid0'] = fid
            st['rmax'] = max(st['rmax'], r); st['nmvs'].append(nmv); st['E'] = E
            D = self.dims()
            if alpha > 0:
                A, lam, B, k0, dd = I1.mixer_split(th, self.n0[j], D[j], self.n0[j + 2], D[j + 2], chi, self.cutoff,
                                                   alpha, ts, ts.pack_dense(th))
            else:
                A, lam, B, k0, dd = svd_split(th, self.n0[j], D[j], self.n0[j + 2], D[j + 2], chi, self.cutoff)
            st['disc'] = max(st['disc'], dd)
            self.n0[j + 1] = k0
            return A, lam, B, k0

        def right(j):
            A, lam, B, k0 = upd(j)
            self.M[j] = A
            self.M[j + 1] = B * X.to_mod(lam, xm)[:, None, None]
            Lenv[j + 1] = X.env_left(Lenv[j], A, Wd, self.n0[j], k0, self.mpo)

        def left(j, keep_lam=False):
            A, lam, B, k0 = upd(j)
            self.M[j + 1] = B
            Renv[j + 1] = X.env_right(Renv[j + 2], B, Wd, k0, self.n0[j + 2], self.mpo)
            self.M[j] = A if keep_lam else A * X.to_mod(lam, xm)[None, None, :]
            return lam, k0

        for _ in range(nfull):
            for j in range(n - 1):
                right(j)
            for j in range(n - 2, -1, -1):
                left(j)
        # right sweep over the whole window
        for j in range(n - 1):
            right(j)
        # left sweep back to the middle; the split of the pair (q-1, q) gives Lam_t at bond q
        Sb = {}
        for j in range(n - 2, q - 2, -1):
            lam, k0 = left(j, keep_lam=(j == q - 1))
            Sb[j + 1] = entropy_lam(lam)
            if j == q - 1:
                lam_t, k0_t = lam, k0
        E = st['E']
        rmax, disc, nmvs, fid0 = st['rmax'], st['disc'], st['nmvs'], st['fid0']
        LPn = X.env_left(Lenv[q - 1], self.M[q - 1], Wd, self.n0[q - 1], self.n0[q], self.mpo)
        RPn = Renv[q]
        Etot = float(E)
        S_cell = [Sb[b] for b in range(q, n)]                   # bonds q .. n-1: one cell of inequivalent bonds
        # keep what the measurements need (cells, Lam_t, Lam_{t-1}, labels)
        self.last = dict(A=[self.M[i] for i in range(q)], B=[self.M[i] for i in range(q, n)], lam=lam_t,
                         lam_prev=self.lam_prev, n0=list(self.n0), D=self.dims())
        dens, Ncell, varN = self.cell_number()
        # absorb, shift the completed channel of LP by the window energy
        Dq = LPn.shape[0]
        LPn[:, 2, :] -= Etot * xm.eye(Dq, dtype=LPn.dtype)
        self.LP, self.RP = LPn, RPn
        self._predict(lam_t)
        self.t += 1
        e_site = Etot / n + self.e_off
        self.Ehist.append(e_site)
        return dict(t=self.t, chi=int(chi), E=Etot, e_site=e_site, S=Sb[q], S_cell=S_cell, Sbar=float(np.mean(S_cell)),
                    disc=disc, alpha=alpha, res_max=rmax, mv_mean=float(np.mean(nmvs)), mv_max=int(np.max(nmvs)),
                    r_pred=r_pred, fid0=fid0, D=int(self.dims()[0]), k0=int(k0_t), Ncell=Ncell, varN=varN,
                    dens=dens, secs=time.time() - t0)

    def _predict(self, lam_t):
        """Next window Lam_t B[q..n-1] Lam_{t-1}^{-1} A[0..q-1] Lam_t, right-canonical with the centre at site 0."""
        xm = self.xm
        n, q = self.n, self.q
        lp = self.lam_prev
        linv = np.where(lp > 1e-10 * lp.max(), 1.0 / np.maximum(lp, 1e-300), 0.0)
        lt = X.to_mod(lam_t, xm)
        Bc = [self.M[i] for i in range(q, n)]
        Ac = [self.M[i] for i in range(q)]
        n0B = [self.n0[i] for i in range(q, n + 1)]          # bonds q .. n
        n0A = [self.n0[i] for i in range(0, q + 1)]          # bonds 0 .. q
        assert n0B[-1] == n0A[0]
        Bc[0] = Bc[0] * lt[:, None, None]
        Bc[-1] = Bc[-1] * X.to_mod(linv, xm)[None, None, :]
        Ac[-1] = Ac[-1] * lt[None, None, :]
        newM = Bc + Ac
        newn0 = n0B + n0A[1:]
        seg = _Seg(newM, newn0)
        for j in range(n - 1, 0, -1):
            X._rq_move_left(seg, j)
        nrm = float(xm.linalg.norm(newM[0]))
        newM[0] = newM[0] / nrm
        self.M, self.n0 = newM, newn0
        self.lam_prev = np.asarray(lam_t)
        self.pred_norm = nrm

    # ------------------------------------------------------------------------------------------ measurements
    def window_density(self):
        """<n_j> on the current window (right-canonical, centre at site 0), host."""
        xm = self.xm
        out = []
        Rs = [None] * (self.n + 1)
        Dn = self.M[-1].shape[2]
        Rs[self.n] = xm.eye(Dn, dtype=np.complex128)
        for j in range(self.n - 1, 0, -1):
            B = self.M[j]
            Rs[j] = B[:, 0, :] @ Rs[j + 1] @ B[:, 0, :].conj().T + B[:, 1, :] @ Rs[j + 1] @ B[:, 1, :].conj().T
        D0 = self.M[0].shape[0]
        Lx = xm.eye(D0, dtype=np.complex128)
        for j in range(self.n):
            M = self.M[j]
            # T[bra, ket] = M_s^dag L M_s; R[ket, bra] = sum_s B_s R B_s^dag (Hermitian, so R.conj() = R.T)
            T1 = M[:, 1, :].conj().T @ Lx @ M[:, 1, :]
            T0 = M[:, 0, :].conj().T @ Lx @ M[:, 0, :]
            R = Rs[j + 1] if j + 1 < self.n else xm.eye(M.shape[2], dtype=np.complex128)
            num = complex(xm.sum(T1 * R.conj()))
            den = complex(xm.sum((T0 + T1) * R.conj()))
            out.append(float((num / den).real))
            Lx = T0 + T1
        return out

    def cell_number(self):
        """Density profile, mean and variance of the particle number of one cell, from the left cell
        A[0..q-1] of the last iteration closed with rho = Lam_t^2 at bond q (identity at bond 0)."""
        xm = self.xm
        L = self.last
        A = L['A']
        lam = X.to_mod(L['lam'] ** 2, xm)
        D0 = A[0].shape[0]
        Id = xm.eye(D0, dtype=np.complex128)
        # right objects R_j = sum_s A_j[s] R_{j+1} A_j[s]^dag with R_q = Lam_t^2 (ket, bra), left objects from Id
        Rr = [None] * (self.q + 1)
        Rr[self.q] = xm.diag(lam).astype(np.complex128)
        for j in range(self.q - 1, -1, -1):
            Aj = A[j]
            Rr[j] = Aj[:, 0, :] @ Rr[j + 1] @ Aj[:, 0, :].conj().T + Aj[:, 1, :] @ Rr[j + 1] @ Aj[:, 1, :].conj().T
        dens = []
        Lx = Id
        for j in range(self.q):
            Aj = A[j]
            T1 = Aj[:, 1, :].conj().T @ Lx @ Aj[:, 1, :]
            T0 = Aj[:, 0, :].conj().T @ Lx @ Aj[:, 0, :]
            dens.append(float(complex(xm.sum(T1 * Rr[j + 1].T)).real))
            Lx = T0 + T1
        # generating environments E0 (1), E1 (N), E2 (N^2)
        E0, E1, E2 = Id, xm.zeros_like(Id), xm.zeros_like(Id)
        for i in range(self.q):
            Ai = A[i]
            T = lambda Xm, s: Ai[:, s, :].conj().T @ Xm @ Ai[:, s, :]
            n0_ = T(E0, 1)
            nE1 = T(E1, 0) + T(E1, 1)
            nE2 = T(E2, 0) + T(E2, 1) + 2 * T(E1, 1) + n0_
            nE1 = nE1 + n0_
            E0 = T(E0, 0) + T(E0, 1)
            E1, E2 = nE1, nE2
        nrm = float(complex(xm.sum(xm.diagonal(E0) * lam)).real)
        N1 = float(complex(xm.sum(xm.diagonal(E1) * lam)).real) / nrm
        N2 = float(complex(xm.sum(xm.diagonal(E2) * lam)).real) / nrm
        return dens, N1, N2 - N1 * N1

    def transfer_spectrum(self, k=6, m=60, tol=1e-10, seed=1):
        """Leading |eigenvalues| of the two-cell transfer map X -> M^dag X M, M = A[0..q-1] Lam_t B[q..n-1]
        Lam_{t-1}^{-1}, per sector r = q_bra xor q_ket (r = 0 even, r = 1 odd). xi_r = -2q / ln(|lam_r| / lam_0)
        in sites (PROOF 3)."""
        xm = self.xm
        L = self.last
        n0, D = L['n0'], L['D']
        q, n = self.q, self.n
        A, B = L['A'], L['B']
        lt = X.to_mod(L['lam'], xm)
        lp = L['lam_prev']
        linv = X.to_mod(np.where(lp > 1e-10 * lp.max(), 1.0 / np.maximum(lp, 1e-300), 0.0), xm)
        a0, Dd = n0[0], D[0]

        def apply(Xf):
            for Ai in A:
                Xf = Ai[:, 0, :].conj().T @ Xf @ Ai[:, 0, :] + Ai[:, 1, :].conj().T @ Xf @ Ai[:, 1, :]
            Xf = lt[:, None] * Xf * lt[None, :]
            for Bi in B:
                Xf = Bi[:, 0, :].conj().T @ Xf @ Bi[:, 0, :] + Bi[:, 1, :].conj().T @ Xf @ Bi[:, 1, :]
            return linv[:, None] * Xf * linv[None, :]

        out = {}
        for r in (0, 1):
            blocks = [(qq, qq ^ r) for qq in (0, 1) if bdim(a0, Dd, qq) and bdim(a0, Dd, qq ^ r)]
            shapes = [(bdim(a0, Dd, a), bdim(a0, Dd, b)) for a, b in blocks]
            sizes = [a * b for a, b in shapes]
            ntot = int(sum(sizes))
            if ntot == 0:
                out[r] = dict(abs=[], arg=[], nmv=0, conv=True, res=[])
                continue

            def embed(v):
                Xf = xm.zeros((Dd, Dd), dtype=np.complex128)
                off = 0
                for (a, b), shp, sz in zip(blocks, shapes, sizes):
                    Xf[bsl(a0, Dd, a), bsl(a0, Dd, b)] = v[off:off + sz].reshape(shp)
                    off += sz
                return Xf

            def extract(Xf):
                return xm.concatenate([Xf[bsl(a0, Dd, a), bsl(a0, Dd, b)].reshape(-1) for a, b in blocks])

            mv = lambda v: extract(apply(embed(v)))
            rng = np.random.default_rng(seed + r)
            v0 = X.to_mod(rng.standard_normal(ntot) + 1j * rng.standard_normal(ntot), xm)
            kk = min(k, max(1, ntot - 1))
            if ntot <= 400:
                Tm = xm.stack([mv(e) for e in xm.eye(ntot, dtype=np.complex128)]).T
                lam = np.linalg.eigvals(X.asnp(Tm))
                lam = lam[np.argsort(-np.abs(lam))][:kk]
                out[r] = dict(abs=np.abs(lam).tolist(), arg=np.angle(lam).tolist(), nmv=ntot, conv=True,
                              res=[0.0] * len(lam))
                continue
            ev, res, nmv, conv = arnoldi_top(mv, v0, k=kk, m=min(m, ntot - 1), tol=tol)
            ev, res = ev[:kk], res[:kk]
            out[r] = dict(abs=np.abs(ev).tolist(), arg=np.angle(ev).tolist(), nmv=int(nmv), conv=bool(conv),
                          res=[float(x) for x in res])
        lam0 = out[0]['abs'][0]
        n2 = 2 * q
        ev_even = out[0]['abs'][1] if len(out[0]['abs']) > 1 else 0.0
        ev_odd = out[1]['abs'][0] if out[1]['abs'] else 0.0
        xe = float(-n2 / np.log(ev_even / lam0)) if 0 < ev_even < lam0 else (np.inf if ev_even >= lam0 else 0.0)
        xo = float(-n2 / np.log(ev_odd / lam0)) if 0 < ev_odd < lam0 else (np.inf if ev_odd >= lam0 else 0.0)
        return dict(lam0=lam0, xi_even=xe, xi_odd=xo, xi=max(xe, xo), spec=out)

    # ------------------------------------------------------------------------------------------ checkpoint
    def save(self, fn, meta):
        d = {f'M{j}': X.asnp(m) for j, m in enumerate(self.M)}
        d.update(LP=X.asnp(self.LP), RP=X.asnp(self.RP), lam_prev=np.asarray(self.lam_prev), n0=np.array(self.n0),
                 ints=np.array([self.t, self.q]), Ehist=np.array(self.Ehist), meta=np.array(json.dumps(meta)))
        tmp = fn + '.tmp.npz'
        np.savez(tmp, **d)
        os.replace(tmp, fn)

    def load(self, fn):
        z = np.load(fn, allow_pickle=False)
        self.xm = np
        self.M = [z[f'M{j}'] for j in range(self.n)]
        self.LP, self.RP, self.lam_prev = z['LP'], z['RP'], z['lam_prev']
        self.n0 = [int(x) for x in z['n0']]
        self.t = int(z['ints'][0])
        assert int(z['ints'][1]) == self.q
        self.Ehist = [float(x) for x in z['Ehist']]
        return json.loads(str(z['meta']))


# ----------------------------------------------------------------------------------------------------------
# Start states for the first (open) window of n = 2q sites
# ----------------------------------------------------------------------------------------------------------
def eff_window(n, N, w_eff, chis=(8, 16, 64, 64, 64), seed=0, **kw):
    """Ground state of the effective model (C173 PROOF eq. 3.2, C176's fibdmrg, gauge 0, units |t2| = 1) on an open
    chain of n sites with exactly N particles, as a Z2-labelled qubit MPS (particle = flipped spin s = 1). Returns
    (Ms, n0, E_eff). The labels are converted as in PROOF 1: right-physical label (N_right, g) at bond i becomes
    the Z2 label N_left mod 2 = (N - N_right) mod 2."""
    import fibdmrg as F
    ch = F.Chain(n, N, w_eff, gauge=0, **kw)
    dm = F.DMRG(ch, seed=seed)
    E = None
    for c in chis:
        E, _ = dm.sweep(c)
    # order the labels of each bond: Z2 label first (even blocks first), then the U(1) label
    order, off, n0 = [None] * (n + 1), [None] * (n + 1), [None] * (n + 1)
    for i in range(1, n + 1):
        labs = sorted(dm.dimR[i].keys(), key=lambda lb: ((N - lb[0]) % 2, lb[0], lb[1]))
        o, d = {}, 0
        for lb in labs:
            o[lb] = d
            d += dm.dimR[i][lb]
        order[i], off[i] = labs, (o, d)
        n0[i] = sum(dm.dimR[i][lb] for lb in labs if (N - lb[0]) % 2 == 0)
    n0[0] = 1
    Ms = []
    # site 0 from the centre C (keys ((0,0), s, q2))
    o1, d1 = off[1]
    M0 = np.zeros((1, 2, d1), dtype=np.complex128)
    for (ql, s, q2), m in dm.C.items():
        M0[0, s, o1[q2]:o1[q2] + m.shape[1]] += m.reshape(-1)
    Ms.append(M0)
    for i in range(1, n):
        oL, dL = off[i]
        oR, dR = off[i + 1]
        Mi = np.zeros((dL, 2, dR), dtype=np.complex128)
        for (s, q2), m in dm.B[i].items():
            h = ch.rfuse(i, s, q2)
            Mi[oL[h]:oL[h] + m.shape[0], s, oR[q2]:oR[q2] + m.shape[1]] = m
        Ms.append(Mi)
    assert off[n][1] == 1
    return Ms, n0, float(np.real(E))


def product_window(n, occ):
    """Product state with particles on the sites occ, as a Z2-labelled MPS of bond dimension 1."""
    Ms, n0 = [], [1]
    par = 0
    for j in range(n):
        s = 1 if j in occ else 0
        M = np.zeros((1, 2, 1), dtype=np.complex128)
        M[0, s, 0] = 1.0
        Ms.append(M)
        par ^= s
        n0.append(1 if par == 0 else 0)
    # a bond of label 1 with dimension 1 has n0 = 0
    return Ms, n0


def random_window(n, seed=0, D0=8):
    p = Z0.random_mps(n, 0, D0=D0, seed=seed, dtype=complex)
    return [np.asarray(m, dtype=np.complex128) for m in p.M], list(p.n0)
