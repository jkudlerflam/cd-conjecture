"""C196 (Claude Opus 5.5 agent, 2026-11-10): two-site infinite DMRG (McCulloch 2008) for the nearest-neighbour
qubit chain, built on the C190 solver (z2dmrg_xp.py, copied unchanged): exact Z2 parity labels, the (s1, s2)-block
Kronecker preconditioned Davidson of C190 PROOF section 5, residual target 1e-9, at most 300 Davidson iterations.

State after iteration t (window of 2(t+1) sites, total parity P = 0, so every bond label is both the parity to the
left and the parity to the right of the bond):
    ... A_t  Lam_t  B_t ...      A_t (D_{t-1}, 2, D_t) left-orthonormal, B_t (D_t, 2, D_{t-1}) right-orthonormal.
The next two sites are inserted between A_t and B_t.  McCulloch's guess
    theta_{t+1}(s1, s2) = Lam_t B_t(s1) Lam_{t-1}^{-1} A_t(s2) Lam_t
contracts only matching bases (PROOF 2).  The unit cell is two sites.  Its left transfer matrix
    T(X) = sum_{s1 s2} M^dag X M,   M(s1, s2) = A_t(s1) Lam_t B_t(s2) Lam_{t-1}^{-1}
closes on the left basis of bond t-1 and is block diagonal in r = q_bra xor q_ket: r = 0 (parity-even operators)
and r = 1 (parity-odd operators).  xi_r = -2 / ln(|lambda_r| / |lambda_0|) in sites.

Model energies: e_phys = e_mpo + e_off with e_off = -shift (the constant C190 adds per site).
"""
import os
import sys
import time
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg as Z0          # noqa: E402
import z2dmrg_xp as X        # noqa: E402
import models as Mo          # noqa: E402

bsl, bdim = X.bsl, X.bdim


# ----------------------------------------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------------------------------------
def build_model(cfg):
    """Return (mpo, W_bulk, e_off).  W_bulk is the bulk MPO tensor (site 4 of an 8-site chain), complex128."""
    kind = cfg['model']
    Lm = 8
    if kind in ('p153', 'rat', 'sym40', 'sym80'):
        mpo, gph, shift, info = Mo.candidate(kind, Lm, cfg['hz'], gauge=None)
        e_off = -shift + mpo.const / (Lm - 1)
    elif kind in ('xx1', 'xx2'):
        mpo, _ = Mo.xx_control(Lm, cfg['mu'], legs=1 if kind == 'xx1' else 2)
        e_off = mpo.const / (Lm - 1)
    elif kind == 'tfim':
        mpo = Z0.MPO(Lm, [Z0.pauli_bond(Jx=-1.0)] * (Lm - 1), [-cfg['h'] * Z0.ZZ_] * Lm)
        e_off = mpo.const / (Lm - 1)
    else:
        raise ValueError(kind)
    W = np.ascontiguousarray(np.asarray(mpo.W[4], dtype=np.complex128))
    mpo.dtype = complex
    return mpo, W, float(e_off)


def exact_energy(cfg):
    """Exact ground energy density and density of the controls (None for the candidates)."""
    kind = cfg['model']
    if kind in ('xx1', 'xx2'):
        # each leg (xx2) or the chain (xx1) is a free-fermion chain with eps_k = -2 cos k - mu, per chain site
        mu = cfg['mu']
        kF = np.arccos(-mu / 2)
        return dict(e=float(-(2 * np.sin(kF) + mu * kF) / np.pi), n=float(kF / np.pi))
    if kind == 'tfim':
        h = cfg['h']
        k = (np.arange(200000) + 0.5) * np.pi / 200000
        e = -np.mean(np.sqrt(1 + h * h - 2 * h * np.cos(k)))
        return dict(e=float(e), xi_odd=float(1 / np.log(h)), xi_even=float(1 / (2 * np.log(h))))
    return None


# ----------------------------------------------------------------------------------------------------------
# Split with the Schmidt values kept separately
# ----------------------------------------------------------------------------------------------------------
def mixer_split(th, n0a, Da, n0c, Dc, chi, cutoff, alpha, ts, x):
    """Split with White's density-matrix mixer on both sides (C190's perturbation, z2dmrg TwoSite.mixer_rho_*).
    The kept left and right bases are the leading eigenvectors of the perturbed reduced density matrices, the
    same number per parity block on both sides; then theta is projected and re-diagonalized by an SVD, so Lam
    stays diagonal and A, B orthonormal.  Used only in the warm-up iterations."""
    xm = X.amod(th)
    blk = {r: X._theta_block(th, n0a, Da, n0c, Dc, r) for r in (0, 1)}
    decL, decR, cand = {}, {}, []
    for r in (0, 1):
        Bm = blk[r]
        if Bm.size == 0:
            decL[r] = decR[r] = None
            continue
        rho = Bm @ Bm.conj().T
        pert = ts.mixer_rho_left(x, r)
        tp = float(xm.trace(pert).real)
        rho = rho + (alpha * pert / tp if tp > 1e-300 else 0)
        ev, U = X._eigh((rho + rho.conj().T) / 2)
        decL[r] = (ev[::-1], U[:, ::-1])
        rhoR = Bm.conj().T @ Bm
        pertR = ts.mixer_rho_right(x, r)
        tpR = float(xm.trace(pertR).real)
        rhoR = rhoR + (alpha * pertR / tpR if tpR > 1e-300 else 0)
        evR, V = X._eigh((rhoR + rhoR.conj().T) / 2)
        decR[r] = (evR[::-1], V[:, ::-1])
        nmax = min(Bm.shape)
        for i, e in enumerate(decL[r][0][:nmax]):
            cand.append((max(float(e), 0.0), r, i))
    cand.sort(key=lambda t: -t[0])
    nk = min(chi, len(cand))
    while nk > 1 and cand[nk - 1][0] <= cutoff * cand[0][0]:
        nk -= 1
    kr = {0: 0, 1: 0}
    for c in cand[:nk]:
        kr[c[1]] += 1
    k0, k1 = kr[0], kr[1]
    Dn = k0 + k1
    A = xm.zeros((Da, 2, Dn), dtype=th.dtype)
    B = xm.zeros((Dn, 2, Dc), dtype=th.dtype)
    lam = np.zeros(Dn)
    tot = float(xm.linalg.norm(th)) ** 2
    for r in (0, 1):
        if kr[r] == 0:
            continue
        rr = slice(0, k0) if r == 0 else slice(k0, Dn)
        UL = decL[r][1][:, :kr[r]]
        VR = decR[r][1][:, :kr[r]]
        C = UL.conj().T @ blk[r] @ VR
        u, s, vh = X._svd(C)
        Ar = UL @ u
        Br = vh @ VR.conj().T
        na0, nc0 = bdim(n0a, Da, r), bdim(n0c, Dc, r)
        A[bsl(n0a, Da, r), 0, rr] = Ar[:na0]
        A[bsl(n0a, Da, 1 - r), 1, rr] = Ar[na0:]
        B[rr, 0, bsl(n0c, Dc, r)] = Br[:, :nc0]
        B[rr, 1, bsl(n0c, Dc, 1 - r)] = Br[:, nc0:]
        lam[rr] = s
    kept_w = float(np.sum(lam ** 2))
    disc = max(0.0, 1.0 - kept_w / tot)
    lam = lam / np.sqrt(kept_w)
    return A, lam, B, k0, disc


def svd_split(th, n0a, Da, n0c, Dc, chi, cutoff):
    """theta (Da, 2, 2, Dc) -> A (Da, 2, Dn) left-orthonormal, lam (Dn,) host, B (Dn, 2, Dc) right-orthonormal,
    k0 (number of kept even states), discarded weight.  Same truncation rule as z2dmrg.split (global chi, relative
    cutoff on s^2)."""
    xm = X.amod(th)
    blk = {r: X._theta_block(th, n0a, Da, n0c, Dc, r) for r in (0, 1)}
    dec, cand = {}, []
    for r in (0, 1):
        Bm = blk[r]
        if Bm.size == 0:
            dec[r] = None
            continue
        U, s, Vh = X._svd(Bm)
        dec[r] = (U, s, Vh)
        for i, e in enumerate(s):
            cand.append((float(e * e), r, i))
    cand.sort(key=lambda t: -t[0])
    tot = sum(c[0] for c in cand)
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
    A = xm.zeros((Da, 2, Dn), dtype=th.dtype)
    B = xm.zeros((Dn, 2, Dc), dtype=th.dtype)
    lam = np.zeros(Dn)
    for r in (0, 1):
        if not kept[r]:
            continue
        idx = kept[r]
        rr = slice(0, k0) if r == 0 else slice(k0, Dn)
        U, s, Vh = dec[r]
        idd = idx if xm is np else xm.asarray(idx)
        Uk, Vk = U[:, idd], Vh[idd]
        na0, nc0 = bdim(n0a, Da, r), bdim(n0c, Dc, r)
        A[bsl(n0a, Da, r), 0, rr] = Uk[:na0]
        A[bsl(n0a, Da, 1 - r), 1, rr] = Uk[na0:]
        B[rr, 0, bsl(n0c, Dc, r)] = Vk[:, :nc0]
        B[rr, 1, bsl(n0c, Dc, 1 - r)] = Vk[:, nc0:]
        lam[rr] = s[idx]
    kept_w = float(np.sum(lam ** 2))
    disc = max(0.0, 1.0 - kept_w / tot)
    lam = lam / np.sqrt(kept_w)
    return A, lam, B, k0, disc


# ----------------------------------------------------------------------------------------------------------
# Arnoldi with Krylov-Schur style thick restart (numpy or cupy vectors, small matrices on the host)
# ----------------------------------------------------------------------------------------------------------
def arnoldi_top(mv, v0, k=6, m=60, maxrestart=60, tol=1e-9):
    """Largest-|lambda| eigenvalues of the linear map mv.  Returns (eigenvalues sorted by |.| descending,
    residual estimates, number of matvecs, converged flag)."""
    xm = X.amod(v0)
    n = v0.size
    m = min(m, n)
    k = min(k, max(1, m - 2))
    p = min(2 * k, m - 1)
    V = xm.zeros((m + 1, n), dtype=np.complex128)
    H = np.zeros((m + 1, m), dtype=np.complex128)
    V[0] = v0 / xm.linalg.norm(v0)
    j0 = 0
    nmv = 0
    ev = res = None
    conv = False
    for rs in range(maxrestart):
        mm = m
        for j in range(j0, m):
            w = mv(V[j])
            nmv += 1
            h = V[:j + 1].conj() @ w
            w = w - V[:j + 1].T @ h
            h2 = V[:j + 1].conj() @ w
            w = w - V[:j + 1].T @ h2
            hh = X.asnp(h + h2)
            beta = float(xm.linalg.norm(w))
            H[:j + 1, j] += hh
            H[j + 1, j] = beta
            if beta < 1e-13:
                mm = j + 1
                break
            V[j + 1] = w / beta
        Hm = H[:mm, :mm]
        lam, Y = np.linalg.eig(Hm)
        order = np.argsort(-np.abs(lam))
        lam, Y = lam[order], Y[:, order]
        Y = Y / np.linalg.norm(Y, axis=0)[None, :]
        beta = H[mm, mm - 1] if mm < H.shape[0] and mm == m else 0.0
        res = np.abs(beta) * np.abs(Y[mm - 1, :])
        ev = lam
        scale = max(abs(lam[0]), 1e-300)
        if mm < m or np.all(res[:k] < tol * scale):
            conv = True
            break
        # thick restart on the span of the p leading Ritz vectors
        Q, _ = np.linalg.qr(Y[:, :p])
        Hn = Q.conj().T @ Hm @ Q
        bn = beta * Q[mm - 1, :]
        Qd = X.to_mod(Q, xm)
        Vn = Qd.T @ V[:mm]
        vlast = V[mm].copy()
        V[:] = 0
        V[:p] = Vn
        V[p] = vlast
        H[:] = 0
        H[:p, :p] = Hn
        H[p, :p] = bn
        j0 = p
    return ev, res, nmv, conv


# ----------------------------------------------------------------------------------------------------------
# iDMRG
# ----------------------------------------------------------------------------------------------------------
def warmup_alpha(t, n1=20, n2=40, a1=1e-3, a2=1e-5):
    """Mixer amplitude of global iteration t (warm-up only; zero from iteration n2 on)."""
    return a1 if t < n1 else (a2 if t < n2 else 0.0)


class IDMRG:
    def __init__(self, mpo, W, cutoff=1e-14, seed=0, noise=1e-6):
        self.mpo, self.Wh = mpo, W
        self.cutoff = cutoff
        self.noise = noise
        self.rng = np.random.default_rng(seed)
        w = mpo.w
        Le = np.zeros((1, w, 1), dtype=np.complex128); Le[0, 0, 0] = 1.0
        Re = np.zeros((1, w, 1), dtype=np.complex128); Re[0, 2, 0] = 1.0
        self.Le, self.Re = Le, Re
        self.n0, self.D = 1, 1                       # labels of the current centre bond t (left basis)
        self.n0R, self.DR = 1, 1                     # right basis (equal to the left one after the first step)
        self.nstep = 0                               # iDMRG steps done (warm-up mixer counter)
        self.A = self.B = None
        self.lam = np.ones(1)                        # Lam_t
        self.lam_prev = np.ones(1)                   # Lam_{t-1}
        self.n0_prev, self.D_prev = 1, 1
        self.t = -1
        self.E = []                                  # Davidson eigenvalue per iteration (MPO units)
        self.xm = np

    # -------------------------------------------------------------------------------------- initialization
    def init_finite(self, L0=32, chi0=32, seed=1, log=None):
        """Start from a converged finite chain of L0 sites (C190's dmrg, random Z2 start with P = 0, mixer), cut at
        the middle bond: Le from its left half (left-orthonormal), Re from its right half (right-orthonormal).
        The first iDMRG step then inserts two sites between the halves from a random guess."""
        mpo = self.mpo
        w = mpo.w
        Wc = self.Wh
        mpoL = Z0.MPO.__new__(Z0.MPO)
        mpoL.__dict__.update(mpo.__dict__)
        mpoL.L = L0
        mpoL.W = [Wc] * L0
        mpoL.dtype = complex
        psi = Z0.random_mps(L0, 0, D0=8, seed=seed, dtype=complex)
        sw = [dict(chi=16, alpha=1e-3, res_tol=1e-7, maxit=300), dict(chi=chi0, alpha=1e-4, res_tol=1e-8, maxit=300),
              dict(chi=chi0, alpha=1e-5, res_tol=1e-8, maxit=300), dict(chi=chi0, alpha=0.0, res_tol=1e-9, maxit=300),
              dict(chi=chi0, alpha=0.0, res_tol=1e-9, maxit=300), dict(chi=chi0, alpha=0.0, res_tol=1e-9, maxit=300)]
        psi, hist = Z0.dmrg(psi, mpoL, sw, cutoff=self.cutoff, log=log)
        h = L0 // 2
        for j in range(h):
            Z0._qr_move_right(psi, j)
        Le = np.zeros((1, w, 1), dtype=np.complex128); Le[0, 0, 0] = 1.0
        for j in range(h):
            Le = X.env_left(Le, psi.M[j], Wc, psi.n0[j], psi.n0[j + 1], mpo)
        self._n0_left = psi.n0[h]
        Z0._rq_move_left(psi, h)
        Re = np.zeros((1, w, 1), dtype=np.complex128); Re[0, 2, 0] = 1.0
        for j in range(L0 - 1, h - 1, -1):
            Re = X.env_right(Re, psi.M[j], Wc, psi.n0[j], psi.n0[j + 1], mpo)
        self.Le, self.Re = Le, Re
        self.D = Le.shape[0]
        self.DR = Re.shape[0]
        # left labels of bond h: those of the left-orthonormal M[h-1] (recorded before the RQ move changed n0[h])
        self.n0 = self._n0_left
        self.n0R = psi.n0[h]
        self.A = self.B = None
        self.t = h - 1
        self.E = [float(hist[-1]['E'])]
        self.xm = np
        return hist

    # -------------------------------------------------------------------------------------- backend
    def to(self, xm):
        if xm is self.xm:
            return
        f = lambda a: X.to_mod(a, xm) if a is not None else None
        self.Le, self.Re, self.A, self.B = f(self.Le), f(self.Re), f(self.A), f(self.B)
        self.xm = xm

    def W(self):
        return X.to_mod(self.Wh, self.xm)

    # -------------------------------------------------------------------------------------- one iteration
    def guess(self):
        xm = self.xm
        if self.A is None:
            shp = (self.D, 2, 2, self.DR)
            th = self.rng.standard_normal(shp) + 1j * self.rng.standard_normal(shp)
            return X.to_mod(th.astype(np.complex128), xm)
        lp = self.lam_prev
        linv = np.where(lp > 1e-10 * lp.max(), 1.0 / np.maximum(lp, 1e-300), 0.0)
        T = self.B * X.to_mod(linv, xm)[None, None, :]
        th = xm.tensordot(T, self.A, axes=([2], [0]))                     # (D_t, s1, s2, D_t)
        lt = X.to_mod(self.lam, xm)
        th = th * lt[:, None, None, None] * lt[None, None, None, :]
        return th

    def step(self, chi, res_tol=1e-9, maxit=300, mmax=32, keep=6, sigma=1e-2, alpha=0.0):
        xm = self.xm
        mpo = self.mpo
        W = self.W()
        ts = X.TwoSite(self.Le, W, W, self.Re, self.n0, self.D, self.n0R, self.DR, mpo)
        x0 = ts.pack_dense(self.guess())
        nx = float(xm.linalg.norm(x0))
        if nx == 0.0:
            x0 = X.to_mod(self.rng.standard_normal(ts.size) + 0j, xm)
            nx = float(xm.linalg.norm(x0))
        x0 = x0 / nx
        # residual of the prediction (diagnostic of the fixed point, not a convergence criterion)
        hx = ts.matvec(x0)
        rq = complex(xm.vdot(x0, hx)).real
        r_pred = float(xm.linalg.norm(hx - rq * x0))
        if self.noise > 0:
            # random admixture so that sectors of hidden conservation laws (two-leg control) stay reachable
            z = self.rng.standard_normal(ts.size) + 1j * self.rng.standard_normal(ts.size)
            z = X.to_mod(z / np.linalg.norm(z), xm)
            x0 = x0 + self.noise * z
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
        th = ts.unpack_dense(v)
        # densities of the two new sites
        ab = xm.abs(th) ** 2
        p1 = X.asnp(ab.sum(axis=(0, 2, 3)))
        p2 = X.asnp(ab.sum(axis=(0, 1, 3)))
        tot = float(p1.sum())
        n_s1, n_s2 = float(p1[1] / tot), float(p2[1] / tot)
        if alpha > 0:
            A, lam, B, k0, disc = mixer_split(th, self.n0, self.D, self.n0R, self.DR, chi, self.cutoff, alpha, ts, v)
        else:
            A, lam, B, k0, disc = svd_split(th, self.n0, self.D, self.n0R, self.DR, chi, self.cutoff)
        Dn = lam.size
        self.Le = X.env_left(self.Le, A, W, self.n0, k0, mpo)
        self.Re = X.env_right(self.Re, B, W, k0, self.n0R, mpo)
        self.lam_prev, self.n0_prev, self.D_prev = self.lam, self.n0, self.D
        self.A, self.B, self.lam = A, B, lam
        self.n0, self.D = k0, Dn
        self.n0R, self.DR = k0, Dn
        self.nstep += 1
        self.t += 1
        self.E.append(E)
        return dict(t=self.t, E=E, res=float(r), r_pred=r_pred, nmv=int(nmv), disc=disc, D=Dn, k0=k0,
                    alpha=alpha, S=Z0.entropy(lam ** 2), S_prev=Z0.entropy(self.lam_prev ** 2), n1=n_s1, n2=n_s2)

    # -------------------------------------------------------------------------------------- transfer matrix
    def cell_tensor(self):
        """M[s1, s2] (D_{t-1} x D_{t-1}) = A_t(s1) Lam_t B_t(s2) Lam_{t-1}^{-1}, on the current backend."""
        xm = self.xm
        lt = X.to_mod(self.lam, xm)
        lp = self.lam_prev
        linv = X.to_mod(np.where(lp > 1e-10 * lp.max(), 1.0 / np.maximum(lp, 1e-300), 0.0), xm)
        Ms = []
        for s1 in (0, 1):
            AL = self.A[:, s1, :] * lt[None, :]
            for s2 in (0, 1):
                Ms.append((AL @ self.B[:, s2, :]) * linv[None, :])
        return Ms

    def transfer_spectrum(self, k=6, m=60, tol=1e-10, seed=1):
        """Leading |eigenvalues| of the cell transfer matrix in the even (r = 0) and odd (r = 1) sectors."""
        xm = self.xm
        Ms = self.cell_tensor()
        n0, D = self.n0_prev, self.D_prev
        out = {}
        for r in (0, 1):
            blocks = [(q, q ^ r) for q in (0, 1) if bdim(n0, D, q) and bdim(n0, D, q ^ r)]
            shapes = [(bdim(n0, D, a), bdim(n0, D, b)) for a, b in blocks]
            sizes = [a * b for a, b in shapes]
            ntot = int(sum(sizes))
            if ntot == 0:
                out[r] = None
                continue

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
                Y = xm.zeros_like(Xf)
                for M in Ms:
                    Y += M.conj().T @ Xf @ M
                return extract(Y)

            rng = np.random.default_rng(seed + r)
            v0 = X.to_mod(rng.standard_normal(ntot) + 1j * rng.standard_normal(ntot), xm)
            kk = min(k, max(1, ntot - 1))
            if ntot <= 400:
                # dense for small cases
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
        xi = {}
        xi['lam0'] = lam0
        ev_even = out[0]['abs'][1] if len(out[0]['abs']) > 1 else 0.0
        ev_odd = out[1]['abs'][0] if out.get(1) else 0.0
        xi['xi_even'] = float(-2.0 / np.log(ev_even / lam0)) if ev_even > 0 else 0.0
        xi['xi_odd'] = float(-2.0 / np.log(ev_odd / lam0)) if ev_odd > 0 else 0.0
        xi['xi'] = max(xi['xi_even'], xi['xi_odd'])
        xi['spec'] = out
        return xi

    # -------------------------------------------------------------------------------------- checkpoint
    def save(self, fn, meta):
        d = dict(Le=X.asnp(self.Le), Re=X.asnp(self.Re), A=X.asnp(self.A), B=X.asnp(self.B), lam=self.lam,
                 lam_prev=self.lam_prev,
                 ints=np.array([self.n0, self.D, self.n0_prev, self.D_prev, self.t, self.n0R, self.DR, self.nstep]),
                 E=np.array(self.E), meta=np.array(json.dumps(meta)))
        tmp = fn + '.tmp.npz'
        np.savez(tmp, **d)
        os.replace(tmp, fn)

    def load(self, fn):
        z = np.load(fn, allow_pickle=False)
        self.xm = np
        self.Le, self.Re, self.A, self.B = z['Le'], z['Re'], z['A'], z['B']
        self.lam, self.lam_prev = z['lam'], z['lam_prev']
        self.n0, self.D, self.n0_prev, self.D_prev, self.t, self.n0R, self.DR, self.nstep = [int(x) for x in z['ints']]
        self.E = [float(x) for x in z['E']]
        return json.loads(str(z['meta']))
