"""Two-site finite DMRG with two abelian charges (payload number, length number).

The c = 3 open chain needs bond dimensions of order 10^3, out of reach of the
dense code in verify_C44/dmrg.py. Here every MPS tensor, environment and
effective Hamiltonian is stored in charge blocks. Charges are encoded as one
integer Q = 1000 N_p + N_l (the bond label is the charge to the left of the
bond); the total charge is fixed by the right boundary.

MPO conventions follow verify_C44/dmrg.py (W[wl, wr, s_out, s_in], channel D-1
"nothing applied yet", channel 0 "done"), with every channel carrying a
definite charge transfer delta_w (the bond is split by transfer before the
operator Schmidt decomposition). Environments L_w, R_w map a ket sector Q to
the bra sector Q + delta_w.

The two-site update uses the half-environments
  LW_w2 = sum_w1 L_w1 (x) W[w1, w2],  RW_w2 = sum_w3 W[w2, w3] (x) R_w3,
so that H_eff theta = sum_w2 LW_w2 theta RW_w2^T, block by block.
"""
import numpy as np

QSHIFT = 1000


def local_charges(n=2):
    """local index 4 l + z: (popcount z, l)."""
    return np.array([bin(s & 3).count('1') * QSHIFT + (s >> 2) for s in range(8)])


def mpo_charged(h2, h1, d, qs, tol=1e-13):
    """Lower-triangular MPO with definite channel charge transfers."""
    M = h2.reshape(d, d, d, d)            # s1' s2' s1 s2
    A, B, dl = [], [], []
    viol = 0.0
    deltas = sorted(set(int(qs[a] - qs[b]) for a in range(d) for b in range(d)))
    Mt = M.transpose(0, 2, 1, 3)          # s1' s1 s2' s2
    covered = np.zeros_like(Mt, dtype=bool)
    scale = np.abs(h2).max()
    for dlt in deltas:
        mask = np.zeros((d, d, d, d), bool)
        for s1p in range(d):
            for s1 in range(d):
                if qs[s1p] - qs[s1] != dlt:
                    continue
                for s2p in range(d):
                    for s2 in range(d):
                        if qs[s2p] - qs[s2] == -dlt:
                            mask[s1p, s1, s2p, s2] = True
        covered |= mask
        X = np.where(mask, Mt, 0).reshape(d * d, d * d)
        if np.abs(X).max() == 0:
            continue
        U, S, Vh = np.linalg.svd(X)
        r = int((S > tol * scale).sum())
        for a in range(r):
            A.append((U[:, a] * np.sqrt(S[a])).reshape(d, d))
            B.append((Vh[a] * np.sqrt(S[a])).reshape(d, d))
            dl.append(dlt)
    viol = float(np.abs(np.where(covered, 0, Mt)).max())
    r = len(A)
    D = r + 2
    W = np.zeros((D, D, d, d), complex)
    W[0, 0] = np.eye(d)
    W[D - 1, D - 1] = np.eye(d)
    W[D - 1, 0] = h1
    for a in range(r):
        W[D - 1, a + 1] = A[a]
        W[a + 1, 0] = B[a]
    delta = np.array([0] + dl + [0])
    rec = sum(np.kron(A[a], B[a]) for a in range(r))
    err = float(np.abs(rec - h2).max())
    return W, delta, err, viol


class BlockMPS:
    def __init__(self, L, qs, Qtot, rng, width=2):
        self.L, self.qs, self.d, self.Qtot = L, qs, len(qs), Qtot
        npt, nlt = divmod(Qtot, QSHIFT)
        # reachable sectors at each bond, near the average filling
        self.bond = []
        for i in range(L + 1):
            secs = {}
            for np_ in range(0, 2 * i + 1):
                for nl in range(0, i + 1):
                    if not (0 <= npt - np_ <= 2 * (L - i) and 0 <= nlt - nl <= L - i):
                        continue
                    if abs(np_ - npt * i / L) > width or abs(nl - nlt * i / L) > width:
                        continue
                    secs[np_ * QSHIFT + nl] = 1
            self.bond.append(secs)
        self.bond[0] = {0: 1}
        self.bond[L] = {Qtot: 1}
        self.A = []
        for i in range(L):
            blk = {}
            for Ql, ml in self.bond[i].items():
                for s in range(self.d):
                    Qr = Ql + int(qs[s])
                    if Qr in self.bond[i + 1]:
                        blk[(Ql, s)] = rng.standard_normal((ml, self.bond[i + 1][Qr])) + 0j
            self.A.append(blk)
        self.prune()
        self.right_canonicalize()

    def prune(self):
        """Drop sectors that cannot connect to both boundaries."""
        for _ in range(2):
            for i in range(1, self.L):
                left_ok = {Ql + int(self.qs[s]) for (Ql, s) in self.A[i - 1]}
                right_ok = {Ql for (Ql, s) in self.A[i]}
                keep = set(self.bond[i]) & left_ok & right_ok
                self.bond[i] = {Q: m for Q, m in self.bond[i].items() if Q in keep}
                self.A[i - 1] = {k: v for k, v in self.A[i - 1].items() if k[0] + int(self.qs[k[1]]) in keep}
                self.A[i] = {k: v for k, v in self.A[i].items() if k[0] in keep}

    def right_layout(self, i):
        """For site i: per left sector Ql, the list of (s, Qr, offset, m)."""
        lay = {}
        for (Ql, s), blk in self.A[i].items():
            lay.setdefault(Ql, [])
        for Ql in lay:
            off = 0
            for s in range(self.d):
                if (Ql, s) in self.A[i]:
                    m = self.A[i][(Ql, s)].shape[1]
                    lay[Ql].append((s, Ql + int(self.qs[s]), off, m))
                    off += m
        return lay

    def right_canonicalize(self):
        for i in range(self.L - 1, 0, -1):
            lay = self.right_layout(i)
            newA, newdims, R = {}, {}, {}
            for Ql, items in lay.items():
                M = np.hstack([self.A[i][(Ql, s)] for s, Qr, off, m in items])
                Q, Rm = np.linalg.qr(M.conj().T)          # M = Rm^dag Q^dag
                k = Q.shape[1]
                Bm = Q.conj().T                            # k x cols
                for s, Qr, off, m in items:
                    newA[(Ql, s)] = Bm[:, off:off + m]
                newdims[Ql] = k
                R[Ql] = Rm.conj().T                        # ml x k
            self.A[i] = newA
            self.bond[i] = newdims
            self.A[i - 1] = {(Qa, s): blk @ R[Qa + int(self.qs[s])] for (Qa, s), blk in self.A[i - 1].items()
                             if Qa + int(self.qs[s]) in R}
        nrm = np.sqrt(sum(np.linalg.norm(b) ** 2 for b in self.A[0].values()))
        for k in self.A[0]:
            self.A[0][k] /= nrm


def _left_combined(mps, i):
    """Rows of theta at bond i+1: for each Qm, list of (Qa, s1, off, m)."""
    lay = {}
    for Qa, ma in mps.bond[i].items():
        for s in range(mps.d):
            Qm = Qa + int(mps.qs[s])
            lay.setdefault(Qm, []).append((Qa, s, ma))
    out = {}
    for Qm, items in lay.items():
        off, lst = 0, []
        for Qa, s, ma in items:
            lst.append((Qa, s, off, ma))
            off += ma
        out[Qm] = (lst, off)
    return out


def _right_combined(mps, j):
    """Columns of theta for site j (right site): for each Qm (bond j), list of (s2, Qb, off, m)."""
    lay = {}
    for Qb, mb in mps.bond[j + 1].items():
        for s in range(mps.d):
            Qm = Qb - int(mps.qs[s])
            lay.setdefault(Qm, []).append((s, Qb, mb))
    out = {}
    for Qm, items in lay.items():
        items.sort()
        off, lst = 0, []
        for s, Qb, mb in items:
            lst.append((s, Qb, off, mb))
            off += mb
        out[Qm] = (lst, off)
    return out


def build_LW(Lenv, W, delta, lc, D):
    """LW[w2][Qm] = dense (rows: bra combined Qm+delta_w2, cols: ket combined Qm)."""
    LW = []
    nz = [[(w1, np.argwhere(np.abs(W[w1, w2]) > 0)) for w1 in range(D) if np.any(W[w1, w2])] for w2 in range(D)]
    for w2 in range(D):
        blocks = {}
        for Qm, (ket, nk) in lc.items():
            Qb = Qm + int(delta[w2])
            if Qb not in lc:
                continue
            bra, nb = lc[Qb]
            braidx = {(Qa, s): (off, m) for Qa, s, off, m in bra}
            X = None
            for w1, pairs in nz[w2]:
                Lw = Lenv[w1]
                if not Lw:
                    continue
                for Qa, s1, off, ma in ket:
                    if Qa not in Lw:
                        continue
                    Lb = Lw[Qa]
                    Qa_b = Qa + int(delta[w1])
                    for s1p, s1k in pairs:
                        if s1k != s1:
                            continue
                        key = (Qa_b, int(s1p))
                        if key not in braidx:
                            continue
                        boff, bm = braidx[key]
                        if X is None:
                            X = np.zeros((nb, nk), complex)
                        X[boff:boff + bm, off:off + ma] += W[w1, w2, s1p, s1k] * Lb
            if X is not None:
                blocks[Qm] = X
        LW.append(blocks)
    return LW


def build_RW(Renv, W, delta, rc, D):
    """RW[w2][Qm] = dense (rows: bra right-combined Qm+delta_w2, cols: ket right-combined Qm)."""
    RW = []
    nz = [[(w3, np.argwhere(np.abs(W[w2, w3]) > 0)) for w3 in range(D) if np.any(W[w2, w3])] for w2 in range(D)]
    for w2 in range(D):
        blocks = {}
        for Qm, (ket, nk) in rc.items():
            Qb = Qm + int(delta[w2])
            if Qb not in rc:
                continue
            bra, nb = rc[Qb]
            braidx = {(s, Qr): (off, m) for s, Qr, off, m in bra}
            X = None
            for w3, pairs in nz[w2]:
                Rw = Renv[w3]
                if not Rw:
                    continue
                for s2, Qr, off, mr in ket:
                    if Qr not in Rw:
                        continue
                    Rb = Rw[Qr]
                    Qr_b = Qr + int(delta[w3])
                    for s2p, s2k in pairs:
                        if s2k != s2:
                            continue
                        key = (int(s2p), Qr_b)
                        if key not in braidx:
                            continue
                        boff, bm = braidx[key]
                        if X is None:
                            X = np.zeros((nb, nk), complex)
                        X[boff:boff + bm, off:off + mr] += W[w2, w3, s2p, s2k] * Rb
            if X is not None:
                blocks[Qm] = X
        RW.append(blocks)
    return RW


def lanczos(matvec, v0, k=60, tol=1e-13):
    v = v0 / np.linalg.norm(v0)
    V = [v]
    alph, beta = [], []
    w = matvec(v)
    e_old = None
    for j in range(k):
        a = np.vdot(V[j], w).real
        alph.append(a)
        w = w - a * V[j] - (beta[-1] * V[j - 1] if j > 0 else 0)
        for u in V:
            w = w - np.vdot(u, w) * u
        b = np.linalg.norm(w)
        T = np.diag(alph) + np.diag(beta, 1) + np.diag(beta, -1)
        ev = np.linalg.eigvalsh(T)
        if e_old is not None and abs(ev[0] - e_old) < tol * max(1, abs(ev[0])):
            break
        e_old = ev[0]
        if b < 1e-14:
            break
        beta.append(b)
        V.append(w / b)
        w = matvec(V[-1])
    T = np.diag(alph) + np.diag(beta[:len(alph) - 1], 1) + np.diag(beta[:len(alph) - 1], -1)
    ev, evec = np.linalg.eigh(T)
    x = sum(evec[i, 0] * V[i] for i in range(len(alph)))
    return ev[0], x / np.linalg.norm(x), len(alph)


def finite_dmrg(W, delta, L, Qtot, chi_schedule, qs, seed=1, lanczos_k=60, diag_ops=None, verbose=False,
                cut=1e-14, width=2):
    d = len(qs)
    D = W.shape[0]
    rng = np.random.default_rng(seed)
    mps = BlockMPS(L, qs, Qtot, rng, width=width)
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    Lenv[0] = [dict() for _ in range(D)]
    Lenv[0][D - 1] = {0: np.ones((1, 1), complex)}
    Renv[L] = [dict() for _ in range(D)]
    Renv[L][0] = {Qtot: np.ones((1, 1), complex)}

    def right_env(j):
        """R at bond j from R at bond j+1 and the right-canonical site j."""
        rc = _right_combined(mps, j)
        RW = build_RW(Renv[j + 1], W, delta, rc, D)
        Bm = {}
        for Qm, (lst, n) in rc.items():
            if Qm not in mps.bond[j]:
                continue
            M = np.zeros((mps.bond[j][Qm], n), complex)
            for s, Qb, off, mb in lst:
                if (Qm, s) in mps.A[j]:
                    M[:, off:off + mb] = mps.A[j][(Qm, s)]
            Bm[Qm] = M
        R = [dict() for _ in range(D)]
        for w in range(D):
            for Qm, X in RW[w].items():
                Qb = Qm + int(delta[w])
                if Qm in Bm and Qb in Bm:
                    R[w][Qm] = Bm[Qb].conj() @ X @ Bm[Qm].T
        return R

    for j in range(L - 1, 0, -1):
        Renv[j] = right_env(j)

    hist = []
    S_prof = np.zeros(L - 1)
    sector_info = [None] * (L - 1)
    dens = {k: np.zeros(L) for k in (diag_ops or {})}

    def update(i, chi, direction):
        lc = _left_combined(mps, i)
        rc = _right_combined(mps, i + 1)
        keys = [Qm for Qm in lc if Qm in rc]
        # theta
        th = {}
        for Qm in keys:
            rows, nr = lc[Qm]
            cols, ncol = rc[Qm]
            X = np.zeros((nr, ncol), complex)
            # middle bond sectors of the current MPS
            if Qm in mps.bond[i + 1]:
                left = np.zeros((nr, mps.bond[i + 1][Qm]), complex)
                for Qa, s, off, ma in rows:
                    if (Qa, s) in mps.A[i]:
                        left[off:off + ma] = mps.A[i][(Qa, s)]
                right = np.zeros((mps.bond[i + 1][Qm], ncol), complex)
                for s, Qb, off, mb in cols:
                    if (Qm, s) in mps.A[i + 1]:
                        right[:, off:off + mb] = mps.A[i + 1][(Qm, s)]
                X = left @ right
            th[Qm] = X
        shapes = [(Qm, th[Qm].shape) for Qm in keys]
        LW = build_LW(Lenv[i], W, delta, lc, D)
        RW = build_RW(Renv[i + 2], W, delta, rc, D)
        terms = []
        for w2 in range(D):
            for Qm in keys:
                Qb = Qm + int(delta[w2])
                if Qm in LW[w2] and Qm in RW[w2] and Qb in th:
                    terms.append((Qm, Qb, LW[w2][Qm], RW[w2][Qm].T))
        sizes = [int(np.prod(s)) for _, s in shapes]
        offs = np.concatenate([[0], np.cumsum(sizes)])
        index = {Qm: k for k, (Qm, _) in enumerate(shapes)}

        def unflat(v):
            return {Qm: v[offs[k]:offs[k + 1]].reshape(s) for k, (Qm, s) in enumerate(shapes)}

        def mv(v):
            t = unflat(v)
            out = np.zeros_like(v)
            for Qm, Qb, A_, B_ in terms:
                k = index[Qb]
                out[offs[k]:offs[k + 1]] += (A_ @ t[Qm] @ B_).ravel()
            return out
        v0 = np.concatenate([th[Qm].ravel() for Qm, _ in shapes])
        if np.linalg.norm(v0) < 1e-12:
            v0 = rng.standard_normal(len(v0)) + 0j
        E, x, nit = lanczos(mv, v0, k=lanczos_k)
        t = unflat(x)
        # truncated SVD across sectors
        svd = {}
        allS = []
        for Qm in keys:
            if min(t[Qm].shape) == 0:
                continue
            U, S, Vh = np.linalg.svd(t[Qm], full_matrices=False)
            svd[Qm] = (U, S, Vh)
            allS += [(sv, Qm) for sv in S]
        allS.sort(key=lambda a: -a[0])
        tot = sum(a[0] ** 2 for a in allS)
        smax = allS[0][0]
        kept = [a for a in allS[:chi] if a[0] > cut * smax]
        trunc = 1 - sum(a[0] ** 2 for a in kept) / tot
        cnt = {}
        for sv, Qm in kept:
            cnt[Qm] = cnt.get(Qm, 0) + 1
        nrm = np.sqrt(sum(a[0] ** 2 for a in kept))
        newA, newB, newdims = {}, {}, {}
        Svals = []
        for Qm, k in cnt.items():
            U, S, Vh = svd[Qm]
            U, S, Vh = U[:, :k], S[:k] / nrm, Vh[:k]
            Svals += list(S)
            if direction == 'R':
                Lmat, Rmat = U, S[:, None] * Vh
            else:
                Lmat, Rmat = U * S[None, :], Vh
            rows, _ = lc[Qm]
            for Qa, s, off, ma in rows:
                blk = Lmat[off:off + ma]
                if np.any(blk):
                    newA[(Qa, s)] = blk
            cols, _ = rc[Qm]
            for s, Qb, off, mb in cols:
                blk = Rmat[:, off:off + mb]
                if np.any(blk):
                    newB[(Qm, s)] = blk
            newdims[Qm] = k
        mps.A[i], mps.A[i + 1], mps.bond[i + 1] = newA, newB, newdims
        Sv = np.array(Svals)
        p = Sv ** 2
        p = p[p > 1e-30]
        ent = float(-(p * np.log(p)).sum())
        # local densities of diagonal operators from the normalized two-site state
        loc = {}
        if diag_ops:
            for name, diag in diag_ops.items():
                v1 = v2 = 0.0
                for Qm in keys:
                    X = t[Qm]
                    rows, _ = lc[Qm]
                    cols, _ = rc[Qm]
                    for Qa, s, off, ma in rows:
                        v1 += diag[s] * np.linalg.norm(X[off:off + ma]) ** 2
                    for s, Qb, off, mb in cols:
                        v2 += diag[s] * np.linalg.norm(X[:, off:off + mb]) ** 2
                loc[name] = (v1, v2)
        # environments
        if direction == 'R':
            Lnew = [dict() for _ in range(D)]
            Am = {}
            for Qm in cnt:
                rows, nr = lc[Qm]
                M = np.zeros((nr, cnt[Qm]), complex)
                for Qa, s, off, ma in rows:
                    if (Qa, s) in newA:
                        M[off:off + ma] = newA[(Qa, s)]
                Am[Qm] = M
            for w in range(D):
                for Qm, X in LW[w].items():
                    Qb = Qm + int(delta[w])
                    if Qm in Am and Qb in Am:
                        Lnew[w][Qm] = Am[Qb].conj().T @ X @ Am[Qm]
            Lenv[i + 1] = Lnew
        else:
            Rnew = [dict() for _ in range(D)]
            Bm = {}
            for Qm in cnt:
                cols, ncol = rc[Qm]
                M = np.zeros((cnt[Qm], ncol), complex)
                for s, Qb, off, mb in cols:
                    if (Qm, s) in newB:
                        M[:, off:off + mb] = newB[(Qm, s)]
                Bm[Qm] = M
            for w in range(D):
                for Qm, X in RW[w].items():
                    Qb = Qm + int(delta[w])
                    if Qm in Bm and Qb in Bm:
                        Rnew[w][Qm] = Bm[Qb].conj() @ X @ Bm[Qm].T
            Renv[i + 1] = Rnew
        return E, ent, trunc, nit, loc, cnt

    for sweep, chi in enumerate(chi_schedule):
        maxtr, maxit = 0.0, 0
        for i in range(L - 1):
            E, ent, tr, nit, loc, cnt = update(i, chi, 'R')
            maxtr, maxit = max(maxtr, tr), max(maxit, nit)
            S_prof[i] = ent
        for i in range(L - 2, -1, -1):
            E, ent, tr, nit, loc, cnt = update(i, chi, 'L')
            maxtr, maxit = max(maxtr, tr), max(maxit, nit)
            S_prof[i] = ent
            sector_info[i] = {str(k): v for k, v in cnt.items()}
            for name, (v1, v2) in loc.items():
                dens[name][i + 1] = v2
                if i == 0:
                    dens[name][0] = v1
        hist.append(dict(sweep=sweep, chi=chi, E=float(E), max_trunc=maxtr, max_lanczos=maxit,
                         max_bond=max(sum(b.values()) for b in mps.bond)))
        if verbose:
            print(hist[-1], flush=True)
    return dict(E=float(E), S=S_prof.copy(), hist=hist, meas={k: list(v) for k, v in dens.items()},
                sectors_mid=sector_info[L // 2 - 1])
