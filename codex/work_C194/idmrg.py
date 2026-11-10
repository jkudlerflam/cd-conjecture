"""C189: blockade-exact, U(1) infinite DMRG (McCulloch 2008) for C173's effective model (PROOF eq. 3.2) at fixed
density n = p/q, with a window of two q-site cells.

Hamiltonian: the MPO of C176's fibdmrg.Chain (imported, so the model is the same object). fermionic=True puts the
string (-1)^n on the middle site of the NNN hop (free fermions when the blockade is off, control ff2v).

Labels (PROOF section 6). Charges are shifted and scaled to integers. A left-physical label at a bond is (a, f) with
a = q N_left - p X (X the global bond position) and f the occupation of the site left of the bond; a right-physical
label is (b, g) with b = q N_right - p m (m the number of sites to the right) and g the occupation of the site right
of the bond. One site shifts a or b by q s - p, so the label sets are the same in every iteration and the
environments of the half-infinite chains can be reused cell after cell. Charge conservation at density p/q:
a + q(s1 + s2) + b = 2p across two sites, a + b = 0 across a bond. The blockade and N are structural zeros, and the
label scheme is C176's (left-physical labels on left-canonical tensors, right-physical labels on right-canonical
ones, mixed labels on the centre), which keeps both bases orthonormal (C176 PROOF 1.1, 1.2).

One iteration: guess the 2q-site window from McCulloch's prediction psi = C_t B_t C_{t-1}^+ A_t C_t, sweep the
window right (2q - 1 two-site updates) and back to the middle (q - 1 updates), split the centre into the bond
matrix C_{t+1}, absorb the left cell into LP and the right cell into RP.

Correlation lengths: the two-cell transfer matrix T = T[A_t] M[C_t] T[B_t] M[C_{t-1}^+] closes on the bond basis
of LP_{t-1}; it is resolved by the charge dQ carried between bra and ket.
"""
import numpy as np
import scipy.sparse.linalg as spla
import time, json, sys, os, pickle
from fibdmrg import Chain, DMRG, ID, BB, BD, NN, PAR, S_, B1, Bd1, B2, Bd2, N1, N2, N3, E_, NCH


def _nz(op):
    return [(s, t, op[s, t]) for s in range(2) for t in range(2) if abs(op[s, t]) > 0]


class Model:
    def __init__(self, p, q, w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, blockade=True, gauge=0, fermionic=False):
        ch = Chain(8, 2, w, t2=t2, V2=V2, V3=V3, blockade=blockade, gauge=gauge)
        ent = []
        for (a, b, nzl) in ch.ent:
            if fermionic and ((a, b) == (B1, B2) or (a, b) == (Bd1, Bd2)):
                nzl = _nz(PAR)
            ent.append((a, b, nzl))
        self.ent = ent
        self.frm = [[[] for _ in range(2)] for _ in range(NCH)]
        for (a, b, nzl) in ent:
            for (so, si, v) in nzl:
                self.frm[a][si].append((b, so, v))
        self.dtype = ch.dtype
        self.p, self.q, self.bk, self.fermionic = p, q, bool(blockade), fermionic
        self.kw = dict(w=w, t2=t2, V2=V2, V3=V3, blockade=blockade, gauge=gauge)

    def lf(self, lab, s):
        a, f = lab
        if self.bk and f and s:
            return None
        return (a + self.q * s - self.p, s if self.bk else 0)

    def rf(self, s, lab):
        b, g = lab
        if self.bk and s and g:
            return None
        return (b + self.q * s - self.p, s if self.bk else 0)


def dadd(d, e, c=1.0):
    for k, v in e.items():
        if k in d:
            d[k] = d[k] + c * v
        else:
            d[k] = c * v
    return d


# ---------------------------------------------------------------- environments
def lupdate(M, Lold, A):
    new = [dict() for _ in range(NCH)]
    for (a, b, nzl) in M.ent:
        La = Lold[a]
        if not La:
            continue
        for (qb, qk), Lm in La.items():
            cache = {}
            for (sp, s, val) in nzl:
                Ab = A.get((qb, sp)); Ak = A.get((qk, s))
                if Ab is None or Ak is None:
                    continue
                if s not in cache:
                    cache[s] = Lm @ Ak
                X = val * (Ab.conj().T @ cache[s])
                key = (M.lf(qb, sp), M.lf(qk, s))
                d = new[b]
                d[key] = d[key] + X if key in d else X
    return new


def rupdate(M, Rold, B):
    new = [dict() for _ in range(NCH)]
    for (a, b, nzl) in M.ent:
        Rb = Rold[b]
        if not Rb:
            continue
        for (qb2, qk2), Rm in Rb.items():
            cache = {}
            for (sp, s, val) in nzl:
                Bb = B.get((sp, qb2)); Bk = B.get((s, qk2))
                if Bb is None or Bk is None:
                    continue
                if s not in cache:
                    cache[s] = Rm @ Bk.T
                X = val * (Bb.conj() @ cache[s])
                key = (M.rf(sp, qb2), M.rf(s, qk2))
                d = new[a]
                d[key] = d[key] + X if key in d else X
    return new


def ket_index(env):
    out = []
    for d in env:
        ix = {}
        for (qb, qk), m in d.items():
            ix.setdefault(qk, []).append((qb, m))
        out.append(ix)
    return out


def entropy_from_s(s):
    p = s ** 2
    p = p[p > 1e-30]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


class Win:
    """The 2q-site window between LP and RP. Positions 0..n-1, bonds 0..n."""

    def __init__(self, M):
        self.M = M
        self.n = 2 * M.q
        n = self.n
        self.A = [dict() for _ in range(n)]
        self.B = [dict() for _ in range(n)]
        self.dimL = [dict() for _ in range(n + 1)]
        self.dimR = [dict() for _ in range(n + 1)]
        self.Lenv = [None] * (n + 1)
        self.Renv = [None] * (n + 1)
        self.C = None
        self.rng = np.random.default_rng(1)
        self.noise = 0.0

    def __getstate__(self):
        d = dict(self.__dict__)
        d['Lenv'] = [None] * (self.n + 1)
        d['Renv'] = [None] * (self.n + 1)
        return d

    # ---- two-site problem
    def ts_keys(self, j):
        M = self.M
        keys = []
        for ql, dl in self.dimL[j].items():
            for s1 in (0, 1):
                if M.lf(ql, s1) is None:
                    continue
                for qr, dr in self.dimR[j + 2].items():
                    for s2 in (0, 1):
                        if M.bk and s1 and s2:
                            continue
                        if ql[0] + M.q * (s1 + s2) + qr[0] != 2 * M.p:
                            continue
                        if M.rf(s2, qr) is None:
                            continue
                        keys.append(((ql, s1, s2, qr), (dl, dr)))
        return keys

    def solve(self, j, theta0, tol):
        M = self.M
        keys = self.ts_keys(j)
        kl = [k for k, _ in keys]; shapes = [sh for _, sh in keys]
        offs = np.cumsum([0] + [a * b for a, b in shapes]); size = int(offs[-1])
        kpos = {k: i for i, k in enumerate(kl)}
        Lidx = ket_index(self.Lenv[j]); Ridx = ket_index(self.Renv[j + 2])
        frm, bk, dt = M.frm, M.bk, M.dtype

        def unpack(v):
            return {k: v[offs[i]:offs[i + 1]].reshape(shapes[i]) for i, k in enumerate(kl)}

        def mv(v):
            th = unpack(v)
            X = {}
            for (ql, s1, s2, qr), T in th.items():
                for a in range(NCH):
                    lst = Lidx[a].get(ql)
                    if not lst:
                        continue
                    fr = frm[a][s1]
                    if not fr:
                        continue
                    for (qbl, Lm) in lst:
                        LT = None
                        for (b, s1p, val) in fr:
                            if bk and qbl[1] and s1p:
                                continue
                            if LT is None:
                                LT = Lm @ T
                            key = (b, qbl, s1p, s2, qr)
                            if key in X:
                                X[key] += val * LT
                            else:
                                X[key] = val * LT
            out = np.zeros(size, dtype=dt)
            for (b, qbl, s1p, s2, qr), Xm in X.items():
                for (c, s2p, val) in frm[b][s2]:
                    if bk and s1p and s2p:
                        continue
                    lst = Ridx[c].get(qr)
                    if not lst:
                        continue
                    for (qbr, Rm) in lst:
                        nn = kpos.get((qbl, s1p, s2p, qbr))
                        if nn is None:
                            continue
                        out[offs[nn]:offs[nn + 1]] += (val * (Xm @ Rm.T)).reshape(-1)
            return out
        v0 = np.zeros(size, dtype=dt)
        for k, m in theta0.items():
            i = kpos.get(k)
            if i is not None and m.shape == shapes[i]:
                v0[offs[i]:offs[i + 1]] = m.reshape(-1)
        if size <= 64 + (10 ** 9 if False else 0):
            H = np.array([mv(e) for e in np.eye(size, dtype=dt)]).T
            H = 0.5 * (H + H.conj().T)
            e, U = np.linalg.eigh(H)
            return e[0], unpack(U[:, 0])
        nv0 = np.linalg.norm(v0)
        if nv0 < 1e-12:
            v0 = self.rng.standard_normal(size).astype(dt)
        elif self.noise > 0:
            # a small random admixture: the two-leg control has a hidden U(1) x U(1) (leg numbers), and a guess in
            # the wrong leg sector would otherwise be an exact eigenvector that Lanczos never leaves
            v0 = v0 + self.noise * nv0 / np.sqrt(size) * self.rng.standard_normal(size).astype(dt)
        op = spla.LinearOperator((size, size), matvec=mv, dtype=dt)
        e, U = spla.eigsh(op, k=1, which='SA', v0=v0, tol=tol, ncv=min(size - 1, 20), maxiter=1000)
        return e[0], unpack(U[:, 0])

    def theta_right(self, j):
        """theta(j, j+1) from the centre at site j (keys (ql, s1, qmid)) and B[j+1] (keys (s2, qr))."""
        M = self.M
        th = {}
        byleft = {}
        for (s2, qr), m in self.B[j + 1].items():
            byleft.setdefault(M.rf(s2, qr), []).append((s2, qr, m))
        for (ql, s1, qmid), Cm in self.C.items():
            for (s2, qr, m) in byleft.get(qmid, ()):
                if M.bk and s1 and s2:
                    continue
                k = (ql, s1, s2, qr)
                th[k] = th[k] + Cm @ m if k in th else Cm @ m
        return th

    def theta_left(self, j):
        """theta(j, j+1) from A[j] (keys (ql, s1)) and the centre at site j+1 (keys (qm, s2, qr))."""
        M = self.M
        th = {}
        byleft = {}
        for (qm, s2, qr), m in self.C.items():
            byleft.setdefault(qm, []).append((s2, qr, m))
        for (ql, s1), Am in self.A[j].items():
            qm = M.lf(ql, s1)
            for (s2, qr, m) in byleft.get(qm, ()):
                if M.bk and s1 and s2:
                    continue
                k = (ql, s1, s2, qr)
                th[k] = th[k] + Am @ m if k in th else Am @ m
        return th

    def split(self, j, th, chi, direction, eps=1e-14):
        """Label-respecting truncated split of theta(j, j+1) at bond j+1 (C176's rule)."""
        M = self.M
        groups = {}
        for k in th:
            ql, s1, s2, qr = k
            g = M.lf(ql, s1) if direction == 'right' else M.rf(s2, qr)
            groups.setdefault(g, []).append(k)
        dec = {}
        allw = []
        for g, ks in groups.items():
            rows = sorted({(k[0], k[1]) for k in ks}); cols = sorted({(k[2], k[3]) for k in ks})
            rdim = {r: self.dimL[j][r[0]] for r in rows}
            cdim = {c: self.dimR[j + 2][c[1]] for c in cols}
            ro = {}; o = 0
            for r in rows:
                ro[r] = o; o += rdim[r]
            co = {}; o2 = 0
            for c in cols:
                co[c] = o2; o2 += cdim[c]
            Mx = np.zeros((o, o2), dtype=M.dtype)
            for k in ks:
                r = (k[0], k[1]); c = (k[2], k[3])
                Mx[ro[r]:ro[r] + rdim[r], co[c]:co[c] + cdim[c]] = th[k]
            try:
                U, s, Vh = np.linalg.svd(Mx, full_matrices=False)
            except np.linalg.LinAlgError:
                U, s, Vh = np.linalg.svd(Mx + 1e-15 * self.rng.standard_normal(Mx.shape), full_matrices=False)
            dec[g] = (U, s, Vh, rows, cols, ro, co, rdim, cdim)
            for x in s:
                allw.append((x * x, g))
        tot = sum(w for w, _ in allw)
        allw.sort(key=lambda t: -t[0])
        kept = [t for t in allw[:chi] if t[0] > eps * tot]
        if not kept:
            kept = allw[:1]
        kcount = {}
        for _, g in kept:
            kcount[g] = kcount.get(g, 0) + 1
        trunc = float(max(0.0, 1.0 - sum(w for w, _ in kept) / tot))
        nk = np.sqrt(sum(w for w, _ in kept))
        newC = {}
        newdim = {}
        if direction == 'right':
            A = {}
            for g, (U, s, Vh, rows, cols, ro, co, rdim, cdim) in dec.items():
                k = kcount.get(g, 0)
                if k == 0:
                    continue
                newdim[g] = k
                for r in rows:
                    A[r] = U[ro[r]:ro[r] + rdim[r], :k].copy()
                SV = (s[:k, None] / nk) * Vh[:k]
                for c in cols:
                    newC[(g, c[0], c[1])] = SV[:, co[c]:co[c] + cdim[c]].copy()
            self.A[j] = A
            self.dimL[j + 1] = newdim
        else:
            Bn = {}
            for g, (U, s, Vh, rows, cols, ro, co, rdim, cdim) in dec.items():
                k = kcount.get(g, 0)
                if k == 0:
                    continue
                newdim[g] = k
                for c in cols:
                    Bn[c] = Vh[:k, co[c]:co[c] + cdim[c]].copy()
                US = U[:, :k] * (s[None, :k] / nk)
                for r in rows:
                    newC[(r[0], r[1], g)] = US[ro[r]:ro[r] + rdim[r], :].copy()
            self.B[j + 1] = Bn
            self.dimR[j + 1] = newdim
        self.C = newC
        return trunc, len(kept)


# ---------------------------------------------------------------- bond-matrix helpers
def bond_svd(M, centre, dimL_left, dimR_right):
    """Split a one-site centre (keys (ql, s, qr)) at its LEFT bond: centre = Cb . B with B right-canonical
    (right-physical labels) and Cb a mixed-label bond matrix (keys (ql, h)). No truncation beyond 1e-14."""
    groups = {}
    for (ql, s, qr), m in centre.items():
        h = M.rf(s, qr)
        groups.setdefault(h, []).append((ql, s, qr))
    Cb, B, dimh = {}, {}, {}
    for h, ks in groups.items():
        rows = sorted({k[0] for k in ks}); cols = sorted({(k[1], k[2]) for k in ks})
        ro = {}; o = 0
        for r in rows:
            ro[r] = o; o += dimL_left[r]
        co = {}; o2 = 0
        for c in cols:
            co[c] = o2; o2 += dimR_right[c[1]]
        Mx = np.zeros((o, o2), dtype=M.dtype)
        for (ql, s, qr) in ks:
            m = centre[(ql, s, qr)]
            Mx[ro[ql]:ro[ql] + m.shape[0], co[(s, qr)]:co[(s, qr)] + m.shape[1]] = m
        U, sv, Vh = np.linalg.svd(Mx, full_matrices=False)
        k = int(np.sum(sv > 1e-14 * max(sv.max(), 1e-300)))
        if k == 0:
            continue
        dimh[h] = k
        US = U[:, :k] * sv[None, :k]
        for r in rows:
            Cb[(r, h)] = US[ro[r]:ro[r] + dimL_left[r]].copy()
        for c in cols:
            B[c] = Vh[:k, co[c]:co[c] + dimR_right[c[1]]].copy()
    return Cb, B, dimh


def bond_pinv(Cb, rcond=1e-10):
    """pseudo-inverse of a mixed bond matrix (keys (ql, qr)), per charge sector; returns keys (qr, ql)."""
    sect = {}
    for (ql, qr), m in Cb.items():
        sect.setdefault(ql[0], []).append((ql, qr, m))
    out = {}
    for a, lst in sect.items():
        ls = sorted({x[0] for x in lst}); rs = sorted({x[1] for x in lst})
        dl = {l: [x[2].shape[0] for x in lst if x[0] == l][0] for l in ls}
        dr = {r: [x[2].shape[1] for x in lst if x[1] == r][0] for r in rs}
        ro = {}; o = 0
        for l in ls:
            ro[l] = o; o += dl[l]
        co = {}; o2 = 0
        for r in rs:
            co[r] = o2; o2 += dr[r]
        Mx = np.zeros((o, o2), dtype=lst[0][2].dtype)
        for (l, r, m) in lst:
            Mx[ro[l]:ro[l] + dl[l], co[r]:co[r] + dr[r]] = m
        P = np.linalg.pinv(Mx, rcond=rcond)
        for r in rs:
            for l in ls:
                out[(r, l)] = P[co[r]:co[r] + dr[r], ro[l]:ro[l] + dl[l]].copy()
    return out


def bond_spectrum(Cb):
    sect = {}
    for (ql, qr), m in Cb.items():
        sect.setdefault(ql[0], []).append((ql, qr, m))
    svs = []
    for a, lst in sect.items():
        ls = sorted({x[0] for x in lst}); rs = sorted({x[1] for x in lst})
        dl = {x[0]: x[2].shape[0] for x in lst}; dr = {x[1]: x[2].shape[1] for x in lst}
        ro = {}; o = 0
        for l in ls:
            ro[l] = o; o += dl[l]
        co = {}; o2 = 0
        for r in rs:
            co[r] = o2; o2 += dr[r]
        Mx = np.zeros((o, o2), dtype=lst[0][2].dtype)
        for (l, r, m) in lst:
            Mx[ro[l]:ro[l] + dl[l], co[r]:co[r] + dr[r]] = m
        sv = np.linalg.svd(Mx, compute_uv=False)
        svs.append(sv)
    s = np.concatenate(svs)
    return np.sort(s)[::-1]


# ---------------------------------------------------------------- the iDMRG driver
class IDMRG:
    def __init__(self, M, chi0=16, seed=0, verbose=True):
        self.M = M
        self.W = Win(M)
        self.verbose = verbose
        self.t = 0
        self.Eshift = 0.0
        self.hist = []
        self._init_from_finite(chi0, seed)

    def _init_from_finite(self, chi0, seed):
        """First window: C176's finite DMRG on 2q sites with N = 2p, labels converted."""
        M, W = self.M, self.W
        n, q, p = W.n, M.q, M.p
        kw = dict(M.kw)
        ch = Chain(n, 2 * p, kw['w'], t2=kw['t2'], V2=kw['V2'], V3=kw['V3'], blockade=kw['blockade'], gauge=kw['gauge'])
        if M.fermionic:
            # rebuild the MPO with the string on the NNN hop
            ent = []
            for (a, b, nzl) in ch.ent:
                if (a, b) == (B1, B2) or (a, b) == (Bd1, Bd2):
                    nzl = _nz(PAR)
                ent.append((a, b, nzl))
            ch.ent = ent
            ch.frm = [[[] for _ in range(2)] for _ in range(NCH)]
            for (a, b, nzl) in ent:
                for (so, si, v) in nzl:
                    ch.frm[a][si].append((b, so, v))
        dm = DMRG(ch, seed=seed)
        for chi in (8, 16, chi0, chi0):
            dm.sweep(chi)
        # dm: centre C at site 0 keys ((0,0), s, q2), B[1..n-1] keys (s, q2) with right labels (N_right, g)

        def cl(i, lab):
            return (q * lab[0] - p * i, lab[1])

        def cr(i, lab):
            return (q * lab[0] - p * (n - i), lab[1])
        W.C = {(cl(0, ql), s, cr(1, q2)): m.copy() for (ql, s, q2), m in dm.C.items()}
        for i in range(1, n):
            W.B[i] = {(s, cr(i + 1, q2)): m.copy() for (s, q2), m in dm.B[i].items()}
        for i in range(1, n + 1):
            W.dimR[i] = {cr(i, lab): d for lab, d in dm.dimR[i].items()} if i < n else {(0, 0): 1}
        W.dimL[0] = {(0, 0): 1}
        self.LP = [dict() for _ in range(NCH)]
        self.LP[S_][((0, 0), (0, 0))] = np.ones((1, 1), dtype=M.dtype)
        self.RP = [dict() for _ in range(NCH)]
        self.RP[E_][((0, 0), (0, 0))] = np.ones((1, 1), dtype=M.dtype)
        self.Cprev = {((0, 0), (0, 0)): np.ones((1, 1), dtype=M.dtype)}
        self.dimLprev = {(0, 0): 1}
        self.dimRprev = {(0, 0): 1}

    def _build_renv(self):
        W, M = self.W, self.M
        n = W.n
        W.Renv[n] = self.RP
        for i in range(n - 1, 0, -1):
            W.Renv[i] = rupdate(M, W.Renv[i + 1], W.B[i])
        W.Lenv[0] = self.LP

    def sweep_window(self, chi, nsw=1, tol=1e-10):
        """full right and left sweeps of the window without absorbing (centre back at site 0); used to converge
        the first window and for the finite-chain validation."""
        W, M = self.W, self.M
        n = W.n
        E = None
        for _ in range(nsw):
            self._build_renv()
            for j in range(n - 1):
                th = W.theta_right(j)
                E, th = W.solve(j, th, tol)
                W.split(j, th, chi, 'right')
                W.Lenv[j + 1] = lupdate(M, W.Lenv[j], W.A[j])
            for j in range(n - 2, -1, -1):
                th = W.theta_left(j)
                E, th = W.solve(j, th, tol)
                W.split(j, th, chi, 'left')
                W.Renv[j + 1] = rupdate(M, W.Renv[j + 2], W.B[j + 1])
        return float(np.real(E))

    def iterate(self, chi, tol=1e-10):
        """one iDMRG iteration on the current window (already prepared: centre at site 0, B[1..n-1])."""
        W, M = self.W, self.M
        n, q = W.n, M.q
        self._build_renv()
        trmax = 0.0
        E = None
        for j in range(n - 1):
            th = W.theta_right(j)
            E, th = W.solve(j, th, tol)
            tr, k = W.split(j, th, chi, 'right')
            trmax = max(trmax, tr)
            W.Lenv[j + 1] = lupdate(M, W.Lenv[j], W.A[j])
        for j in range(n - 2, q - 1, -1):
            th = W.theta_left(j)
            E, th = W.solve(j, th, tol)
            tr, k = W.split(j, th, chi, 'left')
            trmax = max(trmax, tr)
            W.Renv[j + 1] = rupdate(M, W.Renv[j + 2], W.B[j + 1])
        # centre at site q (keys (ql at bond q, s, g at bond q+1)): split at bond q
        Cb, Bq, dimh = bond_svd(M, W.C, W.dimL[q], W.dimR[q + 1])
        W.B[q] = Bq
        W.dimR[q] = dimh
        Renv_q = rupdate(M, W.Renv[q + 1], W.B[q])
        # cells and bond matrices of this iteration
        Acell = [W.A[j] for j in range(q)]
        Bcell = [W.B[j] for j in range(q, n)]
        dimLcell = [W.dimL[j] for j in range(q + 1)]      # bond 0 = previous LP basis, bond q = new
        dimRcell = [W.dimR[j] for j in range(q, n + 1)]   # bond q = new mid, bond n = previous RP basis
        spec = bond_spectrum(Cb)
        S = entropy_from_s(spec)
        Etot = float(np.real(E))
        e_site = (Etot - 0.0) / n   # after the shift below, E is the energy of the newest 2q sites
        rec = dict(t=self.t, chi=int(chi), E=Etot, e_site=e_site, S=S, trunc=trmax,
                   chi_mid=int(sum(W.dimL[q].values())), nsv=len(spec))
        self.last = dict(Acell=Acell, Bcell=Bcell, Cb=Cb, Cprev=self.Cprev, dimLcell=dimLcell, dimRcell=dimRcell,
                         spec=spec)
        # absorb
        self.LP = W.Lenv[q]
        self.RP = Renv_q
        # shift the energy so the next eigenvalue is the energy of the 2q new sites
        for (bra, ket), m in self.LP[E_].items():
            if bra == ket:
                self.LP[E_][(bra, ket)] = m - Etot * np.eye(m.shape[0], dtype=m.dtype)
        # prediction for the next window: psi = C_t B_t C_{t-1}^+ A_t C_t
        self._predict(Acell, Bcell, Cb, self.Cprev, dimLcell, dimRcell, chi)
        self.Cprev = Cb
        self.t += 1
        return rec

    def _predict(self, Acell, Bcell, Cb, Cprev, dimLcell, dimRcell, chi):
        """Build the next window, right-canonical with the centre at site 0, by a right-to-left pass of SVDs
        (truncated to 2 chi, no optimization) through C_t B_t C_{t-1}^+ A_t C_t."""
        W, M = self.W, self.M
        n, q = W.n, M.q
        Cp = bond_pinv(Cprev)              # keys (r at old mid right basis, l at old mid left basis)
        # site tensors of the new window in order, as generic (leftlab, s, rightlab) -> mat
        # sites 0..q-1: B_t cell (right labels both sides); sites q..2q-1: A_t cell (left labels both sides)
        # bond matrices: left boundary Cb (left->right), junction Cp (right->left), right boundary Cb (left->right)
        # start at the right end: X = A_t[q-1] . Cb  -> keys (ql, s, qr)
        Bnew = [dict() for _ in range(n)]
        dimRnew = [dict() for _ in range(n + 1)]
        dimRnew[n] = dict(self._dims_right_of(Cb))
        # current right object: R[(lab_left_bond, h)] = matrix, where lab_left_bond is the label type of the
        # tensor to the left; begin with Cb (left labels -> right labels at the new RP basis)
        Rmat = {(l, r): m for (l, r), m in Cb.items()}       # left label at bond n (A_t's right), right label
        Rdims_right = dimRnew[n]
        for pos in range(n - 1, -1, -1):
            if pos >= q:
                T = Acell[pos - q]                             # keys (ql, s) -> lf(ql, s)
                cen = {}
                for (ql, s), Am in T.items():
                    t = M.lf(ql, s)
                    for (l, r), Rm in Rmat.items():
                        if l != t or M.rf(s, r) is None:
                            continue
                        k = (ql, s, r)
                        cen[k] = cen[k] + Am @ Rm if k in cen else Am @ Rm
                dl = dimLcell[pos - q]
            else:
                T = Bcell[pos]                                 # keys (s, qr) -> rf(s, qr) on the left
                cen = {}
                for (s, qr), Bm in T.items():
                    hin = M.rf(s, qr)
                    for (l, r), Rm in Rmat.items():
                        if l != qr or M.rf(s, r) is None:
                            continue
                        k = (hin, s, r)
                        cen[k] = cen[k] + Bm @ Rm if k in cen else Bm @ Rm
                dl = dimRcell[pos]
            if pos == 0:
                # left boundary: multiply Cb (keys (lnew, rlabel)) on the left, giving the centre at site 0
                C0 = {}
                for (hin, s, r), m in cen.items():
                    for (l, rr), Cm in Cb.items():
                        if rr != hin:
                            continue
                        k = (l, s, r)
                        C0[k] = C0[k] + Cm @ m if k in C0 else Cm @ m
                nrm = np.sqrt(sum(np.sum(np.abs(m) ** 2) for m in C0.values()))
                W.C = {k: m / nrm for k, m in C0.items()}
                break
            # split at the left bond of pos, truncated
            Rmat, Bp, dh = self._lsplit(cen, dl, Rdims_right, 2 * chi)
            Bnew[pos] = Bp
            dimRnew[pos] = dh
            Rdims_right = dh
            if pos == q:
                # junction: Rmat keys (l at old mid left basis, h) -> multiply Cp (keys (r_old, l_old)) on the left
                R2 = {}
                for (r_old, l_old), Pm in Cp.items():
                    for (l, h), Rm in Rmat.items():
                        if l != l_old:
                            continue
                        k = (r_old, h)
                        R2[k] = R2[k] + Pm @ Rm if k in R2 else Pm @ Rm
                Rmat = R2
        for i in range(1, n):
            W.B[i] = Bnew[i]
        for i in range(1, n + 1):
            W.dimR[i] = dimRnew[i]
        W.dimL[0] = {l: m.shape[0] for (l, r), m in Cb.items()}

    @staticmethod
    def _dims_right_of(Cb):
        d = {}
        for (l, r), m in Cb.items():
            d[r] = m.shape[1]
        return d

    def _lsplit(self, cen, dl, dr, chi):
        """cen keys (lab, s, r): split at the left bond, grouping by h = rf(s, r); keep at most chi states."""
        M = self.M
        groups = {}
        for (lab, s, r) in cen:
            groups.setdefault(M.rf(s, r), []).append((lab, s, r))
        dec = {}; allw = []
        for h, ks in groups.items():
            rows = sorted({k[0] for k in ks}); cols = sorted({(k[1], k[2]) for k in ks})
            ro = {}; o = 0
            for r in rows:
                ro[r] = o; o += dl[r]
            co = {}; o2 = 0
            for c in cols:
                co[c] = o2; o2 += dr[c[1]]
            Mx = np.zeros((o, o2), dtype=M.dtype)
            for k in ks:
                m = cen[k]
                Mx[ro[k[0]]:ro[k[0]] + m.shape[0], co[(k[1], k[2])]:co[(k[1], k[2])] + m.shape[1]] = m
            U, sv, Vh = np.linalg.svd(Mx, full_matrices=False)
            dec[h] = (U, sv, Vh, rows, cols, ro, co)
            allw += [(x * x, h) for x in sv]
        tot = sum(w for w, _ in allw)
        allw.sort(key=lambda t: -t[0])
        kept = [t for t in allw[:chi] if t[0] > 1e-14 * tot] or allw[:1]
        kc = {}
        for _, h in kept:
            kc[h] = kc.get(h, 0) + 1
        Rmat, B, dh = {}, {}, {}
        for h, (U, sv, Vh, rows, cols, ro, co) in dec.items():
            k = kc.get(h, 0)
            if k == 0:
                continue
            dh[h] = k
            US = U[:, :k] * sv[None, :k]
            for r in rows:
                Rmat[(r, h)] = US[ro[r]:ro[r] + dl[r]].copy()
            for c in cols:
                B[c] = Vh[:k, co[c]:co[c] + dr[c[1]]].copy()
        return Rmat, B, dh

    # ------------------------------------------------------------ measurements on the last iteration
    def transfer_eigs(self, dQ, string=False, k=6, tol=1e-10):
        """Spectrum of T = T[A_t] M[C_t] T[B_t] M[C_{t-1}^+] (two cells, 2q sites) on (bra, ket) objects at the
        LP_{t-1} bond with ket charge = bra charge + q dQ."""
        M = self.M
        L = self.last
        Acell, Bcell, Cb, Cprev = L['Acell'], L['Bcell'], L['Cb'], L['Cprev']
        Cp = bond_pinv(Cprev)
        dims0 = L['dimLcell'][0]
        op = _nz(PAR) if string else _nz(ID)
        q = M.q
        shapes = {}
        for l1, d1 in dims0.items():
            for l2, d2 in dims0.items():
                if l2[0] == l1[0] + q * dQ:
                    shapes[(l1, l2)] = (d1, d2)
        keys = list(shapes); sh = [shapes[k] for k in keys]
        offs = np.cumsum([0] + [a * b for a, b in sh]); nv = int(offs[-1])
        if nv == 0:
            return np.array([])
        Bin = []
        for Bt in Bcell:
            ix = {}
            for (s, r), m in Bt.items():
                ix.setdefault((s, M.rf(s, r)), []).append((r, m))
            Bin.append(ix)

        def bondmap(X, Cmat):
            # X keys (l1, l2) in the row labels of Cmat; returns keys (r1, r2): sum C^dag X C
            byrow = {}
            for (l, r), m in Cmat.items():
                byrow.setdefault(l, []).append((r, m))
            out = {}
            for (l1, l2), Xm in X.items():
                for (r1, C1) in byrow.get(l1, ()):
                    Y = C1.conj().T @ Xm
                    for (r2, C2) in byrow.get(l2, ()):
                        kk = (r1, r2)
                        Z = Y @ C2
                        out[kk] = out[kk] + Z if kk in out else Z
            return out

        def btransfer(X, ix):
            out = {}
            for (b1, k1), Xm in X.items():
                for (sp, s, v) in op:
                    for (bo, Bb) in ix.get((sp, b1), ()):
                        Y = Bb.conj().T @ Xm
                        for (ko, Bk) in ix.get((s, k1), ()):
                            kk = (bo, ko)
                            Z = v * (Y @ Bk)
                            out[kk] = out[kk] + Z if kk in out else Z
            return out

        def mv(v):
            X = {keys[i]: v[offs[i]:offs[i + 1]].reshape(sh[i]) for i in range(len(keys))}
            for A in Acell:
                X = lupdate1(M, X, op, A)
            X = bondmap(X, Cb)
            for ix in Bin:
                X = btransfer(X, ix)
            X = bondmap(X, Cp)
            out = np.zeros(nv, dtype=complex)
            for i, kk in enumerate(keys):
                if kk in X and X[kk].shape == sh[i]:
                    out[offs[i]:offs[i + 1]] = X[kk].reshape(-1)
            return out
        if nv <= 300:
            T = np.array([mv(e) for e in np.eye(nv, dtype=complex)]).T
            ev = np.linalg.eigvals(T)
        else:
            opx = spla.LinearOperator((nv, nv), matvec=mv, dtype=complex)
            ev = spla.eigs(opx, k=min(k, nv - 2), which='LM', tol=tol, return_eigenvectors=False, maxiter=3000)
        ev = ev[np.argsort(-np.abs(ev))]
        return ev[:k]

    def window_density(self):
        """<n_j> on the window from the last iteration's A cell and the mid bond (left cell only)."""
        M = self.M
        L = self.last
        Acell, Cb = L['Acell'], L['Cb']
        q = M.q
        # rho at bond q (left labels): Cb Cb^dag, in (bra, ket) order conj(C) C^T
        byr = {}
        for (l, r), m in Cb.items():
            byr.setdefault(r, []).append((l, m))
        rho = {}
        for r, lst in byr.items():
            for (l1, m1) in lst:
                for (l2, m2) in lst:
                    kk = (l1, l2)
                    Y = m1.conj() @ m2.T
                    rho[kk] = rho[kk] + Y if kk in rho else Y
        # transfer from each site to the mid bond: need left objects; propagate identity from bond 0
        dims0 = L['dimLcell'][0]
        dens = []
        for j in range(q):
            X = {(l, l): np.eye(d, dtype=M.dtype) for l, d in dims0.items()}
            for i in range(q):
                X = lupdate1(M, X, _nz(NN) if i == j else _nz(ID), Acell[i])
            dens.append(float(np.real(sum(np.sum(X[kk] * rho[kk]) for kk in X if kk in rho))))
        return dens


def lupdate1(M, X, nzl, A):
    out = {}
    for (qb, qk), Xm in X.items():
        for (sp, s, v) in nzl:
            Ab = A.get((qb, sp)); Ak = A.get((qk, s))
            if Ab is None or Ak is None:
                continue
            Y = v * (Ab.conj().T @ (Xm @ Ak))
            kk = (M.lf(qb, sp), M.lf(qk, s))
            out[kk] = out[kk] + Y if kk in out else Y
    return out
