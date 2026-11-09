"""C176: blockade-exact, N-conserving two-site DMRG for the Fibonacci-cusp model (C173 PROOF eq. 3.2).

Model (open chain, units |t2| = 1), hard-core particles with no two adjacent (when blockade=True):
    H = P sum_j [ (w c+_{j+1} c_j + h.c.) + (t2 c+_{j+2} c_j + h.c.) + V2 n_j n_{j+2} + V3 n_j n_{j+3} ] P
with t2 = -1, V2 = -2, V3 = -2/3 by default. Particles never pass, so fermions = hard-core bosons.
The chemical potential is not in the MPO: N is fixed exactly, and mu enters only through the choice of N
(grand potential E(N) - mu N).

gauge=1 runs the model in the frame d_j = (-i)^j c_j (C173 PROOF 5.3): w -> -i w, t2 -> -t2. At phi = pi/2 the
Hamiltonian is then real. Observables that depend on the frame are converted back to the c frame in measure().

Labels. Left-canonical tensors A[i] carry left-physical labels (N_left, n_{i-1}) on bond i; right-canonical
tensors B[i] carry right-physical labels (N_right, n_i) on bond i. Both bases are orthonormal because the labels
are physical quantum numbers of the block they describe. The blockade (n_i n_{i+1} = 0) and N are therefore
structural zeros of every tensor: they hold exactly, not through a penalty. The truncation keeps the label
structure (it is done per label, which costs some bond dimension relative to an unconstrained SVD). The
entanglement at each cut is computed from the full mixed-label centre matrix, so it is the exact von Neumann
entropy of the represented state.

blockade=False gives ordinary U(1) DMRG (labels (N, 0)); with w = 0 this is the c = 2 control of C173
(two decoupled XX legs, the NNN hop passes freely as for hard-core bosons).
"""
import numpy as np
import scipy.sparse.linalg as spla
import time, json, sys, os, pickle

ID = np.eye(2)
BB = np.array([[0., 1.], [0., 0.]])   # b : |1> -> |0>
BD = BB.T.copy()                      # b+
NN = np.diag([0., 1.])
PAR = np.diag([1., -1.])              # (-1)^n, fermion string
S_, B1, Bd1, B2, Bd2, N1, N2, N3, E_ = range(9)
NCH = 9


def _nz(op):
    return [(s, t, op[s, t]) for s in range(2) for t in range(2) if abs(op[s, t]) > 0]


class Chain:
    def __init__(self, L, N, w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, blockade=True, gauge=0):
        self.L, self.N, self.bk, self.gauge = L, N, bool(blockade), gauge
        self.w0, self.t20 = w, t2
        if gauge:
            w, t2 = -1j * w, -t2
        w = complex(w); t2 = complex(t2)
        if abs(w.imag) < 1e-15:
            w = w.real
        if abs(t2.imag) < 1e-15:
            t2 = t2.real
        self.w, self.t2, self.V2, self.V3 = w, t2, V2, V3
        self.dtype = complex if (isinstance(w, complex) or isinstance(t2, complex)) else float
        ent = [(S_, S_, ID), (E_, E_, ID)]
        if w != 0:
            ent += [(S_, B1, BB), (S_, Bd1, BD), (B1, E_, w * BD), (Bd1, E_, np.conj(w) * BB)]
        if t2 != 0:
            if w == 0:
                ent += [(S_, B1, BB), (S_, Bd1, BD)]
            ent += [(B1, B2, ID), (Bd1, Bd2, ID), (B2, E_, t2 * BD), (Bd2, E_, np.conj(t2) * BB)]
        if V2 != 0 or V3 != 0:
            ent += [(S_, N1, NN), (N1, N2, ID)]
            if V2 != 0:
                ent += [(N2, E_, V2 * NN)]
            if V3 != 0:
                ent += [(N2, N3, ID), (N3, E_, V3 * NN)]
        self.ent = [(a, b, _nz(op)) for (a, b, op) in ent]
        # from[a][s_in] -> list of (b, s_out, val)
        self.frm = [[[] for _ in range(2)] for _ in range(NCH)]
        for (a, b, nzl) in self.ent:
            for (so, si, v) in nzl:
                self.frm[a][si].append((b, so, v))
        self._feas()

    # ---------- label bookkeeping ----------
    def _maxp(self, nsites):
        if nsites <= 0:
            return 0
        return (nsites + 1) // 2 if self.bk else nsites

    def _feas(self):
        L, N, bk = self.L, self.N, self.bk
        F = [set() for _ in range(L + 1)]
        F[0] = {(0, 0)}
        for i in range(L):
            for (n, f) in F[i]:
                for s in (0, 1):
                    if bk and f and s:
                        continue
                    if n + s <= N:
                        F[i + 1].add((n + s, s if bk else 0))
        self.Lfeas = []
        for i in range(L + 1):
            ok = set()
            for (n, f) in F[i]:
                rem = N - n
                if 0 <= rem <= self._maxp(L - i - (f if bk else 0)):
                    ok.add((n, f))
            self.Lfeas.append(ok)
        G = [set() for _ in range(L + 1)]
        G[L] = {(0, 0)}
        for i in range(L - 1, -1, -1):
            for (m, g) in G[i + 1]:
                for s in (0, 1):
                    if bk and s and g:
                        continue
                    if m + s <= N:
                        G[i].add((m + s, s if bk else 0))
        self.Rfeas = []
        for i in range(L + 1):
            ok = set()
            for (m, g) in G[i]:
                rem = N - m
                if 0 <= rem <= self._maxp(i - (g if bk else 0)):
                    ok.add((m, g))
            self.Rfeas.append(ok)
        assert (0, 0) in self.Lfeas[0] and (0, 0) in self.Rfeas[L], "N not reachable"

    def lfuse(self, i, q, s):
        n, f = q
        if self.bk and f and s:
            return None
        q2 = (n + s, s if self.bk else 0)
        return q2 if q2 in self.Lfeas[i + 1] else None

    def rfuse(self, i, s, q):
        m, g = q
        if self.bk and s and g:
            return None
        q2 = (m + s, s if self.bk else 0)
        return q2 if q2 in self.Rfeas[i] else None


# ---------- environments ----------
def lupdate(ch, Lold, A, i):
    new = [dict() for _ in range(NCH)]
    for (a, b, nzl) in ch.ent:
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
                key = (ch.lfuse(i, qb, sp), ch.lfuse(i, qk, s))
                d = new[b]
                if key in d:
                    d[key] = d[key] + X
                else:
                    d[key] = X
    return new


def rupdate(ch, Rold, B, i):
    new = [dict() for _ in range(NCH)]
    for (a, b, nzl) in ch.ent:
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
                key = (ch.rfuse(i, sp, qb2), ch.rfuse(i, s, qk2))
                d = new[a]
                if key in d:
                    d[key] = d[key] + X
                else:
                    d[key] = X
    return new


def ket_index(env):
    """env[a][(qb,qk)] -> idx[a][qk] = list of (qb, mat)."""
    out = []
    for d in env:
        ix = {}
        for (qb, qk), m in d.items():
            ix.setdefault(qk, []).append((qb, m))
        out.append(ix)
    return out


# ---------- two-site problem ----------
class TwoSite:
    def __init__(self, ch, j, dimL, dimR):
        self.ch, self.j = ch, j
        keys = []
        for ql, dl in dimL.items():
            for s1 in (0, 1):
                qm = ch.lfuse(j, ql, s1)
                if qm is None:
                    continue
                for qr, dr in dimR.items():
                    for s2 in (0, 1):
                        if ch.bk and s1 and s2:
                            continue
                        if ql[0] + s1 + s2 + qr[0] != ch.N:
                            continue
                        qmr = ch.rfuse(j + 1, s2, qr)
                        if qmr is None:
                            continue
                        keys.append(((ql, s1, s2, qr), (dl, dr)))
        self.keys = [k for k, _ in keys]
        self.shapes = [sh for _, sh in keys]
        self.offs = np.cumsum([0] + [a * b for a, b in self.shapes])
        self.kpos = {k: n for n, k in enumerate(self.keys)}
        self.size = int(self.offs[-1])

    def unpack(self, v):
        return {k: v[self.offs[n]:self.offs[n + 1]].reshape(self.shapes[n]) for n, k in enumerate(self.keys)}

    def pack(self, blocks, dtype):
        v = np.zeros(self.size, dtype=dtype)
        for k, m in blocks.items():
            n = self.kpos.get(k)
            if n is not None:
                v[self.offs[n]:self.offs[n + 1]] = m.reshape(-1)
        return v


def make_matvec(ch, ts, Lidx, Ridx, dtype):
    frm = ch.frm
    bk = ch.bk
    kpos, offs, shapes = ts.kpos, ts.offs, ts.shapes

    def mv(v):
        th = ts.unpack(v)
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
        out = np.zeros(ts.size, dtype=dtype)
        for (b, qbl, s1p, s2, qr), Xm in X.items():
            for (c, s2p, val) in frm[b][s2]:
                if bk and s1p and s2p:
                    continue
                lst = Ridx[c].get(qr)
                if not lst:
                    continue
                for (qbr, Rm) in lst:
                    n = kpos.get((qbl, s1p, s2p, qbr))
                    if n is None:
                        continue
                    out[offs[n]:offs[n + 1]] += (val * (Xm @ Rm.T)).reshape(-1)
        return out
    return mv


def entropy_from_s(s):
    p = s ** 2
    p = p[p > 1e-30]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


class DMRG:
    def __init__(self, ch, seed=0, window=3):
        self.ch = ch
        L, N = ch.L, ch.N
        rng = np.random.default_rng(seed)
        dt = ch.dtype
        # initial right-physical labels within a window around the linear profile
        self.dimR = [dict() for _ in range(L + 1)]
        self.dimL = [dict() for _ in range(L + 1)]
        for i in range(L + 1):
            tgt = N * (L - i) / L
            for q in ch.Rfeas[i]:
                if abs(q[0] - tgt) <= window + 1 or i in (0, L):
                    self.dimR[i][q] = 1
        self.dimL[0] = {(0, 0): 1}
        self.B = [dict() for _ in range(L)]
        for i in range(L - 1, 0, -1):
            for q2 in self.dimR[i + 1]:
                for s in (0, 1):
                    q = ch.rfuse(i, s, q2)
                    if q is not None and q in self.dimR[i]:
                        self.B[i][(s, q2)] = rng.standard_normal((1, 1)).astype(dt)
        C = {}
        for q2 in self.dimR[1]:
            for s in (0, 1):
                q = ch.rfuse(0, s, q2)
                if q is not None and q[0] == N:
                    C[((0, 0), s, q2)] = rng.standard_normal((1, 1)).astype(dt)
        self.C = C
        self.A = [dict() for _ in range(L)]
        self._right_normalize()
        self.Lenv = [None] * (L + 1)
        self.Renv = [None] * (L + 1)
        self.Lenv[0] = [dict() for _ in range(NCH)]
        self.Lenv[0][S_][((0, 0), (0, 0))] = np.ones((1, 1), dtype=dt)
        self.Renv[L] = [dict() for _ in range(NCH)]
        self.Renv[L][E_][((0, 0), (0, 0))] = np.ones((1, 1), dtype=dt)
        for i in range(L - 1, 0, -1):
            self.Renv[i] = rupdate(ch, self.Renv[i + 1], self.B[i], i)
        self.S = np.zeros(L - 1)
        self.trunc = np.zeros(L - 1)
        self.chi_bond = np.zeros(L - 1, dtype=int)

    def _right_normalize(self):
        ch = self.ch
        L = ch.L
        for i in range(L - 1, 0, -1):
            Bi = self.B[i]
            groups = {}
            for (s, q2), m in Bi.items():
                q = ch.rfuse(i, s, q2)
                groups.setdefault(q, []).append((s, q2))
            newdim = {}
            for q, cols in groups.items():
                X = np.hstack([Bi[c] for c in cols])
                Q, R = np.linalg.qr(X.conj().T)
                k = Q.shape[1]
                Qh = Q.conj().T
                o = 0
                for c in cols:
                    w_ = Bi[c].shape[1]
                    Bi[c] = Qh[:, o:o + w_].copy(); o += w_
                Rh = R.conj().T  # (dim q, k)
                newdim[q] = k
                if i - 1 >= 1:
                    for (s, q2) in list(self.B[i - 1].keys()):
                        if q2 == q:
                            self.B[i - 1][(s, q2)] = self.B[i - 1][(s, q2)] @ Rh
                else:
                    for key in list(self.C.keys()):
                        if key[2] == q:
                            self.C[key] = self.C[key] @ Rh
            # drop labels with no continuation
            self.dimR[i] = newdim
            if i - 1 >= 1:
                for key in list(self.B[i - 1].keys()):
                    if key[1] not in newdim:
                        del self.B[i - 1][key]
            else:
                for key in list(self.C.keys()):
                    if key[2] not in newdim:
                        del self.C[key]
        nrm = np.sqrt(sum(np.sum(np.abs(m) ** 2) for m in self.C.values()))
        for k in self.C:
            self.C[k] = self.C[k] / nrm

    # ----- core step -----
    def _solve(self, j, theta0, tol, maxiter=400):
        ch = self.ch
        ts = TwoSite(ch, j, self.dimL[j], self.dimR[j + 2])
        Lidx = ket_index(self.Lenv[j]); Ridx = ket_index(self.Renv[j + 2])
        mv = make_matvec(ch, ts, Lidx, Ridx, ch.dtype)
        v0 = ts.pack(theta0, ch.dtype)
        n = ts.size
        if n <= 64:
            Hd = np.array([mv(e) for e in np.eye(n, dtype=ch.dtype)]).T
            Hd = 0.5 * (Hd + Hd.conj().T)
            e, V = np.linalg.eigh(Hd)
            return e[0], ts.unpack(V[:, 0]), ts
        if np.linalg.norm(v0) < 1e-12:
            v0 = np.random.default_rng(j).standard_normal(n).astype(ch.dtype)
        op = spla.LinearOperator((n, n), matvec=mv, dtype=ch.dtype)
        e, V = spla.eigsh(op, k=1, which='SA', v0=v0, tol=tol, ncv=min(n - 1, 20), maxiter=maxiter)
        return e[0], ts.unpack(V[:, 0]), ts

    def _theta_right(self, j):
        """theta(j, j+1) from centre C[j] (keys (ql,s1,qmid)) and B[j+1] (keys (s2,qr))."""
        ch = self.ch
        th = {}
        Bn = self.B[j + 1]
        for (ql, s1, qmid), Cm in self.C.items():
            for s2 in (0, 1):
                pass
        # index B by left label
        byleft = {}
        for (s2, qr), m in Bn.items():
            byleft.setdefault(ch.rfuse(j + 1, s2, qr), []).append((s2, qr, m))
        for (ql, s1, qmid), Cm in self.C.items():
            for (s2, qr, m) in byleft.get(qmid, []):
                if ch.bk and s1 and s2:
                    continue
                k = (ql, s1, s2, qr)
                th[k] = th.get(k, 0) + Cm @ m
        return th

    def _theta_left(self, j):
        """theta(j, j+1) from A[j] (keys (ql,s1)) and centre C[j+1] (keys (qm,s2,qr))."""
        ch = self.ch
        th = {}
        byleft = {}
        for (qm, s2, qr), m in self.C.items():
            byleft.setdefault(qm, []).append((s2, qr, m))
        for (ql, s1), Am in self.A[j].items():
            qm = ch.lfuse(j, ql, s1)
            for (s2, qr, m) in byleft.get(qm, []):
                if ch.bk and s1 and s2:
                    continue
                k = (ql, s1, s2, qr)
                th[k] = th.get(k, 0) + Am @ m
        return th

    def _split(self, j, th, chi, direction, eps=1e-14):
        """Truncated, label-respecting split of theta at bond j+1. Returns (S, trunc, chi_kept)."""
        ch = self.ch
        groups = {}
        if direction == 'right':
            for k, m in th.items():
                ql, s1, s2, qr = k
                groups.setdefault(ch.lfuse(j, ql, s1), []).append(k)
        else:
            for k, m in th.items():
                ql, s1, s2, qr = k
                groups.setdefault(ch.rfuse(j + 1, s2, qr), []).append(k)
        dec = {}
        allw = []
        for g, ks in groups.items():
            if direction == 'right':
                rows = sorted({(k[0], k[1]) for k in ks}); cols = sorted({(k[2], k[3]) for k in ks})
            else:
                rows = sorted({(k[0], k[1]) for k in ks}); cols = sorted({(k[2], k[3]) for k in ks})
            rdim = {r: (self.dimL[j][r[0]]) for r in rows}
            cdim = {c: (self.dimR[j + 2][c[1]]) for c in cols}
            ro = {}; o = 0
            for r in rows:
                ro[r] = o; o += rdim[r]
            nr = o
            co = {}; o = 0
            for c in cols:
                co[c] = o; o += cdim[c]
            nc = o
            M = np.zeros((nr, nc), dtype=ch.dtype)
            for k in ks:
                r = (k[0], k[1]); c = (k[2], k[3])
                M[ro[r]:ro[r] + rdim[r], co[c]:co[c] + cdim[c]] = th[k]
            try:
                U, s, Vh = np.linalg.svd(M, full_matrices=False)
            except np.linalg.LinAlgError:
                U, s, Vh = np.linalg.svd(M + 1e-15 * np.random.default_rng(0).standard_normal(M.shape), full_matrices=False)
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
        if direction == 'right':
            A = {}
            newdim = {}
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
            self.C = newC
            S = self._mixed_entropy_right(newC)
        else:
            Bn = {}
            newdim = {}
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
            S = self._mixed_entropy_left(newC, j)
        return S, trunc, len(kept)

    def _mixed_entropy_right(self, C):
        """C keys (qm, s2, qr): rows = left-physical labels qm (orthonormal), cols = (s2, qr) orthonormal.
        Exact Schmidt values per N sector of the left block."""
        bysec = {}
        for (qm, s2, qr), m in C.items():
            bysec.setdefault(qm[0], []).append(((qm,), (s2, qr), m))
        return self._sector_entropy(bysec)

    def _mixed_entropy_left(self, C, j):
        """C keys (ql, s1, qmr): rows (ql, s1) orthonormal, cols qmr right-physical labels (orthonormal)."""
        bysec = {}
        for (ql, s1, qmr), m in C.items():
            bysec.setdefault(ql[0] + s1, []).append(((ql, s1), (qmr,), m))
        return self._sector_entropy(bysec)

    def _sector_entropy(self, bysec):
        svals = []
        for sec, items in bysec.items():
            rows = {}; cols = {}
            for r, c, m in items:
                rows.setdefault(r, m.shape[0]); cols.setdefault(c, m.shape[1])
            ro = {}; o = 0
            for r in sorted(rows):
                ro[r] = o; o += rows[r]
            nr = o
            co = {}; o = 0
            for c in sorted(cols):
                co[c] = o; o += cols[c]
            nc = o
            M = np.zeros((nr, nc), dtype=self.ch.dtype)
            for r, c, m in items:
                M[ro[r]:ro[r] + m.shape[0], co[c]:co[c] + m.shape[1]] += m
            svals.append(np.linalg.svd(M, compute_uv=False))
        s = np.concatenate(svals)
        self.last_spec = s
        return entropy_from_s(s)

    def sweep(self, chi, tol=1e-10):
        ch = self.ch
        L = ch.L
        trunc_max = 0.0
        E = None
        # left to right
        for j in range(L - 1):
            th = self._theta_right(j)
            E, th, ts = self._solve(j, th, tol)
            S, tr, k = self._split(j, th, chi, 'right')
            trunc_max = max(trunc_max, tr)
            self.Lenv[j + 1] = lupdate(ch, self.Lenv[j], self.A[j], j)
            if j < L - 2:
                # move centre: C[j+1] already set with keys (qm, s2, qr); qr labels bond j+2
                pass
        # now C is centre at site L-1 with keys (qm, s, qr=(0,0))
        for j in range(L - 2, -1, -1):
            th = self._theta_left(j)
            E, th, ts = self._solve(j, th, tol)
            S, tr, k = self._split(j, th, chi, 'left')
            trunc_max = max(trunc_max, tr)
            self.S[j] = S; self.trunc[j] = tr; self.chi_bond[j] = k
            self.Renv[j + 1] = rupdate(ch, self.Renv[j + 2], self.B[j + 1], j + 1)
        return float(np.real(E)), trunc_max

    # ---------- measurements (MPS right-canonical, centre C at site 0) ----------
    def _site_tensors(self):
        """List of site tensors in B format keyed (s, q_right); site 0 from the centre."""
        T0 = {}
        for (ql, s, q2), m in self.C.items():
            T0[(s, q2)] = m
        return [T0] + [self.B[i] for i in range(1, self.ch.L)]

    def _lfuse_meas(self, i, s, q2):
        if i == 0:
            return ('L', s)  # left boundary pseudo-label
        return self.ch.rfuse(i, s, q2)

    def _transfer(self, Tl, E, i, op):
        """E keyed (qb, qk) at bond i (right-physical labels; bond 0 uses ('L', s)); returns env at bond i+1."""
        new = {}
        T = Tl[i]
        bys = {}
        for (s, q2), m in T.items():
            bys.setdefault(self._lfuse_meas(i, s, q2), []).append((s, q2, m))
        for (qb, qk), Em in E.items():
            lb = bys.get(qb); lk = bys.get(qk)
            if not lb or not lk:
                continue
            for (t, qk2, mk) in lk:
                EK = None
                for (s, qb2, mb) in lb:
                    v = op[s, t]
                    if v == 0:
                        continue
                    if EK is None:
                        EK = Em @ mk
                    X = v * (mb.conj().T @ EK)
                    key = (qb2, qk2)
                    if key in new:
                        new[key] = new[key] + X
                    else:
                        new[key] = X
        return new

    @staticmethod
    def _close(E):
        tot = 0.0
        for (qb, qk), m in E.items():
            if qb == qk:
                tot += np.trace(m)
        return tot

    def measure(self, gwin=None, refs=None):
        ch = self.ch
        L = ch.L
        Tl = self._site_tensors()
        E0 = {(('L', s), ('L', t)): np.ones((1, 1)) for s in (0, 1) for t in (0, 1)}
        Eid = [E0]
        for i in range(L - 1):
            Eid.append(self._transfer(Tl, Eid[-1], i, ID))
        norm = self._close(self._transfer(Tl, Eid[L - 1], L - 1, ID))
        dens = np.array([np.real(self._close(self._transfer(Tl, Eid[i], i, NN))) for i in range(L)]) / np.real(norm)
        nnadj = np.array([np.real(self._close(self._transfer(Tl, self._transfer(Tl, Eid[i], i, NN), i + 1, NN))) for i in range(L - 1)])
        out = dict(norm=float(np.real(norm)), dens=dens.tolist(), nn_adjacent_max=float(np.max(np.abs(nnadj))))
        # one-body density matrices on a central window, boson (no string) and fermion (string (-1)^n)
        if gwin is None:
            gwin = (L // 4, 3 * L // 4)
        a0, a1 = gwin
        W = a1 - a0
        Gb = np.zeros((W, W), dtype=complex); Gf = np.zeros((W, W), dtype=complex)
        for ii in range(a0, a1):
            Gb[ii - a0, ii - a0] = dens[ii]; Gf[ii - a0, ii - a0] = dens[ii]
            Tb = self._transfer(Tl, Eid[ii], ii, BD)
            Tf = Tb
            for jj in range(ii + 1, a1):
                Gb[ii - a0, jj - a0] = self._close(self._transfer(Tl, Tb, jj, BB))
                Gf[ii - a0, jj - a0] = self._close(self._transfer(Tl, Tf, jj, BB))
                if jj + 1 < a1:
                    Tb = self._transfer(Tl, Tb, jj, ID)
                    Tf = self._transfer(Tl, Tf, jj, PAR)
            for jj in range(ii + 1, a1):
                Gb[jj - a0, ii - a0] = np.conj(Gb[ii - a0, jj - a0]); Gf[jj - a0, ii - a0] = np.conj(Gf[ii - a0, jj - a0])
        if ch.gauge:
            x = np.arange(a0, a1)
            ph = (1j) ** ((x[None, :] - x[:, None]) % 4)   # <c+_i c_j> = i^(j-i) <d+_i d_j>
            Gb = Gb * ph; Gf = Gf * ph
        out['gwin'] = [a0, a1]
        out['Gb_re'] = np.real(Gb).tolist(); out['Gb_im'] = np.imag(Gb).tolist()
        out['Gf_re'] = np.real(Gf).tolist(); out['Gf_im'] = np.imag(Gf).tolist()
        # density-density and valley-polarization correlations from reference sites
        if refs is None:
            refs = [L // 2 - 3, L // 2 - 2, L // 2 - 1, L // 2]
        # valley polarization operator P_j = c+_{j+1} c_j + h.c. (c frame); in the d frame (gauge 1)
        # P_j = -i (d+_{j+1} d_j - d+_j d_{j+1}).  As two-site products (op at j, op at j+1):
        if ch.gauge:
            Pterms = [(-1j, BB, BD), (1j, BD, BB)]
        else:
            Pterms = [(1.0, BB, BD), (1.0, BD, BB)]
        Pexp = np.zeros(L - 1)
        for i in range(L - 1):
            v = 0
            for (al, o1, o2) in Pterms:
                v += al * self._close(self._transfer(Tl, self._transfer(Tl, Eid[i], i, o1), i + 1, o2))
            Pexp[i] = np.real(v)
        out['P_exp'] = Pexp.tolist()
        cnn = {}; cpp = {}
        for i0 in refs:
            # n n
            T = self._transfer(Tl, Eid[i0], i0, NN)
            row = []
            for jj in range(i0 + 1, L):
                row.append(float(np.real(self._close(self._transfer(Tl, T, jj, NN)))) - dens[i0] * dens[jj])
                T = self._transfer(Tl, T, jj, ID)
            cnn[i0] = row
            # P P (j >= i0 + 2)
            T = None
            for (al, o1, o2) in Pterms:
                X = self._transfer(Tl, self._transfer(Tl, Eid[i0], i0, o1), i0 + 1, o2)
                X = {k: al * m for k, m in X.items()}
                if T is None:
                    T = X
                else:
                    for k, m in X.items():
                        T[k] = T[k] + m if k in T else m
            row = []
            for jj in range(i0 + 2, L - 1):
                v = 0
                for (al, o1, o2) in Pterms:
                    v += al * self._close(self._transfer(Tl, self._transfer(Tl, T, jj, o1), jj + 1, o2))
                row.append(float(np.real(v)) - Pexp[i0] * Pexp[jj])
                T = self._transfer(Tl, T, jj, ID)
            cpp[i0] = row
        # pair correlations D(i0, j) = < (c_{i0} c_{i0+d})^dag (c_j c_{j+d}) >, d = 2, 3 (c frame), j >= i0 + d + 1.
        # With the blockade the sites strictly between the two members of a pair are empty, so
        # c_j c_{j+d} = -b_j b_{j+d} on the constrained space and the strings cancel between the two pairs.
        # Frame: c_j c_{j+d} = i^(2j+d) d_j d_{j+d}, so D_c(i0, j) = (-1)^(j-i0) D_d(i0, j) in gauge 1.
        dpair = {}
        for dd in (2, 3):
            for i0 in refs:
                T = self._transfer(Tl, Eid[i0], i0, BD)
                for k in range(1, dd):
                    T = self._transfer(Tl, T, i0 + k, ID)
                T = self._transfer(Tl, T, i0 + dd, BD)
                row = []
                for jj in range(i0 + dd + 1, L - dd):
                    X = self._transfer(Tl, T, jj, BB)
                    for k in range(1, dd):
                        X = self._transfer(Tl, X, jj + k, ID)
                    v = self._close(self._transfer(Tl, X, jj + dd, BB))
                    if ch.gauge:
                        v = v * (-1) ** ((jj - i0) % 2)
                    row.append([float(np.real(v)), float(np.imag(v))])
                    T = self._transfer(Tl, T, jj, ID)
                dpair[f"{dd}_{i0}"] = row
        out['Dpair'] = dpair
        out['Cnn'] = {str(k): v for k, v in cnn.items()}
        out['Cpp'] = {str(k): v for k, v in cpp.items()}
        return out

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(dict(C=self.C, B=self.B, dimR=self.dimR, dimL=self.dimL), f)

    def load(self, path):
        with open(path, 'rb') as f:
            d = pickle.load(f)
        self.C, self.B, self.dimR = d['C'], d['B'], d['dimR']
        L = self.ch.L
        for i in range(L - 1, 0, -1):
            self.Renv[i] = rupdate(self.ch, self.Renv[i + 1], self.B[i], i)


def ground(L, N, w, chis, ch_kw=None, seed=0, verbose=False, tol=1e-10, measure_last=False, ckpt=None, record=None):
    ch = Chain(L, N, w, **(ch_kw or {}))
    dm = DMRG(ch, seed=seed)
    if ckpt and os.path.exists(ckpt):
        dm.load(ckpt)
    hist = []
    for chi in chis:
        t0 = time.time()
        E, tr = dm.sweep(chi, tol=tol)
        hist.append((int(chi), E, tr, time.time() - t0, int(dm.chi_bond.max())))
        if verbose:
            print(f"   L={L} N={N} chi={chi} E={E:.12f} trunc={tr:.2e} maxchi={dm.chi_bond.max()} t={time.time()-t0:.1f}s", flush=True)
        if record is not None:
            record(dm, hist)
        if ckpt:
            dm.save(ckpt)
    return dm, hist


if __name__ == "__main__":
    pass
