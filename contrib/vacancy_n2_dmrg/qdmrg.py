#!/usr/bin/env python3
"""Block-sparse two-site finite DMRG on open chains with abelian U(1)^k charges (numpy/scipy only).

Every local state s carries an integer charge vector, encoded linearly as one integer key
(key = sum_i q_i K^(k-1-i), K = 4096), so that sums of keys are keys of sums.  Bond bases are
labelled by the charge accumulated to their left.  The Hamiltonian is a uniform nearest-neighbour
sum H = sum_j h_{j,j+1} with h a (d^2, d^2) matrix (left site the most significant kron factor)
that conserves every charge.  h is decomposed by charge transfer as h = sum_k A_k (x) B_k, where A_k
changes the left-site charge by delta_k and B_k changes the right-site charge by -delta_k.

Left environment of the bond i (sites 0..i-1): HL (block Hamiltonian, dict q -> matrix) and
OL[k] (A_k on site i-1, dict q_ket -> matrix of shape (dim(q_ket + delta_k), dim(q_ket))).
Right environment of bond i (sites i..L-1): HR and OR[k] (B_k on site i, same convention).

The two-site effective Hamiltonian is split as H_eff = F_left + F_right, where F_left (block
Hamiltonian, the coupling to site 1, and h on sites 1, 2) conserves the right bond charge and
F_right (the coupling of site 2 to the right block, and the right block Hamiltonian) conserves the
left bond charge.  Both are assembled as dense matrices per conserved charge, so one matrix-vector
product costs two dense products per bond sector.

Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-28.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
import time
import numpy as np

KEY = 4096


def encode(qs):
    qs = np.atleast_2d(np.asarray(qs, dtype=np.int64))
    out = np.zeros(qs.shape[0], dtype=np.int64)
    for i in range(qs.shape[1]):
        out = out * KEY + qs[:, i]
    return out


def decode(key, k):
    out = []
    key = int(key)
    for _ in range(k):
        r = ((key + KEY // 2) % KEY) - KEY // 2
        out.append(r)
        key = (key - r) // KEY
    return out[::-1]


def decompose_bond(h, q, tol=1e-12):
    """h (d^2, d^2) -> list of (delta, A (d,d), B (d,d)) with h = sum A (x) B, A[out,in]."""
    d = len(q)
    T = h.reshape(d, d, d, d)                      # s1 s2 t1 t2 (out out in in)
    X = T.transpose(0, 2, 1, 3).reshape(d * d, d * d)   # (s1 t1), (s2 t2)
    dl = (q[:, None] - q[None, :]).ravel()          # delta of left factor for (s1, t1)
    dr = (q[:, None] - q[None, :]).ravel()
    terms = []
    for de in np.unique(dl):
        rows = np.flatnonzero(dl == de); cols = np.flatnonzero(dr == -de)
        blk = X[np.ix_(rows, cols)]
        if np.abs(blk).max(initial=0) < tol:
            continue
        u, s, v = np.linalg.svd(blk)
        for j in range(len(s)):
            if s[j] < tol * max(1.0, s[0]):
                continue
            A = np.zeros(d * d); B = np.zeros(d * d)
            A[rows] = u[:, j] * np.sqrt(s[j]); B[cols] = v[j] * np.sqrt(s[j])
            terms.append((int(de), A.reshape(d, d), B.reshape(d, d)))
    # check the decomposition exactly reproduces h (charge-violating entries must be zero)
    R = sum(np.kron(A, B) for _, A, B in terms)
    assert np.abs(R - h).max() < 1e-10, 'bond does not conserve the charges'
    return terms


def lanczos(mv, x0, tol, mmax=40, restarts=8):
    x = x0 / np.linalg.norm(x0)
    e = None
    for _ in range(restarts):
        V = [x]; a = []; bb = []
        w = mv(x)
        for m in range(mmax):
            al = V[-1] @ w; a.append(al)
            w = w - al * V[-1] - (bb[-1] * V[-2] if m > 0 else 0)
            Vm = np.array(V)
            w -= Vm.T @ (Vm @ w)
            be = np.linalg.norm(w)
            if m % 3 == 2 or be < 1e-13 or m == mmax - 1:
                T = np.diag(a) + np.diag(bb, 1) + np.diag(bb, -1)
                ee, y = np.linalg.eigh(T)
                res = abs(be * y[-1, 0])
                if res < tol or be < 1e-13 or m == mmax - 1:
                    break
            bb.append(be); V.append(w / be)
            w = mv(V[-1])
        x = np.array(V[:len(y)]).T @ y[:, 0]
        x /= np.linalg.norm(x)
        e = ee[0]
        if res < tol or be < 1e-13:
            break
    return e, x, res


class QDMRG:
    def __init__(self, h, qloc, L, config):
        """h: bond matrix; qloc: (d, k) integer charges of the local states; config: product-state
        local states (length L) fixing the charge sector."""
        self.h = np.asarray(h, float)
        self.qv = np.asarray(qloc, dtype=np.int64)
        self.k = self.qv.shape[1]
        self.q = encode(self.qv)
        self.d = d = len(self.q)
        self.L = L
        self.terms = decompose_bond(self.h, self.q)
        self.dk = [t[0] for t in self.terms]
        self.Ak = [t[1] for t in self.terms]
        self.Bk = [t[2] for t in self.terms]
        self.nzA = [[(t, s, A[t, s]) for t in range(d) for s in range(d) if A[t, s] != 0] for A in self.Ak]
        self.nzB = [[(t, s, B[t, s]) for t in range(d) for s in range(d) if B[t, s] != 0] for B in self.Bk]
        # pair combos
        self.combos = {}
        for s1 in range(d):
            for s2 in range(d):
                self.combos.setdefault(int(self.q[s1] + self.q[s2]), []).append((s1, s2))
        self.pidx = {c: {p: i for i, p in enumerate(P)} for c, P in self.combos.items()}
        hh = self.h.reshape(d, d, d, d)
        self.hc = {c: np.array([[hh[t1, t2, s1, s2] for (s1, s2) in P] for (t1, t2) in P]) for c, P in self.combos.items()}
        # PB[k][c] : (np(c - delta), np(c)) acting on s1 with B_k ; PA[k][c]: (np(c + delta), np(c)) on s2 with A_k
        self.PB = []; self.PA = []
        for kk in range(len(self.terms)):
            de = self.dk[kk]
            pb = {}; pa = {}
            for c, P in self.combos.items():
                c2 = c - de
                if c2 in self.combos:
                    M = np.zeros((len(self.combos[c2]), len(P)))
                    for j, (s1, s2) in enumerate(P):
                        for t1 in range(d):
                            if self.Bk[kk][t1, s1] != 0 and (t1, s2) in self.pidx[c2]:
                                M[self.pidx[c2][(t1, s2)], j] += self.Bk[kk][t1, s1]
                    if np.any(M):
                        pb[c] = M
                c3 = c + de
                if c3 in self.combos:
                    M = np.zeros((len(self.combos[c3]), len(P)))
                    for j, (s1, s2) in enumerate(P):
                        for t2 in range(d):
                            if self.Ak[kk][t2, s2] != 0 and (s1, t2) in self.pidx[c3]:
                                M[self.pidx[c3][(s1, t2)], j] += self.Ak[kk][t2, s2]
                    if np.any(M):
                        pa[c] = M
            self.PB.append(pb); self.PA.append(pa)
        # product state
        config = list(config)
        assert len(config) == L
        self.bonds = [{0: 1}]
        acc = 0
        for s in config:
            acc += int(self.q[s])
            self.bonds.append({acc: 1})
        self.Qtot = acc
        self.A = []
        for i, s in enumerate(config):
            qa = next(iter(self.bonds[i]))
            self.A.append({(qa, s): np.ones((1, 1))})
        nk = len(self.terms)
        self.HL = [None] * (L + 1); self.OL = [None] * (L + 1)
        self.HR = [None] * (L + 1); self.OR = [None] * (L + 1)
        self.HL[0] = {0: np.zeros((1, 1))}; self.OL[0] = [dict() for _ in range(nk)]
        self.HR[L] = {self.Qtot: np.zeros((1, 1))}; self.OR[L] = [dict() for _ in range(nk)]
        for i in range(L - 1, 1, -1):
            self._update_R(i)
        self.Sbond = [0.0] * (L + 1)
        self.spec = [None] * (L + 1)
        self.trunc = 0.0
        self.dens = np.zeros((L, self.k))

    # ------------------------------------------------------------------ environments
    def _update_L(self, i):
        A = self.A[i]; HL = self.HL[i]; OL = self.OL[i]; q = self.q
        H = {}; O = [dict() for _ in self.terms]

        def add(D, key, M):
            if key in D:
                D[key] += M
            else:
                D[key] = M.copy()
        for (qa, s), M in A.items():
            qb = qa + int(q[s])
            add(H, qb, M.T @ HL[qa] @ M)
        for kk in range(len(self.terms)):
            de = self.dk[kk]
            Ok = OL[kk]
            for (t, s, cf) in self.nzB[kk]:
                for qa, X in Ok.items():
                    Ms = A.get((qa, s)); Mt = A.get((qa + de, t))
                    if Ms is None or Mt is None:
                        continue
                    add(H, qa + int(q[s]), cf * (Mt.T @ (X @ Ms)))
            for (t, s, cf) in self.nzA[kk]:
                for (qa, s_), Ms in A.items():
                    if s_ != s:
                        continue
                    Mt = A.get((qa, t))
                    if Mt is None:
                        continue
                    add(O[kk], qa + int(q[s]), cf * (Mt.T @ Ms))
        self.HL[i + 1] = H; self.OL[i + 1] = O

    def _update_R(self, i):
        B = self.A[i]; HR = self.HR[i + 1]; OR = self.OR[i + 1]; q = self.q
        H = {}; O = [dict() for _ in self.terms]

        def add(D, key, M):
            if key in D:
                D[key] += M
            else:
                D[key] = M.copy()
        for (qa, s), M in B.items():
            qb = qa + int(q[s])
            add(H, qa, M @ HR[qb] @ M.T)
        for kk in range(len(self.terms)):
            de = self.dk[kk]
            Ok = OR[kk]
            for (t, s, cf) in self.nzA[kk]:
                for (qa, s_), Ms in B.items():
                    if s_ != s:
                        continue
                    X = Ok.get(qa + int(q[s])); Mt = B.get((qa, t))
                    if X is None or Mt is None:
                        continue
                    add(H, qa, cf * (Mt @ (X @ Ms.T)))
            for (t, s, cf) in self.nzB[kk]:
                for (qa, s_), Ms in B.items():
                    if s_ != s:
                        continue
                    Mt = B.get((qa + de, t))
                    if Mt is None:
                        continue
                    add(O[kk], qa, cf * (Mt @ Ms.T))
        self.HR[i] = H; self.OR[i] = O

    # ------------------------------------------------------------------ two-site problem
    def _plan(self, i):
        Ls = self.bonds[i]; Rs = self.bonds[i + 2]
        HL, OL, HR, OR = self.HL[i], self.OL[i], self.HR[i + 2], self.OR[i + 2]
        blocks = {}          # (qb, c) -> (offset, np, da, db)
        byqb = {}            # qb -> (start, R, [(c, rowoff)])
        off = 0
        for qb, db in Rs.items():
            start = off; rows = 0; lst = []
            for c, P in self.combos.items():
                qa = qb - c
                if qa in Ls:
                    da = Ls[qa]
                    blocks[(qb, c)] = (off, len(P), da, db)
                    lst.append((c, rows))
                    off += len(P) * da * db
                    rows += len(P) * da
            if rows:
                byqb[qb] = (start, rows, lst)
        n = off
        # left operators per qb
        Fl = {}
        for qb, (start, R, lst) in byqb.items():
            F = np.zeros((R, R))
            rowoff = dict(lst)
            for c, ro in lst:
                qa = qb - c; da = Ls[qa]; npc = len(self.combos[c])
                sl = slice(ro, ro + npc * da)
                F4 = F[sl, sl].reshape(npc, da, npc, da)
                ap_, aa_ = np.arange(npc), np.arange(da)
                F4[ap_, :, ap_, :] += HL[qa][None]
                F4[:, aa_, :, aa_] += self.hc[c][None]
                F[sl, sl] = F4.reshape(npc * da, npc * da)
                for kk in range(len(self.terms)):
                    pb = self.PB[kk].get(c)
                    if pb is None:
                        continue
                    X = OL[kk].get(qa)
                    c2 = c - self.dk[kk]
                    if X is None or c2 not in rowoff:
                        continue
                    ro2 = rowoff[c2]
                    n2, d2 = pb.shape[0], X.shape[0]
                    F[ro2:ro2 + n2 * d2, sl] += (pb[:, None, :, None] * X[None, :, None, :]).reshape(n2 * d2, npc * da)
            Fl[qb] = F
        # right operators per qa, with gather indices
        Fr = {}; gidx = {}
        for qa, da in Ls.items():
            cols = []; coloff = {}; co = 0
            for c, P in self.combos.items():
                qb = qa + c
                if (qb, c) in blocks:
                    o, npc, da_, db = blocks[(qb, c)]
                    coloff[c] = co
                    ix = o + (np.arange(npc)[None, :, None] * da + np.arange(da)[:, None, None]) * db + np.arange(db)[None, None, :]
                    cols.append(ix.reshape(da, npc * db))
                    co += npc * db
            if not cols:
                continue
            C = co
            F = np.zeros((C, C))
            for c, cofs in coloff.items():
                qb = qa + c; npc = len(self.combos[c]); db = Rs[qb]
                sl = slice(cofs, cofs + npc * db)
                F4 = F[sl, sl].reshape(npc, db, npc, db)
                ap_ = np.arange(npc)
                F4[ap_, :, ap_, :] += HR[qb][None]
                F[sl, sl] = F4.reshape(npc * db, npc * db)
                for kk in range(len(self.terms)):
                    pa = self.PA[kk].get(c)
                    if pa is None:
                        continue
                    X = OR[kk].get(qb)
                    c2 = c + self.dk[kk]
                    if X is None or c2 not in coloff:
                        continue
                    co2 = coloff[c2]
                    n2, d2 = pa.shape[0], X.shape[0]
                    F[co2:co2 + n2 * d2, sl] += (pa[:, None, :, None] * X[None, :, None, :]).reshape(n2 * d2, npc * db)
            Fr[qa] = F.T.copy()
            gidx[qa] = np.concatenate(cols, axis=1)
        return dict(n=n, blocks=blocks, byqb=byqb, Fl=Fl, Fr=Fr, gidx=gidx, Ls=Ls, Rs=Rs)

    def _matvec(self, P, x):
        y = np.empty_like(x)
        for qb, (start, R, lst) in P['byqb'].items():
            db = P['Rs'][qb]
            y[start:start + R * db] = (P['Fl'][qb] @ x[start:start + R * db].reshape(R, db)).ravel()
        for qa, ix in P['gidx'].items():
            y[ix] += x[ix] @ P['Fr'][qa]
        return y

    def _theta(self, i, P):
        x = np.zeros(P['n'])
        A1, A2 = self.A[i], self.A[i + 1]
        q = self.q
        for (qb, c), (o, npc, da, db) in P['blocks'].items():
            qa = qb - c
            blk = np.zeros((npc, da, db))
            nzero = False
            for j, (s1, s2) in enumerate(self.combos[c]):
                M1 = A1.get((qa, s1))
                if M1 is None:
                    continue
                M2 = A2.get((qa + int(q[s1]), s2))
                if M2 is None:
                    continue
                blk[j] = M1 @ M2; nzero = True
            if nzero:
                x[o:o + npc * da * db] = blk.ravel()
        return x

    def _split(self, i, P, x, chi, cutoff, right):
        q = self.q
        Ls, Rs = P['Ls'], P['Rs']
        # middle sectors
        rows = {}; cols = {}
        for qa, da in Ls.items():
            for s1 in range(self.d):
                rows.setdefault(qa + int(q[s1]), []).append((qa, s1, da))
        for qb, db in Rs.items():
            for s2 in range(self.d):
                cols.setdefault(qb - int(q[s2]), []).append((s2, qb, db))
        res = []
        for qm in rows:
            if qm not in cols:
                continue
            R = sum(r[2] for r in rows[qm]); C = sum(c[2] for c in cols[qm])
            M = np.zeros((R, C))
            r0 = 0
            for (qa, s1, da) in rows[qm]:
                c0 = 0
                for (s2, qb, db) in cols[qm]:
                    c = int(q[s1] + q[s2])
                    o, npc, _, _ = P['blocks'][(qb, c)]
                    p = self.pidx[c][(s1, s2)]
                    M[r0:r0 + da, c0:c0 + db] = x[o + p * da * db: o + (p + 1) * da * db].reshape(da, db)
                    c0 += db
                r0 += da
            if not np.any(M):
                continue
            try:
                u, s, v = np.linalg.svd(M, full_matrices=False)
            except np.linalg.LinAlgError:
                import scipy.linalg as sl
                u, s, v = sl.svd(M, full_matrices=False, lapack_driver='gesvd')
            res.append((qm, u, s, v))
        allS = np.concatenate([r[2] for r in res])
        nrm2 = np.sum(allS ** 2)
        order = np.sort(allS)[::-1]
        w = np.cumsum(order ** 2) / nrm2
        keep = min(chi, int(np.sum(order > 1e-15 * order[0])))
        k2 = int(np.searchsorted(w, 1 - cutoff)) + 1
        keep = max(1, min(keep, k2))
        thr = order[keep - 1]
        trunc = float(1 - w[keep - 1])
        self.trunc = max(self.trunc, trunc)
        # ties: keep exactly `keep` values
        newA1, newA2, bond = {}, {}, {}
        kept_s = []; kept_q = []
        budget = keep
        for (qm, u, s, v) in sorted(res, key=lambda r: -r[2][0]):
            m = int(np.sum(s > thr)) + int(np.sum(s == thr))
            m = min(m, budget)
            if m <= 0:
                continue
            budget -= m
            bond[qm] = m
            kept_s.append(s[:m]); kept_q += [qm] * m
            U = u[:, :m]; Vh = v[:m]; S = s[:m]
            r0 = 0
            for (qa, s1, da) in rows[qm]:
                blk = U[r0:r0 + da]
                newA1[(qa, s1)] = blk * S[None, :] if not right else blk.copy()
                r0 += da
            c0 = 0
            for (s2, qb, db) in cols[qm]:
                blk = Vh[:, c0:c0 + db]
                newA2[(qm, s2)] = S[:, None] * blk if right else blk.copy()
                c0 += db
        ks = np.concatenate(kept_s); nk = np.sqrt(np.sum(ks ** 2))
        for key in newA1:
            if not right:
                newA1[key] /= nk
        for key in newA2:
            if right:
                newA2[key] /= nk
        bond = dict(sorted(bond.items()))
        self.bonds[i + 1] = bond
        self.A[i] = newA1; self.A[i + 1] = newA2
        p = (ks / nk) ** 2
        p = p[p > 1e-300]
        self.Sbond[i + 1] = float(-np.sum(p * np.log(p)))
        self.spec[i + 1] = (np.array(kept_q), ks / nk)
        return trunc

    def _densities(self, i, P, x):
        """<charge vector> on sites i and i+1 from the normalized two-site state."""
        acc1 = np.zeros(self.k); acc2 = np.zeros(self.k)
        for (qb, c), (o, npc, da, db) in P['blocks'].items():
            w = (x[o:o + npc * da * db].reshape(npc, -1) ** 2).sum(1)
            for j, (s1, s2) in enumerate(self.combos[c]):
                acc1 += w[j] * self.qv[s1]; acc2 += w[j] * self.qv[s2]
        return acc1, acc2

    def step(self, i, chi, cutoff, tol, right, noise, rng):
        P = self._plan(i)
        x0 = self._theta(i, P)
        if not np.any(x0):
            x0 = rng.normal(size=P['n'])
        if P['n'] <= 64:
            H = np.column_stack([self._matvec(P, e) for e in np.eye(P['n'])])
            ee, vv = np.linalg.eigh((H + H.T) / 2); E, x = ee[0], vv[:, 0]
        else:
            E, x, _ = lanczos(lambda v: self._matvec(P, v), x0, tol)
        if noise > 0:      # random admixture of relative norm `noise` (sector mixing)
            g = rng.normal(size=x.shape)
            x = x + noise * g / np.linalg.norm(g)
            x /= np.linalg.norm(x)
        d1, d2 = self._densities(i, P, x)
        self.dens[i] = d1; self.dens[i + 1] = d2
        self._split(i, P, x, chi, cutoff, right)
        return E

    def sweep(self, chi, cutoff=1e-12, tol=1e-8, noise=0.0, rng=None):
        rng = rng or np.random.default_rng(0)
        L = self.L
        self.trunc = 0.0
        for i in range(0, L - 1):
            E = self.step(i, chi, cutoff, tol, True, noise, rng)
            self._update_L(i)
        for i in range(L - 2, -1, -1):
            E = self.step(i, chi, cutoff, tol, False, noise, rng)
            self._update_R(i + 1)
        self.E = E
        return E

    def get_state(self):
        return dict(A=self.A, bonds=self.bonds)

    def set_state(self, st):
        """Load a right-canonical state (centre on site 0), as left by sweep(), and rebuild the right
        environments."""
        self.A = st['A']; self.bonds = st['bonds']
        for i in range(self.L - 1, 1, -1):
            self._update_R(i)

    def chis(self):
        return [sum(b.values()) for b in self.bonds]


def run(h, qloc, L, config, chis, sweeps=2, maxsweeps=8, etol=1e-8, noise0=1e-3, log=print, on_chi=None, seed=0,
        noise_sweeps=2, state=None):
    dm = QDMRG(h, qloc, L, config)
    if state is not None:
        dm.set_state(state)
    rng = np.random.default_rng(seed)
    hist = []
    t0 = time.time()
    for ci, chi in enumerate(chis):
        Eold = None
        for sw in range(maxsweeps):
            nz = noise0 * (0.1 ** sw) if sw < noise_sweeps else 0.0     # mixing at every new chi
            E = dm.sweep(chi, noise=nz, rng=rng)
            log(f'chi={chi} sw={sw} E={E:.12f} trunc={dm.trunc:.2e} maxchi={max(dm.chis())} '
                f'S_mid={dm.Sbond[L // 2]:.6f} t={time.time() - t0:.0f}s')
            if sw + 1 >= sweeps and Eold is not None and abs(E - Eold) < etol * L and nz == 0:
                break
            Eold = E
        rec = dict(chi=chi, E=float(E), trunc=float(dm.trunc), sweeps=sw + 1, maxchi=int(max(dm.chis())),
                   S=list(dm.Sbond), dens=dm.dens.tolist(), time_s=time.time() - t0)
        hist.append(rec)
        if on_chi:
            on_chi(dm, rec)
    return dm, hist
