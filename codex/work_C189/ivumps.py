"""C189: blockade-exact, U(1) VUMPS for C173's effective model (PROOF eq. 3.2) on the infinite chain at fixed density n = p/q.

Hamiltonian: the MPO of C176's fibdmrg.Chain (imported, so the model is the same object), units |t2| = 1.
Option fermionic=True puts the string (-1)^n on the middle site of the NNN hop, which turns the blockade-free
model into free fermions (control ff2v); with the blockade the middle site is empty and the option is moot.

Labels (PROOF section 6). The uniform MPS carries shifted integer charges. A left-physical label at the bond left of
site x is (a, f) with a = q N_left - p x and f = n_{x-1}; a right-physical label is (b, g) with b = q N_right - p m
(m = number of sites to the right) and g = n_x. One step changes a by q s - p, so a mod q runs through all q residues
and labels the position in the q-site unit cell: a single dict of blocks stores the whole cell. Total charge
conservation at density p/q reads a + b = 0 across a bond and a + q s + b = p across a site. The blockade and N are
structural zeros, as in C176. A_L carries left labels on both bonds, A_R right labels, C and A_C mixed labels.

Algorithm: parallel VUMPS (Zauner-Stauber et al. 2018) with per-residue eigensolves (H_AC and H_C are block
diagonal in the residue, i.e. in the cell position), per-block polar decompositions, environments from GMRES, and
two-site subspace expansion to raise the bond dimension.
"""
import numpy as np
import scipy.sparse.linalg as spla
import time, json, sys, os, pickle
from fibdmrg import Chain, ID, BB, BD, NN, PAR, S_, B1, Bd1, B2, Bd2, N1, N2, N3, E_, NCH


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
        self.p, self.q, self.bk = p, q, bool(blockade)
        self.w, self.t2, self.V2, self.V3 = ch.w, ch.t2, V2, V3

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

    def fl(self):
        return (0, 1) if self.bk else (0,)


# ------------------------------------------------------------------ transfers
def ltransfer(M, Ld, nzl, AL):
    out = {}
    for (bra, ket), X in Ld.items():
        for (sp, s, v) in nzl:
            Ab = AL.get((bra, sp)); Ak = AL.get((ket, s))
            if Ab is None or Ak is None:
                continue
            Y = v * (Ab.conj().T @ (X @ Ak))
            k = (M.lf(bra, sp), M.lf(ket, s))
            if k in out:
                out[k] += Y
            else:
                out[k] = Y
    return out


def rtransfer(M, Rd, nzl, AR):
    out = {}
    for (bra, ket), X in Rd.items():
        for (sp, s, v) in nzl:
            Bb = AR.get((sp, bra)); Bk = AR.get((s, ket))
            if Bb is None or Bk is None:
                continue
            Y = v * (Bb.conj() @ (X @ Bk.T))
            k = (M.rf(sp, bra), M.rf(s, ket))
            if k in out:
                out[k] += Y
            else:
                out[k] = Y
    return out


def dadd(d, e, c=1.0):
    for k, v in e.items():
        if k in d:
            d[k] = d[k] + c * v
        else:
            d[k] = c * v
    return d


def pair(x, y):
    t = 0.0
    for k, v in x.items():
        if k in y:
            t += np.sum(v * y[k])
    return t


class Vec:
    """vectorize a dict of blocks over a fixed key list"""
    def __init__(self, shapes):
        self.keys = list(shapes.keys())
        self.shapes = [shapes[k] for k in self.keys]
        self.offs = np.cumsum([0] + [a * b for a, b in self.shapes])
        self.n = int(self.offs[-1])
        self.kpos = {k: i for i, k in enumerate(self.keys)}

    def pack(self, d, dtype):
        v = np.zeros(self.n, dtype=dtype)
        for i, k in enumerate(self.keys):
            if k in d:
                v[self.offs[i]:self.offs[i + 1]] = d[k].reshape(-1)
        return v

    def unpack(self, v):
        return {k: v[self.offs[i]:self.offs[i + 1]].reshape(self.shapes[i]) for i, k in enumerate(self.keys)}


# ------------------------------------------------------------------ the state
class UMPS:
    def __init__(self, M, window=1, seed=0):
        self.M = M
        q = M.q
        rng = np.random.default_rng(seed)
        self.rng = rng
        self.dimL = {}; self.dimR = {}
        for r in range(q):
            for k in range(-window, window + 1):
                a = r + q * k - (q if r > q // 2 else 0)
                for f in M.fl():
                    self.dimL[(a, f)] = 1
                    self.dimR[(-a, f)] = 1
        self.prune()
        self.AC = self.random_AC(); self.C = self.random_C()
        self.normalize()
        self.AL, self.AR = {}, {}
        self.polar()

    def prune(self):
        M = self.M
        changed = True
        while changed:
            changed = False
            for t in list(self.dimL):
                inc = any(M.lf(l, s) == t for l in self.dimL for s in (0, 1))
                out = any(M.lf(t, s) in self.dimL for s in (0, 1))
                cc = any((-t[0], g) in self.dimR and not (M.bk and t[1] and g) for g in M.fl())
                if not (inc and out and cc):
                    del self.dimL[t]; changed = True
            for h in list(self.dimR):
                inc = any(M.rf(s, r) == h for r in self.dimR for s in (0, 1))
                out = any(M.rf(s, h) in self.dimR for s in (0, 1))
                cc = any((-h[0], f) in self.dimL and not (M.bk and h[1] and f) for f in M.fl())
                if not (inc and out and cc):
                    del self.dimR[h]; changed = True

    # key enumerations
    def ac_keys(self):
        M = self.M
        out = {}
        for l, dl in self.dimL.items():
            for s in (0, 1):
                if M.bk and l[1] and s:
                    continue
                b = M.p - l[0] - M.q * s
                for g in M.fl():
                    if M.bk and s and g:
                        continue
                    r = (b, g)
                    if r in self.dimR:
                        out[(l, s, r)] = (dl, self.dimR[r])
        return out

    def c_keys(self):
        M = self.M
        out = {}
        for l, dl in self.dimL.items():
            for g in M.fl():
                if M.bk and l[1] and g:
                    continue
                r = (-l[0], g)
                if r in self.dimR:
                    out[(l, r)] = (dl, self.dimR[r])
        return out

    def random_AC(self):
        dt = self.M.dtype
        return {k: self.rng.standard_normal(sh).astype(dt) for k, sh in self.ac_keys().items()}

    def random_C(self):
        dt = self.M.dtype
        return {k: self.rng.standard_normal(sh).astype(dt) for k, sh in self.c_keys().items()}

    def res(self, a):
        return a % self.M.q

    def normalize(self):
        q = self.M.q
        nc = np.zeros(q); na = np.zeros(q)
        for (l, r), m in self.C.items():
            nc[self.res(l[0])] += np.sum(np.abs(m) ** 2)
        for (l, s, r), m in self.AC.items():
            na[self.res(l[0])] += np.sum(np.abs(m) ** 2)
        for k in self.C:
            self.C[k] = self.C[k] / np.sqrt(nc[self.res(k[0][0])])
        for k in self.AC:
            self.AC[k] = self.AC[k] / np.sqrt(na[self.res(k[0][0])])

    # polar decompositions per block
    @staticmethod
    def _polar(X):
        U, s, Vh = np.linalg.svd(X, full_matrices=False)
        return U @ Vh

    def polar(self):
        M = self.M
        # A_L: target label t = lf(l, s); rows (l) with lf(l,s) = t; X = sum_r AC[l,s,r] C[t,r]^dag
        groups = {}
        for (l, s, r), m in self.AC.items():
            t = M.lf(l, s)
            if t is None or t not in self.dimL:
                continue
            groups.setdefault(t, {}).setdefault((l, s), []).append(r)
        AL = {}
        for t, rows in groups.items():
            rl = sorted(rows)
            blocks = []
            for (l, s) in rl:
                X = np.zeros((self.dimL[l], self.dimL[t]), dtype=M.dtype)
                for r in rows[(l, s)]:
                    Cm = self.C.get((t, r))
                    if Cm is not None:
                        X += self.AC[(l, s, r)] @ Cm.conj().T
                blocks.append(X)
            U = self._polar(np.vstack(blocks))
            o = 0
            for (l, s), X in zip(rl, blocks):
                AL[(l, s)] = U[o:o + X.shape[0]].copy(); o += X.shape[0]
        self.AL = AL
        # A_R: source label h = rf(s, r); cols (r) with rf(s, r) = h; Y = sum_l C[l,h]^dag AC[l,s,r]
        groups = {}
        for (l, s, r), m in self.AC.items():
            h = M.rf(s, r)
            if h is None or h not in self.dimR:
                continue
            groups.setdefault(h, {}).setdefault((s, r), []).append(l)
        AR = {}
        for h, cols in groups.items():
            cl = sorted(cols)
            blocks = []
            for (s, r) in cl:
                Y = np.zeros((self.dimR[h], self.dimR[r]), dtype=M.dtype)
                for l in cols[(s, r)]:
                    Cm = self.C.get((l, h))
                    if Cm is not None:
                        Y += Cm.conj().T @ self.AC[(l, s, r)]
                blocks.append(Y)
            U = self._polar(np.hstack(blocks))
            o = 0
            for (s, r), Y in zip(cl, blocks):
                AR[(s, r)] = U[:, o:o + Y.shape[1]].copy(); o += Y.shape[1]
        self.AR = AR

    def errors(self):
        """||A_C - A_L C|| and ||A_C - C A_R|| (summed over the cell, squared then rooted)."""
        M = self.M
        eL = 0.0; eR = 0.0
        for (l, s, r), m in self.AC.items():
            t = M.lf(l, s); h = M.rf(s, r)
            X = np.zeros_like(m)
            if (l, s) in self.AL and (t, r) in self.C:
                X = self.AL[(l, s)] @ self.C[(t, r)]
            eL += np.sum(np.abs(m - X) ** 2)
            Y = np.zeros_like(m)
            if (s, r) in self.AR and (l, h) in self.C:
                Y = self.C[(l, h)] @ self.AR[(s, r)]
            eR += np.sum(np.abs(m - Y) ** 2)
        return float(np.sqrt(eL)), float(np.sqrt(eR))

    # ------------------------------------------------------------ fixed points and environments
    def rho_L(self):
        """right fixed point of the A_L transfer, in (bra, ket) left labels: conj(C) C^T"""
        out = {}
        for (l1, r), m1 in self.C.items():
            for g in self.M.fl():
                pass
        byr = {}
        for (l, r), m in self.C.items():
            byr.setdefault(r, []).append((l, m))
        for r, lst in byr.items():
            for (l1, m1) in lst:
                for (l2, m2) in lst:
                    if l1[0] != l2[0]:
                        continue
                    k = (l1, l2)
                    Y = m1.conj() @ m2.T
                    out[k] = out[k] + Y if k in out else Y
        return out

    def rho_R(self):
        """left fixed point of the A_R transfer, in (bra, ket) right labels: C^dag C"""
        out = {}
        byl = {}
        for (l, r), m in self.C.items():
            byl.setdefault(l, []).append((r, m))
        for l, lst in byl.items():
            for (r1, m1) in lst:
                for (r2, m2) in lst:
                    if r1[0] != r2[0]:
                        continue
                    k = (r1, r2)
                    Y = m1.conj().T @ m2
                    out[k] = out[k] + Y if k in out else Y
        return out

    def _solve_env(self, Y, rho, ident, transfer, x0, tol):
        """solve x - T(x) + (x|rho) 1 = Y - e 1, e = (Y|rho)/(1|rho)"""
        M = self.M
        nrm = pair(ident, rho)
        e = pair(Y, rho) / nrm
        b = dadd(dict(Y), ident, -e)
        dims = {}
        for (k1, k2), m in ident.items():
            dims[k1] = m.shape[0]
        bych = {}
        for lab in dims:
            bych.setdefault(lab[0], []).append(lab)
        shapes = {}
        for c, labs in bych.items():
            for l1 in labs:
                for l2 in labs:
                    shapes[(l1, l2)] = (dims[l1], dims[l2])
        V = Vec(shapes)
        dt = np.result_type(M.dtype, float)

        def mv(v):
            x = V.unpack(v)
            Tx = transfer(x)
            out = dict(x)
            dadd(out, Tx, -1.0)
            dadd(out, ident, pair(x, rho))
            return V.pack(out, dt)
        op = spla.LinearOperator((V.n, V.n), matvec=mv, dtype=dt)
        bv = V.pack(b, dt)
        v0 = V.pack(x0, dt) if x0 is not None else None
        sol, info = spla.gmres(op, bv, x0=v0, rtol=tol, atol=0.0, restart=40, maxiter=40)
        res = np.linalg.norm(mv(sol) - bv) / max(np.linalg.norm(bv), 1e-300)
        return V.unpack(sol), float(np.real(e)), res, info

    def environments(self, tol=1e-10):
        M = self.M
        AL, AR = self.AL, self.AR
        idL = {(l, l): np.eye(d, dtype=M.dtype) for l, d in self.dimL.items()}
        idR = {(r, r): np.eye(d, dtype=M.dtype) for r, d in self.dimR.items()}
        idnz = _nz(ID)
        # left
        Lw = [dict() for _ in range(NCH)]
        Lw[S_] = idL
        YL = {}
        for b in range(1, NCH):
            acc = {}
            for (a, bb, nzl) in M.ent:
                if bb != b or a == b:
                    continue
                dadd(acc, ltransfer(M, Lw[a], nzl, AL))
            if b == E_:
                YL = acc
            else:
                Lw[b] = acc
        x0 = getattr(self, '_LE', None)
        Lw[E_], eL, resL, infoL = self._solve_env(YL, self.rho_L(), idL, lambda x: ltransfer(M, x, idnz, AL), x0, tol)
        self._LE = Lw[E_]
        # right
        Rw = [dict() for _ in range(NCH)]
        Rw[E_] = idR
        YR = {}
        for a in range(NCH - 2, -1, -1):
            acc = {}
            for (aa, b, nzl) in M.ent:
                if aa != a or a == b:
                    continue
                dadd(acc, rtransfer(M, Rw[b], nzl, AR))
            if a == S_:
                YR = acc
            else:
                Rw[a] = acc
        x0 = getattr(self, '_RS', None)
        Rw[S_], eR, resR, infoR = self._solve_env(YR, self.rho_R(), idR, lambda x: rtransfer(M, x, idnz, AR), x0, tol)
        self._RS = Rw[S_]
        self.Lw, self.Rw = Lw, Rw
        q = M.q
        # (Y|rho)/(1|rho) is the energy per site averaged over the cell
        self.e_site = eL
        return eL, eR, resL, resR

    # ------------------------------------------------------------ effective Hamiltonians
    @staticmethod
    def _idx(env):
        out = []
        for d in env:
            ix = {}
            for (bra, ket), m in d.items():
                ix.setdefault(ket, []).append((bra, m))
            out.append(ix)
        return out

    def hac_matvec(self, V):
        M = self.M
        Lidx, Ridx = self._idx(self.Lw), self._idx(self.Rw)
        frm, bk = M.frm, M.bk
        kpos, offs = V.kpos, V.offs
        dt = M.dtype

        def mv(v):
            blocks = V.unpack(v)
            X = {}
            for (l, s, r), T in blocks.items():
                for a in range(NCH):
                    fr = frm[a][s]
                    if not fr:
                        continue
                    for (lb, Lm) in Lidx[a].get(l, ()):
                        LT = None
                        for (b, so, val) in fr:
                            if bk and lb[1] and so:
                                continue
                            if LT is None:
                                LT = Lm @ T
                            key = (b, lb, so, r)
                            if key in X:
                                X[key] += val * LT
                            else:
                                X[key] = val * LT
            out = np.zeros(V.n, dtype=dt)
            for (b, lb, so, r), Xm in X.items():
                for (rb, Rm) in Ridx[b].get(r, ()):
                    if bk and so and rb[1]:
                        continue
                    n = kpos.get((lb, so, rb))
                    if n is None:
                        continue
                    out[offs[n]:offs[n + 1]] += (Xm @ Rm.T).reshape(-1)
            return out
        return mv

    def hc_matvec(self, V):
        M = self.M
        Lidx, Ridx = self._idx(self.Lw), self._idx(self.Rw)
        bk = M.bk
        kpos, offs = V.kpos, V.offs
        dt = M.dtype

        def mv(v):
            blocks = V.unpack(v)
            out = np.zeros(V.n, dtype=dt)
            for (l, r), Cm in blocks.items():
                for a in range(NCH):
                    Ll = Lidx[a].get(l); Rl = Ridx[a].get(r)
                    if not Ll or not Rl:
                        continue
                    for (lb, Lm) in Ll:
                        LC = Lm @ Cm
                        for (rb, Rm) in Rl:
                            if bk and lb[1] and rb[1]:
                                continue
                            n = kpos.get((lb, rb))
                            if n is None:
                                continue
                            out[offs[n]:offs[n + 1]] += (LC @ Rm.T).reshape(-1)
            return out
        return mv

    def _eig(self, mv, n, v0, tol):
        dt = self.M.dtype
        if n <= 200:
            H = np.array([mv(e) for e in np.eye(n, dtype=dt)]).T
            H = 0.5 * (H + H.conj().T)
            e, U = np.linalg.eigh(H)
            return e[0], U[:, 0]
        op = spla.LinearOperator((n, n), matvec=mv, dtype=dt)
        if np.linalg.norm(v0) < 1e-14:
            v0 = self.rng.standard_normal(n).astype(dt)
        e, U = spla.eigsh(op, k=1, which='SA', v0=v0, tol=tol, ncv=min(n - 1, 24), maxiter=2000)
        return e[0], U[:, 0]

    def solve_local(self, tol):
        """per-residue ground states of H_AC and H_C"""
        M = self.M
        acs = self.ac_keys(); cks = self.c_keys()
        newAC, newC = {}, {}
        eAC = np.zeros(M.q); eC = np.zeros(M.q)
        for r in range(M.q):
            sh = {k: v for k, v in acs.items() if k[0][0] % M.q == r}
            if sh:
                V = Vec(sh)
                e, v = self._eig(self.hac_matvec(V), V.n, V.pack(self.AC, M.dtype), tol)
                eAC[r] = e
                newAC.update(V.unpack(v))
            sh = {k: v for k, v in cks.items() if k[0][0] % M.q == r}
            if sh:
                V = Vec(sh)
                e, v = self._eig(self.hc_matvec(V), V.n, V.pack(self.C, M.dtype), tol)
                eC[r] = e
                newC.update(V.unpack(v))
        # fix the sign gauge of C and AC to the previous iterate (keeps the polar factors smooth)
        for r in range(M.q):
            ov = sum(np.sum(np.conj(self.C[k]) * newC[k]) for k in newC if k in self.C and k[0][0] % M.q == r)
            if abs(ov) > 1e-14:
                ph = np.conj(ov) / abs(ov)
                for k in newC:
                    if k[0][0] % M.q == r:
                        newC[k] = newC[k] * ph
            ov = sum(np.sum(np.conj(self.AC[k]) * newAC[k]) for k in newAC if k in self.AC and k[0][0] % M.q == r)
            if abs(ov) > 1e-14:
                ph = np.conj(ov) / abs(ov)
                for k in newAC:
                    if k[0][0] % M.q == r:
                        newAC[k] = newAC[k] * ph
        self.AC, self.C = newAC, newC
        self.normalize()
        return eAC, eC

    def step(self, tol_eig=1e-8, tol_env=1e-10):
        eL, eR, rL, rR = self.environments(tol_env)
        eAC, eC = self.solve_local(tol_eig)
        self.polar()
        errL, errR = self.errors()
        return dict(e=eL, eR=eR, envres=max(rL, rR), err=max(errL, errR))

    # ------------------------------------------------------------ two-site subspace expansion
    def expand(self, chi_target):
        """Raise every bond (residue) toward chi_target total states using the two-site residual
        N_L^dag H2 (A_L C A_R) N_R^dag (Zauner-Stauber et al. 2018, appendix)."""
        M = self.M
        q, p, bk, dt = M.q, M.p, M.bk, M.dtype
        if not hasattr(self, 'Lw'):
            self.environments()
        # current two-site state theta[(l, s1, s2, r)] = AL[l,s1] C[t,h] AR[s2,r], t = lf(l,s1), h = rf(s2,r)
        ARby = {}
        for (s2, r), m in self.AR.items():
            ARby.setdefault(M.rf(s2, r), []).append((s2, r, m))
        theta = {}
        for (l, s1), A in self.AL.items():
            t = M.lf(l, s1)
            for (tt, h), Cm in self.C.items():
                if tt != t:
                    continue
                AC = A @ Cm
                for (s2, r, Bm) in ARby.get(h, ()):
                    if bk and s1 and s2:
                        continue
                    k = (l, s1, s2, r)
                    Y = AC @ Bm
                    theta[k] = theta[k] + Y if k in theta else Y
        # all valid two-site keys over the current label sets (the target space)
        keys = {}
        for l, dl in self.dimL.items():
            for s1 in (0, 1):
                if bk and l[1] and s1:
                    continue
                for s2 in (0, 1):
                    if bk and s1 and s2:
                        continue
                    b = 2 * p - l[0] - q * (s1 + s2)
                    for g in M.fl():
                        if bk and s2 and g:
                            continue
                        r = (b, g)
                        if r in self.dimR:
                            keys[(l, s1, s2, r)] = (dl, self.dimR[r])
        V = Vec(keys)
        Lidx, Ridx = self._idx(self.Lw), self._idx(self.Rw)
        frm = M.frm
        # H2 theta
        X = {}
        for (l, s1, s2, r), T in theta.items():
            for a in range(NCH):
                fr = frm[a][s1]
                if not fr:
                    continue
                for (lb, Lm) in Lidx[a].get(l, ()):
                    LT = Lm @ T
                    for (b, s1p, val) in fr:
                        if bk and lb[1] and s1p:
                            continue
                        key = (b, lb, s1p, s2, r)
                        X[key] = X[key] + val * LT if key in X else val * LT
        H2 = {}
        for (b, lb, s1p, s2, r), Xm in X.items():
            for (c, s2p, val) in frm[b][s2]:
                if bk and s1p and s2p:
                    continue
                for (rb, Rm) in Ridx[c].get(r, ()):
                    if bk and s2p and rb[1]:
                        continue
                    k = (lb, s1p, s2p, rb)
                    if k not in keys:
                        continue
                    Y = val * (Xm @ Rm.T)
                    H2[k] = H2[k] + Y if k in H2 else Y
        # project out the current left span (A_L columns) and right span (A_R rows)
        # left projection per left-fused label t: rows (l, s1)
        rowsby = {}
        for (l, s1, s2, r) in H2:
            rowsby.setdefault(M.lf(l, s1), set()).add((l, s1))
        for t, rs in rowsby.items():
            rs = sorted(rs)
            if t in self.dimL and all((l, s1) in self.AL for (l, s1) in rs):
                pass
            # build projector action: H2 <- H2 - A (A^dag H2) for blocks with fused label t
            cols = {}
            for (l, s1, s2, r) in H2:
                if M.lf(l, s1) == t:
                    cols.setdefault((s2, r), []).append((l, s1))
            for (s2, r), rl in cols.items():
                acc = None
                for (l, s1) in rl:
                    A = self.AL.get((l, s1))
                    if A is None:
                        continue
                    Y = A.conj().T @ H2[(l, s1, s2, r)]
                    acc = Y if acc is None else acc + Y
                if acc is None:
                    continue
                for (l, s1) in rl:
                    A = self.AL.get((l, s1))
                    if A is not None:
                        H2[(l, s1, s2, r)] = H2[(l, s1, s2, r)] - A @ acc
        colsby = {}
        for (l, s1, s2, r) in H2:
            colsby.setdefault(M.rf(s2, r), {}).setdefault((l, s1), []).append((s2, r))
        for h, rows in colsby.items():
            for (l, s1), cl in rows.items():
                acc = None
                for (s2, r) in cl:
                    B = self.AR.get((s2, r))
                    if B is None:
                        continue
                    Y = H2[(l, s1, s2, r)] @ B.conj().T
                    acc = Y if acc is None else acc + Y
                if acc is None:
                    continue
                for (s2, r) in cl:
                    B = self.AR.get((s2, r))
                    if B is not None:
                        H2[(l, s1, s2, r)] = H2[(l, s1, s2, r)] - acc @ B
        # candidate new left vectors per label t and right vectors per label h
        def tot(dim, rr):
            return sum(d for lab, d in dim.items() if lab[0] % q == rr)
        newL, newR = [], []
        # left: matrix rows (l,s1) of group t, cols all (s2, r)
        for t in rowsby:
            rl = sorted({(l, s1) for (l, s1, s2, r) in H2 if M.lf(l, s1) == t})
            cl = sorted({(s2, r) for (l, s1, s2, r) in H2 if M.lf(l, s1) == t})
            ro = {}; o = 0
            for (l, s1) in rl:
                ro[(l, s1)] = o; o += self.dimL[l]
            co = {}; o2 = 0
            for (s2, r) in cl:
                co[(s2, r)] = o2; o2 += self.dimR[r]
            Mx = np.zeros((o, o2), dtype=dt)
            for (l, s1, s2, r), Y in H2.items():
                if M.lf(l, s1) == t:
                    Mx[ro[(l, s1)]:ro[(l, s1)] + Y.shape[0], co[(s2, r)]:co[(s2, r)] + Y.shape[1]] = Y
            avail = o - self.dimL.get(t, 0)
            if avail <= 0:
                continue
            U, sv, _ = np.linalg.svd(Mx, full_matrices=False)
            for i in range(min(avail, len(sv))):
                newL.append((sv[i], t, i, U[:, i], rl, ro))
        for h in colsby:
            rl = sorted(colsby[h].keys())
            cl = sorted({x for v in colsby[h].values() for x in v})
            ro = {}; o = 0
            for (l, s1) in rl:
                ro[(l, s1)] = o; o += self.dimL[l]
            co = {}; o2 = 0
            for (s2, r) in cl:
                co[(s2, r)] = o2; o2 += self.dimR[r]
            Mx = np.zeros((o, o2), dtype=dt)
            for (l, s1), cs in colsby[h].items():
                for (s2, r) in cs:
                    Y = H2[(l, s1, s2, r)]
                    Mx[ro[(l, s1)]:ro[(l, s1)] + Y.shape[0], co[(s2, r)]:co[(s2, r)] + Y.shape[1]] = Y
            avail = o2 - self.dimR.get(h, 0)
            if avail <= 0:
                continue
            _, sv, Vh = np.linalg.svd(Mx, full_matrices=False)
            for i in range(min(avail, len(sv))):
                newR.append((sv[i], h, i, Vh[i], cl, co))
        # choose per residue of the bond
        addL = {}; addR = {}
        for rr in range(q):
            nl = chi_target - tot(self.dimL, rr)
            cand = sorted([x for x in newL if x[1][0] % q == rr], key=lambda x: -x[0])
            for x in cand[:max(0, nl)]:
                if x[0] > 1e-12:
                    addL.setdefault(x[1], []).append(x)
            nr = chi_target - tot(self.dimR, rr)
            cand = sorted([x for x in newR if x[1][0] % q == rr], key=lambda x: -x[0])
            for x in cand[:max(0, nr)]:
                if x[0] > 1e-12:
                    addR.setdefault(x[1], []).append(x)
        # apply: new columns of A_L blocks (l, s1) -> t; new rows of A_R blocks (s2, r) <- h
        for t, lst in addL.items():
            k = len(lst)
            old = self.dimL.get(t, 0)
            rl, ro = lst[0][4], lst[0][5]
            for (l, s1) in rl:
                A = self.AL.get((l, s1))
                if A is None:
                    A = np.zeros((self.dimL[l], old), dtype=dt)
                add = np.stack([x[3][ro[(l, s1)]:ro[(l, s1)] + self.dimL[l]] for x in lst], axis=1)
                self.AL[(l, s1)] = np.hstack([A, add.astype(dt)])
            self.dimL[t] = old + k
        for h, lst in addR.items():
            k = len(lst)
            old = self.dimR.get(h, 0)
            cl, co = lst[0][4], lst[0][5]
            for (s2, r) in cl:
                B = self.AR.get((s2, r))
                if B is None:
                    B = np.zeros((old, self.dimR[r]), dtype=dt)
                add = np.stack([x[3][co[(s2, r)]:co[(s2, r)] + self.dimR[r]] for x in lst], axis=0)
                self.AR[(s2, r)] = np.vstack([B, add.astype(dt)])
            self.dimR[h] = old + k
        self._pad_all()
        return sum(len(v) for v in addL.values()), sum(len(v) for v in addR.values())

    def _pad_all(self):
        M, dt = self.M, self.M.dtype

        def pad(m, r, c):
            if m.shape == (r, c):
                return m
            out = np.zeros((r, c), dtype=dt)
            out[:m.shape[0], :m.shape[1]] = m
            return out
        for (l, s), m in list(self.AL.items()):
            t = M.lf(l, s)
            self.AL[(l, s)] = pad(m, self.dimL[l], self.dimL[t])
        for (s, r), m in list(self.AR.items()):
            h = M.rf(s, r)
            self.AR[(s, r)] = pad(m, self.dimR[h], self.dimR[r])
        newC = {}
        for k, sh in self.c_keys().items():
            newC[k] = pad(self.C[k], *sh) if k in self.C else np.zeros(sh, dtype=dt)
        self.C = newC
        newAC = {}
        for k, sh in self.ac_keys().items():
            newAC[k] = pad(self.AC[k], *sh) if k in self.AC else np.zeros(sh, dtype=dt)
        self.AC = newAC
        for attr in ('_LE', '_RS'):
            if hasattr(self, attr):
                delattr(self, attr)

    def chi(self):
        q = self.M.q
        return [sum(d for lab, d in self.dimL.items() if lab[0] % q == r) for r in range(q)], \
               [sum(d for lab, d in self.dimR.items() if lab[0] % q == r) for r in range(q)]

    # ------------------------------------------------------------ measurements
    def entropies(self):
        """Schmidt spectrum and von Neumann entropy at every bond (residue r of the left label)"""
        q = self.M.q
        out = []
        specs = []
        for rr in range(q):
            bych = {}
            for (l, r), m in self.C.items():
                if l[0] % q == rr:
                    bych.setdefault(l[0], []).append((l, r, m))
            sv = []
            for a, lst in bych.items():
                ls = sorted({x[0] for x in lst}); rs = sorted({x[1] for x in lst})
                ro = {}; o = 0
                for l in ls:
                    ro[l] = o; o += self.dimL[l]
                co = {}; o2 = 0
                for r in rs:
                    co[r] = o2; o2 += self.dimR[r]
                Mx = np.zeros((o, o2), dtype=self.M.dtype)
                for (l, r, m) in lst:
                    Mx[ro[l]:ro[l] + m.shape[0], co[r]:co[r] + m.shape[1]] = m
                sv.append(np.linalg.svd(Mx, compute_uv=False))
            s = np.concatenate(sv)
            pr = s ** 2; pr = pr[pr > 1e-30]; pr = pr / pr.sum()
            out.append(float(-np.sum(pr * np.log(pr))))
            specs.append(np.sort(s)[::-1])
        return out, specs

    def density(self):
        """<n_x> at every cell position, from A_L and the right fixed point"""
        M = self.M
        rho = self.rho_L()
        idL = {(l, l): np.eye(d, dtype=M.dtype) for l, d in self.dimL.items()}
        T = ltransfer(M, idL, _nz(NN), self.AL)    # env at bond x+1 after n_x
        q = M.q
        nres = np.zeros(q); norm = np.zeros(q)
        for k, m in T.items():
            if k in rho:
                nres[k[0][0] % q] += np.real(np.sum(m * rho[k]))
        for k, m in rho.items():
            if k[0] == k[1]:
                norm[k[0][0] % q] += np.real(np.trace(m))
        return nres / norm, norm

    def transfer_spectrum(self, dQ, string=False, k=6, r0=0, tol=1e-9):
        """eigenvalues of the q-site transfer matrix (product over one cell) on the (bra, ket) blocks with
        ket charge = bra charge + q dQ, starting at residue r0, with or without the fermion string."""
        M = self.M
        q = M.q
        op = _nz(PAR) if string else _nz(ID)
        shapes = {}
        for l1, d1 in self.dimL.items():
            if l1[0] % q != r0:
                continue
            for l2, d2 in self.dimL.items():
                if l2[0] == l1[0] + q * dQ:
                    shapes[(l1, l2)] = (d1, d2)
        V = Vec(shapes)
        dt = complex

        def mv(v):
            x = V.unpack(v)
            for _ in range(q):
                x = ltransfer(M, x, op, self.AL)
            return V.pack(x, dt)
        n = V.n
        if n == 0:
            return np.array([])
        if n <= 400:
            T = np.array([mv(e) for e in np.eye(n, dtype=dt)]).T
            ev = np.linalg.eigvals(T)
        else:
            opx = spla.LinearOperator((n, n), matvec=mv, dtype=dt)
            ev = spla.eigs(opx, k=min(k, n - 2), which='LM', tol=tol, return_eigenvectors=False, maxiter=5000)
        ev = ev[np.argsort(-np.abs(ev))]
        return ev[:k]

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(dict(dimL=self.dimL, dimR=self.dimR, AL=self.AL, AR=self.AR, C=self.C, AC=self.AC), f)

    def load(self, path):
        with open(path, 'rb') as f:
            d = pickle.load(f)
        for k, v in d.items():
            setattr(self, k, v)
        for attr in ('_LE', '_RS', 'Lw', 'Rw'):
            if hasattr(self, attr):
                delattr(self, attr)
