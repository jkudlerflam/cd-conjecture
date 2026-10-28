#!/usr/bin/env python3
"""Exact diagonalization of the spinful mobile selector on a ring, in squeezed coordinates.

Microscopic chain (d = 6, one-site translation invariant, nearest neighbour).  Letters
    L, R, #up, #dn, Xup, Xdn       with arrows  L=(0,0), R=(1,1), #=(1,0), X=(0,1)
(the C8 section 15 convention).  The bond is
    h = K Q + (mu/2)(n_# x 1 + 1 x n_#)
        - t  sum_s (|L Xs><Xs R| + h.c.)          head hop   (C7)
        - u  sum_s (|#s L><R #s| + h.c.)          marker hop (C8 section 13)
        + (J/2)(S+ S- + S- S+) + Jz Sz Sz         on the adjacent wall pairs (#,X) and (X,#)
Q projects on pairs whose touching arrows disagree.  Every Q_i commutes with H (C7, C8), so
the legal sector (all Q_i = 0) is an exact invariant subspace; the code works inside it.

Legal sector = binary arrows b_i on the bonds, i.e. W domain walls (W even) on sites, walls
alternating between X (0->1, a "head") and # (1->0, a "marker"), each wall carrying a spin.
Basis: (wall positions, orientation o = species of the lowest-position wall, 0 = X, 1 = #)
tensor (spins of the W walls listed in order of position) = the squeezed spin chain.
A wall hopping across the bond (L-1, 0) cyclically relabels the squeezed chain; that is the
only place where the two factors talk apart from the exchange, which acts on squeezed
neighbours k, k+1 only when they are physically adjacent.

numpy + scipy only.  Written by Claude (Opus 5.5) for contrib/mobile_cells, 2026-10-27.
"""
import itertools
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, eigsh
from math import comb


def popcount(a):
    a = np.asarray(a, dtype=np.int64)
    c = np.zeros_like(a)
    x = a.copy()
    while np.any(x):
        c += x & 1
        x >>= 1
    return c


class Ring:
    def __init__(self, L, W, nup, t=1.0, u=1.0, mu=0.0, J=1.0, Jz=0.0, dJ=0.0, cells=None):
        """cells: optional tuple of cell lengths (l_1, ..., l_M), sum L, each >= 2.  Restricts to
        configurations with markers at 0, l_1, l_1 + l_2, ... and one head inside each cell, an
        exact invariant subspace when u = 0 (the static selector)."""
        assert W % 2 == 0 and 2 <= W <= L
        self.L, self.W, self.nup = L, W, nup
        self.t, self.u, self.mu, self.J, self.Jz, self.dJ = t, u, mu, J, Jz, dJ
        if cells is None:
            # ---- charge basis: wall masks with W bits, times orientation
            pos = np.array(list(itertools.combinations(range(L), W)), dtype=np.int64)
            masks = (np.int64(1) << pos).sum(axis=1)
            cm = np.concatenate([masks, masks]); co = np.concatenate([np.zeros(len(masks), np.int64), np.ones(len(masks), np.int64)])
            cp = np.concatenate([pos, pos])
        else:
            assert sum(cells) == L and len(cells) == W // 2 and min(cells) >= 2
            x = np.concatenate([[0], np.cumsum(cells)[:-1]])
            heads = itertools.product(*[range(xj + 1, xj + lj) for xj, lj in zip(x, cells)])
            pos = np.array([sorted(list(x) + list(hh)) for hh in heads], dtype=np.int64)
            cm = (np.int64(1) << pos).sum(axis=1)
            co = np.ones(len(cm), np.int64)          # the wall at site 0 is a marker
            cp = pos
        self.cm, self.co, self.cp = cm, co, cp
        self.nc = len(cm)
        code = cm * 2 + co
        order = np.argsort(code)
        self._ccode = code[order]; self._corder = order
        # ---- spin basis: W-bit integers with nup ones
        sb = np.array([sum(1 << p for p in c) for c in itertools.combinations(range(W), nup)], dtype=np.int64)
        sb.sort()
        self.sb = sb
        self.ns = len(sb)
        self.dim = self.nc * self.ns
        self._build()

    # ------------------------------------------------------------------ indexing
    def cidx(self, mask, o):
        code = np.asarray(mask, np.int64) * 2 + np.asarray(o, np.int64)
        j = np.searchsorted(self._ccode, code)
        j = np.clip(j, 0, len(self._ccode) - 1)
        ok = self._ccode[j] == code
        return np.where(ok, self._corder[j], -1)

    def sidx(self, s):
        j = np.searchsorted(self.sb, s)
        j = np.clip(j, 0, self.ns - 1)
        return np.where(self.sb[j] == s, j, -1)

    # ------------------------------------------------------------------ Hamiltonian pieces
    def _build(self):
        L, W = self.L, self.W
        cm, co, cp = self.cm, self.co, self.cp
        species = (co[:, None] + np.arange(W)[None, :]) % 2        # 0 = X (hop t), 1 = # (hop u)
        amp = np.where(species == 0, -self.t, -self.u)
        rows, cols, vals = [], [], []
        # bulk hops (do not cross the bond L-1 | 0)
        for k in range(W):
            x = cp[:, k]
            for dx in (+1, -1):
                y = x + dx
                ok = (y >= 0) & (y <= L - 1)
                yy = np.where(ok, y, 0)
                ok &= ((cm >> yy) & 1) == 0
                newm = cm - (1 << x) + (1 << yy)
                j = self.cidx(newm, co)
                ok &= j >= 0
                i = np.nonzero(ok)[0]
                rows.append(j[i]); cols.append(i); vals.append(amp[i, k])
        self.Hb = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                    shape=(self.nc, self.nc))
        # crossing hops: last wall L-1 -> 0 (becomes first; squeezed chain rotates right)
        ok = (cp[:, W - 1] == L - 1) & ((cm & 1) == 0)
        i = np.nonzero(ok)[0]
        newm = cm[i] - (1 << (L - 1)) + 1
        newo = species[i, W - 1]
        j = self.cidx(newm, newo)
        assert np.all(j >= 0) or self._restricted()
        g = j >= 0
        self.HxR = sparse.csr_matrix((amp[i[g], W - 1], (j[g], i[g])), shape=(self.nc, self.nc))
        # first wall 0 -> L-1 (becomes last; squeezed chain rotates left)
        ok = (cp[:, 0] == 0) & (((cm >> (L - 1)) & 1) == 0)
        i = np.nonzero(ok)[0]
        newm = cm[i] - 1 + (1 << (L - 1))
        newo = species[i, 1]
        j = self.cidx(newm, newo)
        g = j >= 0
        self.HxL = sparse.csr_matrix((amp[i[g], 0], (j[g], i[g])), shape=(self.nc, self.nc))
        # spin permutations
        sb = self.sb
        full = (1 << W) - 1
        rotR = ((sb << 1) & full) | (sb >> (W - 1))     # new bit k = old bit k-1, new bit 0 = old bit W-1
        rotL = (sb >> 1) | ((sb & 1) << (W - 1))        # new bit k = old bit k+1
        self.rotRidx = self.sidx(rotR); self.rotLidx = self.sidx(rotL)
        # exchange on squeezed bonds k, k+1 (k = W-1 is the wrap bond W-1, 0)
        self.ex = []
        for k in range(W):
            k2 = (k + 1) % W
            if k < W - 1:
                adj = cp[:, k + 1] == cp[:, k] + 1
            else:
                adj = (cp[:, W - 1] == L - 1) & (cp[:, 0] == 0)
            r = np.nonzero(adj)[0]
            if len(r) == 0:
                continue
            bk = (sb >> k) & 1; bk2 = (sb >> k2) & 1
            anti = np.nonzero(bk != bk2)[0]
            tgt = self.sidx(sb[anti] ^ ((1 << k) | (1 << k2)))
            zz = (bk - 0.5) * (bk2 - 0.5)
            # (#, X) bonds carry (1 + dJ), (X, #) bonds (1 - dJ); dJ != 0 breaks the reflection
            f = 1.0 + self.dJ * np.where((co[r] + k) % 2 == 1, 1.0, -1.0)
            self.ex.append((r, anti, tgt, zz, f))
        self.const = self.mu * W / 2

    def _restricted(self):
        return self.nc < 2 * comb(self.L, self.W)

    def matvec(self, v):
        Psi = v.reshape(self.nc, self.ns)
        out = self.Hb @ Psi
        if self.HxR.nnz:
            out[:, self.rotRidx] += self.HxR @ Psi
        if self.HxL.nnz:
            out[:, self.rotLidx] += self.HxL @ Psi
        hJ = 0.5 * self.J
        for r, anti, tgt, zz, f in self.ex:
            P = Psi[r]
            blk = np.zeros_like(P)
            if hJ != 0:
                blk[:, tgt] += hJ * P[:, anti]
            if self.Jz != 0:
                blk += self.Jz * zz[None, :] * P
            out[r] += f[:, None] * blk if self.dJ else blk
        out += self.const * Psi
        return out.ravel()

    def op(self):
        return LinearOperator((self.dim, self.dim), matvec=self.matvec, dtype=float)

    def sparse(self):
        """The same Hamiltonian as an explicit sparse matrix on (charge x spin), row-major."""
        ns = self.ns
        I = sparse.identity(ns, format='csr')
        def perm(idx):        # |s> -> |idx[s]>
            return sparse.csr_matrix((np.ones(ns), (idx, np.arange(ns))), shape=(ns, ns))
        H = sparse.kron(self.Hb, I, format='csr')
        if self.HxR.nnz:
            H = H + sparse.kron(self.HxR, perm(self.rotRidx), format='csr')
        if self.HxL.nnz:
            H = H + sparse.kron(self.HxL, perm(self.rotLidx), format='csr')
        for r, anti, tgt, zz, f in self.ex:
            D = sparse.csr_matrix((f, (r, r)), shape=(self.nc, self.nc))
            X = sparse.csr_matrix((0.5 * self.J * np.ones(len(anti)), (tgt, anti)), shape=(ns, ns))
            if self.Jz:
                X = X + sparse.diags(self.Jz * zz)
            H = H + sparse.kron(D, X, format='csr')
        H = H + self.const * sparse.identity(self.dim, format='csr')
        return H.tocsr()

    def lowest(self, k=10, tol=1e-10, v0=None, ncv=None, use_sparse=6e6):
        if self.dim <= 2500:
            H = self.sparse().toarray()
            w, V = np.linalg.eigh(0.5 * (H + H.T))
            return w[:k], V[:, :k]
        k = min(k, self.dim - 2)
        A = self.sparse() if self.dim <= use_sparse else self.op()
        w, V = eigsh(A, k=k, which='SA', tol=tol, v0=v0, ncv=ncv or max(2 * k + 1, 40))
        o = np.argsort(w)
        return w[o], V[:, o]

    # ------------------------------------------------------------------ symmetries
    def _tmaps(self):
        if not hasattr(self, '_tm'):
            L, W = self.L, self.W
            cm, co, cp = self.cm, self.co, self.cp
            wrap = cp[:, W - 1] == L - 1                       # T: last wall L-1 -> 0
            newm = np.where(wrap, ((cm << 1) & ((1 << L) - 1)) | 1, cm << 1)
            j = self.cidx(newm, np.where(wrap, 1 - co, co))
            wrapi = cp[:, 0] == 0                              # T^-1: first wall 0 -> L-1
            newmi = np.where(wrapi, (cm >> 1) | (1 << (L - 1)), cm >> 1)
            ji = self.cidx(newmi, np.where(wrapi, 1 - co, co))
            assert np.all(j >= 0) and np.all(ji >= 0)
            self._tm = (wrap, j, wrapi, ji)
        return self._tm

    def translate(self, v, inverse=False):
        """T: every wall x -> x+1 (mod L); inverse=True gives T^-1."""
        wrap, j, wrapi, ji = self._tmaps()
        if inverse:
            wrap, j, rot = wrapi, ji, self.rotLidx
        else:
            rot = self.rotRidx
        Psi = v.reshape(self.nc, self.ns)
        out = np.zeros_like(Psi)
        nw = ~wrap
        out[j[nw]] = Psi[nw]
        tmp = np.zeros((wrap.sum(), self.ns), dtype=Psi.dtype)
        tmp[:, rot] = Psi[wrap]
        out[j[wrap]] = tmp
        return out.ravel()

    def lowest_momentum(self, kint, k=3, lam=None, tol=1e-9):
        """Lowest k levels at physical momentum P = 2 pi kint / L: eigsh of
        H + lam (1 - cos(Phat - P)), which leaves momentum P unshifted and lifts the others by
        at least lam (1 - cos(2 pi / L))."""
        L = self.L
        P = 2 * np.pi * kint / L
        if lam is None:
            lam = 3.0 / (1 - np.cos(2 * np.pi / L))
        A = self.sparse() if self.dim <= 6e6 else None
        def mv(v):
            v = np.asarray(v).ravel()
            Hv = (A @ v.real + 1j * (A @ v.imag)) if A is not None else (self.matvec(v.real) + 1j * self.matvec(v.imag))
            Tv = self.translate(v); Tiv = self.translate(v, inverse=True)
            return Hv + lam * (v - 0.5 * (np.exp(-1j * P) * Tv + np.exp(1j * P) * Tiv))
        op = LinearOperator((self.dim, self.dim), matvec=mv, dtype=complex)
        w, V = eigsh(op, k=k, which='SA', tol=tol, ncv=max(4 * k, 30))
        o = np.argsort(w)
        w, V = w[o], V[:, o]
        # check the momentum of each returned vector
        Pm = [np.angle(np.vdot(V[:, a], self.translate(V[:, a]))) * L / (2 * np.pi) for a in range(k)]
        return w, np.array(Pm)

    def momenta(self, w, V, degtol=1e-7):
        """Return momenta P (in units of 2 pi / L) for the eigenvectors, diagonalizing T inside
        clusters of degenerate energies."""
        L = self.L
        P = np.zeros(len(w)); i = 0
        while i < len(w):
            j = i + 1
            while j < len(w) and abs(w[j] - w[i]) < degtol * max(1, abs(w[i])):
                j += 1
            B = V[:, i:j]
            TB = np.column_stack([self.translate(B[:, a]) for a in range(j - i)])
            M = B.T @ TB
            ev = np.linalg.eigvals(M)
            ang = np.sort(np.round(np.angle(ev) / (2 * np.pi) * L, 3) % L)
            P[i:j] = ang
            i = j
        return P

    # ------------------------------------------------------------------ entanglement
    def half_entropy(self, v, ell):
        """Von Neumann entropy of physical sites 0..ell-1 in the d = 6 letter basis."""
        L, W = self.L, self.W
        Psi = v.reshape(self.nc, self.ns)
        cp, co = self.cp, self.co
        nA = (cp < ell).sum(axis=1)                         # walls in A are the first nA walls
        maskA = self.cm & ((1 << ell) - 1)
        maskB = self.cm >> ell
        bin_ = co                                           # arrow entering site 0 (b_{L-1}) = o
        sb = self.sb
        # left key: (b_in, maskA, spinsA); right key: (maskB, spinsB) (its entering arrow is fixed by left)
        nAm = nA[:, None]
        spinsA = sb[None, :] & ((1 << nAm) - 1)
        spinsB = sb[None, :] >> nAm
        lk = ((bin_[:, None] * (1 << ell) + maskA[:, None]) << W) + spinsA
        rk = (((bin_ ^ (nA & 1))[:, None] * (1 << (L - ell)) + maskB[:, None]) << W) + spinsB
        lk = lk.ravel(); rk = rk.ravel(); amp = Psi.ravel()
        nz = np.abs(amp) > 1e-15
        lk, rk, amp = lk[nz], rk[nz], amp[nz]
        lu, li = np.unique(lk, return_inverse=True)
        ru, ri = np.unique(rk, return_inverse=True)
        M = sparse.csr_matrix((amp, (li, ri)), shape=(len(lu), len(ru)))
        # block structure: connected components of the bipartite support graph
        from scipy.sparse.csgraph import connected_components
        nL, nR = len(lu), len(ru)
        G = sparse.bmat([[None, (M != 0).astype(float)], [(M != 0).T.astype(float), None]]).tocsr()
        ncomp, lab = connected_components(G, directed=False)
        labL, labR = lab[:nL], lab[nL:]
        p_all = []
        Mc = M.tocsc()
        orderL = np.argsort(labL); orderR = np.argsort(labR)
        sL = np.searchsorted(labL[orderL], np.arange(ncomp + 1))
        sR = np.searchsorted(labR[orderR], np.arange(ncomp + 1))
        Mr = M.tocsr()
        for c in range(ncomp):
            rl = orderL[sL[c]:sL[c + 1]]; rr = orderR[sR[c]:sR[c + 1]]
            if len(rl) == 0 or len(rr) == 0:
                continue
            B = Mr[rl][:, rr].toarray()
            s = np.linalg.svd(B, compute_uv=False)
            p_all.append(s ** 2)
        p = np.concatenate(p_all); p = p / p.sum(); p = p[p > 1e-16]
        return float(-(p * np.log(p)).sum())


def dense_bond(t=1.0, u=1.0, mu=0.0, J=1.0, Jz=0.0, K=10.0, dJ=0.0):
    """The 36 x 36 microscopic bond h (left site = left kron factor), letters
    0 L, 1 R, 2 #up, 3 #dn, 4 Xup, 5 Xdn."""
    arrows = {0: (0, 0), 1: (1, 1), 2: (1, 0), 3: (1, 0), 4: (0, 1), 5: (0, 1)}
    d = 6
    h = np.zeros((d * d, d * d))
    ix = lambda a, b: a * d + b
    for a in range(d):
        for b in range(d):
            if arrows[a][1] != arrows[b][0]:
                h[ix(a, b), ix(a, b)] += K
            h[ix(a, b), ix(a, b)] += 0.5 * mu * ((a in (2, 3)) + (b in (2, 3)))
    for s in (0, 1):
        X = 4 + s; M = 2 + s
        # head hop |L X><X R|
        h[ix(0, X), ix(X, 1)] += -t; h[ix(X, 1), ix(0, X)] += -t
        # marker hop |# L><R #|
        h[ix(M, 0), ix(1, M)] += -u; h[ix(1, M), ix(M, 0)] += -u
    # exchange on (#,X) and (X,#)
    for (A, B), fac in (((2, 4), 1 + dJ), ((4, 2), 1 - dJ)):
        up_a, dn_a, up_b, dn_b = A, A + 1, B, B + 1
        h[ix(up_a, dn_b), ix(dn_a, up_b)] += fac * J / 2; h[ix(dn_a, up_b), ix(up_a, dn_b)] += fac * J / 2
        for sa, a_ in ((+.5, up_a), (-.5, dn_a)):
            for sb_, b_ in ((+.5, up_b), (-.5, dn_b)):
                h[ix(a_, b_), ix(a_, b_)] += fac * Jz * sa * sb_
    return h


def brute_force_check(L=5, W=2, nup=1, seed=1, **par):
    """Compare the squeezed-coordinate spectrum with full d^L diagonalization restricted to the
    same (legal, W, S^z) sector."""
    rng = np.random.default_rng(seed)
    p = dict(t=rng.uniform(.5, 1.5), u=rng.uniform(.2, 1.5), mu=rng.uniform(-1, 1), J=rng.uniform(.2, 1.5), Jz=rng.uniform(-1, 1),
             dJ=rng.uniform(-.5, .5))
    p.update(par)
    r = Ring(L, W, nup, **p)
    w1, _ = r.lowest(k=r.dim)
    h = dense_bond(K=50.0, **p)
    d = 6
    # direct construction by enumerating basis states (small L only)
    states = np.array(list(itertools.product(range(d), repeat=L)))
    arrows = np.array([(0, 0), (1, 1), (1, 0), (1, 0), (0, 1), (0, 1)])
    legal = np.all(arrows[states][:, :, 1] == arrows[np.roll(states, -1, axis=1)][:, :, 0], axis=1)
    walls = np.isin(states, [2, 3, 4, 5]).sum(axis=1)
    ups = np.isin(states, [2, 4]).sum(axis=1)
    sel = legal & (walls == W) & (ups == nup)
    S = states[sel]
    key = {tuple(s): n for n, s in enumerate(S)}
    Hs = np.zeros((len(S), len(S)))
    for n, s in enumerate(S):
        for i in range(L):
            j = (i + 1) % L
            col = s[i] * d + s[j]
            for row in np.nonzero(h[:, col])[0]:
                a, b = divmod(row, d)
                s2 = s.copy(); s2[i] = a; s2[j] = b
                m = key.get(tuple(s2))
                if m is not None:
                    Hs[m, n] += h[row, col]
    w2 = np.linalg.eigvalsh(Hs)
    return len(S), r.dim, float(np.max(np.abs(np.sort(w1) - np.sort(w2)[:len(w1)])))


if __name__ == '__main__':
    for (L, W, nup) in [(4, 2, 1), (5, 2, 1), (5, 4, 2), (6, 4, 2), (6, 2, 0), (6, 4, 1)]:
        print(L, W, nup, brute_force_check(L, W, nup))
