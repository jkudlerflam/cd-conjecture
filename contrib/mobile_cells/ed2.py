#!/usr/bin/env python3
"""Two-flavour payload (c_p = 2): every wall of the mobile selector carries two qubits, a and b.

Microscopic chain, d = 10: letters L, R, #(sa,sb), X(sa,sb).  Bond as in ed.py, with
    + sum_f (J_f/2)(S+_f S-_f + h.c.)  on adjacent (#,X), (X,#) pairs, f = a, b,
    + on every wall  V S^z_a S^z_b - h_a S^z_a - h_b S^z_b   (half on each bond).
Flavours are exchanged only with the same flavour, so S^z_a and S^z_b are separately conserved;
h_b != 0 gives flavour b a different filling from a (incommensurate flavours), and V couples
the two flavours.  Basis: (charge configuration) x (a bits) x (b bits), as in ed.Ring.
Claude (Opus 5.5) for contrib/mobile_cells, 2026-10-27.
"""
import itertools
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from ed import Ring


class Ring2(Ring):
    def __init__(self, L, W, na, nb, t=1.0, u=1.0, mu=0.0, Ja=0.5, Jb=0.35, V=0.3, ha=0.0, hb=0.05):
        self.na, self.nb, self.Ja, self.Jb, self.V, self.ha, self.hb = na, nb, Ja, Jb, V, ha, hb
        super().__init__(L, W, na, t=t, u=u, mu=mu, J=0.0)

    def _build(self):
        super()._build()                     # charge pieces; exchange rebuilt below
        W = self.W
        A = np.array([sum(1 << p for p in c) for c in itertools.combinations(range(W), self.na)], np.int64)
        B = np.array([sum(1 << p for p in c) for c in itertools.combinations(range(W), self.nb)], np.int64)
        sb = (A[:, None] + (B[None, :] << W)).ravel()
        sb.sort()
        self.sb = sb; self.ns = len(sb); self.dim = self.nc * self.ns
        full = (1 << W) - 1
        def rotR(x):
            a, b = x & full, x >> W
            ra = ((a << 1) & full) | (a >> (W - 1)); rb = ((b << 1) & full) | (b >> (W - 1))
            return ra | (rb << W)
        def rotL(x):
            a, b = x & full, x >> W
            ra = (a >> 1) | ((a & 1) << (W - 1)); rb = (b >> 1) | ((b & 1) << (W - 1))
            return ra | (rb << W)
        self.rotRidx = self.sidx(rotR(sb)); self.rotLidx = self.sidx(rotL(sb))
        cp, L = self.cp, self.L
        self.ex = []
        for k in range(W):
            k2 = (k + 1) % W
            adj = (cp[:, k + 1] == cp[:, k] + 1) if k < W - 1 else ((cp[:, W - 1] == L - 1) & (cp[:, 0] == 0))
            r = np.nonzero(adj)[0]
            if len(r) == 0:
                continue
            for f, Jf in ((0, self.Ja), (1, self.Jb)):
                o = f * W
                b1 = (sb >> (k + o)) & 1; b2 = (sb >> (k2 + o)) & 1
                anti = np.nonzero(b1 != b2)[0]
                tgt = self.sidx(sb[anti] ^ ((1 << (k + o)) | (1 << (k2 + o))))
                self.ex.append((r, anti, tgt, Jf))
        # spin-diagonal on-wall terms (same for every charge configuration)
        a_bits = ((sb[:, None] >> np.arange(W)[None, :]) & 1) - 0.5
        b_bits = ((sb[:, None] >> (W + np.arange(W))[None, :]) & 1) - 0.5
        self.diag_s = self.V * (a_bits * b_bits).sum(1) - self.ha * a_bits.sum(1) - self.hb * b_bits.sum(1)

    def matvec(self, v):
        Psi = v.reshape(self.nc, self.ns)
        out = self.Hb @ Psi
        if self.HxR.nnz:
            out[:, self.rotRidx] += self.HxR @ Psi
        if self.HxL.nnz:
            out[:, self.rotLidx] += self.HxL @ Psi
        for r, anti, tgt, Jf in self.ex:
            blk = np.zeros((len(r), self.ns), dtype=Psi.dtype)
            blk[:, tgt] = 0.5 * Jf * Psi[r][:, anti]
            out[r] += blk
        out += (self.const + self.diag_s[None, :]) * Psi
        return out.ravel()

    def sparse(self):
        ns = self.ns
        I = sparse.identity(ns, format='csr')
        perm = lambda idx: sparse.csr_matrix((np.ones(ns), (idx, np.arange(ns))), shape=(ns, ns))
        H = sparse.kron(self.Hb, I, format='csr')
        if self.HxR.nnz:
            H = H + sparse.kron(self.HxR, perm(self.rotRidx), format='csr')
        if self.HxL.nnz:
            H = H + sparse.kron(self.HxL, perm(self.rotLidx), format='csr')
        for r, anti, tgt, Jf in self.ex:
            D = sparse.csr_matrix((np.ones(len(r)), (r, r)), shape=(self.nc, self.nc))
            X = sparse.csr_matrix((0.5 * Jf * np.ones(len(anti)), (tgt, anti)), shape=(ns, ns))
            H = H + sparse.kron(D, X, format='csr')
        H = H + sparse.kron(sparse.identity(self.nc), sparse.diags(self.const + self.diag_s), format='csr')
        return H.tocsr()

    def half_entropy(self, v, ell):
        raise NotImplementedError


def dense_bond2(t=1.0, u=1.0, mu=0.0, Ja=0.5, Jb=0.35, V=0.3, ha=0.0, hb=0.05, K=10.0):
    """100 x 100 bond.  Letters: 0 L, 1 R, 2 + s #, 6 + s X, with s = 2 sa + sb (1 = up)."""
    d = 10
    ar = lambda x: (0, 0) if x == 0 else (1, 1) if x == 1 else (1, 0) if x < 6 else (0, 1)
    sp = lambda x: ((x - 2) % 4 >> 1, (x - 2) % 4 & 1)          # (sa, sb) of a wall letter
    wall = lambda x: x >= 2
    h = np.zeros((d * d, d * d))
    ix = lambda a, b: a * d + b
    for a in range(d):
        for b in range(d):
            if ar(a)[1] != ar(b)[0]:
                h[ix(a, b), ix(a, b)] += K
            for x in (a, b):
                if 2 <= x < 6:
                    h[ix(a, b), ix(a, b)] += 0.5 * mu
                if wall(x):
                    sa, sb = sp(x)
                    za, zb = sa - 0.5, sb - 0.5
                    h[ix(a, b), ix(a, b)] += 0.5 * (V * za * zb - ha * za - hb * zb)
    for s in range(4):
        X, M = 6 + s, 2 + s
        h[ix(0, X), ix(X, 1)] += -t; h[ix(X, 1), ix(0, X)] += -t
        h[ix(M, 0), ix(1, M)] += -u; h[ix(1, M), ix(M, 0)] += -u
    for base1, base2 in ((2, 6), (6, 2)):
        for s1 in range(4):
            for s2 in range(4):
                a1, b1 = s1 >> 1, s1 & 1
                a2, b2 = s2 >> 1, s2 & 1
                if a1 != a2:      # flip flavour a
                    h[ix(base1 + ((a2 << 1) | b1), base2 + ((a1 << 1) | b2)), ix(base1 + s1, base2 + s2)] += Ja / 2
                if b1 != b2:      # flip flavour b
                    h[ix(base1 + ((a1 << 1) | b2), base2 + ((a2 << 1) | b1)), ix(base1 + s1, base2 + s2)] += Jb / 2
    return h


def brute_force_check2(L=5, W=2, na=1, nb=1, seed=2):
    rng = np.random.default_rng(seed)
    p = dict(t=rng.uniform(.5, 1.5), u=rng.uniform(.2, 1.5), mu=rng.uniform(-1, 1), Ja=rng.uniform(.2, 1),
             Jb=rng.uniform(.2, 1), V=rng.uniform(-1, 1), ha=rng.uniform(-.3, .3), hb=rng.uniform(-.3, .3))
    r = Ring2(L, W, na, nb, **p)
    w1 = np.linalg.eigvalsh(r.sparse().toarray())
    h = dense_bond2(K=50.0, **p)
    d = 10
    states = np.array(list(itertools.product(range(d), repeat=L)))
    arr = np.array([(0, 0), (1, 1)] + [(1, 0)] * 4 + [(0, 1)] * 4)
    legal = np.all(arr[states][:, :, 1] == arr[np.roll(states, -1, axis=1)][:, :, 0], axis=1)
    walls = (states >= 2).sum(1)
    sa = np.where(states >= 2, ((states - 2) % 4) >> 1, 0).sum(1)
    sbb = np.where(states >= 2, ((states - 2) % 4) & 1, 0).sum(1)
    S = states[legal & (walls == W) & (sa == na) & (sbb == nb)]
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
    return len(S), r.dim, float(np.max(np.abs(w1 - w2)))


if __name__ == '__main__':
    for args in [(4, 2, 1, 1), (5, 2, 1, 0), (5, 4, 2, 1), (6, 4, 2, 3), (6, 2, 1, 2)]:
        print(args, brute_force_check2(*args))
