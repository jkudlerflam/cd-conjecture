#!/usr/bin/env python3
"""Exact ring spectra of the d = 11 vacancy chain (vacbond.py) in squeezed coordinates.

A basis state on a ring of L sites with M occupied sites is (positions x_0 < ... < x_{M-1}) x (labels of
the occupied sites in that order).  This is an exact change of basis (../counterexample_campaign/
doping_check/README.md, Section 1).  Because the vacancy swap carries a letter with its type, the cyclic
word w of site types (P/B) read along the occupied sites is conserved.  For a word of length M with
smallest period p, a state is (offset s in Z_p) x (letters), where squeezed site j carries type
w[(j + s) mod M].  Then exactly

  H = -t (Tb (x) 1 + B (x) R + B^T (x) R^T) + J sum_j A_j (x) h_{j,j+1} - mu_v M,

Tb = hops not crossing the cut, B = hop L-1 -> 0 (the last particle becomes the first), R = cyclic shift
of the labels with s -> s - 1, A_j = adjacency of squeezed sites j, j+1 (j = M-1 pairs with 0 across the
cut, A_{M-1} = [x_0 = 0 and x_{M-1} = L-1]), h_{j,j+1} = the demonstrator bond restricted to the two
types.  Nothing is truncated.  Sectors: the word, total logical charge 2Q, and physical momentum
k = 2 pi m / L, imposed by the real projector (1/L) sum_r cos(k r) T^r (it spans +-k, which are
degenerate because H is real).  L = M (no vacancies) is the frozen chain.

Storage: W[alpha, c] with alpha the internal sector basis (all offsets) and c the position configuration.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse, itertools
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '4')
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh, LinearOperator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vacbond as VB  # noqa: E402


def period(word):
    M = len(word)
    for p in range(1, M + 1):
        if M % p == 0 and word == word[p:] + word[:p]:
            return p
    return M


def canonical(word):
    """Lexicographically smallest rotation (cyclic words are conserved up to rotation)."""
    return min(word[i:] + word[:i] for i in range(len(word)))


class Internal:
    """Internal (label) space of a cyclic type word at fixed 2Q, all offsets s = 0..p-1."""

    def __init__(self, word, Q2, hb, letters=None, q2letter=None, allowed=None):
        """allowed: optional boolean matrix over the full letter alphabet; only internal states whose every
        cyclically consecutive pair is allowed are kept (exact when that set is conserved, e.g. the arrow
        grammar of the 23-letter bond)."""
        letters = VB.LETTERS if letters is None else letters
        q2letter = VB.Q2_LETTER if q2letter is None else q2letter
        self.dtype = np.result_type(*[h.dtype for h in hb.values()])
        self.word, self.M, self.p = word, len(word), period(word)
        M = self.M
        self.blocks = []
        off = 0
        for s in range(self.p):
            types = [word[(j + s) % M] for j in range(M)]
            dims = [len(letters[x]) for x in types]
            q2 = [q2letter[letters[x]] for x in types]
            tot = int(np.prod(dims))
            dig = np.array(np.unravel_index(np.arange(tot), dims))          # (M, tot)
            charge = np.zeros(tot, int)
            for j in range(M):
                charge += q2[j][dig[j]]
            ok = charge == Q2
            if allowed is not None:
                for j in range(M):
                    j1 = (j + 1) % M
                    ok &= allowed[letters[types[j]][dig[j]], letters[types[j1]][dig[j1]]]
            keep = np.where(ok)[0]
            pos = -np.ones(tot, np.int64)
            pos[keep] = off + np.arange(len(keep))
            strides = np.array([int(np.prod(dims[j + 1:])) for j in range(M)], np.int64)
            self.blocks.append(dict(types=types, dims=dims, keep=keep, dig=dig[:, keep], pos=pos,
                                    strides=strides, off=off, n=len(keep)))
            off += len(keep)
        self.D = off
        self.ops = [self._bond_op(j, hb) for j in range(M)]
        self.R = self._rot()

    def _bond_op(self, j, hb):
        M = self.M
        j1 = (j + 1) % M
        rows, cols, vals = [], [], []
        for b in self.blocks:
            if b['n'] == 0:
                continue
            ta, tb = b['types'][j], b['types'][j1]
            h = hb[ta + tb]
            da, db = b['dims'][j], b['dims'][j1]
            nz = np.argwhere(np.abs(h) > 1e-15)
            ia, ib = b['dig'][j], b['dig'][j1]
            full = b['keep']
            for o, i in nz:
                oa, ob = divmod(o, db)
                xa, xb = divmod(i, db)
                sel = np.where((ia == xa) & (ib == xb))[0]
                if len(sel) == 0:
                    continue
                newfull = full[sel] + (oa - xa) * b['strides'][j] + (ob - xb) * b['strides'][j1]
                r = b['pos'][newfull]
                assert np.all(r >= 0), 'charge leak'
                rows.append(r); cols.append(b['off'] + sel); vals.append(np.full(len(sel), h[o, i]))
        if not rows:
            return sps.csr_matrix((self.D, self.D), dtype=self.dtype)
        rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
        return sps.csr_matrix((vals.astype(self.dtype), (rows, cols)), shape=(self.D, self.D))

    def _rot(self):
        """R: labels (d_0..d_{M-1}) -> (d_{M-1}, d_0, ..., d_{M-2}), offset s -> s - 1 mod p."""
        rows, cols = [], []
        for s, b in enumerate(self.blocks):
            if b['n'] == 0:
                continue
            nb = self.blocks[(s - 1) % self.p]
            newdig = np.roll(b['dig'], 1, axis=0)
            full = (newdig * nb['strides'][:, None]).sum(0)
            r = nb['pos'][full]
            assert np.all(r >= 0)
            rows.append(r); cols.append(b['off'] + np.arange(b['n']))
        rows = np.concatenate(rows); cols = np.concatenate(cols)
        return sps.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(self.D, self.D))


class Positions:
    def __init__(self, L, M):
        self.L, self.M = L, M
        confs = list(itertools.combinations(range(L), M))
        self.confs = np.array(confs, int).reshape(len(confs), M)
        idx = {c: i for i, c in enumerate(confs)}
        self.C = len(confs)
        r, c = [], []
        br, bc = [], []
        A = np.zeros((M, self.C))
        tmap = np.zeros(self.C, int)
        wrap = np.zeros(self.C, bool)
        for i, x in enumerate(confs):
            occ = set(x)
            for a in range(M):
                for dx in (-1, 1):
                    y = x[a] + dx
                    if 0 <= y < L and y not in occ:
                        nx = list(x); nx[a] = y
                        r.append(idx[tuple(nx)]); c.append(i)
            if M < L and x[-1] == L - 1 and 0 not in occ:
                nx = (0,) + x[:-1]
                br.append(idx[nx]); bc.append(i)
            for j in range(M - 1):
                A[j, i] = float(x[j + 1] == x[j] + 1)
            A[M - 1, i] = float(x[0] == 0 and x[-1] == L - 1) if M < L else 1.0
            if x[-1] == L - 1:
                wrap[i] = True
                tmap[i] = idx[(0,) + tuple(v + 1 for v in x[:-1])]
            else:
                tmap[i] = idx[tuple(v + 1 for v in x)]
        self.Tb = sps.csr_matrix((np.ones(len(r)), (r, c)), shape=(self.C, self.C))
        self.B = sps.csr_matrix((np.ones(len(br)), (br, bc)), shape=(self.C, self.C))
        self.A = A
        self.tmap, self.wrap = tmap, wrap


class Ring:
    def __init__(self, L, word, Q2, t, J, mu_v, lam, alph=None, **kw):
        """alph = None: the demonstrator alphabet.  Otherwise a dict with 'letters' (type -> letter
        indices), 'q2' (integer charge label per letter) and 'hb' (type pair -> restricted bond); a complex
        bond switches to complex arithmetic and the complex momentum projector."""
        self.L, self.word, self.Q2 = L, word, Q2
        self.t, self.J, self.mu_v = t, J, mu_v
        if alph is None:
            hb, _ = VB.type_bonds(lam, **kw)
            self.inn = Internal(word, Q2, hb)
        else:
            self.inn = Internal(word, Q2, alph['hb'], alph['letters'], alph['q2'], alph.get('allowed'))
        self.cplx = np.issubdtype(self.inn.dtype, np.complexfloating)
        self.dtype = complex if self.cplx else float
        self.pos = Positions(L, len(word))
        self.D, self.C = self.inn.D, self.pos.C
        self.N = self.D * self.C
        self.shape = (self.D, self.C)
        self.const = -mu_v * len(word)
        self.BT = self.pos.B.T.tocsr()
        self.RT = self.inn.R.T.tocsr()
        nw = ~self.pos.wrap
        self.nw_src = np.where(nw)[0]; self.nw_dst = self.pos.tmap[nw]
        self.w_src = np.where(self.pos.wrap)[0]; self.w_dst = self.pos.tmap[self.pos.wrap]

    def H(self, v):
        W = v.reshape(self.shape)
        out = self.const * W.astype(np.result_type(W.dtype, self.inn.dtype))
        if self.t != 0 and self.pos.M < self.L:
            hop = (self.pos.Tb @ W.T).T
            hop += self.inn.R @ (self.pos.B @ W.T).T          # (R x B):     W -> R W B^T
            hop += self.RT @ (self.BT @ W.T).T                 # (R^T x B^T): W -> R^T W B
            out -= self.t * hop
        if self.J != 0:
            for j, O in enumerate(self.inn.ops):
                a = self.pos.A[j]
                if not a.any():
                    continue
                out += self.J * (O @ W) * a[None, :]
        return out.ravel()

    def T(self, v):
        W = v.reshape(self.shape)
        out = np.empty_like(W)
        out[:, self.nw_dst] = W[:, self.nw_src]
        if len(self.w_src):
            out[:, self.w_dst] = self.inn.R @ W[:, self.w_src]
        return out.ravel()

    def proj(self, v, m):
        L = self.L
        k = 2 * np.pi * m / L
        out = np.zeros_like(v)
        w = v.copy()
        if self.cplx:                                   # P_k = (1/L) sum_r exp(-i k r) T^r
            for r in range(L):
                out += np.exp(-1j * k * r) * w
                w = self.T(w)
            return out / L
        for r in range(L):
            c = np.cos(k * r)
            if abs(c) > 1e-14:
                out += c * w
            w = self.T(w)
        fac = 1.0 if (m % L == 0 or 2 * m == L) else 2.0
        return fac * out / L

    def lowest(self, m, nev=3, shift=None, tol=1e-11, seed=0, vecs=False):
        rng = np.random.default_rng(seed)
        if shift is None:
            shift = 10.0 + abs(self.const) + 4 * self.L
        v0 = rng.normal(size=self.N) + (1j * rng.normal(size=self.N) if self.cplx else 0)
        v0 = self.proj(v0, m)
        nproj = self.L // 2 + 1
        if self.N // max(nproj, 1) <= 300:
            # small: dense diagonalization of the projected matrix
            E = np.eye(self.N, dtype=self.dtype)
            Pm = np.column_stack([self.proj(e, m) for e in E])
            Hm = np.column_stack([self.H(e) for e in E])
            Hp = Pm.conj().T @ (Hm - shift * np.eye(self.N)) @ Pm
            w, U = np.linalg.eigh((Hp + Hp.conj().T) / 2)
            sel = w < -1e-9
            w, U = w[sel][:nev] + shift, U[:, sel][:, :nev]
            return (w, U) if vecs else w

        def mv(v):
            return self.proj(self.H(v) - shift * v, m)
        op = LinearOperator((self.N, self.N), matvec=mv, dtype=self.dtype)
        k = min(nev, self.N - 2)
        out = eigsh(op, k=k, which='SA', v0=v0, tol=tol, ncv=max(2 * k + 10, 24), maxiter=50000,
                    return_eigenvectors=vecs)
        if vecs:
            w, U = out
            o = np.argsort(w)
            return w[o] + shift, U[:, o]
        return np.sort(out) + shift

    # ------------------------------------------------------------------ observables of a state
    def adjacency(self, v):
        """<A_j> for each squeezed bond role, averaged: returns dict by type pair and the total."""
        W = v.reshape(self.shape)
        prob = (W ** 2).sum(0)                                       # over internal
        prob = prob / prob.sum()
        return float((self.pos.A.sum(0) * prob).sum())


def lowest_free(L, M, theta):
    """Hard-core bosons with hopping -1 on L sites, flux theta on the cut: ground energy via fermions."""
    phi = theta + np.pi * (M - 1)
    best = np.inf
    for l in range(-M - 2, L + 2):
        k = (2 * np.pi * (np.arange(l, l + M)) + phi) / L
        best = min(best, -2 * np.cos(k).sum())
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--word', required=True)
    ap.add_argument('--Q2', type=int, nargs='+', default=None)
    ap.add_argument('--m', type=int, nargs='+', default=None)
    ap.add_argument('--nev', type=int, default=3)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.25)
    ap.add_argument('--mu_v', type=float, default=0.0)
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--tol', type=float, default=1e-11)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    M = len(a.word)
    Q2s = a.Q2 if a.Q2 is not None else ([0, 2] if M % 2 == 0 else [1, 3])
    ms = a.m if a.m is not None else list(range(a.L // 2 + 1))
    res = dict(args=vars(a), levels={}, dims={})
    t0 = time.time()
    for Q2 in Q2s:
        ring = Ring(a.L, a.word, Q2, a.t, a.J, a.mu_v, a.lam)
        res['dims'][str(Q2)] = [ring.D, ring.C]
        for m in ms:
            w = ring.lowest(m, a.nev, tol=a.tol)
            res['levels'][f'{Q2}_{m}'] = [float(x) for x in w]
            print(f'L={a.L} word={a.word} 2Q={Q2} m={m} dim={ring.D}x{ring.C} levels={np.round(w, 10).tolist()} '
                  f't={time.time() - t0:.0f}s', flush=True)
            if a.out:
                json.dump(res, open(a.out + '.tmp', 'w'), indent=1)
                os.replace(a.out + '.tmp', a.out)


if __name__ == '__main__':
    main()
