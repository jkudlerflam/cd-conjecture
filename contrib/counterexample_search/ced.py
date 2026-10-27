#!/usr/bin/env python3
"""Generic symmetry-resolved exact diagonalization for translation-invariant
nearest-neighbour chains on (C^d)^L with periodic boundary conditions.

Structure follows codex/work_C4_D3/d3.py (additive charges, meet-in-the-middle
basis, translation-orbit momentum blocks, charge-resolved Schmidt values), but
is generic in d, in the bond, and in the Cartan charges, and adds global Casimir
labels. Site 0 is the least significant base-d digit. numpy/scipy only.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import numpy as np
from scipy import sparse, linalg
from scipy.sparse.linalg import eigsh


# ----------------------------------------------------------------- local algebra
def spin(s):
    m = np.arange(s, -s - 1, -1)
    z = np.diag(m)
    plus = np.zeros((len(m), len(m)))
    for j in range(1, len(m)):
        plus[j - 1, j] = np.sqrt(s * (s + 1) - m[j] * (m[j] + 1))
    return (plus + plus.T) / 2, (plus - plus.T) / (2j), z


def swap(d):
    return np.eye(d * d).reshape(d, d, d, d).transpose(1, 0, 2, 3).reshape(d * d, d * d)


def two_site_casimir(gens):
    d = gens[0].shape[0]
    one = np.eye(d)
    C = 0
    for t in gens:
        T = np.kron(t, one) + np.kron(one, t)
        C = C + T @ T
    return C


def irrep_projectors(gens, tol=1e-8):
    """Spectral projectors of the two-site quadratic Casimir.
    Returns list of (casimir_value, rank, projector) sorted by value."""
    C = two_site_casimir(gens)
    C = (C + C.conj().T) / 2
    w, v = linalg.eigh(C)
    out = []
    i = 0
    while i < len(w):
        j = i
        while j < len(w) and abs(w[j] - w[i]) < tol:
            j += 1
        V = v[:, i:j]
        P = V @ V.conj().T
        P[np.abs(P) < 1e-14] = 0
        if np.allclose(P.imag, 0, atol=1e-12):
            P = P.real
        out.append((float(np.mean(w[i:j])), j - i, P))
        i = j
    return out


# ----------------------------------------------------------------- basis
def enumerate_basis(d, L, charges, target):
    """All base-d words of length L whose summed integer charge equals target."""
    charges = np.asarray(charges, dtype=np.int64)
    if charges.ndim == 1:
        charges = charges[:, None]
    if charges.shape[1] == 0:
        return np.arange(d ** L, dtype=np.int64)
    ell = L // 2
    other = L - ell

    def half(n):
        code = np.arange(d ** n, dtype=np.int64)
        tmp = code.copy()
        q = np.zeros((len(code), charges.shape[1]), dtype=np.int64)
        for _ in range(n):
            q += charges[tmp % d]
            tmp //= d
        groups = {}
        uq, inv = np.unique(q, axis=0, return_inverse=True)
        inv = inv.ravel()
        order = np.argsort(inv, kind='stable')
        counts = np.bincount(inv)
        starts = np.r_[0, np.cumsum(counts)]
        for i, qi in enumerate(uq):
            groups[tuple(qi)] = code[order[starts[i]:starts[i + 1]]]
        return groups

    left = half(ell)
    right = half(other)
    out = []
    tgt = np.asarray(target, dtype=np.int64)
    for q, ls in left.items():
        rs = right.get(tuple(tgt - np.array(q)))
        if rs is not None:
            out.append((ls[:, None] + (d ** ell) * rs[None, :]).ravel())
    return np.sort(np.concatenate(out)) if out else np.empty(0, dtype=np.int64)


class Sector:
    """Fixed L and additive-charge sector, with translation orbits."""

    def __init__(self, d, L, hs, charges=None, target=None, cache=True):
        self.d, self.L = d, L
        self.hs = [np.asarray(h) for h in hs]
        if charges is None:
            charges = np.zeros((d, 0), dtype=np.int64)
        charges = np.asarray(charges, dtype=np.int64)
        if charges.ndim == 1:
            charges = charges[:, None]
        self.charges = charges
        if target is None:
            target = np.zeros(charges.shape[1], dtype=np.int64)
        self.target = np.asarray(target, dtype=np.int64)
        self.cache = cache
        self.basis = enumerate_basis(d, L, charges, self.target)
        b = self.basis
        rot = b.copy()
        rep = b.copy()
        ashift = np.zeros(len(b), dtype=np.int16)
        period = np.zeros(len(b), dtype=np.int16)
        top = d ** (L - 1)
        for a in range(1, L + 1):
            rot = rot // d + (rot % d) * top
            ix = rot < rep
            rep[ix] = rot[ix]
            ashift[ix] = a
            ix = (rot == b) & (period == 0)
            period[ix] = a
        del rot
        isrep = rep == b
        self.reps = b[isrep]
        self.rid = np.searchsorted(self.reps, rep)
        del rep
        self.period = period
        self.shift = ((period.astype(np.int32) - ashift) % period).astype(np.int16)
        self.rperiod = period[isrep]
        self._blocks = {}

    def check_charge(self):
        """Verify every bond term conserves the declared charges."""
        d = self.d
        for h in self.hs:
            for out in range(d * d):
                for inp in range(d * d):
                    if abs(h[out, inp]) > 1e-12:
                        a, b = divmod(inp, d)
                        c, e = divmod(out, d)
                        if np.any(self.charges[a] + self.charges[b] != self.charges[c] + self.charges[e]):
                            return False
        return True

    def blocks(self, k):
        k = k % self.L
        if k in self._blocks:
            return self._blocks[k]
        L, d = self.L, self.d
        allowed = (k * self.rperiod.astype(np.int64)) % L == 0
        ids = np.flatnonzero(allowed)
        reps = self.reps[ids]
        inverse = np.full(len(self.reps), -1, dtype=np.int64)
        inverse[ids] = np.arange(len(ids))
        n = len(ids)
        mats = []
        for h in self.hs:
            if n == 0:
                mats.append(sparse.csr_matrix((0, 0)))
                continue
            rows, cols, data = [], [], []
            for site in range(L):
                j = (site + 1) % L
                pi, pj = d ** site, d ** j
                ai = (reps // pi) % d
                aj = (reps // pj) % d
                pair = ai * d + aj
                for inp in range(d * d):
                    ci = np.flatnonzero(pair == inp)
                    if not len(ci):
                        continue
                    a, b = divmod(inp, d)
                    for op in np.flatnonzero(np.abs(h[:, inp]) > 1e-13):
                        c, e = divmod(int(op), d)
                        new = reps[ci] + (c - a) * pi + (e - b) * pj
                        bi = np.searchsorted(self.basis, new)
                        bi[bi >= len(self.basis)] = 0
                        if np.any(self.basis[bi] != new):
                            raise ValueError('bond term violates declared charges')
                        ri = inverse[self.rid[bi]]
                        use = ri >= 0
                        v = h[op, inp] * np.sqrt(self.rperiod[ids[ci[use]]] / self.period[bi[use]])
                        if k:
                            v = v * np.exp(2j * np.pi * k * self.shift[bi[use]] / L)
                        rows.append(ri[use])
                        cols.append(ci[use])
                        data.append(v)
            if data:
                H = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                                      shape=(n, n)).tocsr()
            else:
                H = sparse.csr_matrix((n, n))
            H.sum_duplicates()
            H.eliminate_zeros()
            if k == 0 or 2 * k == L:
                if H.nnz and np.max(np.abs(np.imag(H.data))) < 1e-12:
                    H = H.real
            mats.append(H)
        if self.cache:
            self._blocks[k] = (ids, mats)
        return ids, mats

    def solve(self, coeff, k, n=3, tol=1e-10):
        ids, hs = self.blocks(k)
        H = sum(float(c) * h for c, h in zip(coeff, hs) if c != 0)
        if not self.cache:
            self._blocks.clear()
        dim = len(ids)
        if dim == 0:
            return np.empty(0), np.empty((0, 0)), 0.0
        if np.isscalar(H):
            H = sparse.csr_matrix((dim, dim))
        if dim <= max(400, n + 2):
            e, v = linalg.eigh(H.toarray())
            e, v = e[:n], v[:, :n]
        else:
            rng = np.random.default_rng(1234 + 7 * k + 101 * self.L)
            v0 = rng.normal(size=dim)
            if np.iscomplexobj(H.data):
                v0 = v0 + 1j * rng.normal(size=dim)
            e, v = eigsh(H, k=min(n, dim - 2), which='SA', tol=tol, maxiter=20000, v0=v0,
                         ncv=min(dim - 1, max(40, 3 * n + 10)))
            ix = np.argsort(e)
            e, v = e[ix], v[:, ix]
        res = float(np.max(np.linalg.norm(H @ v - v * e, axis=0))) if len(e) else 0.0
        return e, v, res

    def expand(self, k, vec):
        ids = np.flatnonzero((k * self.rperiod.astype(np.int64)) % self.L == 0)
        rv = np.zeros(len(self.reps), dtype=complex)
        rv[ids] = vec
        psi = rv[self.rid] * np.exp(-2j * np.pi * k * self.shift / self.L) / np.sqrt(self.period.astype(float))
        nrm = np.vdot(psi, psi).real
        if abs(nrm - 1) > 1e-8:
            raise AssertionError('expansion normalization %g' % nrm)
        return psi

    def entropies(self, psi, cuts=None, renyi=False):
        """Von Neumann (and optionally Renyi-2) entropies of contiguous blocks
        of sites 0..ell-1, for ell in cuts (default 1..L//2)."""
        d, L = self.d, self.L
        cuts = range(1, L // 2 + 1) if cuts is None else cuts
        S, S2 = [], []
        for ell in cuts:
            left = self.basis % d ** ell
            right = self.basis // d ** ell
            if self.charges.shape[1]:
                lu = np.unique(left)
                q = np.zeros((len(lu), self.charges.shape[1]), dtype=np.int64)
                t = lu.copy()
                for _ in range(ell):
                    q += self.charges[t % d]
                    t //= d
                _, group = np.unique(q, axis=0, return_inverse=True)
                group = group.ravel()
                gl = group[np.searchsorted(lu, left)]
            else:
                gl = np.zeros(len(left), dtype=np.int64)
            w = []
            order = np.argsort(gl, kind='stable')
            counts = np.bincount(gl)
            starts = np.r_[0, np.cumsum(counts)]
            for g in range(len(counts)):
                sel = order[starts[g]:starts[g + 1]]
                ls, rs = left[sel], right[sel]
                us, li = np.unique(ls, return_inverse=True)
                vs, ri = np.unique(rs, return_inverse=True)
                mat = np.zeros((len(us), len(vs)), dtype=complex)
                mat[li.ravel(), ri.ravel()] = psi[sel]
                if min(mat.shape) == 0:
                    continue
                sv = linalg.svdvals(mat, check_finite=False)
                w.extend((sv * sv).tolist())
            w = np.array(w)
            w = w[w > 1e-16]
            if abs(w.sum() - 1) > 1e-8:
                raise AssertionError('Schmidt normalization')
            S.append(float(-np.dot(w, np.log(w))))
            S2.append(float(-np.log(np.dot(w, w))))
        return (S, S2) if renyi else S

    def casimir(self, gens, psis):
        """<C2> and variance ||(C2-<C2>)psi||^2 for columns of psis (sector vectors).
        gens: Hermitian d x d generators, normalized by the caller."""
        d, L = self.d, self.L
        psis = np.atleast_2d(psis.T).T
        ans = np.zeros(psis.shape[1])
        applied = np.zeros_like(psis, dtype=complex)
        key = id(gens)
        if getattr(self, '_cas_key', None) == key:
            ops = self._cas_ops
        else:
            ops = None
        if ops is not None:
            for op in ops:
                tr = op @ psis
                ans += np.sum(np.abs(tr) ** 2, axis=0)
                applied += op.getH() @ tr
            var = np.sum(np.abs(applied - psis * ans) ** 2, axis=0)
            return ans, var
        ops = []
        for g in gens:
            rows, cols, data = [], [], []
            for site in range(L):
                p = d ** site
                digit = self.basis // p % d
                for a in range(d):
                    ix = np.flatnonzero(digit == a)
                    for b in np.flatnonzero(np.abs(g[:, a]) > 1e-13):
                        rows.append(self.basis[ix] + (int(b) - a) * p)
                        cols.append(ix)
                        data.append(np.full(len(ix), g[b, a], dtype=complex))
            targets = np.concatenate(rows)
            ci = np.concatenate(cols)
            vals = np.concatenate(data)
            uniq, ri = np.unique(targets, return_inverse=True)
            op = sparse.coo_matrix((vals, (ri.ravel(), ci)), shape=(len(uniq), len(self.basis))).tocsr()
            ops.append(op)
            tr = op @ psis
            ans += np.sum(np.abs(tr) ** 2, axis=0)
            applied += op.getH() @ tr
        if self.cache:
            self._cas_key, self._cas_ops = key, ops
        var = np.sum(np.abs(applied - psis * ans) ** 2, axis=0)
        return ans, var


# ----------------------------------------------------------------- measurement
def chord(ell, L):
    return (L / np.pi) * np.sin(np.pi * np.asarray(ell) / L)


def entropy_fit(S, L, lmin=2, lmax=None):
    ell = np.arange(1, len(S) + 1)
    lmax = len(S) if lmax is None else lmax
    m = (ell >= lmin) & (ell <= lmax)
    x = np.log(chord(ell[m], L))
    y = np.asarray(S)[m]
    if len(y) < 2:
        return None
    A = np.stack([x, np.ones_like(x)], 1)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(3 * coef[0])


def measure(sec, coeff, nlev=3, momenta=None, entropy=True, gens=None, ncas=0):
    """Low spectrum at all momenta 0..L//2, ground-state entropies and optional
    Casimir labels of the lowest ncas states in each momentum block."""
    L = sec.L
    ks = range(L // 2 + 1) if momenta is None else momenta
    spec, vecs, res = {}, {}, []
    for k in ks:
        e, v, r = sec.solve(coeff, k, nlev)
        if len(e):
            spec[k] = e
            vecs[k] = v
            res.append(r)
    k0 = min(spec, key=lambda k: spec[k][0])
    e0 = spec[k0][0]
    out = {'L': L, 'k0': int(k0), 'E0': float(e0), 'dim': int(len(sec.basis)),
           'spec': {int(k): [float(x) for x in e] for k, e in spec.items()},
           'residual': float(max(res))}
    deg = [k for k in spec if abs(spec[k][0] - e0) < 1e-9]
    out['ground_degenerate_momenta'] = [int(k) for k in deg]
    if entropy:
        psi = sec.expand(k0, vecs[k0][:, 0])
        S, S2 = sec.entropies(psi, renyi=True)
        out['S'] = S
        out['S2'] = S2
    if gens is not None and ncas:
        cas = {}
        for k, v in vecs.items():
            m = min(ncas, v.shape[1])
            psis = np.stack([sec.expand(k, v[:, i]) for i in range(m)], 1)
            c, var = sec.casimir(gens, psis)
            cas[int(k)] = [[float(a), float(b)] for a, b in zip(c, var)]
        out['casimir'] = cas
    return out
