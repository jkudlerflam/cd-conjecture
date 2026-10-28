#!/usr/bin/env python3
"""Momentum-resolved exact diagonalization of H = sum_j h_{j,j+1} on a periodic ring of
L sites with local dimension d (numpy/scipy only).

Site j of basis state s is the base-d digit (s // d**j) % d.  h is (d*d, d*d) with the
left site the more significant kron factor, as in twobody.py and vumps.py.

Functions
  sector_levels(h, L, k, nev)    lowest nev levels at momentum 2 pi k / L
  spectrum(h, L, nev)            all momenta
  ground_state_full(h, L, k)     lowest state at momentum k as a full d**L vector
  entropies(psi, L, d)           von Neumann entropy of the first l sites, l = 1..L-1
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

_cache = {}


def orbit_tables(L, d):
    key = (L, d)
    if key in _cache:
        return _cache[key]
    N = d ** L
    s = np.arange(N, dtype=np.int64)
    rep = s.copy(); shift = np.zeros(N, np.int64)
    cur = s.copy()
    period = np.full(N, L, np.int64)
    top = d ** (L - 1)
    for l in range(1, L):
        # translate by one site: site j -> site j+1 (cyclic)
        cur = (cur % top) * d + cur // top
        better = cur < rep
        rep = np.where(better, cur, rep); shift = np.where(better, l, shift)
        same = (cur == s) & (period == L)
        period = np.where(same, l, period)
    _cache[key] = (rep, shift, period)
    return rep, shift, period


def sector_basis(L, d, k):
    rep, shift, period = orbit_tables(L, d)
    s = np.arange(d ** L, dtype=np.int64)
    isrep = rep == s
    r = s[isrep]; p = period[isrep]
    ok = (k * p) % L == 0           # e^{iKp} = 1
    return r[ok], p[ok]


def sector_matrix(h, L, k, d=2, tol=1e-14):
    rep, shift, period = orbit_tables(L, d)
    reps, per = sector_basis(L, d, k)
    index = -np.ones(d ** L, np.int64); index[reps] = np.arange(len(reps))
    K = 2 * np.pi * k / L
    rows, cols, vals = [], [], []
    h = np.asarray(h, complex)
    for j in range(L):
        jr = (j + 1) % L
        a = (reps // d ** j) % d; b = (reps // d ** jr) % d
        bin_ = a * d + b
        base = reps - a * d ** j - b * d ** jr
        for ib in range(d * d):
            m = bin_ == ib
            if not m.any():
                continue
            src = np.where(m)[0]
            for ob in range(d * d):
                amp = h[ob, ib]
                if abs(amp) < tol:
                    continue
                na, nb = divmod(ob, d)
                new = base[m] + na * d ** j + nb * d ** jr
                r2 = rep[new]; l2 = shift[new]
                tgt = index[r2]
                ok = tgt >= 0
                # |K,r> = (sqrt(p_r)/L) sum_{l<L} e^{-iKl} T^l |r>, and H|r> = amp T^{l2}|r'>, so
                # <K,r'|H|K,r> = amp e^{iK l2} sqrt(p_r/p_r')
                norm = np.sqrt(per[src[ok]] / per[tgt[ok]])
                rows.append(tgt[ok]); cols.append(src[ok])
                vals.append(amp * np.exp(1j * K * l2[ok]) * norm)
    n = len(reps)
    H = sparse.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    return H, reps, per


def sector_levels(h, L, k, nev=4, d=2, return_vecs=False):
    H, reps, per = sector_matrix(h, L, k, d)
    n = H.shape[0]
    if n <= 600:
        Hd = H.toarray(); Hd = (Hd + Hd.conj().T) / 2
        w, v = np.linalg.eigh(Hd)
        w, v = w[:nev], v[:, :nev]
    else:
        # explicit random start: ARPACK's own start vector lies in the range of H and misses
        # exact zero modes such as a decoupled vacuum
        rng = np.random.default_rng(12345)
        v0 = rng.normal(size=n) + 1j * rng.normal(size=n)
        # and ARPACK also drops exact null vectors of the operator itself, so shift H to be
        # positive definite (Lanczos is shift invariant otherwise)
        sh = 1.0 + float(abs(H).sum(axis=1).max())
        Hs = H + sh * sparse.identity(n, format='csr')
        w, v = eigsh(Hs, k=min(nev, n - 2), which='SA', tol=1e-11, ncv=max(2 * nev + 10, 30), v0=v0)
        w = w - sh
        o = np.argsort(w); w, v = w[o], v[:, o]
    if return_vecs:
        return w, v, reps, per
    return w


def spectrum(h, L, nev=4, d=2):
    return {k: sector_levels(h, L, k, nev, d) for k in range(L)}


def ground_state_full(h, L, k=0, d=2, which=0):
    w, v, reps, per = sector_levels(h, L, k, which + 1, d, return_vecs=True)
    K = 2 * np.pi * k / L
    psi = np.zeros(d ** L, complex)
    top = d ** (L - 1)
    cur = reps.copy()
    # accumulating l < L visits each orbit element L/p times, so weight v_r sqrt(p_r)
    coef = v[:, which] * np.sqrt(per)
    for l in range(L):
        np.add.at(psi, cur, coef * np.exp(-1j * K * l))
        cur = (cur % top) * d + cur // top
    psi /= np.linalg.norm(psi)
    return w[which], psi


def entropies(psi, L, d=2, alpha=1):
    """S(l) for the block of sites 0..l-1 (the least significant digits)."""
    out = []
    for l in range(1, L):
        M = psi.reshape(d ** (L - l), d ** l)
        s = np.linalg.svd(M, compute_uv=False) ** 2
        s = s[s > 1e-15]
        out.append(float(-np.sum(s * np.log(s))) if alpha == 1 else float(np.log(np.sum(s ** alpha)) / (1 - alpha)))
    return np.array(out)
