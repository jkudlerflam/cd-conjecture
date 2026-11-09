"""C175: momentum-resolved ring ED of model (3.2), hard-core bosons (= the spin chain), periodic ring.
Blockade (no two adjacent, cyclically), NN hop w (amplitude to move +1), NNN hop -1 over an empty site,
V2 = -2 at distance 2, V3 = -2/3 at distance 3.  Units |t2| = 1.  Optional flux theta: hop +d gets e^{i theta d / L}.
Usage as a module: levels(L, N, w, m, k) -> lowest k eigenvalues at momentum K = 2 pi m / L."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, itertools
from functools import lru_cache

def configs(L, N):
    out = []
    def rec(start, left, acc, first):
        if left == 0:
            if first is None or True: out.append(acc)
            return
        for x in range(start, L):
            if first is not None and x == L-1 and first == 0: continue   # cyclic blockade
            rec(x+2, left-1, acc | (1 << x), x if first is None else first)
    rec(0, N, 0, None)
    return np.array(out, dtype=np.int64)

def rot(b, s, L):
    mask = (1 << L)-1
    return ((b << s) | (b >> (L-s))) & mask if s else b

@lru_cache(maxsize=None)
def reps(L, N):
    b = configs(L, N)
    allr = np.stack([rot(b, s, L) for s in range(L)])   # (L, n)
    rep = allr.min(axis=0)
    shift = allr.argmin(axis=0)                          # rot(b, shift) = rep
    R = np.unique(rep)
    # period of each rep
    per = np.array([next(s for s in range(1, L+1) if rot(int(r), s % L, L) == int(r)) for r in R])
    return b, rep, shift, R, per

def hamiltonian(L, N, w, m, t2=-1.0, V2=-2.0, V3=-2/3, theta=0.0):
    b, rep, shift, R, per = reps(L, N)
    K = 2*np.pi*m/L
    # momentum states exist only if K * per is a multiple of 2 pi
    good = np.abs(np.exp(1j*K*per)-1) < 1e-9
    Rg = R[good]; pg = per[good]
    index = {int(r): i for i, r in enumerate(Rg)}
    n = len(Rg)
    norm = np.sqrt(L/pg)       # |K,r> = (1/sqrt(per..)) sum_s e^{-iKs} T^s r ; use standard normalization factor
    bits = lambda arr, j: (arr >> (j % L)) & 1
    # diagonal
    diag = np.zeros(n)
    for j in range(L):
        diag += V2*(bits(Rg, j) & bits(Rg, j+2)) + V3*(bits(Rg, j) & bits(Rg, j+3))
    rows = [np.arange(n)]; cols = [np.arange(n)]; vals = [diag.astype(complex)]
    lookup_b = {int(x): (int(rr), int(ss)) for x, rr, ss in zip(b, rep, shift)}
    for d, amp in ((1, w), (-1, np.conj(w)), (2, t2), (-2, t2)):
        for j in range(L):
            a = j; c = (j+d) % L
            ok = (bits(Rg, a) == 1) & (bits(Rg, c) == 0) & (bits(Rg, c+1 if d > 0 else c-1) == 0)
            if abs(d) == 2:
                ok &= bits(Rg, j+d//2) == 0
            # blockade of target on the other side
            src = np.nonzero(ok)[0]
            for i in src:
                t = (int(Rg[i]) & ~(1 << a)) | (1 << c)
                if t not in lookup_b: continue
                rr, ss = lookup_b[t]
                if rr not in index: continue
                k = index[rr]
                # H|K,r_i> contains T^{-?}... phase: t = T^{-ss} rr  (rot(t, ss) = rr)
                ph = np.exp(-1j*K*ss) * norm[k]/norm[i] * np.exp(1j*theta*d/L)
                rows.append(np.array([k])); cols.append(np.array([i])); vals.append(np.array([amp*ph]))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    return H, Rg

def levels(L, N, w, m, k=4, **kw):
    H, Rg = hamiltonian(L, N, w, m, **kw)
    err = abs(H-H.getH()).max() if H.shape[0] else 0
    assert err < 1e-9, err
    n = H.shape[0]
    if n == 0: return np.array([])
    if n <= 400:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(spla.eigsh(H, k=min(k, n-2), which='SA', tol=1e-10)[0])
