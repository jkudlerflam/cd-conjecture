"""Vectorized real-space ring ED of (3.2) as hard-core bosons (the qubit chain, periodic), with flux Phi:
a hop by d sites forward carries e^{i Phi d / L}. Validated against ringk.py (validate_ringv.out)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from itertools import combinations

def basis(L, N):
    out = []
    def rec(pos, left, cur):
        if left == 0:
            out.append(cur); return
        if L - pos < 2 * left - 1:
            return
        rec(pos + 2, left - 1, cur | (1 << pos)); rec(pos + 1, left, cur)
    rec(0, N, 0)
    b = np.array(sorted(out), dtype=np.int64)
    return b[~((((b >> (L - 1)) & 1) == 1) & ((b & 1) == 1))]

def build(L, N, w, Phi=0.0, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, b=None):
    if b is None:
        b = basis(L, N)
    n = len(b); bit = lambda j: (b >> (j % L)) & 1
    diag = np.zeros(n)
    for j in range(L):
        diag += V2 * (bit(j) & bit(j + 2)) + V3 * (bit(j) & bit(j + 3))
    rows = [np.arange(n)]; cols = [np.arange(n)]; vals = [diag.astype(complex)]
    for d, amp in ((1, w), (2, t2)):
        for j in range(L):
            a, c = j, (j + d) % L
            ok = (bit(a) == 1) & (bit(c) == 0) & (bit(c + 1) == 0)
            if d == 1:
                pass   # left neighbour of c is a, vacated
            else:
                ok &= (bit(j + 1) == 0)
            src = np.nonzero(ok)[0]
            if len(src) == 0:
                continue
            tgt = np.searchsorted(b, (b[src] & ~(1 << a)) | (1 << c))
            f = amp * np.exp(1j * Phi * d / L)
            rows += [tgt, src]; cols += [src, tgt]; vals += [np.full(len(src), f), np.full(len(src), np.conj(f))]
    return b, sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))

def e0(L, N, w, Phi=0.0, k=1, b=None):
    b, H = build(L, N, w, Phi, b=b)
    if H.shape[0] < 800:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(spla.eigsh(H, k=k, which='SA', tol=1e-11)[0])
