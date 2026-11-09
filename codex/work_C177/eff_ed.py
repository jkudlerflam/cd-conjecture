"""Vectorized ED of the blockaded second-order model (PROOF.md section 3), open chain, fixed N.

Units |t2| = 1.  H = P sum_j [ (w c+_{j+1} c_j + h.c.) + t2 (c+_{j+2} c_j + h.c.) + V2 n_j n_{j+2} + V3 n_j n_{j+3} ] P
with t2 = -1, V2 = 2 t2 = -2, V3 = (2/3) t2 (the values derived from the pairing), P = no two adjacent particles.
On an open chain particles never pass each other, so all fermion signs are +1.
Options: blockade off and V2 = V3 = 0 gives free fermions with a two-harmonic band (c = 2 control when four
Fermi points are present).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from math import comb


def basis_fixed(L, N, blockade=True):
    out = []

    def rec(pos, left, cur):
        if left == 0:
            out.append(cur); return
        if pos >= L:
            return
        # need room
        if blockade and L - pos < 2 * left - 1:
            return
        if not blockade and L - pos < left:
            return
        rec(pos + (2 if blockade else 1), left - 1, cur | (1 << pos))
        rec(pos + 1, left, cur)
    rec(0, N, 0)
    return np.array(sorted(out), dtype=np.int64)


def build(L, N, w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, blockade=True):
    b = basis_fixed(L, N, blockade)
    n = len(b)
    bit = lambda j: (b >> j) & 1
    diag = np.zeros(n)
    for j in range(L - 2):
        diag += V2 * (bit(j) & bit(j + 2))
    for j in range(L - 3):
        diag += V3 * (bit(j) & bit(j + 3))
    rows = [np.arange(n)]; cols = [np.arange(n)]; vals = [diag.astype(complex)]
    for (dist, amp) in ((1, w), (2, t2)):
        for j in range(L - dist):
            a, c = j, j + dist  # move a -> c
            ok = (bit(a) == 1) & (bit(c) == 0)
            if dist == 2 and blockade:
                ok &= (bit(j + 1) == 0)  # automatic under the blockade
            if blockade:
                if c + 1 < L:
                    ok &= (bit(c + 1) == 0)
                # left neighbour of c is c-1; for dist 1 it is a (now empty); for dist 2 it is j+1 (empty)
            src = np.nonzero(ok)[0]
            if len(src) == 0:
                continue
            tgt_states = (b[src] & ~(1 << a)) | (1 << c)
            tgt = np.searchsorted(b, tgt_states)
            # sign: number of particles strictly between a and c
            if dist == 2 and not blockade:
                sgn = 1.0 - 2.0 * bit(j + 1)[src]  # free fermions: passing an occupied site gives -1
            else:
                sgn = np.ones(len(src))
            # forward move a -> c carries amplitude amp for c+_{c} c_a (c = a + dist)
            rows += [tgt, src]; cols += [src, tgt]
            vals += [amp * sgn, np.conj(amp) * sgn]
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    return b, H


def ground(L, N, w, k=1, **kw):
    b, H = build(L, N, w, **kw)
    if H.shape[0] <= 400:
        e, v = np.linalg.eigh(H.toarray())
        return b, e[:k], v[:, :k]
    e, v = spla.eigsh(H, k=k, which='SA', tol=1e-10)
    o = np.argsort(e)
    return b, e[o], v[:, o]


def entropies(b, psi, L):
    """Bipartite von Neumann entropy for every cut l = 1..L-1 (left block = sites < l)."""
    S = []
    for l in range(1, L):
        left = b & ((1 << l) - 1)
        right = b >> l
        ul, il = np.unique(left, return_inverse=True)
        ur, ir = np.unique(right, return_inverse=True)
        M = np.zeros((len(ul), len(ur)), dtype=complex)
        M[il, ir] = psi
        s = np.linalg.svd(M, compute_uv=False)
        p = s ** 2; p = p[p > 1e-14]
        S.append(float(-(p * np.log(p)).sum()))
    return np.array(S)


def cc_fit(S, L, lo=None, hi=None):
    """Fit S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + a over the central window, also averaged over l parity."""
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    lo = lo or max(2, L // 6); hi = hi or L - lo
    sel = (l >= lo) & (l <= hi)
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    c6, a = np.linalg.lstsq(A, S[sel], rcond=None)[0]
    # parity-smoothed: average neighbouring cuts
    Sm = 0.5 * (S[:-1] + S[1:]); xm = 0.5 * (x[:-1] + x[1:]); lm = 0.5 * (l[:-1] + l[1:])
    sm = (lm >= lo) & (lm <= hi)
    A2 = np.vstack([xm[sm], np.ones(sm.sum())]).T
    c62 = np.linalg.lstsq(A2, Sm[sm], rcond=None)[0][0]
    return 6 * c6, 6 * c62
