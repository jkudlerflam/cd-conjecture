"""Momentum-resolved ring ED of C173's effective model (PROOF eq. 3.2), written independently of C173's code.

Hard-core bosons (the qubit chain itself, periodic), no two adjacent particles (also across the seam),
    H = sum_j [ w b+_{j+1} b_j + t2 b+_{j+2} b_j + h.c. ] + V2 sum n_j n_{j+2} + V3 sum n_j n_{j+3} - mu N,
with t2 = -1, V2 = -2, V3 = -2/3 (units |t2| = 1). Optional flux Phi: a hop by d sites forward gets e^{i Phi d / L}.
Basis states |K, r> = sqrt(per_r)/L sum_d e^{-iKd} T^d |r>, T = translation by one site.
"""
import numpy as np
import scipy.sparse as sp
from itertools import combinations


def ring_configs(L, N):
    out = []
    for c in combinations(range(L), N):
        s = 0
        ok = True
        for a in c:
            s |= 1 << a
        if s & ((s << 1) | (s >> (L - 1))) & ((1 << L) - 1):
            continue
        out.append(s)
    return np.array(out, dtype=np.int64)


def rot(s, d, L):
    d %= L
    m = (1 << L) - 1
    return ((s << d) | (s >> (L - d))) & m


def reps(L, N):
    cf = ring_configs(L, N)
    # representative = min over rotations; shift l with config = T^l rep, T: site j -> j+1 is rot by 1
    allrot = np.stack([rot(cf, d, L) for d in range(L)])  # allrot[d] = T^d cf
    rep = allrot.min(axis=0)
    # l such that cf = T^l rep  <=> rep = T^{-l} cf = allrot[(L-l)%L]
    dmin = allrot.argmin(axis=0)          # rep = T^{dmin} cf  -> cf = T^{-dmin} rep
    shift = (-dmin) % L
    per = np.array([next(p for p in range(1, L + 1) if rot(r, p, L) == r) for r in np.unique(rep)])
    ureps = np.unique(rep)
    return cf, rep, shift, ureps, per


def build_k(L, N, K, w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, mu=0.0, Phi=0.0, data=None):
    """Hamiltonian in the momentum sector K (K = 2 pi m / L). Returns sparse matrix and the rep list kept."""
    if data is None:
        data = reps(L, N)
    cf, rep, shift, ureps, per = data
    keep = np.array([abs(np.exp(-1j * K * p) - 1) < 1e-9 for p in per])
    R = ureps[keep]; P = per[keep]
    idx = {int(r): i for i, r in enumerate(R)}
    cmap = {int(c): (int(r), int(l)) for c, r, l in zip(cf, rep, shift)}
    n = len(R)
    rows, cols, vals = [], [], []
    bit = lambda s, j: (s >> (j % L)) & 1
    for i, r in enumerate(R):
        r = int(r)
        diag = -mu * N
        for j in range(L):
            if bit(r, j):
                diag += V2 * bit(r, j + 2) + V3 * bit(r, j + 3)
        rows.append(i); cols.append(i); vals.append(diag)
        for j in range(L):
            if not bit(r, j):
                continue
            for d, amp in ((1, w), (-1, np.conj(w)), (2, t2), (-2, t2)):
                tgt = (j + d) % L
                if bit(r, tgt):
                    continue
                s = r & ~(1 << j) | (1 << tgt)
                if s not in cmap:   # violates the blockade
                    continue
                rs, l = cmap[s]
                if rs not in idx:
                    continue
                ph = np.exp(1j * Phi * d / L)
                rows.append(idx[rs]); cols.append(i)
                vals.append(amp * ph * np.exp(1j * K * l) * np.sqrt(P[i] / P[idx[rs]]))
    H = sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)
    return H, R, P


def rstat(e, frac=0.1):
    e = np.sort(e)
    n = len(e); e = e[int(frac * n): int((1 - frac) * n)]
    s = np.diff(e)
    s = s[s > 1e-10 * max(1.0, np.abs(e).max())] if False else s
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    r = r[np.isfinite(r)]
    return r.mean(), len(r), int((s < 1e-9).sum())
