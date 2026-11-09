"""Momentum-resolved ring ED of C173's Fibonacci-cusp model (PROOF eq. 3.2) in the blockaded space.

Hard-core bosons (spin language, so the ring is the periodic spin chain), no two adjacent when blockade=True.
    H = sum_j [ w b+_{j+1} b_j + conj(w) b+_j b_{j+1} + t2 (b+_{j+2} b_j + h.c.) + V2 n_j n_{j+2} + V3 n_j n_{j+3} ]
Units |t2| = 1 (t2 = -1, V2 = -2, V3 = -2/3 for the cusp model). A uniform twist theta gives each hop by d the
phase exp(i theta d / L) (flux theta through the ring). theta = pi reproduces C173's fermion-periodic rings at even N.

Basis: representatives r (minimal integer over rotations), |r,K> = N_r^{-1/2} sum_l e^{-iKl} T^l |r>, T: j -> j+1.
If y = T^m s with s a representative, <s,K|H|r,K> gets h * e^{iKm} * sqrt(P_r/P_s) (P = period).

Controls:
  hard rods (c = 1):     blockade=True, t2=0, V2=V3=0
  two XX legs (c = 2):   blockade=False, w=0, t2=-1, V2=V3=0 (NNN hop ignores the middle site)
"""
import numpy as np
import numba as nb
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@nb.njit(cache=True)
def rotr(x, p, L, mask):
    if p == 0:
        return x
    return ((x >> p) | (x << (L - p))) & mask


@nb.njit(cache=True)
def popc(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=True)
def rep_of(y, L, mask):
    """Return (rep, m) with y = T^m rep, i.e. rep = rotr(y, m). Rep = minimal integer over rotations.
    The minimum always has bit 0 set, so only rotations bringing a particle to site 0 are tried."""
    best = -1
    bm = 0
    for p in range(L):
        if (y >> p) & 1:
            z = rotr(y, p, L, mask)
            if best < 0 or z < best:
                best = z
                bm = p
    return best, bm


@nb.njit(cache=True)
def period(x, L, mask):
    for p in range(1, L + 1):
        if rotr(x, p, L, mask) == x:
            return p
    return L


@nb.njit(cache=True)
def _gen(L, N, gap, out, count_only):
    """Enumerate configurations with bit 0 set, N particles, cyclic gaps >= gap; keep the representatives."""
    mask = (np.int64(1) << L) - 1
    pos = np.zeros(N, np.int64)
    cnt = 0
    # iterative odometer over positions p_1 < ... < p_{N-1}, p_0 = 0
    if N == 0:
        return 0
    k = 1
    pos[0] = 0
    if N == 1:
        x = np.int64(1)
        r, m = rep_of(x, L, mask)
        if r == x:
            if not count_only:
                out[cnt] = x
            cnt += 1
        return cnt
    pos[1] = gap - 1  # will be incremented
    while k >= 1:
        pos[k] += 1
        # last allowed position for particle k: leave room for remaining and cyclic gap to particle 0
        lim = L - gap - (N - 1 - k) * gap
        if pos[k] > lim:
            k -= 1
            continue
        if k == N - 1:
            x = np.int64(0)
            for i in range(N):
                x |= np.int64(1) << pos[i]
            r, m = rep_of(x, L, mask)
            if r == x:
                if not count_only:
                    out[cnt] = x
                cnt += 1
        else:
            k += 1
            pos[k] = pos[k - 1] + gap - 1
    return cnt


def basis(L, N, blockade=True):
    gap = 2 if blockade else 1
    n = _gen(L, N, gap, np.zeros(1, np.int64), True)
    out = np.zeros(n, np.int64)
    _gen(L, N, gap, out, False)
    out.sort()
    mask = (1 << L) - 1
    per = np.array([period(int(x), L, mask) for x in out], np.int64) if n < 2000 else _periods(out, L)
    return out, per


@nb.njit(cache=True, parallel=True)
def _periods(reps, L):
    mask = (np.int64(1) << L) - 1
    P = np.zeros(len(reps), np.int64)
    for i in nb.prange(len(reps)):
        P[i] = period(reps[i], L, mask)
    return P


@nb.njit(cache=True)
def _lookup(reps, y):
    i = np.searchsorted(reps, y)
    if i < len(reps) and reps[i] == y:
        return i
    return -1


@nb.njit(cache=True)
def _ok_block(y, t, L):
    # after moving a particle to t, sites t-1 and t+1 must be empty
    return (((y >> ((t + 1) % L)) & 1) == 0) and (((y >> ((t - 1) % L)) & 1) == 0)


@nb.njit(cache=True)
def _row_terms(x, L, mask, blockade, ds, amps, V2, V3, reps, P, K, theta, rcol, tgt, val, ptr0):
    """Append the matrix elements of H|r> for rep x (column rcol). Returns new pointer."""
    p = ptr0
    # diagonal
    d = 0.0
    if V2 != 0.0:
        d += V2 * popc(x & rotr(x, 2, L, mask))
    if V3 != 0.0:
        d += V3 * popc(x & rotr(x, 3, L, mask))
    if d != 0.0:
        tgt[p] = rcol
        val[p] = d
        p += 1
    for j in range(L):
        if ((x >> j) & 1) == 0:
            continue
        x0 = x & ~(np.int64(1) << j)
        for a in range(len(ds)):
            dd = ds[a]
            t = (j + dd) % L
            if (x0 >> t) & 1:
                continue
            if blockade and not _ok_block(x0, t, L):
                continue
            y = x0 | (np.int64(1) << t)
            s, m = rep_of(y, L, mask)
            idx = _lookup(reps, s)
            if idx < 0:
                continue
            ph = amps[a] * np.exp(1j * (theta * dd / L + K * m)) * np.sqrt(P[rcol] / P[idx])
            tgt[p] = idx
            val[p] = ph
            p += 1
    return p


@nb.njit(cache=True, parallel=True)
def _build(reps, P, L, blockade, ds, amps, V2, V3, K, theta, maxper):
    n = len(reps)
    mask = (np.int64(1) << L) - 1
    tgt = np.zeros(n * maxper, np.int64)
    val = np.zeros(n * maxper, np.complex128)
    cnt = np.zeros(n, np.int64)
    for i in nb.prange(n):
        x = reps[i]
        # skip states incompatible with K: e^{-iK P} != 1
        ok = abs(np.exp(-1j * K * P[i]) - 1.0) < 1e-9
        if not ok:
            cnt[i] = 0
            continue
        p = _row_terms(x, L, mask, blockade, ds, amps, V2, V3, reps, P, K, theta, i, tgt[i * maxper:(i + 1) * maxper],
                       val[i * maxper:(i + 1) * maxper], 0)
        cnt[i] = p
    return tgt, val, cnt


def hamiltonian(L, N, K, w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, blockade=True, theta=0.0, B=None):
    """Return (H as csr in the compatible-rep subspace, reps, periods, index of kept reps)."""
    reps, P = B if B is not None else basis(L, N, blockade)
    ds, amps = [], []
    if w != 0:
        ds += [1, -1]; amps += [w, np.conj(w)]
    if t2 != 0:
        ds += [2, -2]; amps += [t2, t2]
    ds = np.array(ds, np.int64); amps = np.array(amps, np.complex128)
    maxper = 1 + N * len(ds)
    tgt, val, cnt = _build(reps, P, L, blockade, ds, amps, float(V2), float(V3), float(K), float(theta), maxper)
    n = len(reps)
    keep = np.nonzero(np.abs(np.exp(-1j * K * P) - 1.0) < 1e-9)[0]
    newidx = -np.ones(n, np.int64); newidx[keep] = np.arange(len(keep))
    rows, cols, vals = [], [], []
    starts = np.arange(n) * maxper
    sel = np.concatenate([starts[i] + np.arange(cnt[i]) for i in keep]) if len(keep) else np.zeros(0, np.int64)
    colsrc = np.repeat(keep, cnt[keep])
    r = newidx[tgt[sel]]; c = newidx[colsrc]; v = val[sel]
    good = r >= 0
    H = sp.csr_matrix((v[good], (r[good], c[good])), shape=(len(keep), len(keep)))
    H.sum_duplicates()
    return H, reps, P, keep


def lowest(H, k=6):
    n = H.shape[0]
    if n == 0:
        return np.zeros(0), np.zeros((0, 0))
    if n <= 400:
        e, v = np.linalg.eigh(H.toarray())
        return e[:k], v[:, :k]
    k = min(k, n - 2)
    e, v = spla.eigsh(H, k=k, which='SA', tol=1e-11, ncv=max(2 * k + 10, 30))
    o = np.argsort(e)
    return e[o], v[:, o]


# ---------------- observables ----------------

@nb.njit(cache=True, parallel=True)
def _expect_hops(reps, P, keep, newidx, L, K, blockade, fermion, psi, bra, rs):
    """g[r] = <psi| sum_j b+_{j+r} [(-1)^{# particles strictly between j and j+r going +}] b_j |psi> (projected)."""
    mask = (np.int64(1) << L) - 1
    n = len(keep)
    g = np.zeros(L, np.complex128)
    for r in rs:
        acc = 0.0 + 0.0j
        for ci in nb.prange(n):
            i = keep[ci]
            x = reps[i]
            loc = 0.0 + 0.0j
            for j in range(L):
                if ((x >> j) & 1) == 0:
                    continue
                t = (j + r) % L
                x0 = x & ~(np.int64(1) << j)
                if (x0 >> t) & 1:
                    continue
                if blockade and not _ok_block(x0, t, L):
                    continue
                sg = 1.0
                if fermion and r > 1:
                    mk = rotr((np.int64(1) << (r - 1)) - 1, (L - (j + 1) % L) % L, L, mask)
                    if popc(x0 & mk) % 2 == 1:
                        sg = -1.0
                y = x0 | (np.int64(1) << t)
                s, m = rep_of(y, L, mask)
                idx = _lookup(reps, s)
                if idx < 0:
                    continue
                k2 = newidx[idx]
                if k2 < 0:
                    continue
                loc += sg * np.exp(1j * K * m) * np.sqrt(P[i] / P[idx]) * np.conj(bra[k2]) * psi[ci]
            acc += loc
        g[r] = acc
    return g


def momentum_distribution(L, N, K, psi, reps, P, keep, blockade=True, fermion=True):
    """n(k) of JW fermions (fermion=True; boundary condition e^{ikL} = (-1)^{N-1} for periodic bosons) or bosons."""
    newidx = -np.ones(len(reps), np.int64); newidx[keep] = np.arange(len(keep))
    psi = np.ascontiguousarray(psi, dtype=np.complex128)
    g = _expect_hops(reps, P, keep, newidx, L, float(K), blockade, fermion, psi, psi, np.arange(1, L)) / L
    g[0] = N / L
    sh = 0.0 if (not fermion or (N - 1) % 2 == 0) else 0.5
    ks = 2 * np.pi * (np.arange(L) + sh) / L
    rr = np.arange(L)
    nk = np.array([(np.exp(1j * k * rr) * g).sum().real for k in ks])
    return ks, nk


def valley_imbalance(ks, nk):
    c = np.cos(ks)
    return nk[c > 1e-12].sum() - nk[c < -1e-12].sum()


@nb.njit(cache=True)
def _inv(x, L):
    y = np.int64(0)
    for j in range(L):
        if (x >> j) & 1:
            y |= np.int64(1) << ((L - j) % L)
    return y


def S_matrix(L, N, K, reps, P, keep, oddmask=None):
    """S = inversion x (-1)^{n on odd sites}, valid at twist 0, L even, in a self-conjugate sector 2K = N pi mod 2pi.
    S|r,K> = sigma(r) e^{-i(K+N pi) m} sqrt(P_s/P_r)... returns sparse matrix in the kept basis."""
    mask = (1 << L) - 1
    odd = sum(1 << j for j in range(1, L, 2))
    newidx = -np.ones(len(reps), np.int64); newidx[keep] = np.arange(len(keep))
    rows, cols, vals = [], [], []
    for ci, i in enumerate(keep):
        x = int(reps[i])
        sg = -1.0 if bin(x & odd).count('1') % 2 else 1.0
        y = int(_inv(np.int64(x), L))
        s, m = rep_of(np.int64(y), L, np.int64(mask))
        idx = np.searchsorted(reps, s)
        k2 = newidx[idx]
        # S|r,K> = sigma e^{-i(K+N pi) m} sqrt(N_s/N_r) |s,K'>,  N_r = L^2/P_r -> sqrt(P_r/P_s)
        rows.append(k2); cols.append(ci)
        vals.append(sg * np.exp(-1j * (K + N * np.pi) * m) * np.sqrt(P[i] / P[idx]))
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(keep), len(keep)))


@nb.njit(cache=True, parallel=True)
def csr_matvec(indptr, indices, data, x):
    n = len(indptr) - 1
    y = np.zeros(n, np.complex128)
    for i in nb.prange(n):
        acc = 0.0 + 0.0j
        for p in range(indptr[i], indptr[i + 1]):
            acc += data[p] * x[indices[p]]
        y[i] = acc
    return y


def lowest_par(H, k=6, tol=1e-10):
    """Lanczos with a numba-parallel matvec (scipy's csr matvec is single-threaded)."""
    n = H.shape[0]
    if n <= 1500:
        import scipy.linalg as sla
        e, v = sla.eigh(H.toarray())
        return e[:k], v[:, :k]
    ip, ix, dt = H.indptr.astype(np.int64), H.indices.astype(np.int64), H.data.astype(np.complex128)
    op = spla.LinearOperator((n, n), matvec=lambda v: csr_matvec(ip, ix, dt, np.ascontiguousarray(v, dtype=np.complex128).ravel()), dtype=np.complex128)
    k = min(k, n - 2)
    e, v = spla.eigsh(op, k=k, which='SA', tol=tol, ncv=max(2 * k + 12, 32))
    o = np.argsort(e)
    return e[o], v[:, o]


def imbalance_matrix(L, N, K, V, reps, P, keep, blockade=True, full=True):
    """Matrix of the valley-imbalance operator M = sum_k sign(cos k) n_F(k) in the span of the columns of V."""
    newidx = -np.ones(len(reps), np.int64); newidx[keep] = np.arange(len(keep))
    sh = 0.0 if (N - 1) % 2 == 0 else 0.5
    ks = 2 * np.pi * (np.arange(L) + sh) / L
    sgn = np.sign(np.cos(ks)); sgn[np.abs(np.cos(ks)) < 1e-12] = 0.0
    rr = np.arange(L)
    if full:
        coef = np.array([(sgn * np.exp(1j * ks * r)).sum() for r in rr])  # M = sum_r coef[r] g[r]
        rs = np.arange(1, L)
    else:  # C1 = sum_k cos(k) n_k, needs only r = 1 and r = L - 1
        coef = np.array([(np.cos(ks) * np.exp(1j * ks * r)).sum() for r in rr])
        rs = np.array([1, L - 1])
    d = V.shape[1]
    Mm = np.zeros((d, d), complex)
    Vc = [np.ascontiguousarray(V[:, a], dtype=np.complex128) for a in range(d)]
    for a in range(d):
        for b in range(d):
            g = _expect_hops(reps, P, keep, newidx, L, float(K), blockade, True, Vc[b], Vc[a], rs) / L
            g[0] = (N / L) if a == b else 0.0
            Mm[a, b] = (coef * g).sum()
    return Mm


@nb.njit(cache=True, parallel=True)
def _count(reps, P, L, blockade, ds, K, V2, V3):
    n = len(reps)
    mask = (np.int64(1) << L) - 1
    cnt = np.zeros(n, np.int64)
    for i in nb.prange(n):
        if abs(np.exp(-1j * K * P[i]) - 1.0) > 1e-9:
            continue
        x = reps[i]
        c = 1
        for j in range(L):
            if ((x >> j) & 1) == 0:
                continue
            x0 = x & ~(np.int64(1) << j)
            for a in range(len(ds)):
                t = (j + ds[a]) % L
                if (x0 >> t) & 1:
                    continue
                if blockade and not _ok_block(x0, t, L):
                    continue
                c += 1
        cnt[i] = c
    return cnt


@nb.njit(cache=True, parallel=True)
def _fill(reps, P, L, blockade, ds, amps, V2, V3, K, theta, keep, newidx, indptr, indices, data):
    """Row ci of H (Hermitian) = conj of column ci: H[ci, s] = conj(<s|H|r>)."""
    mask = (np.int64(1) << L) - 1
    for ci in nb.prange(len(keep)):
        i = keep[ci]
        x = reps[i]
        p = indptr[ci]
        d = 0.0
        if V2 != 0.0:
            d += V2 * popc(x & rotr(x, 2, L, mask))
        if V3 != 0.0:
            d += V3 * popc(x & rotr(x, 3, L, mask))
        indices[p] = ci
        data[p] = d
        p += 1
        for j in range(L):
            if ((x >> j) & 1) == 0:
                continue
            x0 = x & ~(np.int64(1) << j)
            for a in range(len(ds)):
                dd = ds[a]
                t = (j + dd) % L
                if (x0 >> t) & 1:
                    continue
                if blockade and not _ok_block(x0, t, L):
                    continue
                y = x0 | (np.int64(1) << t)
                s, m = rep_of(y, L, mask)
                idx = _lookup(reps, s)
                k2 = newidx[idx]
                if idx < 0 or k2 < 0:
                    indices[p] = ci
                    data[p] = 0.0
                else:
                    val = amps[a] * np.exp(1j * (theta * dd / L + K * m)) * np.sqrt(P[i] / P[idx])
                    indices[p] = k2
                    data[p] = np.conj(val)
                p += 1


def hamiltonian2(L, N, K, w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, blockade=True, theta=0.0, B=None):
    """Same as hamiltonian() but builds CSR directly in numba (memory-lean, for large sectors)."""
    reps, P = B if B is not None else basis(L, N, blockade)
    ds, amps = [], []
    if w != 0:
        ds += [1, -1]; amps += [w, np.conj(w)]
    if t2 != 0:
        ds += [2, -2]; amps += [t2, t2]
    ds = np.array(ds, np.int64); amps = np.array(amps, np.complex128)
    cnt = _count(reps, P, L, blockade, ds, float(K), float(V2), float(V3))
    keep = np.nonzero(cnt > 0)[0]
    newidx = -np.ones(len(reps), np.int64); newidx[keep] = np.arange(len(keep))
    indptr = np.zeros(len(keep) + 1, np.int64); indptr[1:] = np.cumsum(cnt[keep])
    indices = np.zeros(indptr[-1], np.int64); data = np.zeros(indptr[-1], np.complex128)
    _fill(reps, P, L, blockade, ds, amps, float(V2), float(V3), float(K), float(theta), keep, newidx, indptr, indices, data)
    H = sp.csr_matrix((data, indices, indptr), shape=(len(keep), len(keep)))
    return H, reps, P, keep


@nb.njit(cache=True, parallel=True)
def _S_build(reps, P, keep, newidx, L, N, K, oddmask):
    mask = (np.int64(1) << L) - 1
    n = len(keep)
    rows = np.zeros(n, np.int64)
    vals = np.zeros(n, np.complex128)
    for ci in nb.prange(n):
        i = keep[ci]
        x = reps[i]
        sg = -1.0 if popc(x & oddmask) % 2 == 1 else 1.0
        y = _inv(x, L)
        s, m = rep_of(y, L, mask)
        idx = _lookup(reps, s)
        rows[ci] = newidx[idx]
        vals[ci] = sg * np.exp(-1j * (K + N * np.pi) * m) * np.sqrt(P[i] / P[idx])
    return rows, vals


def S_matrix_fast(L, N, K, reps, P, keep):
    odd = np.int64(sum(1 << j for j in range(1, L, 2)))
    newidx = -np.ones(len(reps), np.int64); newidx[keep] = np.arange(len(keep))
    rows, vals = _S_build(reps, P, keep, newidx, L, N, float(K), odd)
    return sp.csr_matrix((vals, (rows, np.arange(len(keep)))), shape=(len(keep), len(keep)))
