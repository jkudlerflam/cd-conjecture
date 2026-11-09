"""C186: exact momentum-resolved ring ED of the mixed second-order cusp model, all particle numbers (no truncation).

    H = sum_j P[ w b+_{j+1} b_j + h.c. - (b+_{j+2} b_j + h.c.) + V2 n_j n_{j+2} + V3 n_j n_{j+3} - mu n_j
                 + lam (b+_j + b_j) + ( r b+_j (n_{j+2} - n_{j-2}) + h.c. ) ]P
Units |t2| = 1, spin language, ring of L sites, P = no two adjacent particles (also across the seam).
Bloch states |K;s> = N^-1/2 sum_x e^{iKx} T^x |s>, T shifts every particle by +1.
"""
import numpy as np, numba as nb, scipy.sparse as sp, scipy.sparse.linalg as spla, sys, json, time


def pconfigs(L):
    A0 = np.array([0], np.int64); A1 = np.array([1], np.int64)     # length-1 strings ending in 0 / 1
    for n in range(1, L):
        A0, A1 = np.concatenate([A0, A1]), A0 | (np.int64(1) << n)
    S = np.concatenate([A0, A1])
    ok = ~(((S >> (L - 1)) & 1).astype(bool) & ((S & 1).astype(bool)))
    return np.sort(S[ok])


@nb.njit(cache=True)
def rotl(s, l, L):
    l = l % L
    mask = (np.int64(1) << L) - 1
    return ((s << l) | (s >> (L - l))) & mask


@nb.njit(cache=True)
def rep_of(s, L):
    best = s; shift = 0
    t = s
    for l in range(1, L):
        t = rotl(s, l, L)
        if t < best:
            best = t; shift = l
    # best = rotl(s, shift) = T^shift s  ->  s = T^{-shift} best
    return best, shift


def pconfigs_rot(S, L):
    """empty helper (reps_and_periods only needs each orbit's representative to be present)"""
    return np.zeros(0, np.int64)


@nb.njit(cache=True)
def period_of(s, L):
    for l in range(1, L + 1):
        if rotl(s, l, L) == s:
            return l
    return L


@nb.njit(cache=True)
def reps_and_periods(S, L):
    n = S.shape[0]
    isrep = np.zeros(n, np.bool_)
    for i in range(n):
        r, sh = rep_of(S[i], L)
        isrep[i] = (r == S[i])
    R = S[isrep]
    per = np.empty(R.shape[0], np.int64)
    for i in range(R.shape[0]):
        per[i] = period_of(R[i], L)
    return R, per


@nb.njit(cache=True)
def occ(s, j, L):
    return (s >> (j % L)) & 1


@nb.njit(cache=True)
def build_coo(R, per, L, K, wr, wi, V2, V3, mu, rr, lamr, maxnnz, Nmax):
    n = R.shape[0]
    rows = np.empty(maxnnz, np.int64); cols = np.empty(maxnnz, np.int64)
    vr = np.empty(maxnnz, np.float64); vi = np.empty(maxnnz, np.float64)
    nz = 0
    Nsite = np.empty(n, np.float64)
    for i in range(n):
        s = R[i]
        N = 0; dia = 0.0
        for j in range(L):
            if occ(s, j, L):
                N += 1
                if occ(s, j + 2, L): dia += V2
                if occ(s, j + 3, L): dia += V3
        Nsite[i] = N
        dia -= mu * N
        rows[nz] = i; cols[nz] = i; vr[nz] = dia; vi[nz] = 0.0; nz += 1
        # off-diagonal moves: list of (target config, amplitude re, im)
        for j in range(L):
            if occ(s, j, L):
                base = s & ~(np.int64(1) << j)
                # NN hop j -> j+1 (amp w), j -> j-1 (amp conj w)
                for (dj, ar, ai) in ((1, wr, wi), (-1, wr, -wi)):
                    tj = (j + dj) % L
                    if occ(s, j + 2 * dj, L) == 0 and L > 4:
                        t = base | (np.int64(1) << tj)
                        r_, sh = rep_of(t, L)
                        k = np.searchsorted(R, r_)
                        if k < n and R[k] == r_:
                            # t = T^{-sh} r_  -> l = -sh ; phase e^{-iK l} = e^{iK sh}
                            ph = K * sh
                            fac = np.sqrt(per[i] / per[k])
                            c, sn = np.cos(ph), np.sin(ph)
                            rows[nz] = k; cols[nz] = i
                            vr[nz] = fac * (ar * c - ai * sn); vi[nz] = fac * (ar * sn + ai * c); nz += 1
                # NNN hops
                for dj in (2, -2):
                    tj = (j + dj) % L
                    sg = 1 if dj > 0 else -1
                    if occ(s, j + dj, L) == 0 and occ(s, j + dj + sg, L) == 0 and L > 6:
                        t = base | (np.int64(1) << tj)
                        r_, sh = rep_of(t, L)
                        k = np.searchsorted(R, r_)
                        if k < n and R[k] == r_:
                            ph = K * sh
                            fac = np.sqrt(per[i] / per[k])
                            rows[nz] = k; cols[nz] = i
                            vr[nz] = -fac * np.cos(ph); vi[nz] = -fac * np.sin(ph); nz += 1
                # annihilation at j: amp conj(lam + r (n_{j+2} - n_{j-2})), r and lam real
                amp = lamr + rr * (occ(s, j + 2, L) - occ(s, j - 2, L))
                if amp != 0.0:
                    t = base
                    r_, sh = rep_of(t, L)
                    k = np.searchsorted(R, r_)
                    if k < n and R[k] == r_:
                        ph = K * sh
                        fac = np.sqrt(per[i] / per[k])
                        rows[nz] = k; cols[nz] = i
                        vr[nz] = fac * amp * np.cos(ph); vi[nz] = fac * amp * np.sin(ph); nz += 1
            else:
                # creation at j (needs j-1, j+1 empty), only below the particle-number cap
                if N < Nmax and occ(s, j - 1, L) == 0 and occ(s, j + 1, L) == 0:
                    amp = lamr + rr * (occ(s, j + 2, L) - occ(s, j - 2, L))
                    if amp != 0.0:
                        t = s | (np.int64(1) << j)
                        r_, sh = rep_of(t, L)
                        k = np.searchsorted(R, r_)
                        if k < n and R[k] == r_:
                            ph = K * sh
                            fac = np.sqrt(per[i] / per[k])
                            rows[nz] = k; cols[nz] = i
                            vr[nz] = fac * amp * np.cos(ph); vi[nz] = fac * amp * np.sin(ph); nz += 1
    return rows[:nz], cols[:nz], vr[:nz], vi[:nz], Nsite


def anchored(L, Nmax):
    """all P configs on the ring with site 0 occupied and 1 <= N <= Nmax, as int64 masks"""
    out = []
    def rec(mask, last, n):
        out.append(mask)
        if n == Nmax:
            return
        for x in range(last + 2, L - 1):      # x <= L-2 keeps the seam gap (x, 0) >= 2
            rec(mask | (1 << x), x, n + 1)
    rec(1, 0, 1)
    return np.array(out, np.int64)


class Ring:
    def __init__(self, L, Nmax=None):
        self.L = L
        self.Nmax = L if Nmax is None else Nmax
        if Nmax is None:
            S = pconfigs(L)
        else:
            assert L <= 62
            S = np.unique(np.concatenate([[0], anchored(L, Nmax)]))
            # canonical representatives of these orbits
            S = np.unique(np.array([rep_of(int(x), L)[0] for x in S], np.int64))
            S = np.sort(np.concatenate([S, pconfigs_rot(S, L)]))
        self.R, self.per = reps_and_periods(S, L)

    def H(self, m, K):
        L = self.L
        # keep only reps compatible with K: K * period = 0 mod 2pi
        kk = int(round(K * L / (2 * np.pi))) % L
        ok = ((kk * self.per) % L) == 0
        R = self.R[ok]; per = self.per[ok]
        maxnnz = R.shape[0] * (1 + 8 * min(self.Nmax, L))
        rows, cols, vr, vi, N = build_coo(R, per, L, 2 * np.pi * kk / L, m['w'].real, m['w'].imag,
                                          m['V2'], m['V3'], m['mu'], m['r'], m.get('lam', 0.0), maxnnz, self.Nmax)
        n = R.shape[0]
        H = sp.csr_matrix((vr + 1j * vi, (rows, cols)), shape=(n, n))
        return H, N

    def lowest(self, m, K, k=2, vecs=False):
        H, N = self.H(m, K)
        n = H.shape[0]
        if n <= 600:
            e, V = np.linalg.eigh(H.toarray())
            e, V = e[:k], V[:, :k]
        else:
            e, V = spla.eigsh(H, k=k, which='SA', tol=1e-10)
            o = np.argsort(e); e, V = e[o], V[:, o]
        Nexp = (np.abs(V) ** 2 * N[:, None]).sum(0)
        return (e, Nexp, V, N) if vecs else (e, Nexp)


def model(u, phi, v, mu, g=2/3, r=None, lam=0.0):
    return dict(w=u * np.exp(1j * phi), V2=-2.0 + v, V3=-g, mu=mu, r=np.sqrt(v) if r is None else r, lam=lam)


if __name__ == "__main__":
    # self-test: hermiticity, and the single-particle band at r = 0
    L = 14
    ring = Ring(L)
    m = model(1.0, np.pi / 2, 0.3, -2.1)
    for kk in range(L):
        H, N = ring.H(m, 2 * np.pi * kk / L)
        assert abs(H - H.getH()).max() < 1e-12, kk
    m0 = model(1.0, 1.3, 0.0, 0.0, r=0.0)
    H, N = ring.H(m0, 2 * np.pi * 3 / L)
    e = np.linalg.eigvalsh(H.toarray()[np.ix_(N == 1, N == 1)])
    k = 2 * np.pi * 3 / L
    print("hermitian ok; one-particle level", e, "eps =", 2 * np.cos(k - 1.3) - 2 * np.cos(2 * k))
