"""C164 Task 2: quantum one-particle band above a uniform gapped ground state, from momentum-resolved ring ED
(est.py tables, exact to 1e-14). Claude Opus agent, 2026-11-08.

band(h, L, nev): lowest nev levels in each momentum sector k = 2 pi m / L (est.py convention), the ground state
(E0, m0) and the one-particle band E1(k) = lowest level above E0 in each sector (second level at k = m0).
interp(E1): trigonometric interpolation of the L samples (exact for an isolated band up to aliasing
e^{-L/(2 xi)}), evaluated on a fine grid; local minima of the interpolant.
two_particle_floor(Ek, k): bottom of the two-particle continuum min_q [E(q) + E(k - q)] from the interpolant.
"""
import numpy as np
import scipy.linalg as sla
import scipy.sparse.linalg as spla
import est


def lowest(Hk, nev):
    n = Hk.shape[0]
    if n <= 600:
        return sla.eigh(Hk.toarray(), subset_by_index=[0, nev - 1], driver='evr')
    w, V = spla.eigsh(Hk, k=nev, which='SA', tol=1e-10, ncv=max(20, 2 * nev + 8))
    o = np.argsort(w)
    return w[o], V[:, o]


def band(h, L, nev=3):
    T = est.table(L)
    lev = np.full((L, nev), np.inf)
    vecs = {}
    for m in range(L):
        Hk = T.sector_matrix(h, m)
        ne = min(nev, Hk.shape[0])
        w, V = lowest(Hk, ne)
        lev[m, :ne] = w
    E0 = lev[:, 0].min(); m0 = int(np.argmin(lev[:, 0]))
    E1 = lev[:, 0].copy(); E1[m0] = lev[m0, 1]
    return dict(L=L, lev=lev, E0=E0, m0=m0, E1=E1 - E0, gap0=float(np.sort(lev.ravel())[1] - E0))


def interp(E1, ngrid=2000):
    L = len(E1)
    c = np.fft.fft(E1) / L
    r = np.fft.fftfreq(L, 1.0 / L)                    # integers
    ks = np.linspace(0, 2 * np.pi, ngrid, endpoint=False)
    # symmetric treatment of the Nyquist harmonic
    Ek = np.zeros(ngrid)
    for ci, ri in zip(c, r):
        if L % 2 == 0 and abs(ri) == L // 2:
            Ek += (ci * np.cos(ri * ks)).real          # cos part only
        else:
            Ek += (ci * np.exp(1j * ri * ks)).real
    return ks, Ek


def local_minima(ks, Ek):
    i = np.arange(len(Ek))
    m = (Ek < np.roll(Ek, 1)) & (Ek <= np.roll(Ek, -1))
    return [(float(ks[j]), float(Ek[j])) for j in i[m]]


def two_particle_floor(ks, Ek):
    n = len(ks)
    # E2(K) = min_q E(q) + E(K - q) on the same grid (cyclic convolution in min-plus algebra)
    E2 = np.full(n, np.inf)
    for s in range(n):
        E2[s] = np.min(Ek + np.roll(Ek[::-1], s + 1))
    return E2


def analyze(h, L=16, nev=3):
    b = band(h, L, nev)
    ks, Ek = interp(b['E1'])
    mins = local_minima(ks, Ek)
    mins.sort(key=lambda t: t[1])
    out = dict(L=L, m0=b['m0'], E1=b['E1'].tolist(), gap=float(b['E1'].min()), mins=mins,
               lev=(b['lev'] - b['E0']).tolist())
    # is the band isolated at its minima? second level in the same sector vs. the lowest
    return out
