"""Thermodynamic helpers: lattice entropy from spectra, SU(N)_1 torus entropy, free fermions,
and extraction of the conformal window T*_eps.

Conventions: entropy per site s(T); CFT line s_CFT(T) = pi c T/(3 v); the window edge
T*_eps = sup{T : |s(T')/s_ref(T') - 1| <= eps for all T' in [T_lo, T]}, where s_ref is either the
thermodynamic CFT line or the finite-ring CFT entropy, and T_lo is the lowest temperature used.
"""
import glob
import itertools
import numpy as np
from uls_schur_weyl import gl_mult


def entropy_from_levels(E, g, Ts):
    """Entropy (total, not per site) of levels E with degeneracies g at temperatures Ts."""
    E = np.asarray(E, float); g = np.asarray(g, float)
    keep = g > 0
    E, g = E[keep], g[keep]
    E0 = E.min()
    out = []
    for T in Ts:
        x = -(E - E0) / T
        xm = x.max()
        w = g * np.exp(x - xm)
        Z = w.sum()
        Eavg = (w * (E - E0)).sum() / Z
        out.append(np.log(Z) + xm + Eavg / T)
    return np.array(out)


def load_uls(L, path="data"):
    """Merge the per-partition spectra of the ULS ring of length L."""
    files = sorted(glob.glob(f"{path}/uls_L{L}.npz")) or sorted(glob.glob(f"{path}/uls_L{L}_r*_part*.npz"))
    spec = {}
    for f in files:
        z = np.load(f)
        for k in z.files:
            if k.startswith("E_"):
                key = k[2:]
                lam = tuple(int(x) for x in key.split("_"))
                spec[lam] = (z["E_" + key], z["w_" + key], int(z["exact_" + key]))
    return spec


def uls_levels(spec, d):
    E, g = [], []
    for lam, (e, w, _) in spec.items():
        m = gl_mult(lam, d)
        if m > 0:
            E.append(e); g.append(w * m)
    return np.concatenate(E), np.concatenate(g)


def maxrows_available(spec):
    return max(len(l) for l in spec)


# ---------- SU(N)_1 WZW on a torus (diagonal invariant, ring length L = 0 mod N) ----------

def sun1_levels(N, Dmax):
    """Scaling dimensions Delta = h + hbar <= Dmax with multiplicities, for sum_j |chi_j|^2."""
    K = int(np.ceil(np.sqrt(2 * Dmax))) + 1
    nmax = int(Dmax) + 1
    p = np.zeros(nmax + 1); p[0] = 1.0  # coefficients of prod (1-q^n)^{-(N-1)}
    for _ in range(N - 1):
        for n in range(1, nmax + 1):
            for m in range(n, nmax + 1):
                p[m] += p[m - n]
    grid = np.array(np.meshgrid(*([np.arange(-K, K + 1)] * (N - 1)), indexing="ij")).reshape(N - 1, -1).T
    chis = []
    for j in range(N):
        last = j - grid.sum(1)  # fix the N-th entry so that sum(y) = j
        ok = np.abs(last) <= K + j
        y2 = (grid[ok] ** 2).sum(1) + last[ok] ** 2
        h = 0.5 * (y2 - j * j / N)
        h = h[h <= Dmax + 1e-9]
        hv, hm = np.unique(np.round(h, 9), return_counts=True)
        lev = {}
        for h0, m0 in zip(hv, hm):
            for n in range(0, nmax + 1):
                hh = round(h0 + n, 9)
                if hh <= Dmax:
                    lev[hh] = lev.get(hh, 0) + m0 * p[n]
        chis.append(lev)
    D = {}
    for lev in chis:
        for h1, m1 in lev.items():
            for h2, m2 in lev.items():
                Dl = round(h1 + h2, 9)
                if Dl <= Dmax:
                    D[Dl] = D.get(Dl, 0) + m1 * m2
    Ds = np.array(sorted(D)); ms = np.array([D[x] for x in Ds])
    return Ds, ms


_LEVELS = {}


def sun1_torus_entropy(N, L, v, Ts, Dmax=14.0):
    """Entropy per site of SU(N)_1 (diagonal invariant) on a ring of length L, velocity v.
    Uses Z(x) = Z(1/x), x = v/(L T), and the channel with x >= 1, where Delta <= Dmax converges."""
    key = (N, Dmax)
    if key not in _LEVELS:
        _LEVELS[key] = sun1_levels(N, Dmax)
    Ds, ms = _LEVELS[key]
    c = N - 1
    out = []
    for T in np.atleast_1d(Ts):
        x = v / (L * T)
        y = x if x >= 1 else 1 / x
        a = -2 * np.pi * y * (Ds - c / 12.0)
        am = a.max()
        w = ms * np.exp(a - am)
        G = w.sum()
        mean = (w * (Ds - c / 12.0)).sum() / G
        # S = log Z - x d(log Z)/dx; direct channel: +2 pi y <Delta - c/12>, crossed channel: -2 pi y <...>
        out.append(np.log(G) + am + (1 if x >= 1 else -1) * 2 * np.pi * y * mean)
    return np.array(out) / L, Dmax


# ---------- free fermions in the thermodynamic limit ----------

def ff_entropy(eps_k, Ts):
    """Entropy per site of free fermions with single-particle energies eps_k on a uniform k grid."""
    out = []
    for T in Ts:
        x = np.abs(eps_k) / T
        s = np.log1p(np.exp(-x)) + x / (np.exp(x) + 1.0)
        out.append(s.mean())
    return np.array(out)


def window(Ts, ratio, eps):
    """Largest T in the grid with |ratio - 1| <= eps on all grid points up to T (None if the first fails)."""
    ok = np.abs(np.asarray(ratio) - 1.0) <= eps
    if not ok[0]:
        return None
    bad = np.nonzero(~ok)[0]
    return Ts[-1] if len(bad) == 0 else Ts[bad[0] - 1]
