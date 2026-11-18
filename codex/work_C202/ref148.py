"""C202 referee code for C148 (easy-axis U(1) class). Written from scratch, no C148 code reused.
ED (sector, arbitrary hopping phase), circle Bethe solver, restricted continuum by Gauss-Legendre Nystrom."""
import numpy as np
from math import comb, pi
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.optimize import brentq

# ---------------------------------------------------------------- ED
def sector_states(L, n):
    """bitmasks with n set bits (set bit = down spin), sorted"""
    from itertools import combinations
    if n == 0:
        return np.array([0], dtype=np.int64)
    arr = np.fromiter((sum(1 << i for i in c) for c in combinations(range(L), n)), dtype=np.int64, count=comb(L, n))
    arr.sort()
    return arr

def sector_energy(L, n, Jz, hop=-1.0 + 0j):
    """ground energy of sum over oriented bonds x -> x+1 (both oriented bonds at L = 2) of
    (1/2)(hop S^+_x S^-_{x+1} + conj(hop) S^-_x S^+_{x+1}) + Jz S^z_x S^z_{x+1}, n down spins.
    hop = J_perp + i D. hop = -1 is the stoquastic gauge (flux pi per bond)."""
    st = sector_states(L, n)
    dim = len(st)
    diag = np.zeros(dim)
    rows, cols, vals = [], [], []
    for x in range(L):
        y = (x + 1) % L
        bx = (st >> x) & 1
        by = (st >> y) & 1
        diag += Jz * np.where(bx == by, 0.25, -0.25)
        # S^+_x S^-_y: x down -> up (bit 1 -> 0), y up -> down (bit 0 -> 1)
        sel = (bx == 1) & (by == 0)
        src = st[sel]
        dst = src ^ ((1 << x) | (1 << y))
        idx = np.searchsorted(st, dst)
        rows.append(idx); cols.append(np.nonzero(sel)[0]); vals.append(np.full(len(idx), 0.5 * hop))
        rows.append(np.nonzero(sel)[0]); cols.append(idx); vals.append(np.full(len(idx), 0.5 * np.conj(hop)))
    cplx = abs(np.imag(hop)) > 1e-15
    r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
    if not cplx:
        v = v.real
    H = sp.coo_matrix((v, (r, c)), shape=(dim, dim)).tocsr() + sp.diags(diag)
    if dim <= 600:
        return float(np.linalg.eigvalsh(H.toarray())[0])
    return float(eigsh(H, k=1, which='SA', tol=1e-13, ncv=min(dim - 1, 40))[0][0])

def full_energy(L, Jperp, D, Jz, B):
    """E_L of the ring with bond h = Jperp(SxSx+SySy) + D(SxSy-SySx) + Jz SzSz - (B/2)(Sz1 + 1Sz).
    S^x S^x + S^y S^y = (S^+S^- + S^-S^+)/2, S^x S^y - S^y S^x = (i/2)(S^+ S^- - S^- S^+),
    so the coefficient of S^+_x S^-_{x+1} is (Jperp + iD)/2. Field: -B sum_x S^z_x = -B (L/2 - n)."""
    best = None
    for n in range(L + 1):
        e = sector_energy(L, n, Jz, hop=Jperp + 1j * D) - B * (L / 2 - n)
        best = e if best is None else min(best, e)
    return best

def half_spread(J, Jz, B):
    ev = [Jz / 4 + B / 2, Jz / 4 - B / 2, -Jz / 4 + J / 2, -Jz / 4 - J / 2]
    return (max(ev) - min(ev)) / 2

# ---------------------------------------------------------------- circle kernels
def a1(x, lam):
    return np.sinh(lam) / (pi * (np.cosh(lam) - np.cos(2 * x)))

def Kk(x, lam):
    return np.sinh(2 * lam) / (pi * (np.cosh(2 * lam) - np.cos(2 * x)))

def ebare(x, lam):
    return -np.sinh(lam) ** 2 / (np.cosh(lam) - np.cos(2 * x))

def r0(x, lam):
    """(I + K*)^{-1} a1 on the circle, by its Fourier series"""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    kmax = int(60 / lam) + 30
    k = np.arange(1, kmax + 1)
    return (0.5 + (np.cos(2 * np.outer(x, k)) / np.cosh(k * lam)).sum(1)) / pi

def Jhole(x, lam):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    kmax = int(40 / lam) + 30
    k = np.arange(1, kmax + 1)
    t = np.exp(-2 * k * lam)
    return (0.5 + 2 * (np.cos(2 * np.outer(x, k)) * (t / (1 + t))).sum(1)) / pi

# ---------------------------------------------------------------- Bethe
def _G(x, L, lam, I):
    d = x[:, None] - x[None, :]
    th = 2 * np.arctan(np.tan(d) / np.tanh(lam)) + 2 * pi * np.round(d / pi)
    return L * 2 * np.arctan(np.tan(x) / np.tanh(lam / 2)) - th.sum(1) - 2 * pi * I

def _Jac(x, L, lam):
    d = x[:, None] - x[None, :]
    Km = Kk(d, lam)
    return 2 * pi * (np.diag(L * a1(x, lam) - Km.sum(1)) + Km)

def newton(x, L, lam, I, tol=None, itmax=200):
    if tol is None:
        tol = max(1e-12, 2e-14 * L)   # rounding floor of L p(x) grows with L
    g = _G(x, L, lam, I); ng = np.abs(g).max()
    for it in range(itmax):
        if ng < tol:
            break
        dx = np.linalg.solve(_Jac(x, L, lam), -g)
        s = 1.0
        while True:
            xn = x + s * dx
            ok = np.all(np.diff(xn) > 0) and xn[-1] - xn[0] < pi and abs(xn).max() < pi / 2
            if ok:
                gn = _G(xn, L, lam, I); nn = np.abs(gn).max()
                if nn < ng:
                    break
                if ng < 1e-9 and s < 1e-3:   # stagnation at rounding level
                    return x, ng
                if s < 1e-6:
                    break
            s /= 2
            if s < 1e-12:
                raise RuntimeError('line search failed')
        x, g, ng = xn, gn, nn
    return x, ng

def bethe_roots(L, n, lam, lam_start=6.0, ratio=0.8):
    """ground-state roots, continued from the Ising-side guess x_j = pi I_j/(L - n) down to lam"""
    I = np.arange(n) - (n - 1) / 2
    x = pi * I / (L - n)
    l = max(lam, lam_start)
    x, res = newton(x, L, l, I)
    r = ratio
    while l > lam * (1 + 1e-13):
        lt = max(lam, l * r)
        try:
            # Euler predictor along the branch: dx/dlam = -Jac^{-1} dG/dlam
            hh = 1e-6 * l
            dG = (_G(x, L, l + hh, I) - _G(x, L, l - hh, I)) / (2 * hh)
            guess = x + np.linalg.solve(_Jac(x, L, l), -dG) * (lt - l)
            if not (np.all(np.diff(guess) > 0) and abs(guess).max() < pi / 2):
                guess = x.copy()
            xt, res = newton(guess, L, lt, I)
            if res > 1e-9:
                raise RuntimeError('residual')
        except RuntimeError:
            r = np.sqrt(r)          # adaptive: shorten the continuation step
            if 1 - r < 1e-7:
                raise RuntimeError(f'no convergence L={L} n={n} lam={lt}')
            continue
        x, l = xt, lt
        r = max(ratio, r ** 1.5)
    return x

def bethe_energy(L, n, lam, x=None):
    if n == 0:
        return np.cosh(lam) * L / 4, np.array([])
    if x is None:
        x = bethe_roots(L, n, lam)
    # pairwise symmetric summation for accuracy
    return np.cosh(lam) * L / 4 + float(np.sum(np.sort(ebare(x, lam)))), x

def sigma_L(xg, x, L, lam):
    xg = np.atleast_1d(xg)
    return a1(xg, lam) - Kk(xg[:, None] - x[None, :], lam).sum(1) / L

def z_L(xg, x, L, lam):
    xg = np.atleast_1d(xg)
    d = xg[:, None] - x[None, :]
    th = 2 * np.arctan(np.tan(d) / np.tanh(lam)) + 2 * pi * np.round(d / pi)
    return (2 * np.arctan(np.tan(xg) / np.tanh(lam / 2)) - th.sum(1) / L) / (2 * pi)

# ---------------------------------------------------------------- continuum
_gl_cache = {}
def gl(N):
    if N not in _gl_cache:
        _gl_cache[N] = np.polynomial.legendre.leggauss(N)
    return _gl_cache[N]

def continuum(Q, lam, N=None):
    """restricted continuum on [-Q, Q]: returns m(Q), E(Q) (per site), r_Q(Q)"""
    if N is None:
        N = int(min(5000, 40 * Q / lam + 200))
    t, w = gl(N)
    y = Q * t; wy = Q * w
    A = np.eye(N) + Kk(y[:, None] - y[None, :], lam) * wy[None, :]
    r = np.linalg.solve(A, a1(y, lam))
    m = float(wy @ r)
    E = np.cosh(lam) / 4 + float(wy @ (ebare(y, lam) * r))
    return m, E

def ebulk_half(lam):
    k = np.arange(1, int(60 / lam) + 30)
    return np.cosh(lam) / 4 - np.sinh(lam) * (0.5 + 2 * np.sum(1 / (np.exp(2 * k * lam) + 1)))

def ebulk(nu, lam):
    if nu <= 0:
        return np.cosh(lam) / 4
    mh, Eh = continuum(pi / 2, lam)
    if nu >= mh - 1e-13:
        return Eh
    Q = brentq(lambda q: continuum(q, lam)[0] - nu, 1e-9, pi / 2, xtol=1e-15, rtol=1e-15)
    return continuum(Q, lam)[1]

def M0(lam):
    return (1 + np.cosh(lam)) / 4
