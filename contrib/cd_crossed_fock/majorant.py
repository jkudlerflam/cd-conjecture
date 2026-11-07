"""Bosonic majorant and the XX exact crossed spectrum (m = infinity).

Majorant:  log Maj(L) = K + S F(L/(V' beta)),  F(x) = -sum_{k>=1} log(1 - exp(-2 pi k x)).
For a crossed log-partition function g = log Z~(L) the constraint g <= K + S F(x) with x = L/(V' beta)
gives V' >= L/(beta F^{-1}((g-K)/S)) whenever g > K.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def F(x):
    """-sum log(1-q^k), q = exp(-2 pi x); uses the modular transform for small x."""
    x = np.asarray(x, float)
    out = np.empty_like(x)
    small = x < 1.0
    xs = x[small]
    if xs.size:
        # F(x) = pi/(12x) + (1/2) log x - pi x/12 + F(1/x)
        out[small] = np.pi / (12 * xs) + 0.5 * np.log(xs) - np.pi * xs / 12 + _Fdirect(1 / xs)
    xl = x[~small]
    if xl.size:
        out[~small] = _Fdirect(xl)
    return out


def _Fdirect(x):
    k = np.arange(1, 60)[:, None]
    return -np.log1p(-np.exp(-2 * np.pi * k * x[None, :])).sum(0)


_xg = np.geomspace(1e-6, 50, 20000)
_Fg = F(_xg)


def Finv(y):
    """inverse of the decreasing function F on (0, inf); y > 0."""
    y = np.asarray(y, float)
    return np.interp(-y, -_Fg, _xg, left=1e-6, right=50.)


def vmin_needed(g, L, beta, S, K):
    """smallest V' such that g <= K + S F(L/(V' beta)); 0 if g <= K."""
    g = np.asarray(g, float)
    out = np.zeros_like(g)
    mask = g > K
    x = Finv((g[mask] - K) / S)
    out[mask] = (np.broadcast_to(L, g.shape)[mask] / (np.broadcast_to(beta, g.shape)[mask] * x))
    return out


def K_needed(g, L, beta, S, V):
    return np.max(np.asarray(g) - S * F(np.asarray(L, float) / (V * np.asarray(beta, float))))


# ---------------- XX exact crossed spectrum (C125 PROOF section 3) -----------------
# h = SxSx + SySy, M = 1/2, v = 1, A = 1.

def xx_delta_beta(beta):
    """delta_beta = p_c - p_s = (1/2pi) int log coth(beta|cos k|/2) dk."""
    if beta < 80:
        val = 2 / np.pi * quad(lambda k: -np.log(np.tanh(beta * np.cos(k) / 2)), 0, np.pi / 2,
                               epsabs=1e-13, limit=200)[0]
    else:
        val = 2 / (np.pi * beta) * quad(lambda y: -np.log(np.tanh(y / 2)) / np.sqrt(1 - (y / beta) ** 2),
                                        0, min(40, beta * 0.999999), epsabs=1e-13, limit=200)[0]
    return val


def xx_logZt(L, beta, nmax=None):
    """exact log Z~(L) = log Z_L + L eps_0 for even (or real > 1) L, from formula (6) of C125."""
    Lmin = float(np.min(L))
    if nmax is None:
        nmax = int(max(2000, 40 * beta)) if Lmin > 3 else int(max(2e6, 400 * beta))
    n = np.arange(nmax)
    ko = np.arcsinh(np.pi * (2 * n + 1) / beta)
    ke = np.arcsinh(2 * np.pi * (n + 1) / beta)
    L = np.atleast_1d(np.asarray(L, float))
    out = []
    db = xx_delta_beta(beta)
    for l in L:
        a = 4 * np.log1p(np.exp(-l * ko)).sum()
        bm = 4 * np.log1p(-np.exp(-l * ko)).sum()
        first = np.logaddexp(a, bm) - np.log(2)
        second = np.log(2) - l * db + 4 * np.log1p(np.exp(-l * ke)).sum()
        out.append(np.logaddexp(first, second))
    return np.array(out)
