"""Delta = cosh(lam) > 1: Bethe roots on the circle x in (-pi/2, pi/2], J=1.
Kernels a1, K=a2 (period pi), r0 = periodized sech, hole kernel J = a1*r0.
Restricted continuum on [-Q,Q] (Nystrom), finite-L diagnostics. Diagnostics only, floating point."""
import numpy as np
from scipy.optimize import brentq
import bethe_scan as bs
np.seterr(over="ignore", invalid="ignore")


class Circ:
    def __init__(s, delta):
        s.d = delta; s.l = l = float(np.arccosh(delta))
        s.a = np.pi*np.sinh(l); s.b = np.pi/l; s.M0 = (1+delta)/4
    def a1(s, x): return np.sinh(s.l)/(np.pi*(np.cosh(s.l)-np.cos(2*x)))
    def a1p(s, x): return -2*np.sinh(s.l)*np.sin(2*x)/(np.pi*(np.cosh(s.l)-np.cos(2*x))**2)
    def K(s, x): return np.sinh(2*s.l)/(np.pi*(np.cosh(2*s.l)-np.cos(2*x)))
    def Kp(s, x): return -2*np.sinh(2*s.l)*np.sin(2*x)/(np.pi*(np.cosh(2*s.l)-np.cos(2*x))**2)
    def e(s, x): return -s.a*s.a1(x)
    def ep(s, x): return -s.a*s.a1p(x)
    def r0(s, x):
        x = np.asarray(x, float)
        if s.l < 1.5:
            m = np.arange(-6, 7)
            return (1/(2*s.l*np.cosh(s.b*(x[..., None]+np.pi*m)))).sum(-1)
        k = np.arange(1, int(min(199, 700/s.l))+1)
        return (0.5+(np.cos(2*k*x[..., None])/np.cosh(k*s.l)).sum(-1))/np.pi
    def Z0(s, x):
        """int_0^x r0"""
        x = np.asarray(x, float)
        if s.l < 1.5:
            m = np.arange(-6, 7)
            A = lambda y: np.arctan(np.tanh(s.b*y/2))/np.pi
            return (A(x[..., None]+np.pi*m)-A(np.pi*m)).sum(-1)
        k = np.arange(1, int(min(199, 700/s.l))+1)
        return (x/2+(np.sin(2*k*x[..., None])/(2*k*np.cosh(k*s.l))).sum(-1))/np.pi
    def Jfour(s, k):
        t = np.exp(-2*abs(k)*s.l); return t/(1+t)
    def J(s, x):
        x = np.asarray(x, float)
        k = np.arange(1, int(60/s.l)+50)
        return (0.5+2*(np.cos(2*k*x[..., None])*s.Jfour(k)).sum(-1))/np.pi


def nodes(Q, l, N=8):
    P = int(min(3000, max(6, np.ceil(2*Q/(min(l, 1.0)/3)))))
    t, w = np.polynomial.legendre.leggauss(N)
    ed = np.linspace(-Q, Q, P+1)
    X = ((ed[:-1, None]+ed[1:, None])/2+(ed[1:, None]-ed[:-1, None])/2*t[None, :]).ravel()
    W = ((ed[1:, None]-ed[:-1, None])/2*w[None, :]).ravel()
    return X, W


class Restricted:
    def __init__(s, M, Q, N=8):
        s.M, s.Q = M, Q
        s.X, s.W = nodes(Q, M.l, N)
        A = np.eye(len(s.X))+M.K(s.X[:, None]-s.X[None, :])*s.W[None, :]
        rhs = np.stack([M.a1(s.X), np.ones_like(s.X), M.ep(s.X)], 1)
        sol = np.linalg.solve(A, rhs)
        s.r, s.Z, s.epsp = sol[:, 0], sol[:, 1], sol[:, 2]
        s.m = float(s.W@s.r)
        s.E = M.d/4+float(s.W@(M.e(s.X)*s.r))
    def interp(s, f, g, x):
        x = np.atleast_1d(np.asarray(x, float))
        return g(x)-(s.M.K(x[:, None]-s.X[None, :])*s.W[None, :])@f
    def r_at(s, x): return s.interp(s.r, s.M.a1, x)
    def cum_r(s, x):
        """int_{-Q}^x r, by Nystrom interpolation on a fine grid"""
        g = np.linspace(-s.Q, s.Q, 20001); rg = s.r_at(g)
        c = np.concatenate([[0], np.cumsum((rg[1:]+rg[:-1])/2*np.diff(g))])
        return np.interp(x, g, c)


def Q_of_nu(M, nu):
    if nu >= 0.5: return np.pi/2
    f = lambda Q: Restricted(M, Q).m-nu
    if f(np.pi/2-1e-12) <= 0: return np.pi/2
    return brentq(f, 1e-9, np.pi/2-1e-12, xtol=1e-14, rtol=1e-14)


def e_bulk(M, nu):
    if nu <= 0: return M.d/4
    return Restricted(M, Q_of_nu(M, nu)).E


def e_half(M):
    k = np.arange(1, int(min(1999, 300/M.l))+1)
    return M.d/4-np.sinh(M.l)*(0.5+2*np.sum(1/(np.exp(2*k*M.l)+1)))


def seed_fn(M):
    def seed(I, L, n):
        g = np.linspace(0, np.pi/2, 40001); z = M.Z0(g)
        t = np.abs(I)/L
        x = np.interp(np.minimum(t, z[-1]*(1-1e-9)), z, g)
        return np.sign(I)*x
    return seed


def solve(L, n, delta, **kw):
    """Bethe solve with an r0-quantile seed (patched into bethe_scan)."""
    M = Circ(delta)
    orig = bs.functions
    def patched(d):
        out = list(orig(d)); out[5] = seed_fn(M); return tuple(out)
    bs.functions = patched
    try:
        return bs.solve(L, n, delta, **kw)
    finally:
        bs.functions = orig


def finite_diag(L, n, delta, grid=40001):
    """sigma_L on the circle, counting function, Q_L, sawtooth, G, density error."""
    M = Circ(delta)
    E, dg, x = solve(L, n, delta, return_roots=True)
    p, th, dp, dt, bare, _ = bs.functions(delta)
    nu = n/L
    X = np.linspace(-np.pi/2, np.pi/2, grid)
    sig = M.a1(X)-M.K(X[:, None]-x[None, :]).sum(1)/L
    z = p(X)/(2*np.pi)-th(X[:, None]-x[None, :]).sum(1)/(2*np.pi*L)
    r0 = M.r0(X)
    xn = x[-1]
    inn = (X >= x[0]) & (X <= xn)
    out = ~inn
    res = dict(L=L, n=n, delta=delta, E=float(E),
               min_sig_over_r0_in=float(np.min(sig[inn]/r0[inn])),
               min_sig_over_r0_all=float(np.min(sig/r0)),
               L_r0_xn=float(L*M.r0(np.array([xn]))[0]), xn=float(xn),
               negmass_out=float(np.sum(np.maximum(-sig[out], 0))*(X[1]-X[0])))
    # Q_L: first x >= x_n with z = nu/2
    j0 = np.searchsorted(X, xn)
    idx = np.where(z[j0:] >= nu/2-1e-15)[0]
    QL = np.pi/2 if len(idx) == 0 else float(X[j0+idx[0]])
    res['Q_L'] = QL
    if n >= 1:
        sel = (X >= -QL) & (X <= QL)
        Cmu = np.searchsorted(x, X[sel], side='right')/L-nu/2
        cumR = Cmu-z[sel]
        res['L_sup_cumR'] = float(L*np.max(np.abs(cumR)))
        R = Restricted(M, QL)
        rr = R.r_at(X[sel])
        G = Cmu+nu/2-R.cum_r(X[sel])
        res['L_sup_G'] = float(L*np.max(np.abs(G)))
        res['L_sup_sig_minus_r'] = float(L*np.max(np.abs(sig[sel]-rr)))
        res['L_nu_minus_m'] = float(L*abs(nu-R.m))
        eb = e_bulk(M, nu)
        res['L_F_minus_Le_over_M0'] = float(L*(E-L*eb)/M.M0)
    return res
