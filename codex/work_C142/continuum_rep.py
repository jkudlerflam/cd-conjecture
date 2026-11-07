"""Restricted continuum dressing for 0<Delta<1 (positive kernel). J=1.
r_Q=(I+K_Q)^{-1}a1, Z_Q=(I+K_Q)^{-1}1, eps_Q'=(I+K_Q)^{-1}e' (valid since eps_Q(+-Q)=0).
Nystrom with Gauss-Legendre nodes on [-Q,Q], graded toward the center for narrow kernels.
Diagnostics only, floating point."""
import numpy as np
from scipy.optimize import brentq

class Model:
    def __init__(s, delta):
        s.d = delta; s.g = g = np.arccos(delta)
        s.a = np.pi*np.sin(g); s.b = np.pi/g
        s.q = 1-2*g/np.pi; s.B = np.pi/(2*g); s.qJ = (np.pi-2*g)/(2*(np.pi-g))
    def a1(s, x): return np.sin(s.g)/(np.pi*(np.cosh(2*x)-np.cos(s.g)))
    def a1p(s, x): return -2*np.sin(s.g)*np.sinh(2*x)/(np.pi*(np.cosh(2*x)-np.cos(s.g))**2)
    def K(s, x): return np.sin(2*s.g)/(2*np.pi*(np.sinh(x)**2+np.sin(s.g)**2))
    def Kp(s, x): return -np.sin(2*s.g)*np.sinh(2*x)/(2*np.pi*(np.sinh(x)**2+np.sin(s.g)**2)**2)
    def r0(s, x): return 1/(2*s.g*np.cosh(s.b*x))
    def e(s, x): return -s.a*s.a1(x)
    def ep(s, x): return -s.a*s.a1p(x)

def nodes(Q, g, N):
    # composite Gauss-Legendre with panel width ~ g/4 near kernel scale
    P = max(4, int(np.ceil(2*Q/(g/3))))
    P = min(P, 4000)
    t, w = np.polynomial.legendre.leggauss(N)
    edges = np.linspace(-Q, Q, P+1)
    X = ((edges[:-1, None]+edges[1:, None])/2 + (edges[1:, None]-edges[:-1, None])/2*t[None, :]).ravel()
    W = ((edges[1:, None]-edges[:-1, None])/2*w[None, :]).ravel()
    return X, W

class Restricted:
    def __init__(s, M, Q, N=8):
        s.M, s.Q = M, Q
        s.X, s.W = nodes(Q, M.g, N)
        A = np.eye(len(s.X)) + M.K(s.X[:, None]-s.X[None, :])*s.W[None, :]
        rhs = np.stack([M.a1(s.X), np.ones_like(s.X), M.ep(s.X)], 1)
        sol = np.linalg.solve(A, rhs)
        s.r, s.Z, s.epsp = sol[:, 0], sol[:, 1], sol[:, 2]
        s.m = float(s.W@s.r)
        s.E = M.d/4 + float(s.W@(M.e(s.X)*s.r))
    def interp(s, f, g, x):
        """f = g - K_Q f evaluated at arbitrary x (Nystrom interpolation)."""
        return g(x) - (s.M.K(np.asarray(x)[:, None]-s.X[None, :])*s.W[None, :])@f
    def edge(s):
        Q = np.array([s.Q])
        rQ = s.interp(s.r, s.M.a1, Q)[0]
        ZQ = s.interp(s.Z, lambda x: np.ones_like(x), Q)[0]
        eQ = s.interp(s.epsp, s.M.ep, Q)[0]
        return rQ, ZQ, eQ
    def s_full(s, x):
        """counting density s_Q = a1 - K*(r 1_Q) on the whole line"""
        return s.interp(s.r, s.M.a1, x)

def Q_of_nu(M, nu, N=8):
    if nu >= 0.5: return np.inf
    f = lambda Q: Restricted(M, Q, N).m - nu
    hi = 1.0
    while f(hi) < 0: hi *= 1.6
    return brentq(f, 1e-9, hi, xtol=1e-13, rtol=1e-13)

def e_bulk(M, nu, N=8):
    if nu <= 0: return M.d/4
    if nu >= 0.5:
        # closed form at half filling: Delta/4 + int e r0
        from scipy.integrate import quad
        v, _ = quad(lambda x: 2*M.e(x)*M.r0(x), 0, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)
        return M.d/4 + v
    Q = Q_of_nu(M, nu, N)
    return Restricted(M, Q, N).E
