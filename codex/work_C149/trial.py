"""C149 trial configurations: quantiles of a continuum counting function.
Z(u) = int_0^u rho for an even piecewise-linear rho on a mesh of [0, Qmax]; nodes Z(y_j) = t_j.
Floating point only."""
import numpy as np
from bethe_ferro import Model, residual
from continuum import Kernel, Restricted, grading, bulk_at

class PL:
    """even piecewise-linear density on mesh y (y[0]=0), cumulative Z (odd)."""
    def __init__(self, y, rho):
        self.y = np.asarray(y, float); self.rho = np.asarray(rho, float)
        h = np.diff(self.y)
        self.cum = np.concatenate([[0.0], np.cumsum(h*(self.rho[:-1]+self.rho[1:])/2)])
    def Z(self, u):
        s = np.sign(u); a = np.abs(u)
        i = np.clip(np.searchsorted(self.y, a)-1, 0, len(self.y)-2)
        y0, h = self.y[i], self.y[i+1]-self.y[i]
        r0, r1 = self.rho[i], self.rho[i+1]
        d = a-y0
        return s*(self.cum[i]+r0*d+(r1-r0)*d*d/(2*h))
    def rho_at(self, u):
        return np.interp(np.abs(u), self.y, self.rho)
    def inv(self, t):
        s = np.sign(t); a = np.abs(t)
        i = np.clip(np.searchsorted(self.cum, a)-1, 0, len(self.y)-2)
        y0, h = self.y[i], self.y[i+1]-self.y[i]
        r0, r1 = self.rho[i], self.rho[i+1]
        c = a-self.cum[i]; k = (r1-r0)/(2*h)
        # k d^2 + r0 d - c = 0
        disc = np.sqrt(np.maximum(r0*r0+4*k*c, 0))
        d = np.where(np.abs(k*c) < 1e-14*r0*r0+1e-300, c/r0, 2*c/(r0+disc))
        return s*(y0+d)

def quantiles(Zinv, L, n, shift=0.0):
    I = np.arange(1, n+1)-(n+1)/2
    return Zinv(I/L+shift)

def continuum_trial(eta, L, n, N=900):
    ker = Kernel(eta)
    nu = n/L
    R = solve_Q(ker, nu, N=N)
    return R, PL(R.y, R.r)

def solve_Q(ker, nu, N=700, tol=1e-13):
    """Q with m(Q) = nu by Brent's method in log Q (m is increasing, C136 (2.3))."""
    from scipy.optimize import brentq
    cache = {}
    def f(lq):
        R = Restricted(ker, float(np.exp(lq)), N=N)
        cache[lq] = R
        return R.m-nu
    lq = brentq(f, np.log(1e-4), np.log(60.0), xtol=1e-13, rtol=1e-14)
    R = cache.get(lq) or Restricted(ker, float(np.exp(lq)), N=N)
    return R
