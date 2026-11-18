"""C149 restricted dressing solver on a uniform grid (h <= eta/hfac) with Gregory end corrections
and FFT Toeplitz matvecs; Krylov solve of (I - T_Q W) f = s.  Floating point only, not certified.
Accuracy is controlled by exp(-2 pi eta/h) in the interior and (h/eta)^8 at the two ends."""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.optimize import brentq
from scipy.fft import rfft, irfft, next_fast_len

def gregory_end(K=7):
    """Gregory end weights (units of h) for the left end: trapezoid plus corrections with the
    classical coefficients -1/12, -1/24, -19/720, -3/160, -863/60480, -275/24192 on forward differences."""
    from math import comb
    g = [-1/12, -1/24, -19/720, -3/160, -863/60480, -275/24192][:K-1]
    w = np.zeros(K); w[0] = 0.5; w[1:] = 1.0
    for k, gk in enumerate(g, start=1):
        # sign pattern: left end contributes -gk * (-1)^(k+1)?  use exactness test below to fix sign
        for i in range(k+1):
            w[i] += gk*(-1)**(k-i)*comb(k, i)*(-1)**(k+1)*(-1)
    return w

class Kern:
    def __init__(self, eta):
        self.eta = eta; self.te = np.tan(eta); self.th = np.tan(eta/2)
        self.gamma = np.pi-eta; self.b = np.pi/self.gamma; self.a = np.pi*np.sin(eta)
    def T(self, u):
        t = np.tanh(u); return self.te*(1-t*t)/(np.pi*(t*t+self.te**2))
    def a1(self, x):
        t = np.tanh(x); return self.th*(1-t*t)/(np.pi*(1+(t*self.th)**2))

class Uniform:
    def __init__(self, K, Q, hfac=10.0, Nmin=400, Kg=7):
        self.K, self.Q = K, Q
        N = int(max(Nmin, np.ceil(2*Q*hfac/K.eta)))
        if N % 2: N += 1
        self.N = N; self.h = h = 2*Q/N
        self.x = -Q+h*np.arange(N+1)
        w = np.ones(N+1); c = gregory_end(Kg)
        if len(c) > N//2: raise ValueError('grid too small')
        w[:Kg] = c; w[-Kg:] = c[::-1]
        self.w = h*w
        m = np.arange(-N, N+1)
        t = K.T(m*h)
        L = next_fast_len(3*N+3)
        self.Lfft = L
        self.tf = rfft(np.concatenate([t, np.zeros(L-len(t))]))
    def Tq(self, f):
        g = self.w*f
        c = irfft(self.tf*rfft(g, self.Lfft), self.Lfft)
        return c[self.N:2*self.N+1]
    def solve(self, s, tol=1e-13):
        """CG on the symmetrized system (I - W^(1/2) T W^(1/2)) g = W^(1/2) s, f = W^(-1/2) g."""
        from scipy.sparse.linalg import cg
        n = self.N+1; sw = np.sqrt(self.w)
        def mv(g):
            return g-sw*self.Tq(g/sw)
        A = LinearOperator((n, n), matvec=mv, dtype=float)
        g0 = sw*s
        g, info = cg(A, g0, x0=g0.copy(), rtol=tol, atol=0, maxiter=20000)
        f = g/sw
        res = np.max(np.abs(f-self.Tq(f)-s))
        if res > 1e-9*np.max(np.abs(s)):
            raise RuntimeError(f'cg info={info} res={res}')
        return f

class Restricted2:
    def __init__(self, K, Q, withZ=True, **kw):
        self.K, self.Q = K, Q
        U = Uniform(K, Q, **kw); self.U = U
        self.r = U.solve(K.a1(U.x)); self.Z = U.solve(np.ones(U.N+1)) if withZ else None
        self.m = float(np.sum(U.w*self.r))
        self.x = U.x
    def outside(self, u):
        """continuum counting density s_Q(u) = a1(u) + int_{-Q}^{Q} T(u-v) r(v) dv for |u| >= Q."""
        u = np.atleast_1d(u)
        return self.K.a1(u)+np.array([np.sum(self.U.w*self.K.T(ui-self.x)*self.r) for ui in u])

def solve_Q2(K, nu, Qlo=1e-4, Qhi=40.0, withZ=True, **kw):
    cache = {}
    def f(lq):
        R = Restricted2(K, float(np.exp(lq)), withZ=withZ, **kw); cache[lq] = R; return R.m-nu
    lq = brentq(f, np.log(Qlo), np.log(Qhi), xtol=1e-12, rtol=1e-13)
    return cache.get(lq) or Restricted2(K, float(np.exp(lq)), **kw)

if __name__ == '__main__':
    import time
    print(gregory_end(8))
    for eps in [1e-2, 1e-4, 1e-6]:
        eta = 2*np.arcsin(np.sqrt(eps/2)); K = Kern(eta)
        for Q in [1.0, 8.0]:
            for hf in [6, 10, 14]:
                t = time.time(); R = Restricted2(K, Q, hfac=hf)
                print(eps, Q, hf, R.U.N, R.m, R.r[R.U.N//2], 1/(2*K.gamma), R.r[-1], R.Z[-1], round(time.time()-t, 2))
