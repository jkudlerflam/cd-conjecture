"""C138 continuum solver for the restricted dressing equations near Delta=-1.

J=1, gamma = pi-eta. On [-Q,Q] (even functions, reduced to [0,Q]) solve

    f(x) - int_{-Q}^{Q} T(x-y) f(y) dy = s(x),   T = -theta'/(2 pi) >= 0,

for s = a1 (root density r_Q), s = 1 (dressed charge Z_Q), s = e - h
(dressed energy eps_Q, h = -a r_Q(Q)/Z_Q(Q)).  Piecewise-linear product
integration on a mesh graded at scale eta toward the Fermi point x=Q; the
kernel is integrated exactly up to composite Gauss quadrature that resolves
its width eta.  Floating point, not certified.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

GX, GW = leggauss(16)
G6X, G6W = leggauss(8)


class Kernel:
    def __init__(self, eta):
        self.eta = eta
        self.te = np.tan(eta)
        self.th = np.tan(eta/2)
        self.gamma = np.pi-eta
        self.a = np.pi*np.sin(eta)          # a = pi sin(gamma)
        self.Delta = -np.cos(eta)
        self.eps = 2*np.sin(eta/2)**2

    def T(self, u):
        t = np.tanh(u)
        return self.te*(1-t*t)/(np.pi*(t*t+self.te**2))

    def a1(self, x):
        t = np.tanh(x)
        return self.th*(1-t*t)/(np.pi*(1+(t*self.th)**2))

    def e(self, x):
        return -np.sin(self.eta)**2/(2*(np.sinh(x)**2+np.cos(self.eta/2)**2))


def _moments(K, x, a, b):
    """int_a^b T(x-y) * {1, (y-a)/(b-a)} dy for arrays x,a,b (same shape).

    The y-interval is split at breakpoints x +- eta*2^j so that T is resolved.
    """
    x, a, b = np.broadcast_arrays(np.asarray(x, float), np.asarray(a, float),
                                  np.asarray(b, float))
    h = b-a
    eta = K.eta
    offs = eta*2.0**np.arange(-4, 12)
    bps = np.concatenate([[0.0], offs, -offs])
    out0 = np.zeros(x.shape)
    out1 = np.zeros(x.shape)
    # breakpoints in y: x + bps, clipped to [a,b], plus a and b
    Y = x[..., None]+bps[None, :] if x.ndim == 1 else None
    Y = np.clip(Y, a[:, None], b[:, None])
    Y = np.concatenate([a[:, None], Y, b[:, None]], axis=1)
    Y.sort(axis=1)
    lo, hi = Y[:, :-1], Y[:, 1:]
    mid, rad = (lo+hi)/2, (hi-lo)/2
    for gx, gw in zip(GX, GW):
        y = mid+rad*gx
        w = rad*gw*K.T(x[:, None]-y)
        out0 += w.sum(axis=1)
        out1 += (w*(y-a[:, None])).sum(axis=1)
    return out0, out1/h


def mesh(K, Q, N, c=None):
    """Smooth graded mesh y = Q - Q sinh(c u)/sinh(c), u uniform in [0,1] (u=0 at the Fermi point).

    c is fixed by the base resolution so that Richardson extrapolation in N is clean."""
    u = np.linspace(0, 1, N+1)
    if c is None or c < 1e-8:
        d = Q*u
    else:
        d = Q*np.sinh(c*u)/np.sinh(c)
    return np.sort(Q-d)


def grading(K, Q, N0, edge_fraction=10):
    """c such that the edge spacing at N0 is about eta/edge_fraction."""
    target = K.eta/edge_fraction
    if Q/N0 <= target:
        return 0.0
    lo, hi = 1e-6, 60.0
    for _ in range(200):
        c = (lo+hi)/2
        if Q*c/(N0*np.sinh(c)) > target:
            lo = c
        else:
            hi = c
    return c


def build(K, Q, N, c):
    y = mesh(K, Q, N, c)
    y[0] = 0.0
    N1 = len(y)
    A = np.eye(N1)
    yl, yr = y[:-1], y[1:]
    # direct and reflected images: int_0^Q [T(x-y)+T(x+y)] f(y) dy
    for i in range(N1):
        xi = np.full(N1-1, y[i])
        m0, m1 = _far_or_near(K, xi, yl, yr)
        r0, r1 = _far_or_near(K, -xi, yl, yr)
        c0, c1 = m0+r0, m1+r1
        A[i, :-1] -= c0-c1
        A[i, 1:] -= c1
    return y, A


def _far_or_near(K, x, a, b):
    """Moments with a cheap 8-point Gauss rule on cells far from x."""
    h = b-a
    mid = (a+b)/2
    far = np.abs(x-mid) > 4*h+4*K.eta
    m0 = np.zeros_like(a)
    m1 = np.zeros_like(a)
    if far.any():
        xf, af, hf = x[far], a[far], h[far]
        s0 = np.zeros(far.sum())
        s1 = np.zeros(far.sum())
        for gx, gw in zip(G6X, G6W):
            yy = af+hf*(gx+1)/2
            w = hf*gw/2*K.T(xf-yy)
            s0 += w
            s1 += w*(yy-af)/hf
        m0[far], m1[far] = s0, s1
    nf = ~far
    if nf.any():
        m0[nf], m1[nf] = _moments(K, x[nf], a[nf], b[nf])
    return m0, m1


class Restricted:
    def __init__(self, K, Q, N=600, c=None):
        self.K, self.Q = K, Q
        if c is None:
            c = grading(K, Q, N)
        self.c = c
        self.y, A = build(K, Q, N, c)
        y = self.y
        S = np.column_stack([K.a1(y), np.ones_like(y), K.e(y)])
        X = np.linalg.solve(A, S)
        self.r, self.Z = X[:, 0], X[:, 1]
        self.h = -K.a*self.r[-1]/self.Z[-1]
        self.epsQ = X[:, 2]-self.h*self.Z
        self.m = 2*self.integrate(self.r)
        self.E = K.Delta/4+2*self.integrate(None, prod=(K.e, self.r))

    def integrate(self, f, prod=None):
        """int_0^Q of piecewise-linear f (times a smooth function if prod given)."""
        y = self.y
        tot = 0.0
        a, b = y[:-1], y[1:]
        for gx, gw in zip(G6X, G6W):
            t = (gx+1)/2
            yy = a+(b-a)*t
            if prod is None:
                val = f[:-1]*(1-t)+f[1:]*t
            else:
                fun, g = prod
                val = fun(yy)*(g[:-1]*(1-t)+g[1:]*t)
            tot += np.sum((b-a)*gw/2*val)
        return tot

    def interp(self, f, x):
        return np.interp(np.abs(x), self.y, f)

    def curvature(self):
        """e_bulk''(nu) at nu=m(Q): eps_Q'(Q) / (2 r(Q) Z(Q)^2), one-sided slope."""
        y, e = self.y, self.epsQ
        c = np.polyfit(y[-4:]-y[-1], e[-4:], 3)
        slope = c[2]
        return slope/(2*self.r[-1]*self.Z[-1]**2)


def half_filled_bulk(K):
    """e_bulk(1/2) = Delta/4 + int e(x) r_inf(x) dx, r_inf = 1/(2 gamma cosh(pi x/gamma))."""
    x, w = leggauss(200)
    X = 30*(x+1)/2
    W = 30*w/2
    g = K.gamma
    return K.Delta/4+2*np.sum(W*K.e(X)/(2*g*np.cosh(np.pi*X/g)))


def bulk_at(K, nu, Q0=None, tol=1e-12, **kw):
    """e_bulk(nu) via m(Q)=nu, Newton in Q with m'(Q)=2 r(Q) Z(Q)."""
    Q = Q0 if Q0 is not None else 1.0
    for _ in range(60):
        R = Restricted(K, Q, **kw)
        f = R.m-nu
        if abs(f) < tol:
            return R
        dQ = -f/(2*R.r[-1]*R.Z[-1])
        Q = max(Q/4, min(4*Q, Q+dQ))
    raise RuntimeError('bulk_at did not converge')


if __name__ == '__main__':
    K = Kernel(0.3)
    c = grading(K, 2.0, 200)
    for N in [200, 400, 800]:
        R = Restricted(K, 2.0, N=N, c=c)
        print(N, R.m, R.E, R.Z[0], R.r[-1], R.curvature(), R.h)
    R = Restricted(K, 12.0, N=1200)
    print('Q=12', R.m, R.E, 'half-filled closed form', half_filled_bulk(K))
