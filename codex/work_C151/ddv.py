"""C151: Destri-de Vega contour tools for the attractive XXZ ring near Delta=-1.

J=1, eta = pi - gamma, Delta = -cos(eta). Conventions of C136/C138 (bethe_ferro.py):
    theta(z) = (L/2pi) p(z) + sum_k Theta(z - x_k),   Theta(u) = (1/pi) atan(tanh u / tan eta),
so theta(x_j) = I_j = j - (n+1)/2, theta' = L sigma, sigma = a1 + (1/L) sum_k T(. - x_k).
w(z) = (-1)^n exp(2 pi i theta(z)), ell(z) = Log(1 + w(z)) (principal), analytic for 0 < Im z < min(eta, pi/4).
Contour identity (PROOF Part 1): for f analytic near the rectangle [-Q,Q] x [-eps,eps], real on R,
    sum_j f(x_j) - L int_{-Q}^{Q} f sigma = -(1/pi) Im int_{Gamma_+} f'(z) ell(z) dz,
with Gamma_+ any path from Q to -Q in the closed upper half-strip (theta(+-Q) = +-n/2).
Floating point only, no enclosures.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

GX, GW = leggauss(20)


class DdV:
    def __init__(self, eta, L, roots):
        self.eta = float(eta)
        self.L = L
        self.x = np.asarray(roots, float)
        self.n = len(self.x)
        self.te = np.tan(eta)
        self.th = np.tan(eta/2)
        self.sgn = (-1.0)**self.n
        self.a = np.pi*np.sin(eta)

    # ---- kernels, valid for complex arguments in the strip ----
    def T(self, u):
        e = self.eta
        return np.sin(2*e)/(np.pi*(np.cosh(2*u)-np.cos(2*e)))

    def a1(self, z):
        e = self.eta
        return np.sin(e)/(np.pi*(np.cosh(2*z)+np.cos(e)))

    def Theta(self, u):
        return np.arctan(np.tanh(u)/self.te)/np.pi

    def p(self, z):
        return 2*np.arctan(self.th*np.tanh(z))

    def theta(self, z, chunk=4000):
        z = np.atleast_1d(np.asarray(z, complex))
        out = np.empty(z.shape, complex)
        for i in range(0, len(z), chunk):
            zz = z[i:i+chunk]
            out[i:i+chunk] = self.L*self.p(zz)/(2*np.pi) + self.Theta(zz[:, None]-self.x[None, :]).sum(axis=1)
        return out

    def sigma(self, z, chunk=4000):
        z = np.atleast_1d(np.asarray(z, complex))
        out = np.empty(z.shape, complex)
        for i in range(0, len(z), chunk):
            zz = z[i:i+chunk]
            out[i:i+chunk] = self.a1(zz) + self.T(zz[:, None]-self.x[None, :]).sum(axis=1)/self.L
        return out

    def ell(self, z):
        return np.log(1+self.sgn*np.exp(2j*np.pi*self.theta(z)))

    def dell(self, z):
        """ell'(z) = 2 pi i L sigma w/(1+w)."""
        th = self.theta(z)
        w = self.sgn*np.exp(2j*np.pi*th)
        return 2j*np.pi*self.L*self.sigma(z)*w/(1+w)

    def Psi(self, x, eps):
        """2 pi Im theta(x + i eps) = -log|w(x+i eps)|."""
        return 2*np.pi*np.imag(self.theta(np.asarray(x)+1j*eps))

    def counting_boundary(self):
        """Q_L with theta(Q_L) = n/2 (bisection on the real axis)."""
        lo, hi = self.x[-1], self.x[-1]+1.0
        while np.real(self.theta(hi))[0] < self.n/2:
            hi += 1.0
            if hi > self.x[-1]+200:
                return np.inf
        for _ in range(200):
            mid = (lo+hi)/2
            if np.real(self.theta(mid))[0] < self.n/2:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2


def gauss_nodes(a, b, h):
    """Composite Gauss-Legendre nodes/weights on [a,b] with panels of length <= h (a, b real or complex)."""
    m = max(1, int(np.ceil(abs(b-a)/h)))
    edges = a+(b-a)*np.linspace(0, 1, m+1)
    lo, hi = edges[:-1], edges[1:]
    mid, rad = (lo+hi)/2, (hi-lo)/2
    z = (mid[:, None]+rad[:, None]*GX[None, :]).ravel()
    w = (rad[:, None]*GW[None, :]).ravel()
    return z, w


def path_integral(F, pts, h):
    """int of F along the polygonal path through pts (complex), panels of length <= h."""
    tot = 0j
    for a, b in zip(pts[:-1], pts[1:]):
        z, w = gauss_nodes(a, b, h)
        tot += np.sum(w*F(z))
    return tot


def real_integral(F, a, b, roots, h, eta):
    """int_a^b F(x) dx on the real axis with panels refined near the roots (F has log singularities there)."""
    br = [a, b] + [r for r in roots if a < r < b]
    br = np.sort(np.array(br))
    tot = 0.0
    for lo, hi in zip(br[:-1], br[1:]):
        # graded panels toward both endpoints (roots): geometric near the ends
        L = hi-lo
        k = np.geomspace(1e-12, 1, 30)
        cuts = np.unique(np.concatenate([lo+L/2*k, hi-L/2*k, [lo, hi], lo+np.arange(0, L, h)]))
        for u, v in zip(cuts[:-1], cuts[1:]):
            if v <= u:
                continue
            z, w = gauss_nodes(u, v, max(v-u, 1e-300))
            tot += np.sum(w*F(z))
    return tot
