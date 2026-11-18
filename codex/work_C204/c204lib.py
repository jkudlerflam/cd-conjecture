"""C204 referee: independent Bethe and continuum solvers (written from scratch, no C149/C200 code).

Conventions (C200 PROOF.md): eta = pi - gamma, Delta = -cos(eta), J = 1.
p(x) = 2 atan(tan(eta/2) tanh x), theta(x) = -2 atan(tanh x / tan eta),
a1 = p'/(2pi), T = -theta'/(2pi), Theta = int_0^x T, q = 1 - 2 eta/pi.
Bethe map: F_j = L p(x_j) - sum_k theta(x_j - x_k) - 2 pi I_j, I_j = j - (n+1)/2.
"""
import numpy as np
from numpy.fft import rfft, irfft
from scipy.sparse.linalg import LinearOperator, cg

PI = np.pi


class Kern:
    def __init__(self, eta):
        self.eta = eta
        self.te = np.tan(eta)
        self.th = np.tan(eta / 2)
        self.q = 1 - 2 * eta / PI
        self.s2 = np.sin(2 * eta)
        self.se2 = np.sin(eta) ** 2

    def p(self, x):
        return 2 * np.arctan(self.th * np.tanh(x))

    def a1(self, x):
        # sin(eta)/(pi (cosh 2x + cos eta)), written with sinh^2 to avoid overflow issues
        e = self.eta
        return np.sin(e) / (PI * (2 * np.sinh(np.abs(x)) ** 2 + 1 + np.cos(e)))

    def theta(self, x):
        return -2 * np.arctan(np.tanh(x) / self.te)

    def T(self, x):
        return self.s2 / (2 * PI * (np.sinh(x) ** 2 + self.se2))

    def Theta(self, x):
        return np.arctan(np.tanh(x) / self.te) / PI

    def B(self, x):
        return self.q / 2 - self.Theta(x)


def bethe_solve(eta, L, n, x0=None, tol=1e-13, maxit=200, verbose=False):
    """Newton with exact line search on the strictly convex Yang-Yang action."""
    K = Kern(eta)
    I = np.arange(1, n + 1) - (n + 1) / 2
    if x0 is None:
        x0 = np.linspace(-1, 1, n) * 0.5
    x = np.array(x0, dtype=float)

    def F(x):
        d = x[:, None] - x[None, :]
        return L * K.p(x) - K.theta(d).sum(1) - 2 * PI * I

    def Jac(x):
        d = x[:, None] - x[None, :]
        Tm = K.T(d)
        np.fill_diagonal(Tm, 0.0)
        J = -2 * PI * Tm
        diag = L * 2 * PI * K.a1(x) + 2 * PI * Tm.sum(1)
        J[np.diag_indices_from(J)] = diag
        return J

    f = F(x)
    for it in range(maxit):
        J = Jac(x)
        dx = -np.linalg.solve(J, f)
        # exact line search on phi(s)=S(x+s dx): phi'(s) = F(x+s dx).dx, increasing in s
        s = 1.0
        g1 = F(x + dx) @ dx
        if g1 > 0:
            lo, hi = 0.0, 1.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if F(x + mid * dx) @ dx > 0:
                    hi = mid
                else:
                    lo = mid
            s = 0.5 * (lo + hi)
        x = x + s * dx
        f = F(x)
        if verbose:
            print(it, s, np.abs(f).max(), np.abs(dx).max())
        if (np.abs(s * dx).max() < tol and s == 1.0) or np.abs(f).max() < 1e-12 * max(L, 1):
            break
    x.sort()
    return x, np.abs(F(x)).max(), it


def gregory_weights(m):
    """End corrections for the trapezoid rule, exact for polynomials of degree < 2m on long grids.
    Solved numerically on a symmetric grid of N = 6m intervals; returns w_0..w_{m-1}."""
    N = 6 * m
    xs = np.arange(N + 1, dtype=float)
    c = N / 2
    A = np.zeros((m, m))
    rhs = np.zeros(m)
    for row, k in enumerate(range(0, 2 * m, 2)):
        f = (xs - c) ** k
        exact = 2 * c ** (k + 1) / (k + 1)
        base = f.sum()  # all weights one
        # corrections at both ends, symmetric
        for j in range(m):
            A[row, j] = 2 * f[j]
        rhs[row] = exact - base
    corr = np.linalg.solve(A, rhs)
    return 1.0 + corr


class Continuum:
    """r = a1 + T_Q r on [-Q, Q], trapezoid with Gregory end corrections, PCG with FFT."""

    def __init__(self, eta, Q, hfac=10, m=6):
        self.K = K = Kern(eta)
        self.eta, self.Q = eta, Q
        N = max(int(np.ceil(2 * Q / (eta / hfac))), 8 * m)
        self.N = N
        self.h = h = 2 * Q / N
        self.v = -Q + h * np.arange(N + 1)
        w = np.ones(N + 1)
        g = gregory_weights(m)
        w[:m] = g
        w[-m:] = g[::-1]
        self.w = w
        assert w.min() > 0
        # circulant embedding of Toeplitz T(v_i - v_j)
        M = N + 1
        P = 1
        while P < 2 * M:
            P *= 2
        self.P = P
        col = np.zeros(P)
        k = np.arange(M)
        col[:M] = K.T(k * h)
        col[P - M + 1:] = K.T((M - 1 - np.arange(M - 1)) * h)
        self.Tf = rfft(col)
        # preconditioner: inverse of (I - h T) on whole line (symbol 1 - h That > 1 - q > 0)
        sym = 1 - h * np.real(self.Tf)
        self.prec = 1.0 / np.maximum(sym, 1e-3)
        self.sw = np.sqrt(w)
        self.r = self.solve()

    def Tmul(self, y):
        z = np.zeros(self.P)
        z[: self.N + 1] = y
        return irfft(rfft(z) * self.Tf, self.P)[: self.N + 1]

    def solve(self, src=None, tol=1e-13):
        if src is None:
            src = self.K.a1(self.v)
        h, sw = self.h, self.sw
        n1 = self.N + 1
        A = LinearOperator((n1, n1), matvec=lambda y: y - h * sw * self.Tmul(sw * y), dtype=float)

        def pm(y):
            z = np.zeros(self.P)
            z[:n1] = y
            return irfft(rfft(z) * self.prec, self.P)[:n1]

        Mop = LinearOperator((n1, n1), matvec=pm, dtype=float)
        y, info = cg(A, sw * src, rtol=tol, atol=0, maxiter=5000, M=Mop)
        assert info == 0, info
        return y / sw

    # Nystrom extension and cumulative functions
    def s(self, u):
        u = np.atleast_1d(u)
        out = self.K.a1(u).copy()
        wr = self.h * self.w * self.r
        for i0 in range(0, len(u), 256):
            uu = u[i0:i0 + 256]
            out[i0:i0 + 256] += self.K.T(uu[:, None] - self.v[None, :]) @ wr
        return out

    def N0(self, u):
        """int_{-Q}^{min(u,Q)} r, via the Nystrom interpolant (exact kernel antiderivatives)."""
        u = np.atleast_1d(u)
        uc = np.clip(u, -self.Q, self.Q)
        K = self.K
        # int a1 from -Q to uc: (p(uc) - p(-Q))/(2pi)
        out = (K.p(uc) - K.p(-self.Q)) / (2 * PI)
        wr = self.h * self.w * self.r
        for i0 in range(0, len(u), 256):
            uu = uc[i0:i0 + 256]
            out[i0:i0 + 256] += (K.Theta(uu[:, None] - self.v[None, :]) - K.Theta(-self.Q - self.v[None, :])) @ wr
        return out

    def mass(self):
        return float(self.N0(np.array([self.Q]))[0])

    def Z0(self, u):
        """continuum counting function p/(2pi) + T*N0 - q nu/2, computed as
        p/(2pi) + int Theta(u-v) r(v) dv (equal, see PROOF 4.1)."""
        u = np.atleast_1d(u)
        out = self.K.p(u) / (2 * PI)
        wr = self.h * self.w * self.r
        for i0 in range(0, len(u), 256):
            uu = u[i0:i0 + 256]
            out[i0:i0 + 256] += self.K.Theta(uu[:, None] - self.v[None, :]) @ wr
        return out


def continuum_for_nu(eta, nu, Qguess=1.0, hfac=10, tol=1e-12):
    """Find Q with m(Q) = nu by secant in log Q."""
    def m(Q):
        return Continuum(eta, Q, hfac=hfac)

    Q0 = Qguess
    c0 = m(Q0)
    f0 = c0.mass() - nu
    Q1 = Q0 * (1.3 if f0 < 0 else 0.75)
    c1 = m(Q1)
    f1 = c1.mass() - nu
    for _ in range(40):
        if abs(f1) < tol:
            break
        lq = np.log(Q1) - f1 * (np.log(Q1) - np.log(Q0)) / (f1 - f0)
        lq = np.clip(lq, np.log(Q1) - 1.0, np.log(Q1) + 1.0)
        Q0, f0 = Q1, f1
        Q1 = float(np.exp(lq))
        c1 = m(Q1)
        f1 = c1.mass() - nu
    return c1
