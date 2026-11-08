"""C162 (Claude Opus agent, 2026-11-08). General K-stratum blow-up functional (C160 eq. 2.2 with D = z):

    E[n] = int dx [ |n'|^2/2 - (n x n')_z + W(n) ],   W = s * ( a n_x^2/2 + b n_y^2/2 - hx n_x - hy n_y ),

gauge S_zz = 0 (|n| = 1). Z2: n_z -> -n_z.  For hy != 0 the field is off the in-plane principal axes, which breaks the
extra classical symmetry x -> -x, n_y -> -n_y of C161's family, so solitons need not be left-right symmetric.
Helpers: in-plane background, tail masses, planar soliton (DM is topological for planar textures), out-of-plane
Hessian, full 3-component minimization on an interval and on a ring.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.linalg import eigh_tridiagonal


class KFun:
    def __init__(self, a, b, hx, hy, s=1.0):
        self.a, self.b, self.hx, self.hy, self.s = a, b, hx, hy, s

    def W(self, n):
        return self.s * (0.5 * self.a * n[..., 0] ** 2 + 0.5 * self.b * n[..., 1] ** 2 - self.hx * n[..., 0] - self.hy * n[..., 1])

    def dW(self, n):
        g = np.zeros_like(n)
        g[..., 0] = self.s * (self.a * n[..., 0] - self.hx)
        g[..., 1] = self.s * (self.b * n[..., 1] - self.hy)
        return g

    def Veq(self, p):
        return self.W(np.stack([np.cos(p), np.sin(p), 0 * p], -1))

    def background(self):
        ps = np.linspace(0, 2 * np.pi, 2001)
        p0 = ps[np.argmin(self.Veq(ps))]
        r = minimize_scalar(lambda p: self.Veq(np.array(p)), bracket=(p0 - 0.01, p0, p0 + 0.01))
        pb = float(r.x) % (2 * np.pi)
        # global check on S^2
        th = np.linspace(0, np.pi, 181)[:, None]; ph = np.linspace(0, 2 * np.pi, 361)[None, :]
        n = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th) + 0 * ph], -1)
        glob = self.W(n).min() >= self.Veq(np.array(pb)) - 1e-9
        return pb, glob

    def masses(self, pb):
        c, sn = np.cos(pb), np.sin(pb)
        Q = 0.5 * (self.a * c * c + self.b * sn * sn); L = self.hx * c + self.hy * sn
        mw2 = self.s * (-2 * Q + L)
        eps = 1e-4
        mu2 = (self.Veq(np.array(pb + eps)) - 2 * self.Veq(np.array(pb)) + self.Veq(np.array(pb - eps))) / eps ** 2
        return mu2, mw2

    def d2theta(self, p):
        c, sn = np.cos(p), np.sin(p)
        return self.s * (-(self.a * c * c + self.b * sn * sn) + self.hx * c + self.hy * sn)


def planar(F, pb, X=30.0, dx=0.1):
    """planar 2 pi soliton phi: pb -> pb + 2 pi; returns E_gp (gradient + potential, DM excluded), grid, phi."""
    N = int(round(2 * X / dx)); xs = -X + dx * np.arange(1, N)
    V0 = F.Veq(np.array(pb))
    def f(p):
        ph = np.r_[pb, p, pb + 2 * np.pi]
        d = np.diff(ph) / dx
        E = np.sum(0.5 * d ** 2) * dx + np.sum(F.Veq(p) - V0) * dx
        n = np.stack([np.cos(p), np.sin(p), 0 * p], -1)
        g1 = F.dW(n)
        dV = -g1[:, 0] * np.sin(p) + g1[:, 1] * np.cos(p)
        g = (d[:-1] - d[1:]) + dV * dx
        return E, g
    p0 = pb + np.pi * (1 + np.tanh(xs))
    r = minimize(f, p0, jac=True, method='L-BFGS-B', options=dict(maxiter=50000, gtol=1e-11))
    return r.fun, xs, r.x


def tilt_lowest(F, pb, xs, p):
    dx = xs[1] - xs[0]
    ph = np.r_[pb, p, pb + 2 * np.pi]
    d = np.diff(ph) / dx; dp = 0.5 * (d[:-1] + d[1:])
    U = dp * (2 - dp) + F.d2theta(p)
    w, v = eigh_tridiagonal(2 / dx ** 2 + U, -np.ones(len(p) - 1) / dx ** 2, select='i', select_range=(0, 0))
    return w[0], v[:, 0]


def ring_energy(v, dx, F, nb):
    """discrete energy on a ring, relative to the uniform background nb (per configuration)."""
    V = v.reshape(-1, 3); nr = np.linalg.norm(V, axis=1); n = V / nr[:, None]
    n1 = np.roll(n, -1, 0); d = (n1 - n) / dx
    cz = n[:, 0] * n1[:, 1] - n[:, 1] * n1[:, 0]
    E = np.sum(0.5 * np.sum(d ** 2, 1)) * dx - np.sum(cz) + np.sum(F.W(n) - F.W(nb)) * dx
    G = -d + np.roll(d, 1, 0)
    n0 = np.roll(n, 1, 0)
    G[:, 0] += -n1[:, 1] + n0[:, 1]
    G[:, 1] += n1[:, 0] - n0[:, 0]
    G += F.dW(n) * dx
    G = (G - np.sum(G * n, 1)[:, None] * n) / nr[:, None]
    return E, G.ravel()


def ring_min(F, pb, T, tilts, dx=0.1, shift=None, n_init=None):
    """m = len(tilts) solitons (each +2 pi) on a ring of length m*T, initial tilt signs given."""
    m = len(tilts); L = m * T; N = int(round(L / dx)); dx = L / N
    x = np.arange(N) * dx
    nb = np.array([np.cos(pb), np.sin(pb), 0.0])
    if n_init is None:
        ph = np.full(N, pb); nz = np.zeros(N)
        for j in range(m):
            c = (j + 0.5) * T
            ph += np.pi * (1 + np.tanh(1.0 * (x - c)))
            nz += tilts[j] * np.exp(-(x - c) ** 2 / 4)
        n0 = np.stack([np.cos(ph), np.sin(ph), nz], 1); n0 /= np.linalg.norm(n0, axis=1)[:, None]
    else:
        n0 = n_init
    r = minimize(ring_energy, n0.ravel(), args=(dx, F, nb), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=80000, gtol=1e-11, ftol=1e-16))
    V = r.x.reshape(-1, 3); n = V / np.linalg.norm(V, axis=1)[:, None]
    return r.fun, n, dx
