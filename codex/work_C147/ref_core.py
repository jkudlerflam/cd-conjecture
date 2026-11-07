"""C147 referee: independent solvers (written from scratch, no code reused from C139-C145).
J=1. Hamiltonian at flux pi per bond:
    H = sum_i [ -(1/2)(S+_i S-_(i+1) + h.c.) + Delta Sz_i Sz_(i+1) ],  periodic, L bonds.
Bethe (both signs of Delta): L p(x_j) - sum_k theta(x_j - x_k) = 2 pi I_j, I_j = j-(n+1)/2,
    p(x) = 2 atan(tanh(x) cot(gamma/2)), theta(x) = 2 atan(tanh(x) cot(gamma)), Delta = cos(gamma),
    E = Delta L/4 + sum_j e(x_j),  e(x) = -sin(gamma)^2/(cosh(2x) - cos(gamma)).
Attractive side uses the step-free form L p(x_j) = sum_(k!=j) phi(x_j-x_k), phi(y)=2 atan(tan(eta)/tanh(y)).
"""
import numpy as np
from itertools import combinations
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


# ---------------- exact diagonalization in a fixed S^z sector ----------------
def ed_sector(L, n, Delta, k=1):
    """lowest eigenvalue(s) of H in the sector with n down spins (bits = down)."""
    states = [sum(1 << i for i in c) for c in combinations(range(L), n)]
    idx = {s: i for i, s in enumerate(states)}
    D = len(states)
    rows, cols, vals = [], [], []
    diag = np.zeros(D)
    for a, s in enumerate(states):
        for i in range(L):
            j = (i + 1) % L
            if L == 2 and i == 1:
                j = 0
            bi, bj = (s >> i) & 1, (s >> j) & 1
            diag[a] += Delta * (0.25 if bi == bj else -0.25)
            if bi != bj:
                t = s ^ ((1 << i) | (1 << j))
                rows.append(a); cols.append(idx[t]); vals.append(-0.5)
    H = sp.coo_matrix((vals, (rows, cols)), shape=(D, D)).tocsr() + sp.diags(diag)
    if D <= 400:
        w = np.linalg.eigvalsh(H.toarray())
        return w[:k] if k > 1 else w[0]
    w = eigsh(H, k=max(k, 1), which='SA', return_eigenvectors=False, tol=1e-13)
    w = np.sort(w)
    return w[:k] if k > 1 else w[0]


# ---------------- repulsive Bethe (0 < gamma <= pi/2) ----------------
class Rep:
    def __init__(self, gamma):
        self.g = gamma
        self.D = np.cos(gamma)
        self.b = np.pi / gamma

    def p(self, x):
        return 2 * np.arctan(np.tanh(x) / np.tan(self.g / 2))

    def th(self, x):
        return 2 * np.arctan(np.tanh(x) / np.tan(self.g))

    def a1(self, x):
        return np.sin(self.g) / (np.pi * (np.cosh(2 * x) - self.D))

    def K(self, x):
        return np.sin(2 * self.g) / (2 * np.pi * (np.sinh(x) ** 2 + np.sin(self.g) ** 2))

    def e(self, x):
        return -np.sin(self.g) ** 2 / (np.cosh(2 * x) - self.D)

    def r0(self, x):
        return 1 / (2 * self.g * np.cosh(self.b * x))


def bethe_rep(L, n, gamma, steps=None, tol=1e-14):
    """continue the free solution at gamma=pi/2 down to gamma; returns roots, min Hessian eig along path."""
    I = np.arange(1, n + 1) - (n + 1) / 2
    x = np.arctanh(np.tan(np.pi * I / L))
    if steps is None:
        steps = 60 + int(40 * np.log(np.pi / (2 * gamma) + 1))
    gs = np.pi / 2 - (np.pi / 2 - gamma) * (np.linspace(0, 1, steps + 1)) ** 1.0
    minhess = np.inf
    for g in gs[1:]:
        M = Rep(g)
        for it in range(100):
            dx = x[:, None] - x[None, :]
            F = L * M.p(x) - M.th(dx).sum(1) - 2 * np.pi * I
            Kd = 2 * np.pi * M.K(dx)
            Jm = Kd.copy()
            np.fill_diagonal(Jm, 0)
            Jm = Jm
            diag = L * 2 * np.pi * M.a1(x) - (Kd.sum(1) - np.diag(Kd))
            Jm = Jm + np.diag(diag)
            step = np.linalg.solve(Jm, F)
            x = x - step
            if np.max(np.abs(step)) < tol * max(1, np.max(np.abs(x))):
                break
        else:
            raise RuntimeError('Newton failed at gamma=%g' % g)
        if not np.all(np.diff(x) > 0):
            raise RuntimeError('order lost')
        minhess = min(minhess, np.linalg.eigvalsh(Jm).min())
    return x, minhess


def energy_rep(L, n, gamma, x):
    M = Rep(gamma)
    return M.D * L / 4 + M.e(x).sum()


# ---------------- attractive Bethe (gamma = pi - eta, 0 < eta <= pi/2) ----------------
class Att:
    def __init__(self, eta):
        self.eta = eta
        self.g = np.pi - eta
        self.D = -np.cos(eta)
        self.b = np.pi / self.g
        self.q = 1 - 2 * eta / np.pi
        self.a = np.pi * np.sin(eta)

    def p(self, x):
        return 2 * np.arctan(np.tanh(x) * np.tan(self.eta / 2))

    def a1(self, x):
        return np.sin(self.eta) / (np.pi * (np.cosh(2 * x) + np.cos(self.eta)))

    def T(self, x):
        return np.sin(2 * self.eta) / (2 * np.pi * (np.sinh(x) ** 2 + np.sin(self.eta) ** 2))

    def Tp(self, x):
        s = np.sinh(x)
        return -np.sin(2 * self.eta) / (2 * np.pi) * 2 * s * np.cosh(x) / (s ** 2 + np.sin(self.eta) ** 2) ** 2

    def phi(self, y):
        return 2 * np.arctan(np.tan(self.eta) / np.tanh(y))

    def dphi(self, y):
        # d/dy 2 atan(t/tanh y) = -2 t sech^2 y/(tanh^2 y + t^2) = -2 pi T(y)
        return -2 * np.pi * self.T(y)

    def e(self, x):
        return -np.sin(self.eta) ** 2 / (np.cosh(2 * x) + np.cos(self.eta))

    def r0(self, x):
        return 1 / (2 * self.g * np.cosh(self.b * x))


def bethe_att(L, n, eta, steps=None, tol=1e-14):
    I = np.arange(1, n + 1) - (n + 1) / 2
    x = np.arctanh(np.tan(np.pi * I / L))
    if steps is None:
        steps = 80 + int(60 * np.log(np.pi / (2 * eta) + 1))
    # geometric approach to small eta
    es = np.exp(np.linspace(np.log(np.pi / 2), np.log(eta), steps + 1))
    for et in es[1:]:
        M = Att(et)
        for it in range(200):
            dx = x[:, None] - x[None, :]
            np.fill_diagonal(dx, 1.0)
            ph = M.phi(dx)
            np.fill_diagonal(ph, 0)
            F = L * M.p(x) - ph.sum(1)
            dph = M.dphi(dx)
            np.fill_diagonal(dph, 0)
            Jm = dph.copy()  # dF_j/dx_k = +dphi(x_j-x_k) for k != j ... sign: -d/dx_k phi(x_j-x_k) = +phi'
            diag = L * 2 * np.pi * M.a1(x) - dph.sum(1)
            Jm = Jm + np.diag(diag)
            step = np.linalg.solve(Jm, F)
            lam = 1.0
            xn = x - step
            while not np.all(np.diff(xn) > 0):
                lam /= 2
                xn = x - lam * step
                if lam < 1e-8:
                    raise RuntimeError('ordering')
            x = xn
            if np.max(np.abs(lam * step)) < tol * max(1, np.max(np.abs(x))):
                break
        else:
            raise RuntimeError('Newton failed at eta=%g' % et)
    return x


def energy_att(L, n, eta, x):
    M = Att(eta)
    return M.D * L / 4 + M.e(x).sum()
