"""C138 Bethe solver near the ferromagnet, parametrized by eta = pi - gamma.

J=1. Delta = cos(gamma) = -cos(eta), epsilon = 1+Delta = 2 sin(eta/2)^2.
Same equations and conventions as C129/C134/C136 (bethe_scan.py there):

    L p(x_j) - sum_k theta(x_j - x_k) = 2 pi I_j,   I_j = j-(n+1)/2,

with p(x)=2 atan[tanh x cot(gamma/2)], theta(x)=2 atan[tanh x cot gamma].
For ordered roots, theta(y) = -pi sign(y) + phi(y) with
phi(y) = 2 atan[tan(eta)/tanh(y)], and sum_k -pi sign(x_j-x_k) = -2 pi I_j,
so the equations are equivalent to the step-free form

    L p(x_j) = sum_{k != j} phi(x_j - x_k),                       (*)

which is what is solved here (no O(1) cancellation as eta -> 0).
F_L(n) = Delta L/4 + G_L(n),  G_L(n) = sum_j e(x_j),
e(x) = -sin(eta)^2 / [2 (sinh(x)^2 + cos(eta/2)^2)].
Floating point only, not certified enclosures.
"""
import numpy as np


class Model:
    def __init__(self, eta):
        self.eta = float(eta)
        self.te = np.tan(eta)
        self.th = np.tan(eta/2)
        self.s2 = np.sin(eta)**2
        self.c2 = np.cos(eta/2)**2
        self.eps = 2*np.sin(eta/2)**2
        self.Delta = -np.cos(eta)

    def p(self, x):
        return 2*np.arctan(np.tanh(x)*self.th)

    def dp(self, x):
        u = np.tanh(x)
        return 2*self.th*(1-u*u)/(1+(u*self.th)**2)

    def phi(self, y):
        return 2*np.arctan(self.te/np.tanh(y))

    def dphi(self, y):
        # d/dy 2 atan(te/tanh y) = -2 te sech^2 y / (tanh^2 y + te^2)
        u = np.tanh(y)
        return -2*self.te*(1-u*u)/(u*u+self.te**2)

    def e(self, x):
        return -self.s2/(2*(np.sinh(x)**2+self.c2))


def residual(m, x, L):
    d = x[:, None]-x[None, :]
    np.fill_diagonal(d, 1.0)
    ph = m.phi(d)
    np.fill_diagonal(ph, 0.0)
    return L*m.p(x)-ph.sum(axis=1)


def jacobian(m, x, L):
    d = x[:, None]-x[None, :]
    np.fill_diagonal(d, 1.0)
    W = m.dphi(d)            # negative entries
    np.fill_diagonal(W, 0.0)
    J = W.copy()             # d res_j / d x_k = + phi'(x_j-x_k)
    J[np.diag_indices_from(J)] = L*m.dp(x)-W.sum(axis=1)
    return J


def seed(m, L, n):
    """Continuum-motivated seed: half-filled sech counting, compressed for n<L/2."""
    I = np.arange(1, n+1)-(n+1)/2
    g = np.pi-m.eta
    # z_inf(x) = (1/(2pi)) * 2 atan(tanh(b x/2)) ... use gd inverse of 2 pi I/L
    t = np.clip(2*np.pi*I/L, -np.pi/2+1e-12, np.pi/2-1e-12)
    return (g/np.pi)*np.arcsinh(np.tan(t))


def solve(L, n, eta, x0=None, tol=1e-13, maxit=200):
    m = Model(eta)
    if n == 0:
        return 0.0, np.zeros(0), dict(steps=0, res=0.0)
    x = seed(m, L, n) if x0 is None else np.array(x0, float)
    if n == 1:
        return float(m.e(0.0)), np.zeros(1), dict(steps=0, res=0.0)
    scale = L*m.dp(0.0)  # natural size of the equations
    for step in range(maxit):
        r = residual(m, x, L)
        nr = np.max(np.abs(r))/scale
        if nr < tol:
            break
        J = jacobian(m, x, L)
        dx = np.linalg.solve(J, -r)
        lam = 1.0
        for _ in range(60):
            xt = x+lam*dx
            xt = (xt-xt[::-1])/2
            if np.all(np.diff(xt) > 0):
                rt = np.max(np.abs(residual(m, xt, L)))/scale
                if rt < nr*(1-1e-4*lam) or rt < tol:
                    x = xt
                    break
            lam /= 2
        else:
            raise RuntimeError(f'line search failed L={L} n={n} eta={eta} res={nr}')
    else:
        raise RuntimeError(f'no convergence L={L} n={n} eta={eta} res={nr}')
    G = float(np.sum(m.e(x).astype(np.longdouble)))
    return G, x, dict(steps=step, res=float(nr))


if __name__ == '__main__':
    import time
    for eta in [0.3, 0.05, 0.002]:
        for (L, n) in [(8, 4), (16, 3), (64, 32), (512, 256), (512, 17)]:
            t = time.time()
            G, x, d = solve(L, n, eta)
            print(eta, L, n, G, d, round(time.time()-t, 3), x[-1])
