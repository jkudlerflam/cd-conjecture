"""XX bosonic majorization (PROOF.md section 4).

Theorem. Let X > 1/2, V1 = X/asinh(X), r > 1 + 6/(pi^2 X V1), V' = r V1. For every even L >= 2 and beta > 0,
    log Z~_XX(L) <= K(X, r) + 2 F(L/(V' beta)),
    K(X, r) = sup_{x>0} Phi(x),
    Phi(x) = 2F(x) - 2F(x/r) + log(theta(x)^2 + theta1(x)^2/2) + 4(2X)^(-2) + 1/(pi X V1 x),
theta(x) = sum_Q exp(-pi x Q^2), theta1(x) = sum_Q exp(-pi x (Q^2+Q)), x = L/(V1 beta).
This script evaluates K on a fine grid (with the analytic small- and large-x behaviour checked) and tests the
inequality against the exact crossed product (C125 eq. 6) on a large (L, beta) grid.
"""
import numpy as np
from majorant import F, xx_logZt


def theta(x):
    Q = np.arange(-60, 61)[:, None]
    return np.exp(-np.pi * np.asarray(x)[None, :] * Q ** 2).sum(0)


def theta1(x):
    Q = np.arange(-60, 61)[:, None]
    return np.exp(-np.pi * np.asarray(x)[None, :] * (Q ** 2 + Q)).sum(0)


def theta_small(x):
    """Poisson-resummed forms valid for small x (used below x = 0.05)"""
    Q = np.arange(-20, 21)[:, None]
    x = np.asarray(x)[None, :]
    th = x[0] ** -0.5 * np.exp(-np.pi * Q ** 2 / x).sum(0)
    th1 = np.exp(np.pi * x[0] / 4) * x[0] ** -0.5 * (((-1.0) ** Q) * np.exp(-np.pi * Q ** 2 / x)).sum(0)
    return th, th1


def Phi(x, X, r):
    V1 = X / np.arcsinh(X)
    x = np.asarray(x, float)
    th = np.where(x > 0.05, theta(np.maximum(x, 0.05)), theta_small(np.minimum(x, 0.05))[0])
    th1 = np.where(x > 0.05, theta1(np.maximum(x, 0.05)), theta_small(np.minimum(x, 0.05))[1])
    return 2 * F(x) - 2 * F(x / r) + np.log(th ** 2 + th1 ** 2 / 2) + 4 * (2 * X) ** -2 + 1 / (np.pi * X * V1 * x)


def K_of(X, r, grid=None):
    V1 = X / np.arcsinh(X)
    assert X > 0.5 and r > 1 + 6 / (np.pi ** 2 * X * V1)
    xs = np.geomspace(1e-4, 1e3, 200001) if grid is None else grid
    ph = Phi(xs, X, r)
    i = int(np.argmax(ph))
    return ph[i], xs[i], ph[0], ph[-1]


def check(X, r, K):
    V1 = X / np.arcsinh(X); Vp = r * V1
    worst = -np.inf; where = None
    for beta in np.geomspace(0.05, 4096, 70):
        Ls = np.array([4, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64, 96, 128, 256, 512, 1024, 2048, 4096])
        g = np.concatenate([xx_logZt([2.], beta), xx_logZt(Ls, beta)])
        LL = np.concatenate([[2], Ls])
        margin = K + 2 * F(LL / (Vp * beta)) - g
        j = int(np.argmin(margin))
        if -margin[j] > worst:
            worst = -margin[j]; where = (int(LL[j]), float(beta))
    return worst, where


if __name__ == '__main__':
    for X, r in [(2.0, 1.25), (2.0, 1.5), (2.0, 2.0), (1.5, 1.5), (3.0, 1.3)]:
        V1 = X / np.arcsinh(X)
        K, xm, left, right = K_of(X, r)
        worst, where = check(X, r, K)
        print("X=%.2f r=%.2f: V1=%.4f V'=%.4f (V'/M=%.3f)  K=%.4f at x=%.3g [Phi(1e-4)=%.1f, Phi(1e3)=%.4f]  exact check: max(log Z~ - bound) = %.4f at (L,beta)=%s" % (
            X, r, V1, r * V1, 2 * r * V1, K, xm, left, right, worst, where))
