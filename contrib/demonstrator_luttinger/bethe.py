#!/usr/bin/env python3
"""Bethe-ansatz Luttinger parameter and velocity of H = J sum (Sx Sx + Sy Sy + Delta Sz Sz), 0 <= Delta < 1,
at magnon (down-spin) density n < 1/2, Delta = cos(gamma).

Fermi rapidity B, kernels a_k(x) = sin(k gamma) / [2 pi (cosh x - cos k gamma)]:
    rho + a_2 * rho = a_1,   Z + a_2 * Z = 1,   eps + a_2 * eps = h - 2 pi J sin(gamma) a_1,
all convolutions over [-B, B];  n = int rho,  K = Z(B)^2,  v = eps'(B) / (2 pi rho(B)),  h from eps(B) = 0.
Checks (run as a script): Delta = 0 gives K = 1, v = J sin(pi n); B -> infinity approaches
K = pi / (2 (pi - gamma)), v = (pi/2) J sin(gamma)/gamma.
Claude (Opus 5.5), 2026-10-28.
"""
import numpy as np
from scipy.optimize import brentq


def _a(k, g, x):
    return np.sin(k * g) / (2 * np.pi * (np.cosh(x) - np.cos(k * g)))


def _da(k, g, x):
    return -np.sin(k * g) * np.sinh(x) / (2 * np.pi * (np.cosh(x) - np.cos(k * g)) ** 2)


def solve_B(Delta, B, J=1.0, npts=400):
    g = np.arccos(Delta)
    x, wq = np.polynomial.legendre.leggauss(npts)
    x = B * x; wq = B * wq
    Kmat = _a(2, g, x[:, None] - x[None, :]) * wq[None, :]
    A = np.eye(npts) + Kmat
    src_eps = -2 * np.pi * J * np.sin(g) * _a(1, g, x)
    rho = np.linalg.solve(A, _a(1, g, x))
    Z = np.linalg.solve(A, np.ones(npts))
    eJ = np.linalg.solve(A, src_eps)
    # Nystrom extension to x = B
    ext = lambda f, src, y: src - np.sum(_a(2, g, y - x) * wq * f)
    ZB = ext(Z, 1.0, B)
    eJB = ext(eJ, -2 * np.pi * J * np.sin(g) * _a(1, g, B), B)
    h = -eJB / ZB
    eps = eJ + h * Z
    rhoB = ext(rho, _a(1, g, B), B)
    depsB = -2 * np.pi * J * np.sin(g) * _da(1, g, B) - np.sum(_da(2, g, B - x) * wq * eps)
    n = float(np.sum(wq * rho))
    return dict(n=n, K=float(ZB ** 2), v=float(depsB / (2 * np.pi * rhoB)), h=float(h), B=B)


def at_density(Delta, n, J=1.0, npts=400):
    B = brentq(lambda b: solve_B(Delta, b, J, npts)['n'] - n, 1e-3, 30.0, xtol=1e-13)
    return solve_B(Delta, B, J, npts)


if __name__ == '__main__':
    r = at_density(0.0, 1 / 3)
    print('free n=1/3:', r, 'expected v', np.sin(np.pi / 3))
    for D in (0.03708791208791209, 0.07417582417582419, 0.14835164835164835, 0.5):
        g = np.arccos(D)
        r = solve_B(D, 14.0, npts=800)
        print(f'Delta={D:.5f} B=14: n={r["n"]:.6f} K={r["K"]:.6f} (half filling {np.pi / (2 * (np.pi - g)):.6f}) '
              f'v={r["v"]:.6f} ({np.pi / 2 * np.sin(g) / g:.6f})')
        r = at_density(D, 1 / 3)
        print(f'   n=1/3: K={r["K"]:.8f} v={r["v"]:.8f} h={r["h"]:.6f}  (pert. K = 1 - 2 Delta sin(kF)/pi = {1 - 2 * D * np.sin(np.pi / 3) / np.pi:.6f})')
    # quadrature convergence
    for npts in (100, 200, 400):
        r = at_density(0.14835164835164835, 1 / 3, npts=npts)
        print(npts, r['K'], r['v'])
