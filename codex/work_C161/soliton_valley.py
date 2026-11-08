"""C161 (Claude Opus agent, 2026-11-08). Isolated soliton of the pure-DM-wedge functional (C160 eq. 2.2)
in the K stratum, and its tilt (out-of-plane) instability over a stable, Z2-symmetric uniform background.

    E[n] = int dx [ |n'|^2/2 - (n x n')_z + W(n) ],   W = -h n_x + (kappa/2) n_z^2 + (mu/2) n_y^2.

Z2: n_z -> -n_z (on spins Theta R_pi(z), the K-stratum antiunitary). Background +x is the global minimum of W
iff h + kappa > 0 and h + mu > 0.
For each (kappa, mu, h):
  - planar soliton phi(x) (n_z = 0), winding +2 pi, energy E_pl relative to the background;
  - lowest eigenvalue lam0 of the out-of-plane Hessian  -d^2 + U(x),  U = -phi'^2 + 2 phi' + h cos phi + kappa - mu sin^2 phi;
  - full 3-component minimization from tilted starts: energy E_t, max|n_z|, and the solid angle Omega enclosed by the
    loop (Haldane: soliton crystal momentum P = s Omega = Omega/2 mod 2 pi for s = 1/2; planar -> P = pi).
Lattice bond (est.py coordinates, eps small): v = (-1, -1 + mu eps^2/2, -1 + kappa eps^2/2, 0, 0, -eps, -h eps^2, 0, 0).
"""
import sys
import numpy as np
from scipy.optimize import minimize

X = 20.0
N = 500
dx = 2 * X / N
xs = -X + dx * np.arange(1, N)          # interior points, n fixed to +x at both ends


def planar(h, mu):
    def f(p):
        ph = np.r_[0.0, p, 2 * np.pi]
        d = np.diff(ph) / dx
        E = np.sum(0.5 * d ** 2) * dx - 2 * np.pi + np.sum(h * (1 - np.cos(p)) + 0.5 * mu * np.sin(p) ** 2) * dx
        g = (d[:-1] - d[1:]) + (h * np.sin(p) + mu * np.sin(p) * np.cos(p)) * dx
        return E, g
    p0 = np.pi * (1 + np.tanh(xs))
    r = minimize(f, p0, jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-12))
    return r.fun, r.x


def hz_lowest(h, kappa, mu, p):
    ph = np.r_[0.0, p, 2 * np.pi]
    d = np.diff(ph) / dx
    dp = 0.5 * (d[:-1] + d[1:])
    U = -dp ** 2 + 2 * dp + h * np.cos(p) + kappa - mu * np.sin(p) ** 2
    M = len(p)
    H = np.diag(2 / dx ** 2 + U) - np.diag(np.ones(M - 1) / dx ** 2, 1) - np.diag(np.ones(M - 1) / dx ** 2, -1)
    from scipy.linalg import eigh_tridiagonal
    w, v = eigh_tridiagonal(2 / dx ** 2 + U, -np.ones(M - 1) / dx ** 2, select='i', select_range=(0, 0))
    return w[0], v[:, 0]


def energy3(v, h, kappa, mu):
    V = v.reshape(-1, 3)
    nr = np.linalg.norm(V, axis=1)
    n = V / nr[:, None]
    ex = np.array([[1.0, 0, 0]])
    full = np.vstack([ex, n, ex])
    d = np.diff(full, axis=0) / dx
    a, b = full[:-1], full[1:]
    cz = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    W = -h * n[:, 0] + 0.5 * kappa * n[:, 2] ** 2 + 0.5 * mu * n[:, 1] ** 2 + h
    E = np.sum(0.5 * np.sum(d ** 2, 1)) * dx - np.sum(cz) + np.sum(W) * dx
    G = (d[:-1] - d[1:])                      # d/dn_i of sum |d|^2/2 dx
    # DM: - sum_i (a_i x b_i)_z ; n_i appears as b_{i-1} and a_i
    G[:, 0] += -(-a[:-1, 1]) - (b[1:, 1])
    G[:, 1] += -(a[:-1, 0]) - (-b[1:, 0])
    G[:, 0] += -h * dx
    G[:, 1] += mu * n[:, 1] * dx
    G[:, 2] += kappa * n[:, 2] * dx
    G = (G - np.sum(G * n, 1)[:, None] * n) / nr[:, None]
    return E, G.ravel()


def solid_angle(n):
    """oriented solid angle of the closed loop (+x -> texture -> +x), via the sum of signed spherical triangles
    with apex at the south-pole-avoiding reference -x... use apex -x is on the loop's far side; use apex e_ref."""
    full = np.vstack([[1.0, 0, 0], n, [1.0, 0, 0]])
    ref = np.array([0.0, 0.0, 1.0])
    # choose the apex as the pole away from the loop: compare mean n_z
    if full[:, 2].mean() > 0:
        ref = -ref
    tot = 0.0
    for i in range(len(full) - 1):
        a, b = full[i], full[i + 1]
        num = np.dot(ref, np.cross(a, b))
        den = 1 + np.dot(ref, a) + np.dot(ref, b) + np.dot(a, b)
        tot += 2 * np.arctan2(num, den)
    return tot   # signed area seen from apex 'ref'; planar equator gives +-2 pi


def full_min(h, kappa, mu, p, tilt):
    n0 = np.stack([np.cos(p), np.sin(p), tilt * np.exp(-xs ** 2 / 4)], 1)
    n0 /= np.linalg.norm(n0, axis=1)[:, None]
    r = minimize(energy3, n0.ravel(), args=(h, kappa, mu), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=40000, gtol=1e-10))
    V = r.x.reshape(-1, 3)
    n = V / np.linalg.norm(V, axis=1)[:, None]
    return r.fun, n


def report(h, kappa, mu):
    Epl, p = planar(h, mu)
    lam0, _ = hz_lowest(h, kappa, mu, p)
    best = None
    for tilt in (0.6, -0.6):
        Et, n = full_min(h, kappa, mu, p, tilt)
        Om = solid_angle(n)
        rec = (Et, np.abs(n[:, 2]).max(), np.sign(n[:, 2].sum()), Om)
        if best is None or Et < best[0]:
            best = rec
        out = rec if tilt > 0 else best
    Et, mz, sgn, Om = best
    return Epl, lam0, Et, mz, Om


if __name__ == '__main__':
    print('kappa   mu     h      E_planar   lam0(H_z)   E_full    max|n_z|   Omega/(2pi)  P/pi=Omega/(2pi) mod 2')
    grid = [(0.0, 0.0), (0.0, 2.0), (0.0, 4.0), (0.0, 8.0), (0.5, 8.0), (-0.2, 4.0)]
    for kappa, mu in grid:
        for h in (0.3, 0.6, 1.0, 1.5, 2.5):
            if h + kappa <= 0 or h + mu <= 0:
                continue
            Epl, lam0, Et, mz, Om = report(h, kappa, mu)
            print(f'{kappa:+.2f}  {mu:5.2f}  {h:4.2f}  {Epl:+9.4f}  {lam0:+9.4f}  {Et:+9.4f}  {mz:8.3f}   {Om/(2*np.pi):+8.4f}')
        sys.stdout.flush()
