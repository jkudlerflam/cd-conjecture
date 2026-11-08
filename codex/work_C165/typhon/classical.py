#!/usr/bin/env python3
"""C165 (Claude Opus agent, 2026-11-08): classical ground states of a nearest-neighbour bond
E = sum_j [ n_j^T J n_{j+1} + B . n_j ],  |n_j| = 1  (classical image of the Pauli bond of B_2).

Tools: energy and gradient on rings (pbc) or open chains, multi-start L-BFGS minimisation, Fourier
diagnostics (uniform M, staggered N, largest incommensurate peak), and the tangent-space Hessian.

Model family used throughout (parent axis z, see PROOF section 1):
  J = [[-cos a + g, sin a + gxy, Gxz], [-sin a + gxy, -cos a - g, Gyz], [Gxz', Gyz', Dz]]   (built in fam())
  B = (hx, hy, hz)
"""
import numpy as np
from scipy.optimize import minimize


def fam(a=0.8, Dz=1.0, g=0.0, gxy=0.0, G=0.0, Gp=None, Dx=0.0, Dy=0.0, hx=0.0, hy=0.0, hz=0.0):
    """XXZ + DM_z parent (in-plane FM after the rotating frame, |J_perp_eff| = 1, pitch q = -a) with
    charge-2 harmonic g (XX - YY), gxy (XY + YX), charge-1 symmetric harmonic G (XZ + ZX) [Gp for an
    asymmetric ZX coefficient], charge-1 DM (Dx, Dy) and fields."""
    if Gp is None:
        Gp = G
    J = np.array([[-np.cos(a) + g, np.sin(a) + gxy, G],
                  [-np.sin(a) + gxy, -np.cos(a) - g, Dx],
                  [Gp, -Dx, Dz]], float)
    J[0, 2] += -Dy; J[2, 0] += Dy
    return J, np.array([hx, hy, hz], float)


def unpack(x, L):
    v = x.reshape(L, 3)
    r = np.linalg.norm(v, axis=1)
    return v / r[:, None], v, r


def egrad(x, J, B, L, pbc=True):
    n, v, r = unpack(x, L)
    if pbc:
        nr = np.roll(n, -1, axis=0); nl = np.roll(n, 1, axis=0)
        E = np.einsum("ia,ab,ib->", n, J, nr) + np.sum(n @ B)
        g = nr @ J.T + nl @ J + B
    else:
        E = np.einsum("ia,ab,ib->", n[:-1], J, n[1:]) + np.sum(n @ B)
        g = np.zeros_like(n) + B
        g[:-1] += n[1:] @ J.T
        g[1:] += n[:-1] @ J
    # chain rule through n = v/|v|
    gv = (g - n * np.sum(g * n, axis=1)[:, None]) / r[:, None]
    return E, gv.ravel()


def seeds(L, rng, q0=None, nrand=6):
    out = []
    j = np.arange(L)
    for _ in range(nrand):
        out.append(rng.standard_normal(3 * L))
    qs = [q0] if q0 is not None else []
    qs += list(rng.uniform(0, 2 * np.pi, 3))
    for q in qs:
        for th in (np.pi / 2, 1.2, 0.8, 0.4):
            n = np.stack([np.sin(th) * np.cos(q * j), np.sin(th) * np.sin(q * j), (-1) ** j * np.cos(th)], 1)
            out.append((n + 0.01 * rng.standard_normal(n.shape)).ravel())
            n2 = np.stack([np.sin(th) * np.cos(q * j), np.sin(th) * np.sin(q * j), np.cos(th) + 0 * j], 1)
            out.append((n2 + 0.01 * rng.standard_normal(n.shape)).ravel())
    n = np.stack([0 * j, 0 * j, (-1.0) ** j], 1)
    out.append((n + 0.01 * rng.standard_normal(n.shape)).ravel())
    return out


def ground(J, B, L, pbc=True, seed=0, q0=None, nrand=6, gtol=1e-10, extra=()):
    rng = np.random.default_rng(seed)
    best = None
    for x0 in list(seeds(L, rng, q0, nrand)) + list(extra):
        res = minimize(egrad, x0, args=(J, B, L, pbc), jac=True, method='L-BFGS-B',
                       options=dict(maxiter=20000, gtol=gtol, ftol=1e-15, maxcor=30))
        if best is None or res.fun < best.fun - 1e-12:
            best = res
    n, _, _ = unpack(best.x, L)
    return best.fun / L, n


def diag(n, cut=0):
    """Uniform M, staggered N (vectors), incommensurate peak (k, amplitude) of |FT n|, over the bulk."""
    L = len(n)
    m = n[cut:L - cut] if cut else n
    Lb = len(m); j = np.arange(Lb)
    M = m.mean(0)
    N = ((-1.0) ** j[:, None] * m).mean(0)
    F = np.fft.fft(m - M - ((-1.0) ** j)[:, None] * N, axis=0) / Lb
    amp = np.sqrt(np.sum(np.abs(F) ** 2, axis=1))
    k = 2 * np.pi * np.fft.fftfreq(Lb)
    i = np.argmax(amp)
    return dict(M=M, N=N, absN=float(np.linalg.norm(N)), absM=float(np.linalg.norm(M)),
                kpk=float(abs(k[i])), apk=float(amp[i] * np.sqrt(2)))


def hessian(n, J, B, pbc=True):
    """Tangent-space Hessian (2L x 2L) of E at a stationary configuration n."""
    L = len(n)
    # tangent frame
    e1 = np.zeros_like(n); e2 = np.zeros_like(n)
    for i in range(L):
        a = np.array([1.0, 0, 0]) if abs(n[i, 0]) < 0.9 else np.array([0, 1.0, 0])
        u = a - n[i] * (a @ n[i]); u /= np.linalg.norm(u)
        e1[i] = u; e2[i] = np.cross(n[i], u)
    T = np.stack([e1, e2], 1)  # L x 2 x 3
    nr = np.roll(n, -1, 0) if pbc else np.vstack([n[1:], np.zeros(3)])
    nl = np.roll(n, 1, 0) if pbc else np.vstack([np.zeros(3), n[:-1]])
    h = nr @ J.T + nl @ J + B           # local field (gradient)
    lam = np.sum(h * n, 1)               # Lagrange multiplier
    H = np.zeros((2 * L, 2 * L))
    for i in range(L):
        H[2 * i:2 * i + 2, 2 * i:2 * i + 2] -= lam[i] * np.eye(2)
        jn = (i + 1) % L
        if pbc or i < L - 1:
            blk = T[i] @ J @ T[jn].T
            H[2 * i:2 * i + 2, 2 * jn:2 * jn + 2] += blk
            H[2 * jn:2 * jn + 2, 2 * i:2 * i + 2] += blk.T
    return H, T


if __name__ == '__main__':
    import sys
    # quick check: gradient and the rigid-cone linearity (PROOF Prop. 1)
    rng = np.random.default_rng(1)
    J, B = fam(a=0.8, Dz=1.1, g=0.3, G=0.2, Dx=0.1, hx=0.2, hz=0.3)
    L = 40
    x = rng.standard_normal(3 * L)
    E, g = egrad(x, J, B, L)
    d = rng.standard_normal(3 * L) * 1e-6
    print('grad check', (egrad(x + d, J, B, L)[0] - egrad(x - d, J, B, L)[0]) / 2, g @ d)


def nfloat(eta=0.1, b=1.05):
    """C68 symmetry-free floating family (contrib/d2_floating_ising cut N): S = diag(1+eta, 1-eta, 0.4),
    D = (eta, 2 eta, 0.5) in the est.J_f convention, f = b (3, 4, 5)/sqrt(50)."""
    J = np.diag([1 + eta, 1 - eta, 0.4])
    Dx, Dy, Dz = eta, 2 * eta, 0.5
    J[0, 1] += Dz; J[1, 0] -= Dz; J[1, 2] += Dx; J[2, 1] -= Dx; J[2, 0] += Dy; J[0, 2] -= Dy
    return J, b * np.array([3., 4., 5.]) / np.sqrt(50)
