"""C164 calibration of the lattice classical code (clat.py) against the continuum blow-up results.
(1) C161: K bond v = (-1, -1 + mu e^2/2, -1 + kappa e^2/2, 0, 0, -e, -h e^2, 0, 0), kappa = 0.05, mu = 2.3.
    Lattice planar 2 pi soliton, its energy (lattice E = e * E_resc) and the lowest tangent-Hessian eigenvalue,
    converted with H_lat = e^2 H_resc. Continuum: h_c = 0.0238, lam0 = -0.067.
(2) C162 lemma_check (psi = 0.3, field rotated in plane): rings of length 2T with two solitons, same or opposite tilt.
    Continuum at h = 0.02886: E_sol = +0.0116, U_same(10) = -0.0088, U_alt(10) = +0.0169.
usage: calib_classical.py part eps
"""
import sys, time
import numpy as np
from scipy.optimize import minimize, brentq
import clat

kappa, mu = 0.05, 2.3


def bond(e, h, psi=0.0):
    J = np.diag([-1.0, -1 + mu * e * e / 2, -1 + kappa * e * e / 2])
    J[0, 1] += -e; J[1, 0] -= -e            # Dz = -e
    f = -h * e * e * np.array([np.cos(psi), np.sin(psi), 0.0])
    return J, f


def planar_lattice(J, f, bg, N):
    """planar (n_z = 0) 2 pi soliton on an open chain of N free sites, bg clamped at both ends."""
    pb = np.arctan2(bg[1], bg[0])
    Jxy = J[:2, :2]; fxy = f[:2]
    ref = (N + 1) * (bg[:2] @ Jxy @ bg[:2]) + N * (fxy @ bg[:2])
    def fun(p):
        ph = np.r_[pb, p, pb + 2 * np.pi]
        n = np.stack([np.cos(ph), np.sin(ph)], 1)
        E = np.einsum('ia,ab,ib->', n[:-1], Jxy, n[1:]) + np.sum(n[1:-1] @ fxy) - ref
        G = n[2:] @ Jxy.T + n[:-2] @ Jxy + fxy        # dE/dn_i
        t = np.stack([-np.sin(ph[1:-1]), np.cos(ph[1:-1])], 1)
        return E, np.sum(G * t, 1)
    j = np.arange(1, N + 1)
    best = None
    for w in (2.0, 5.0, 12.0):
        p0 = pb + np.pi * (1 + np.tanh((j - (N + 1) / 2) / w))
        r = minimize(fun, p0, jac=True, method='L-BFGS-B', options=dict(maxiter=50000, gtol=1e-12, ftol=1e-16))
        if best is None or r.fun < best[0]:
            best = (r.fun, r.x)
    p = best[1]
    n = np.stack([np.cos(p), np.sin(p), 0 * p], 1)
    return best[0], n


def part1(e):
    N = int(round(24 / e))
    print(f'# part 1, eps = {e}, N = {N} free sites (continuum length 24)')
    def Epl(h):
        J, f = bond(e, h); bg = np.array([1.0, 0, 0])
        return planar_lattice(J, f, bg, N)[0] / e
    hc = brentq(Epl, 0.002, 0.2, xtol=1e-6)
    for h in (0.0238, hc, 0.03):
        J, f = bond(e, h); bg = np.array([1.0, 0, 0])
        E, n = planar_lattice(J, f, bg, N)
        ch = clat.Chain(J, f, bg, N)
        H, T = ch.hessian_tangent(n)
        w, V = np.linalg.eigh(H)
        # out-of-plane weight of each eigenvector
        zc = np.array([np.sum((T[:, 2, :] * V[:, k].reshape(N, 2)).sum(1) ** 2) for k in range(4)])
        print(f'h={h:.5f}  E_pl/eps={E / e:+.5f}  lowest Hessian eigenvalues / eps^2 = {np.round(w[:3] / e**2, 4)}  '
              f'out-of-plane weight {np.round(zc[:3], 3)}' + ('   <- lattice planar onset' if h == hc else ''))


def ring_min(J, f, bg, L, centres, tilts, width):
    j = np.arange(L)
    ph = np.zeros(L); nz = np.zeros(L)
    for c, t in zip(centres, tilts):
        ph += np.pi * (1 + np.tanh((j - c) / width))
        nz += t * np.exp(-((j - c) / (2 * width)) ** 2)
    pb = np.arctan2(bg[1], bg[0])
    x0 = np.stack([np.cos(ph + pb), np.sin(ph + pb), nz], 1)
    r = minimize(clat.ring_energy_grad, x0.ravel(), args=(J, f), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=100000, gtol=1e-11, ftol=1e-16))
    V = r.x.reshape(-1, 3); n = V / np.linalg.norm(V, axis=1)[:, None]
    e_bg = bg @ J @ bg + f @ bg
    return r.fun - L * e_bg, n


def part2(e, hs):
    psi = 0.3
    for h in hs:
        J, f = bond(e, h, psi)
        bg, ebg = clat.uniform_min(J, f)
        width = 1.5 / e
        L1 = int(round(60 / e))
        E1, n1 = ring_min(J, f, bg, L1, [L1 / 2], [0.6], width)
        print(f'# part 2, eps={e}, psi={psi}, h={h:.5f}: bg={np.round(bg, 4)}, isolated tilted soliton E/eps = {E1 / e:+.6f}, '
              f'max n_z = {n1[:, 2].max():.3f}, min n_z = {n1[:, 2].min():+.2e}')
        print('T      U_same/eps     U_alt/eps    (continuum C162 at h=0.02886: T=10 -0.0088 / +0.0169)')
        for T in (10.0, 15.0, 22.0):
            L = int(round(2 * T / e))
            Es, ns = ring_min(J, f, bg, L, [L / 4, 3 * L / 4], [0.6, 0.6], width)
            Ea, na = ring_min(J, f, bg, L, [L / 4, 3 * L / 4], [0.6, -0.6], width)
            ok = na[:, 2].max() > 0.05 and na[:, 2].min() < -0.05
            print(f'{T:5.1f}  {(Es / 2 - E1) / e:+.6e}  {(Ea / 2 - E1) / e:+.6e}  {"" if ok else "(alt relaxed to same)"}', flush=True)


if __name__ == '__main__':
    t0 = time.time()
    part = sys.argv[1]; e = float(sys.argv[2])
    if part == '1':
        part1(e)
    else:
        hs = [float(x) for x in sys.argv[3:]] or [0.02886]
        part2(e, hs)
    print(f'# {time.time() - t0:.0f} s')
