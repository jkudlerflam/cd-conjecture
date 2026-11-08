"""C160 (Claude Opus agent, 2026-11-08). Classical ground textures of random lattice bonds near h_F are
uniform or singly periodic (PROOF section 3.3, Mather + Poincare-Bendixson in the continuum).

Bond (est.py convention): J = -1 + eps A(D) + eps^2 S', B = eps^2 b, with D, S' (symmetric), b random,
i.e. the doubly rescaled wedge around the pure-DM direction where all terms compete.
Classical energy E = sum_j n_j^T J n_{j+1} + B . n_j on an open chain of N sites, minimized by L-BFGS
over angles from several starts (random, helices of several pitches). In the bulk (middle half) we take
the Fourier spectrum of n_j and test whether all spectral weight sits on harmonics m q1 of one wave vector
(periodic) or needs a second incommensurate q2 (quasi-periodic).
"""
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(7)


def A(D):
    return np.array([[0, D[2], -D[1]], [-D[2], 0, D[0]], [D[1], -D[0], 0]])


def energy_grad(p, J, B, N):
    th, ph = p[:N], p[N:]
    n = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    E = np.sum(np.einsum('ia,ab,ib->i', n[:-1], J, n[1:])) + np.sum(n @ B)
    G = np.zeros_like(n)
    G[:-1] += n[1:] @ J.T
    G[1:] += n[:-1] @ J
    G += B
    dth = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)], 1)
    dph = np.stack([-np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), 0 * th], 1)
    return E, np.r_[np.sum(G * dth, 1), np.sum(G * dph, 1)]


def ground(J, B, N, eps, starts=6):
    best = None
    for s in range(starts):
        if s < 2:
            p0 = np.r_[np.arccos(rng.uniform(-1, 1, N)), rng.uniform(0, 2 * np.pi, N)]
        else:
            q = eps * rng.uniform(0.3, 2.0)
            ax = rng.normal(size=3); ax /= np.linalg.norm(ax)
            # helix about a random axis
            e1 = np.cross(ax, [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0]); e1 /= np.linalg.norm(e1)
            e2 = np.cross(ax, e1)
            j = np.arange(N)
            n = np.outer(np.cos(q * j), e1) + np.outer(np.sin(q * j), e2)
            p0 = np.r_[np.arccos(np.clip(n[:, 2], -1, 1)), np.arctan2(n[:, 1], n[:, 0])]
        r = minimize(energy_grad, p0, args=(J, B, N), jac=True, method='L-BFGS-B',
                     options=dict(maxiter=40000, gtol=1e-12, maxcor=30))
        if best is None or r.fun < best[0]:
            best = (r.fun, r.x)
    th, ph = best[1][:N], best[1][N:]
    n = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    return best[0], n


def periodicity(n):
    N = len(n); b = n[N // 4: 3 * N // 4]
    M = len(b)
    win = np.hanning(M)[:, None]
    F = np.fft.rfft((b - b.mean(0)) * win, n=8 * M, axis=0)
    P = np.sum(np.abs(F) ** 2, 1)
    k = np.fft.rfftfreq(8 * M) * 2 * np.pi
    tot = P.sum()
    if tot < 1e-10 * M:
        return 'uniform', 0.0, 0.0
    q1 = k[np.argmax(P)]
    width = 2 * 2 * np.pi / M
    mask = np.zeros_like(P, bool)
    for m in range(1, 60):
        km = np.mod(m * q1 + np.pi, 2 * np.pi) - np.pi
        mask |= np.abs(k - abs(km)) < width
    return 'textured', q1, 1 - P[mask].sum() / tot


if __name__ == '__main__':
    N = 400
    for trial in range(int(__import__("sys").argv[1]) if len(__import__("sys").argv) > 1 else 10):
        eps = 0.12
        D = rng.normal(size=3); D /= np.linalg.norm(D)
        S = rng.normal(size=(3, 3)); S = (S + S.T) / 2; S -= np.trace(S) / 3 * np.eye(3)
        b = rng.normal(size=3)
        sc = rng.uniform(0.1, 0.8)
        J = -np.eye(3) + eps * A(D) + sc * eps ** 2 * S
        B = sc * eps ** 2 * b
        E, n = ground(J, B, N, eps)
        kind, q1, offw = periodicity(n)
        print(f'trial {trial}: scale {sc:.2f}, E/N = {E / N:+.6f}, {kind}, q1 = {q1:.4f}, weight off the harmonics of q1 = {offw:.2e}')
