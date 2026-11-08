#!/usr/bin/env python3
"""Classical two-sublattice ground state and linear spin-wave theory, rotated frame,
H' = sum_j 4 S_j^T diag(1,1,-u) S_j+1 + 2b sum_j (-1)^j S^z_j, spin 1/2 (S = sigma/2), cell = (A even, B odd).
Prints the classical phase, canting, acoustic velocity and the optical gap.
Written for codex/work_C153 (Claude Opus agent, 2026-11-07)."""
import numpy as np
from scipy.optimize import minimize
S = 0.5


def vec(t, p):
    return np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])


def classical(u, b):
    J = 4 * np.diag([1, 1, -u]); hA = np.array([0, 0, 2 * b]); hB = -hA
    def E(x):
        vA, vB = vec(x[0], 0.0), vec(x[1], x[2])
        return S * S * 2 * vA @ J @ vB + S * (hA @ vA + hB @ vB)
    best = None
    rng = np.random.default_rng(0)
    for _ in range(40):
        x0 = rng.uniform(0, np.pi, 3) * np.array([1, 1, 2])
        res = minimize(E, x0, method='BFGS', options=dict(gtol=1e-12))
        if best is None or res.fun < best.fun - 1e-12:
            best = res
    return best.x, best.fun / 2  # energy per site


def frame(t, p):
    v = vec(t, p)
    ex = np.array([np.cos(t) * np.cos(p), np.cos(t) * np.sin(p), -np.sin(t)])
    ey = np.cross(v, ex)
    return v, ex + 1j * ey


def lswt(u, b, ks):
    x, e = classical(u, b)
    J = 4 * np.diag([1, 1, -u]); h = [np.array([0, 0, 2 * b]), np.array([0, 0, -2 * b])]
    fr = [frame(x[0], 0.0), frame(x[1], x[2])]
    # bonds: (s, s', delta in cells): A->B same cell (0), B->A next cell (+1)
    bonds = [(0, 1, 0), (1, 0, 1)]
    out = []
    for k in ks:
        M = np.zeros((2, 2), complex); D = np.zeros((2, 2), complex)
        for s in range(2):
            M[s, s] += -(h[s] @ fr[s][0])
        for (s, t, dl) in bonds:
            vs, us = fr[s]; vt, ut = fr[t]
            M[s, s] += -S * vs @ J @ vt; M[t, t] += -S * vs @ J @ vt
            tt = (S / 2) * us @ J @ np.conj(ut)
            pp = (S / 2) * us @ J @ ut
            ph = np.exp(1j * k * dl)
            M[s, t] += tt * ph; M[t, s] += np.conj(tt * ph)
            D[s, t] += pp * ph; D[t, s] += pp * np.conj(ph)
        Mm = np.zeros((2, 2), complex)  # M(-k)
        for s in range(2):
            Mm[s, s] += -(h[s] @ fr[s][0])
        for (s, t, dl) in bonds:
            vs, us = fr[s]; vt, ut = fr[t]
            Mm[s, s] += -S * vs @ J @ vt; Mm[t, t] += -S * vs @ J @ vt
            tt = (S / 2) * us @ J @ np.conj(ut); ph = np.exp(-1j * k * dl)
            Mm[s, t] += tt * ph; Mm[t, s] += np.conj(tt * ph)
        H = np.block([[M, D], [np.conj(D.T), Mm.T]])
        eta = np.diag([1, 1, -1, -1])
        w = np.linalg.eigvals(eta @ H)
        w = np.sort(np.real_if_close(w))
        out.append(w)
    return x, e, np.array(out)


if __name__ == '__main__':
    ks = np.concatenate([[0.01], np.linspace(-np.pi, np.pi, 401)])
    print('# u b | thetaA thetaB phi | e_cl | omega_ac(k->0)/k  min optical  max|Im|')
    for b in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]:
        r = np.sqrt(1 + b * b / 4)
        for u in [b / 2 - 1.2, b / 2 - 0.9, r - 0.5, r - 0.2, r - 0.05, r - 0.01, r - 1e-3, r + 0.05]:
            if u <= -1: continue
            x, e, w = lswt(u, b, ks)
            pos = np.sort(np.abs(np.real(w)), axis=1)[:, ::2]  # two branches (each appears as +-omega)
            im = np.max(np.abs(np.imag(w)))
            i0 = 0
            print(f'{u:8.4f} {b:4.1f} | {x[0]:.4f} {x[1]:.4f} {x[2] % (2*np.pi):.4f} | {e:.6f} | '
                  f'{pos[i0,0]/0.01/2:.4f} {pos[:,1].min():.4f} {im:.1e}', flush=True)
