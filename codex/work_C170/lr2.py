#!/usr/bin/env python3
"""C170 (Claude Opus 5.5 referee, independent code): exact second-order (linear-response) relaxation
energy of the classical uniform nu-frame state at the parent switch Dz = 1, hz = 0, for an arbitrary
harmonic (dJ, dB) added to the U(1) parent.  Ring of L sites, pitch a with a*L in 2 pi Z.
sigma-frame configuration of the nu-frame uniform state L:  n_j = Rot_j L,
Rot_j = R_z(q j) diag(1, 1, (-1)^j), q = -a.  Parent energy at the switch is -L per ring for every L.
dE2(L) = -(1/2) F^T K^+ F with K the tangent Hessian of the parent at the uniform state.
"""
import numpy as np

def Rz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])

def parent(a, Dz=1.0):
    return np.array([[-np.cos(a), np.sin(a), 0], [-np.sin(a), -np.cos(a), 0], [0, 0, Dz]])

def frames(L, a):
    q = -a
    return [Rz(q * j) @ np.diag([1, 1, (-1) ** j]) for j in range(L)]

def tangent(Lv):
    Lv = Lv / np.linalg.norm(Lv)
    t = np.array([1., 0, 0]) if abs(Lv[0]) < 0.9 else np.array([0, 1., 0])
    e1 = t - Lv * (t @ Lv); e1 /= np.linalg.norm(e1); e2 = np.cross(Lv, e1)
    return Lv, e1, e2

def energy(n, J, B):
    return np.einsum('ia,ab,ib->', n, J, np.roll(n, -1, 0)) + np.sum(n @ B)

def dE2(Lv, dJ, dB, L=100, a=2*np.pi*13/100, Dz=1.0, check=False):
    J0 = parent(a, Dz); R = frames(L, a)
    Lv, e1, e2 = tangent(np.asarray(Lv, float))
    # sigma-frame tangent vectors at each site
    T = np.array([[R[j] @ e1, R[j] @ e2] for j in range(L)])   # L x 2 x 3
    n0 = np.array([R[j] @ Lv for j in range(L)])
    # first-order check: harmonic energy on the uniform state must vanish
    e1st = energy(n0, dJ, dB)
    # force of the harmonic, projected on the tangent spaces
    g = np.roll(n0, -1, 0) @ dJ.T + np.roll(n0, 1, 0) @ dJ + dB      # dE/dn_j
    F = np.einsum('jkа,ja->jk'.replace('а','a'), T, g).ravel()
    # Hessian of the parent on the sphere at n0 (exact, numerical from second derivatives)
    gp = np.roll(n0, -1, 0) @ J0.T + np.roll(n0, 1, 0) @ J0
    lam = np.sum(gp * n0, 1)       # Lagrange multipliers
    K = np.zeros((2 * L, 2 * L))
    for j in range(L):
        k = (j + 1) % L
        blk = T[j] @ J0 @ T[k].T     # coupling n_j^T J n_k
        K[2*j:2*j+2, 2*k:2*k+2] += blk
        K[2*k:2*k+2, 2*j:2*j+2] += blk.T
        K[2*j:2*j+2, 2*j:2*j+2] -= lam[j] * np.eye(2)
    w, V = np.linalg.eigh(K)
    keep = w > 1e-9
    Fk = V.T @ F
    if check:
        return w[:4], e1st
    return -0.5 * np.sum(Fk[keep] ** 2 / w[keep]) / L, e1st, np.sum(np.abs(Fk[~keep]))

def harm(g=0, gxy=0, G=0, Gp=None, Dx=0, Dy=0, hx=0, hy=0, hz=0):
    if Gp is None: Gp = G
    dJ = np.array([[g, gxy, G], [gxy, -g, Dx], [Gp, -Dx, 0.]])
    dJ[0, 2] += -Dy; dJ[2, 0] += Dy
    return dJ, np.array([hx, hy, hz], float)

if __name__ == '__main__':
    a = 2 * np.pi * 13 / 100; q = -a; L = 100
    print('Hessian low modes, first-order energy:', dE2([0.3, 0.2, 0.9], *harm(g=0.1), check=True))
    us = np.linspace(0, 1, 11)
    cases = {'g': dict(g=0.1), 'gxy': dict(gxy=0.1), 'G (sym)': dict(G=0.1), 'Dx': dict(Dx=0.1),
             'Dy': dict(Dy=0.1), 'G+Dx': dict(G=0.1, Dx=0.1), 'G+Dy': dict(G=0.1, Dy=0.1),
             'G only xz (Gp=0)': dict(G=0.1, Gp=0.0), 'hx': dict(hx=0.1), 'hz': dict(hz=0.1),
             'g+hx+G+Dx': dict(g=0.1, hx=0.1, G=0.1, Dx=0.1)}
    cq = np.cos(q)
    for name, kw in cases.items():
        rows = []
        for u in us:
            for chi in (0.0, 0.7):
                Lv = [np.sqrt(1 - u) * np.cos(chi), np.sqrt(1 - u) * np.sin(chi), np.sqrt(u)]
                e, e1, zm = dE2(Lv, *harm(**kw), L=L, a=a)
                rows.append((u, chi, e, e1, zm))
        r = np.array(rows)
        # chi-independence and polynomial fit in u
        spread = max(abs(r[r[:, 1] == 0, 2] - r[r[:, 1] == 0.7, 2]))
        p = np.polyfit(r[r[:, 1] == 0, 0], r[r[:, 1] == 0, 2], 3)
        print(f'{name:18s} dE2(u) = {p[3]:+.6e} {p[2]:+.6e} u {p[1]:+.6e} u^2 {p[0]:+.1e} u^3 | chi-spread {spread:.1e} | max|E1| {np.max(abs(r[:,3])):.1e}')
    print('predictions (C165): c2 = g^2 cot^2 q/4 =', 0.01 / np.tan(q) ** 2 / 4,
          ' c1 = G^2 (1-cos q)^2/(2(1+cos q)) =', 0.01 * (1 - cq) ** 2 / (2 * (1 + cq)))
    print('referee prediction for DM: u^2 coefficient -Dx^2 sin^2(q/2) =', -0.01 * np.sin(q / 2) ** 2,
          ' and dE2 = -(Dx^2 (1-cos q)/2) (1-u)^2, constant term', -0.01 * (1 - cq) / 2)
