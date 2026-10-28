#!/usr/bin/env python3
"""Velocities from dynamical structure factors on exact rings (sqed.py), as in ../mobile_cells/run_modes.py.

Ground state |0> in momentum sector 0; probe O = sum_x cos(2 pi x / L) o_x with
  o = n   (occupied-site density: couples to the charge boson),
  o = q   (logical U(1) charge density: couples to the payload, glued to a charge twist).
O|0> is projected on the P = +-2 pi/L sector and a Lanczos run from it gives Ritz energies and spectral
weights; the dominant low peak w gives v = L w / (2 pi).  Also reports the twist level of a period-3 word
(physical momentum 2 pi M / (3 L)) and E0.
usage: modes.py --L 10 --word BPPBPP --Q2 0 --J 0.25 --mu_v 0.6 [--k 120] --out f.json
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
import numpy as np
from scipy.linalg import eigh_tridiagonal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sqed  # noqa: E402
import vacbond as VB  # noqa: E402


def lanczos_spectral(ring, v0, m, k, shift, reorth=True):
    """Ritz values and weights |<psi_i|v0>|^2 of the projected Hamiltonian from start vector v0.
    reorth=False keeps two vectors only (spurious copies of converged Ritz values can appear; the
    low-lying weights are unaffected)."""
    def Hp(v):
        return ring.proj(ring.H(v), m)
    nrm0 = np.linalg.norm(v0)
    q = v0 / nrm0
    qs = [q]
    a, b = [], []
    qprev = np.zeros_like(q)
    beta = 0.0
    for j in range(k):
        w = Hp(q) - beta * qprev
        alpha = float(np.vdot(q, w).real)
        w = w - alpha * q
        if reorth:
            for qq in qs:                   # full reorthogonalization
                w = w - np.vdot(qq, w) * qq
        a.append(alpha)
        beta = float(np.linalg.norm(w))
        if beta < 1e-10 or j == k - 1:
            break
        b.append(beta)
        qprev, q = q, w / beta
        if reorth:
            qs.append(q)
    ev, U = eigh_tridiagonal(np.array(a), np.array(b[:len(a) - 1]))
    wts = (U[0] ** 2) * nrm0 ** 2
    return ev, wts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--word', required=True)
    ap.add_argument('--Q2', type=int, default=0)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.25)
    ap.add_argument('--mu_v', type=float, default=0.6)
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--k', type=int, default=120)
    ap.add_argument('--out', default=None)
    ap.add_argument('--noreorth', action='store_true')
    a = ap.parse_args()
    t0 = time.time()
    ring = sqed.Ring(a.L, a.word, a.Q2, a.t, a.J, a.mu_v, a.lam)
    w, U = ring.lowest(0, nev=2, tol=1e-12, vecs=True)
    E0, g = w[0], U[:, 0]
    L, M = a.L, len(a.word)
    x = ring.pos.confs                                   # (C, M) positions
    fx = np.cos(2 * np.pi * x / L)                       # (C, M)
    res = dict(args=vars(a), E0=float(E0), E1_sector0=float(w[1]), dims=[ring.D, ring.C], probes={})
    G = g.reshape(ring.shape)
    # o = n: diagonal in positions, identity on labels
    vn = (G * fx.sum(1)[None, :]).ravel()
    # o = q: logical charge of the label at squeezed site j times cos(2 pi x_j / L)
    qlab = np.zeros((ring.D, M))
    for b in ring.inn.blocks:
        if b['n'] == 0:
            continue
        for j in range(M):
            lets = VB.LETTERS[b['types'][j]][b['dig'][j]]
            qlab[b['off']:b['off'] + b['n'], j] = VB.Q2_LETTER[lets] / 2.0
    vq = np.einsum('aj,cj,ac->ac', qlab, fx, G).ravel()
    for name, v in (('n', vn), ('q', vq)):
        v1 = ring.proj(v, 1)
        if np.linalg.norm(v1) < 1e-12:
            continue
        ev, wt = lanczos_spectral(ring, v1, 1, a.k, 0.0, reorth=not a.noreorth)
        om = ev - E0
        order = np.argsort(om)
        peaks = [dict(omega=float(om[i]), weight=float(wt[i]), v=float(L * om[i] / (2 * np.pi)))
                 for i in order if wt[i] > 1e-3 * wt.sum()][:12]
        top = max(peaks, key=lambda p: p['weight'])
        res['probes'][name] = dict(norm2=float(np.linalg.norm(v1) ** 2), peaks=peaks, dominant=top)
        print(f"probe {name}: dominant omega={top['omega']:.6f} v={top['v']:.5f} weight frac="
              f"{top['weight'] / wt.sum():.3f}; peaks {[(round(p['omega'], 5), round(p['weight'] / wt.sum(), 3)) for p in peaks[:6]]}",
              flush=True)
    if sqed.period(a.word) == 3:
        mt = M // 3
        mt = min(mt, L - mt)
        et = ring.lowest(mt, nev=2, tol=1e-12)
        res['twist'] = dict(m=int(mt), gap=float(et[0] - E0))
        print(f'twist level at m={mt}: gap {et[0] - E0:.6f}', flush=True)
    res['time_s'] = time.time() - t0
    print(f'E0={E0:.12f}  t={res["time_s"]:.0f}s', flush=True)
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
