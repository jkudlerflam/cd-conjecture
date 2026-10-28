#!/usr/bin/env python3
"""Velocities of the two gapless modes from their own excitations at q = 2 pi / L.

The wall-density operator rho_q = sum_i e^{iqi} n^wall_i couples to the length (charge) boson
and S^z_q = sum_i e^{iqi} S^z_i to the payload (spin) boson.  From the ground state |0>, a
Lanczos run started on O_q|0> gives the spectral measure of O_q; its dominant Ritz value
omega_O sets v_O = L omega_O / (2 pi).  z = 1 means v_O tends to a constant.
usage: run_modes.py NAME --u U --mu MU --J J --LW 9:6 10:6 ... [--steps 80]
Writes data/ed/modes_NAME.json.  Claude (Opus 5.5), 2026-10-27."""
import os, json, time, argparse
import numpy as np
from ed import Ring
ap = argparse.ArgumentParser()
ap.add_argument('name')
for k_, v_ in dict(t=1.0, u=1.0, mu=-1.0, J=0.5).items():
    ap.add_argument('--' + k_, type=float, default=v_)
ap.add_argument('--LW', nargs='+', required=True)
ap.add_argument('--steps', type=int, default=80)
a = ap.parse_args()


def lanczos_measure(Hmv, v0, m):
    nrm0 = np.linalg.norm(v0)
    q = v0 / nrm0; qold = np.zeros_like(q); beta = 0.0
    al, be = [], []
    Q = [q]
    for j in range(m):
        w = Hmv(q) - beta * qold
        alpha = np.vdot(q, w).real
        w = w - alpha * q
        for qq in Q:                       # full reorthogonalization (small m)
            w = w - np.vdot(qq, w) * qq
        al.append(alpha)
        beta = np.linalg.norm(w)
        if beta < 1e-12 or j == m - 1:
            break
        be.append(beta)
        qold, q = q, w / beta
        Q.append(q)
    T = np.diag(al) + np.diag(be[:len(al) - 1], 1) + np.diag(be[:len(al) - 1], -1)
    e, U = np.linalg.eigh(T)
    return e, (U[0] ** 2) * nrm0 ** 2


out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ed', f'modes_{a.name}.json')
res = json.load(open(out)) if os.path.exists(out) else []
done = {(r['L'], r['W']) for r in res}
for lw in a.LW:
    L, W = map(int, lw.split(':'))
    if (L, W) in done:
        continue
    t0 = time.time()
    r = Ring(L, W, W // 2, t=a.t, u=a.u, mu=a.mu, J=a.J)
    A = r.sparse() if r.dim <= 6e6 else None
    Hmv = (lambda v: A @ v.real + 1j * (A @ v.imag)) if A is not None else (lambda v: r.matvec(v.real) + 1j * r.matvec(v.imag))
    w, V = r.lowest(k=1)
    E0 = w[0]; g = V[:, 0].reshape(r.nc, r.ns)
    q = 2 * np.pi / L
    ph = np.exp(1j * q * r.cp)                                  # (nc, W)
    rho = (ph.sum(axis=1)[:, None] * g).ravel()
    bits = ((r.sb[None, :] >> np.arange(W)[:, None]) & 1) - 0.5  # (W, ns)
    sz = ((ph @ bits) * g).ravel()
    rec = dict(L=L, W=W, dim=r.dim, E0=float(E0))
    for name, v0 in (('rho', rho), ('Sz', sz)):
        e, wt = lanczos_measure(lambda v: Hmv(v) - E0 * v, v0, a.steps)
        k = int(np.argmax(wt))
        rec[name] = dict(omega=e.tolist(), weight=wt.tolist(), peak=float(e[k]), peak_weight=float(wt[k] / wt.sum()),
                         v=float(L * e[k] / (2 * np.pi)))
    res.append(rec)
    json.dump(res, open(out, 'w'), indent=1)
    print(f"L={L} W={W} dim={r.dim} v_rho={rec['rho']['v']:.4f} (w {rec['rho']['peak_weight']:.2f}) "
          f"v_Sz={rec['Sz']['v']:.4f} (w {rec['Sz']['peak_weight']:.2f}) ({time.time()-t0:.0f}s)", flush=True)
