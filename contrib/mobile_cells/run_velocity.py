#!/usr/bin/env python3
"""Momentum-resolved ring levels for z = 1 and the Casimir test.

For each ring length L and wall number W (S^z = 0): the lowest levels at P = 0 and at
P = 2 pi / L.  v_n(L) = L (E_n(P = 2pi/L) - E_0) / (2 pi) are the velocities of the two
gapless modes (their current excitations sit at P = 2 pi / L); z = 1 means they tend to
L-independent constants.
usage: run_velocity.py NAME --u U --mu MU --J J --LW 9:6 10:6 ... [--k 3]
Writes data/ed/velocity_NAME.json.  Claude (Opus 5.5), 2026-10-27."""
import os, sys, json, time, argparse
import numpy as np
from ed import Ring
ap = argparse.ArgumentParser()
ap.add_argument('name')
for k_, v_ in dict(t=1.0, u=1.0, mu=-1.0, J=0.5).items():
    ap.add_argument('--' + k_, type=float, default=v_)
ap.add_argument('--LW', nargs='+', required=True)
ap.add_argument('--k', type=int, default=3)
a = ap.parse_args()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ed', f'velocity_{a.name}.json')
res = json.load(open(out)) if os.path.exists(out) else []
done = {(r['L'], r['W']) for r in res}
for lw in a.LW:
    L, W = map(int, lw.split(':'))
    if (L, W) in done:
        continue
    t0 = time.time()
    r = Ring(L, W, W // 2, t=a.t, u=a.u, mu=a.mu, J=a.J)
    w0, P0 = r.lowest_momentum(0, k=2)
    w1, P1 = r.lowest_momentum(1, k=a.k)
    v = (L * (w1 - w0[0]) / (2 * np.pi)).tolist()
    rec = dict(L=L, W=W, dim=r.dim, E0=float(w0[0]), E_P0=w0.tolist(), E_P1=w1.tolist(), P1_check=P1.tolist(),
               v=v, time_s=time.time() - t0)
    res.append(rec)
    json.dump(res, open(out, 'w'), indent=1)
    print(f"L={L} W={W} dim={r.dim} E0={w0[0]:.8f} v_n={np.round(v, 4).tolist()} P1={np.round(P1, 3).tolist()} "
          f"({rec['time_s']:.0f}s)", flush=True)
