#!/usr/bin/env python3
"""Ring ground states of the two-flavour payload (ed2.py): scan the sectors (W, n_a, n_b) near the
optimum, report the global ground state, its degeneracy, momentum and the gap to the next level.
usage: run_ed2.py NAME --L 8 9 ... [--u 0.3 --mu -0.4 --Ja 0.5 --Jb 0.35 --V 0.3 --hb 0.06]
Writes data/ed/c2_NAME_L<L>.json per ring.  Claude (Opus 5.5), 2026-10-27."""
import os, json, time, argparse
import numpy as np
from ed2 import Ring2
ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--L', type=int, nargs='+', required=True)
for k_, v_ in dict(t=1.0, u=0.3, mu=-0.4, Ja=0.5, Jb=0.35, V=0.3, hb=0.06, density=0.65).items():
    ap.add_argument('--' + k_, type=float, default=v_)
ap.add_argument('--maxdim', type=float, default=8e6)
a = ap.parse_args()
par = dict(t=a.t, u=a.u, mu=a.mu, Ja=a.Ja, Jb=a.Jb, V=a.V, ha=0.0, hb=a.hb)
HERE = os.path.dirname(os.path.abspath(__file__))
for L in a.L:
    t0 = time.time()
    wc = 2 * round(a.density * L / 2)
    rows = []
    for W in (wc - 2, wc, wc + 2):
        if W < 2 or W > L:
            continue
        for na in (W // 2 - 1, W // 2):
            for nb in (W // 2, W // 2 + 1, W // 2 + 2):
                if not (0 <= na <= W and 0 <= nb <= W):
                    continue
                r = Ring2(L, W, na, nb, **par)
                if r.dim > a.maxdim:
                    rows.append(dict(W=W, na=na, nb=nb, dim=r.dim, skipped=True)); continue
                w, V = r.lowest(k=2)
                P = r.momenta(w[:1], V[:, :1])
                rows.append(dict(W=W, na=na, nb=nb, dim=r.dim, E=w.tolist(), P0=float(P[0])))
    done = [x for x in rows if 'E' in x]
    levels = sorted([(e, x['W'], x['na'], x['nb']) for x in done for e in x['E']])
    E0 = levels[0][0]
    deg = sum(1 for l in levels if abs(l[0] - E0) < 1e-8)
    g = [x for x in done if abs(x['E'][0] - E0) < 1e-12][0]
    gap = levels[deg][0] - E0
    # is the winning sector strictly inside the scanned window?
    rec = dict(L=L, E0=E0, ground=dict(W=g['W'], na=g['na'], nb=g['nb'], P=g['P0']), degeneracy=deg,
               L_gap_over_2pi=L * gap / (2 * np.pi), next=levels[deg][1:], sectors=rows, params=par, time_s=time.time() - t0)
    json.dump(rec, open(os.path.join(HERE, 'data', 'ed', f'c2_{a.name}_L{L}.json'), 'w'), indent=1)
    print(f"L={L} E0={E0:.8f} ground (W,na,nb)=({g['W']},{g['na']},{g['nb']}) P={g['P0']} deg={deg} "
          f"L*gap/2pi={rec['L_gap_over_2pi']:.4f} next={rec['next']} skipped={[ (x['W'],x['na'],x['nb']) for x in rows if 'E' not in x]} "
          f"({rec['time_s']:.0f}s)", flush=True)
