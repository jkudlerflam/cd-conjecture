#!/usr/bin/env python3
"""Characterize the q = 0 first excited ring level against the ground state: Z2 parity P = prod Z,
total magnetization M = sum_j <Z_j>, and s = L dE/2pi.  A moving charge zero mode changes M by about
+-2 (one particle) or +-4 (a pair, same parity); a neutral Ising energy level has dM ~ 0 and P = +1.
usage: levelchar.py TAG [TAG ...] --L 16 18   (TAGs from data/ring_summary.json)
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'd2_classification'))
import numpy as np
import ed, models_fi as MF
ap = argparse.ArgumentParser(); ap.add_argument('tags', nargs='+'); ap.add_argument('--L', type=int, nargs='+', default=[16])
a = ap.parse_args()
S = json.load(open(os.path.join(HERE, 'data', 'ring_summary.json')))
out = {}
for t in a.tags:
    J, f = MF.family(S[t]['fam'], S[t]['p']); h = ed.bond(J, f)
    for L in a.L:
        m0 = S[t][f'L{L}']['m0'] if f'L{L}' in S[t] else None
        if m0 is None:
            E = [ed.lowest(ed.sector_hamiltonian(L, h, m), k=1)[0][0] for m in range(L)]; m0 = int(np.argmin(E))
        e, u = ed.lowest(ed.sector_hamiltonian(L, h, m0), k=3)
        idx = np.arange(1 << L); pc = np.array([bin(i).count('1') for i in idx])
        P = (-1.0) ** pc; M = L - 2 * pc
        res = []
        for n in range(3):
            psi = ed.expand(L, m0, u[:, n]); w = np.abs(psi) ** 2
            res.append(dict(s=float(L * (e[n] - e[0]) / (2 * np.pi)), P=float(w @ P), M=float(w @ M)))
        out[f'{t}_L{L}'] = res
        print(t, L, 'm0', m0, ' '.join(f"[s={r['s']:.3f} P={r['P']:+.2f} M={r['M']:+.2f}]" for r in res), flush=True)
json.dump(out, open(os.path.join(HERE, 'data', 'levelchar.json'), 'w'), indent=1)
