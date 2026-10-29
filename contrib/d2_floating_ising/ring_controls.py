#!/usr/bin/env python3
"""Ring-detector controls: sJ, x0 (any q = 0 level) and x0n (neutral q = 0 level) for XX, XXZ (Delta = 0.5)
and the U(1) anchor-in-field point (g = 0, u = 0.4, D = 0.5, b = 1.05), all exactly one boson, at L = 16, 18, 20.
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29)."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'd2_classification'))
import numpy as np, ed, models_fi as MF
ctl = {'xx': (np.diag([1., 1, 0]), [0, 0, 0]), 'xxz05': (np.diag([1., 1, .5]), [0, 0, 0]),
       'u1_field': MF.family('Z2', dict(g=0, u=0.4, D=0.5, b=1.05))}
out = {}
for nm, (J, f) in ctl.items():
    h = ed.bond(J, f)
    for L in (16, 18, 20):
        E = [ed.lowest(ed.sector_hamiltonian(L, h, m), k=1)[0][0] for m in range(L)]; m0 = int(np.argmin(E))
        e, u = ed.lowest(ed.sector_hamiltonian(L, h, m0), k=4)
        sJ = L * (min(E[(m0 + 1) % L], E[(m0 - 1) % L]) - e[0]) / (2 * np.pi)
        pc = np.array([bin(i).count('1') for i in range(1 << L)]); M = L - 2 * pc
        Ms = [float(np.abs(ed.expand(L, m0, u[:, n])) ** 2 @ M) for n in range(4)]
        s = L * (e - e[0]) / (2 * np.pi)
        neut = [s[n] for n in range(1, 4) if abs(Ms[n] - Ms[0]) < 1]
        out[f'{nm}_L{L}'] = dict(sJ=float(sJ), x0=float(s[1] / sJ), x0n=float(neut[0] / sJ) if neut else None)
        print(nm, L, out[f'{nm}_L{L}'], flush=True)
json.dump(out, open(os.path.join(HERE, 'data', 'ring_controls.json'), 'w'), indent=1)
