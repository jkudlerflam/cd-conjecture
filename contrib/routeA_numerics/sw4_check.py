#!/usr/bin/env python3
"""Convergence check of the linked-cluster expansion: exact four-cell low band of the FULL chain
(48^4 = 5.3e6 states, open ends) against the logical model h2_12 + h2_23 + h2_34 + h3_123 + h3_234.
The difference is the connected four-cell term h4, which the logical chain omits.
usage: sw4_check.py LAM   (reads data/heff_c8_n3_lamLAM.npz; writes data/sw4_LAM.json)
Written by Claude (Opus 5.5) for contrib/routeA_numerics, 2026-10-27."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from cell import load_c8, swap
from sw_cluster import heff

lam = float(sys.argv[1])
t0 = time.time()
cell = load_c8()
z = np.load(os.path.join(HERE, 'data', f'heff_c8_n3_lam{sys.argv[1]}.npz'))
h2, h3 = z[f'h2_{lam}'], z[f'h3_{lam}']
I3 = np.eye(3)
model = (np.kron(np.kron(h2, I3), I3) + np.kron(np.kron(I3, h2), I3) + np.kron(np.kron(I3, I3), h2)
         + np.kron(h3, I3) + np.kron(I3, h3))
H4, info = heff(cell, 4, lam, tol=1e-12)
h4 = H4 - model
Pc = swap(3) - np.eye(9) / 3
b = np.trace(Pc @ h2).real / np.trace(Pc @ Pc).real
h4c = h4 - np.trace(h4) / 81 * np.eye(81)
e_full = np.linalg.eigvalsh(H4); e_mod = np.linalg.eigvalsh(model)
rec = dict(lam=lam, b=b, h4_norm_over_b=float(np.linalg.norm(h4c, 2) / b), h4_trace_per_cell_over_b=float(np.trace(h4).real / 81 / b),
           spec_maxdiff_over_b=float(np.abs((e_full - e_full[0]) - (e_mod - e_mod[0])).max() / b),
           bandwidth_over_b=float((e_full[-1] - e_full[0]) / b), info=info, time_s=time.time() - t0)
print(json.dumps(rec), flush=True)
np.savez(os.path.join(HERE, 'data', f'heff_c8_n4_lam{sys.argv[1]}.npz'), H4=H4, h4=h4)
json.dump(rec, open(os.path.join(HERE, 'data', f'sw4_{sys.argv[1]}.json'), 'w'), indent=1)
