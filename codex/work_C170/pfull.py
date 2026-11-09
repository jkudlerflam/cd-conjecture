#!/usr/bin/env python3
"""C170: open-chain ED of C165's plateau-edge bond (a = 0.9, Dz = 2, g = 0.3, hx = 0.2, Dx = 0.03, hz)
bulk energy per site e ~ [E(L) - E(L-4)]/4, against C165's VUMPS energies (loc/P_hz*.json)."""
import numpy as np, json, sys
sys.path.insert(0, '.')
from qed import ham, low
J = np.array(json.load(open('../work_C165/loc/P_hz1.0.json'))['args']['J']).reshape(3, 3)
print('J =', J.tolist())
for hz in (1.0, 2.2, 2.6):
    B = [0.2, 0.0, hz]
    E = {L: low(ham(J, B, L)[0])[0] for L in (10, 14)}
    ev = json.load(open(f'../work_C165/loc/P_hz{hz}.json'))['results'][-1]['e']
    print(f'hz={hz}: ED bulk e = {(E[14]-E[10])/4:.5f}   C165 VUMPS e = {ev:.5f}', flush=True)
