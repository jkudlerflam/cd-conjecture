#!/usr/bin/env python3
"""Table of the exact logical couplings of the C8 chain (from data/heff_c8_n3_lam*.npz):
swap coefficient b, dipole/quadrupole exchange (logical qutrit = spin 1 in the Cartesian basis),
the equivalent bilinear-biquadratic angle, the uniform on-site field per cell (from h2 and h3),
the rest of the SU(3)-breaking two-body part, and the size of the three-cell term."""
import os, sys, itertools, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from cell import hermitian_basis, swap

T = hermitian_basis(3)
dip, quad = [2, 4, 6], [1, 3, 5, 7, 8]
rows = []
for lam in [0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2]:
    z = np.load(os.path.join(HERE, 'data', f'heff_c8_n3_lam{lam}.npz'))
    h2, h3 = z[f'h2_{lam}'], z[f'h3_{lam}']
    c2 = np.array([[np.trace(np.kron(A, B) @ h2).real for B in T] for A in T])
    gD, gQ = c2[dip, dip].mean(), c2[quad, quad].mean()
    b = np.trace(c2[1:, 1:]) / 8
    theta = np.arctan2(gQ, (gD + gQ) / 2) / np.pi
    c3 = np.zeros((9, 9, 9))
    for i, j, k in itertools.product(range(9), repeat=3):
        if (i == 0) + (j == 0) + (k == 0) >= 2:
            c3[i, j, k] = np.trace(np.kron(np.kron(T[i], T[j]), T[k]) @ h3).real
    f = (c2[1:, 0] + c2[0, 1:] + c3[1:, 0, 0] + c3[0, 1:, 0] + c3[0, 0, 1:]) / b
    F = sum(f[m] * T[m + 1] for m in range(8))
    so3 = c2[1:, 1:] / b - np.diag([gD / b if m + 1 in dip else gQ / b for m in range(8)])
    h3c = h3 - np.trace(h3) / 27 * np.eye(27)
    row = dict(lam=lam, b=b, b_over_lam2=b / lam ** 2, dq=(gD - gQ) / b, theta_over_pi=theta,
               field_eigs=np.linalg.eigvalsh(F).tolist(), other_2body=float(np.linalg.norm(so3, 2)),
               h3=float(np.linalg.norm(h3c, 2) / b))
    rows.append(row)
    print(f"lam={lam:<7} b/lam^2={row['b_over_lam2']:.6f} (gD-gQ)/b={row['dq']:+.4f} theta/pi={theta:.5f} "
          f"field eigs/b={np.round(row['field_eigs'], 4).tolist()} other 2-body/b={row['other_2body']:.4f} |h3|/b={row['h3']:.4f}")
json.dump(rows, open(os.path.join(HERE, 'data', 'sw_table.json'), 'w'), indent=1)
