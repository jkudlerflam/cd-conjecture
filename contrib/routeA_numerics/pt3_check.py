#!/usr/bin/env python3
"""Independent check of the O(lam^3) bond: third-order degenerate perturbation theory for two C8
cells, H3 = P V R V R V P with R = -Q H0^{-1} Q (first order vanishes, so this is the whole third-order
term), compared with the lam^3 coefficient of the exact Schrieffer-Wolff bond (data/heff_c8_n2.npz).
Also prints its dipole / quadrupole / field decomposition (logical spin-1 in the Cartesian basis)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from cell import load_c8, hermitian_basis, swap

c = load_c8()
D, N, e = c['D'], c['N'], c['e']
I = np.eye(D)
H0 = np.kron(np.diag(e), I) + np.kron(I, np.diag(e))
Vb = sum(np.kron(c['FR'][a], c['FL'][a]) for a in range(3))
d0 = np.diag(H0).real
code = np.array([i * D + j for i in range(N) for j in range(N)])
Rm = np.where(d0 > 1e-9, -1.0 / np.where(d0 > 1e-9, d0, 1), 0.0)
VP = Vb[:, code]
H3 = VP.conj().T @ (Rm[:, None] * (Vb @ (Rm[:, None] * VP)))
H2 = VP.conj().T @ (Rm[:, None] * VP)
z = np.load(os.path.join(HERE, 'data', 'heff_c8_n2.npz'))
l1, l2 = 0.0025, 0.005
h1, h2 = z[f'h2_{l1}'], z[f'h2_{l2}']
# h(l) = l^2 H2 + l^3 H3 + l^4 H4 ...
r1 = (h1 - l1 ** 2 * H2) / l1 ** 3
r2 = (h2 - l2 ** 2 * H2) / l2 ** 3
est3 = (l2 * r1 - l1 * r2) / (l2 - l1)          # removes the O(l) (fourth-order) contamination
print('|H2 - C8 formula| =', np.linalg.norm(H2 - (-7177 / 8192 * np.eye(9) + 27 / 8192 * swap(3))))
print('|H3_PT - H3_exactSW (Richardson)| / |H3_PT| =', np.linalg.norm(H3 - est3) / np.linalg.norm(H3))
T = hermitian_basis(3)
cf = np.array([[np.trace(np.kron(A, B) @ H3).real for B in T] for A in T])
dip, quad = [2, 4, 6], [1, 3, 5, 7, 8]
np.set_printoptions(precision=5, suppress=True)
print('H3 exchange: dipole diag', cf[dip, dip], ' quadrupole diag', cf[quad, quad])
print('H3 off-diagonal exchange norm', np.linalg.norm(cf[1:, 1:] - np.diag(np.diag(cf[1:, 1:]))))
print('H3 one-site terms (left cell, right cell):', np.linalg.norm(cf[1:, 0]), np.linalg.norm(cf[0, 1:]))
