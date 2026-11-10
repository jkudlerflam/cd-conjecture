"""L = 36, N = 6, phi = pi/2: the first excited level is only 0.011 above the ground state, so the flux step 0.05 of
manybody.py crosses levels. Redo the flux curvature of the ground state with dPhi = 0.002 and 0.004 (Richardson check),
and print the three lowest levels at each flux. v/K = 0.2035 is taken from manybody.out (E(N +- 2), no flux)."""
import numpy as np
from ringv import e0, basis
L, N = 36, 6; b = basis(L, N)
E = {}
for P in (0.0, 0.002, -0.002, 0.004, -0.004):
    E[P] = e0(L, N, 1j, Phi=P, k=3, b=b); print(f"Phi={P:+.3f}: lowest levels {E[P]}", flush=True)
for d in (0.002, 0.004):
    c = (E[d][0] + E[-d][0] - 2 * E[0.0][0]) / d ** 2
    vK = np.pi * L * c
    print(f"dPhi={d}: vK = {vK:.4f}, K_rho = {np.sqrt(vK/0.2035) if vK>0 else float('nan'):.3f}", flush=True)
