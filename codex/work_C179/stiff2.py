"""Flux stiffness when the ring ground state is a degenerate momentum doublet (E(Phi) then has a cusp at Phi = 0 for the
lowest level alone). Uses the two lowest levels: if they are degenerate at Phi = 0 to 1e-8, the stiffness is taken from
the doublet mean (E_0 + E_1)/2, which is smooth in Phi. Prints the gap E_1 - E_0 at Phi = 0 and both estimates."""
import numpy as np, sys
from ringv import e0, basis
w = 1j
for (L, N) in [tuple(int(x) for x in a.split(',')) for a in sys.argv[1:]]:
    b = basis(L, N); dP = 0.05
    E = {P: e0(L, N, w, Phi=P, k=2, b=b) for P in (0.0, dP, -dP)}
    c0 = (E[dP][0] + E[-dP][0] - 2 * E[0.0][0]) / dP ** 2
    cm = (E[dP].mean() + E[-dP].mean() - 2 * E[0.0].mean()) / dP ** 2
    print(f"L={L} N={N} phi=pi/2: E1-E0 at Phi=0 = {E[0.0][1]-E[0.0][0]:.2e}; vK from lowest level = {np.pi*L*c0:.4f}, from doublet mean = {np.pi*L*cm:.4f}", flush=True)
