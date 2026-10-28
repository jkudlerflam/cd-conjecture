"""C7's finite-volume hypothesis on ULS rings: number of lattice states with
L (E - E_0)/(2 pi v) <= Delta_c = c/4 + 2 against the SU(d)_1 count N_CFT(Delta_c),
and C7's lower bound 1 + (pi e/12) c exp(pi c/6). Exact spectra only (L <= 12)."""
import numpy as np
from thermo import load_uls, uls_levels, sun1_levels
for d, L in [(2, 8), (2, 12), (3, 9), (3, 12), (4, 8), (4, 12), (5, 10), (6, 12)]:
    c, v = d - 1, 2 * np.pi / d
    Dc = c / 4 + 2
    E, g = uls_levels(load_uls(L), d)
    x = L * (E - E.min()) / (2 * np.pi * v)
    Ds, ms = sun1_levels(d, Dc + 0.5)
    ncft = ms[Ds <= Dc + 1e-9].sum()
    nlat = g[x <= Dc + 1e-9].sum()
    lb = 1 + np.pi * np.e / 12 * c * np.exp(np.pi * c / 6)
    print(f"d={d} L={L:2d} Delta_c={Dc:.2f}  N_lat={nlat:8.0f}  N_CFT={ncft:8.0f}  ratio={nlat/ncft:6.3f}  C7 lower bound={lb:6.1f}  d^L={float(d)**L:.3g}")
