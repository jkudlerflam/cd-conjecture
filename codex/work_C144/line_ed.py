"""Ring spectra on the AT plane near the B-critical line (beta = 0.3): lowest levels and L*(E_n - E_0)."""
import sys
import numpy as np
from xzlib import term, ground_energy
be = 0.3
for al in [float(x) for x in sys.argv[1].split(",")]:
    h = term("xz") + al * term("zx") + be * term("yy")
    for L in (8, 12, 16, 20):
        w = ground_energy(h, L, k=8)
        print(f"alpha={al} L={L}: E0/L={w[0]/L:.6f} L*(E-E0)={np.round(L*(w-w[0]),3).tolist()}", flush=True)
