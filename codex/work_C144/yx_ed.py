"""Ring spectra of XZ + eps YX: lowest levels, L*(E_n - E_0), and the Casimir combination from E_L."""
import sys
import numpy as np
from xzlib import term, ground_energy
for eps in [float(x) for x in sys.argv[1].split(",")]:
    h = term("xz") + eps * term("yx")
    E = {}
    for L in (8, 10, 12, 14, 16, 18):
        w = ground_energy(h, L, k=6)
        E[L] = w[0]
        print(f"eps={eps} L={L}: E0/L={w[0]/L:.6f}  L*(E-E0)={np.round(L*(w-w[0]),3).tolist()}", flush=True)
    Ls = sorted(E)
    for i in range(len(Ls) - 2):
        # fit E = e L - pi A/(6 L) over three sizes
        L3 = np.array(Ls[i:i + 3], float); M = np.stack([L3, -np.pi / (6 * L3), np.ones(3)], 1)
        e, A, c0 = np.linalg.solve(M, np.array([E[int(l)] for l in L3]))
        print(f"   A_eff (Casimir c v) from L={Ls[i]},{Ls[i+1]},{Ls[i+2]}: {A:.4f}")
