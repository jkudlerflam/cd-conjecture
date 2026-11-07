"""ED check of the free-fermion dispersion in the paramagnetic region (unique ground state): c=1, g=0.5, B=1.5.
Below a_c the many-body gap should match 2 min_k f over allowed single-particle momenta; at a_c it closes like 1/L^2."""
import numpy as np
from scipy.optimize import brentq
from free_pt import f, fmin, spin_bond
from ed import spectrum, all_levels
c, g, B = 1.0, 0.5, 1.5
ac = brentq(lambda a: fmin(a, c, B, g)[0], 1e-6, 5.0)
print("a_c =", ac)
for L in [8, 10, 12, 14]:
    out = []
    for a in [0.5 * ac, 0.8 * ac, ac]:
        w, k = all_levels(spectrum(spin_bond(a, c, B, g), L)); w = np.sort(w)
        # allowed momenta for one quasiparticle above the even-parity vacuum: k = 2 pi m/L (periodic) or (2m+1) pi/L
        kp = 2 * np.pi * np.arange(L) / L; ka = (2 * np.arange(L) + 1) * np.pi / L
        pred = 2 * min(f(kp, a, c, B, g).min(), f(ka, a, c, B, g).min())
        out.append((w[1] - w[0], pred))
    print(f"L={L}: " + "; ".join(f"gap {x:.5f} (free-fermion lowest 2f on k-grids {p:.5f})" for x, p in out),
          f" L^2*gap(a_c)={L*L*out[2][0]:.3f}")
