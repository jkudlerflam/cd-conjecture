"""C160: fine kappa scan across the CSL -> canted transition at h = 0.45 (csl_tilt.py functional).
For each kappa: best texture over T and several initial tilts, against the planar CSL and the uniform state."""
import numpy as np
from csl_tilt import csl, full_energy, uniform_min
h = 0.45
T0, e0, t, dx = csl(h)
print(f'h={h}: planar CSL T={T0:.4f} e={e0:+.6f}')
for kappa in np.arange(-0.86, -0.905, -0.005):
    best = None
    for tilt0 in (0.05, 0.2, 0.5):
        for T in np.linspace(0.9 * T0, 1.8 * T0, 7):
            ee, nz, p = full_energy(h, kappa, T, tilt0=tilt0)
            if nz.std() > 1e-3 and (best is None or ee < best[0]):     # genuine textures only
                best = (ee, T, nz.mean(), np.abs(nz).max())
    eu = uniform_min(h, kappa)
    s = f'kappa={kappa:+.3f}: uniform {eu:+.6f}, planar CSL {e0:+.6f}'
    if best: s += f', best texture {best[0]:+.6f} (T={best[1]:.2f}, <n_z>={best[2]:+.3f})'
    print(s)
