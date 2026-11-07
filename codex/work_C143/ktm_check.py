"""Kurmann-Thomas-Muller factorizing field for XYZ in a uniform z field, H = sum jx XX + jy YY + jz ZZ + b sum Z,
b_f = 2 sqrt((jx + jz)(jy + jz)). Check frustration-freeness (E_0/L equals the gauge-optimized bond minimum),
ground degeneracy and gap on rings; compare a generic KTM point with points approaching jx = -jy."""
import numpy as np
from ed import bond, spectrum, all_levels, gauge_lmin
def run(jx, jy, jz, Ls=(6, 8, 10, 12)):
    b = 2 * np.sqrt((jx + jz) * (jy + jz))
    h = bond(np.diag([jx, jy, jz]), (0, 0, b))
    lmin, _ = gauge_lmin(h)
    print(f"jx={jx} jy={jy} jz={jz} b_f={b:.4f}")
    for L in Ls:
        w, _ = all_levels(spectrum(h, L)); w = np.sort(w)
        deg = int(np.sum(w - w[0] < 1e-8)); gap = w[deg] - w[0]
        print(f"   L={L}: E0/L - lmin = {w[0]/L - lmin:.2e}, ground deg = {deg}, gap = {gap:.5f}, L^2 gap = {L*L*gap:.3f}")
run(0.6, 0.2, 0.3)
run(0.4, -0.2, 0.3)
run(0.3, -0.25, 0.3)
run(0.3, -0.3, 0.3)
