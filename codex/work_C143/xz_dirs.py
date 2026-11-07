"""Splitting of the 2^(L/2-1) ground manifold at the XZ point along each tangent direction, versus eps."""
import numpy as np
from ed import bond, spectrum, all_levels
lab = "xyz"
J0 = np.zeros((3, 3)); J0[0, 2] = 1.0
dirs = []
for a in range(3):
    for b in range(3):
        if (a, b) == (0, 2): continue
        Jd = np.zeros((3, 3)); Jd[a, b] = 1.0
        dirs.append((lab[a] + lab[b], Jd, np.zeros(3)))
for a in range(3):
    Bd = np.zeros(3); Bd[a] = 1.0
    dirs.append(("B" + lab[a], np.zeros((3, 3)), Bd))
L = 10
D = 2 ** (L // 2 - 1)
print(f"L={L}, manifold size D={D}; width = E_D - E_1 of lowest D states; gap = E_(D+1) - E_D")
for nm, Jd, Bd in dirs:
    row = []
    for eps in [0.02, 0.04, 0.08]:
        w, _ = all_levels(spectrum(bond(J0 + eps * Jd, eps * Bd), L))
        w = np.sort(w)
        row.append((w[D - 1] - w[0], w[D] - w[D - 1]))
    widths = [r[0] for r in row]
    order = np.log(widths[2] / max(widths[0], 1e-300)) / np.log(4) if widths[0] > 1e-12 else float('nan')
    print(f"{nm}: widths {['%.2e' % x for x in widths]}  apparent order in eps {order:.2f};  gaps above {['%.3f' % r[1] for r in row]}")
