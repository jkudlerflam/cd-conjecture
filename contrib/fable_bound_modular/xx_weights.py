"""XX chain: the defect e_0 W - E_0(K_w) for other weight profiles, against the CFT functional
    e_0 W - E_0(K_w) = (A/(24 pi)) int (w')^2 / w dx + O(1)
(A = 1 for the XX chain).  Tent w = min(x, ell-x): integral 2 log ell, same slope 1/(12 pi).
Quartic w = x^2 (ell-x)^2 / ell^3: integral 4/3, so the defect should approach a constant.
Linear ramp w = x on [0, ell] with a hard wall at ell: integral log ell, slope 1/(24 pi), but the wall
adds a boundary energy -f_b * w(ell) that grows linearly in ell."""
import numpy as np
from xx_modular import kmat
e0 = -1 / np.pi
def defect(wts):
    k = kmat(wts); e = np.linalg.eigvalsh(k); return e0 * wts.sum() - e[e < 0].sum()
prev = {}
for ell in [128, 256, 512, 1024, 2048]:
    x = np.arange(ell - 1) + 1.0
    prof = {'tent': np.minimum(x, ell - x), 'quartic': x**2 * (ell - x)**2 / ell**3, 'ramp+wall': x}
    line = f"ell={ell:5d}"
    for name, w in prof.items():
        d = defect(w)
        if name in prev:
            line += f"   {name}: defect={d:.5f} slope={(d - prev[name]) / np.log(2):.5f}"
        else:
            line += f"   {name}: defect={d:.5f}"
        prev[name] = d
    print(line)
print("CFT slopes: tent 1/(12 pi) = %.5f, quartic 0, ramp+wall 1/(24 pi) = %.5f plus wall term (-f_b) ell = %.4f ell" % (1/(12*np.pi), 1/(24*np.pi), 2/np.pi - 0.5))
