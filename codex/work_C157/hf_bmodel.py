"""Uniform Hartree-Fock (hfk.py) of the second-order B chain of XZ + alpha ZX + beta YY (+ dimerized-A variants):
    h_B = lam sum b_j - 2 sum_{r>=2} K(r) b_j b_{j+r},   lam = t / beta^2,  i.e. t = lam, U_r = -2 K(r),
with K the YY kernel of kernel_q.py on a ring of 400 (r <= 200). Reports the number of Majorana zeros of the
mean-field band and the Lifshitz indicator v = sum_r r h_r (slope at k = 0; c = 3/2 requires it to change sign).
Also: the same with K replaced by its r = 2 term only (pure RZFA, g = 2 K(2) > 0)."""
import numpy as np
import hfk
from kernel_q import kernel_rows

L = 400
for d in (0.0, 0.05, 0.2):
    hop = 1.0 + d * (-1.0) ** np.arange(L)
    K, _ = kernel_rows(L, hop)
    Kav = 0.5 * (K[0] + K[1])
    U = {r: -2 * Kav[r] for r in range(2, 199)}
    print(f"A dimerization d = {d}: K(2) = {Kav[2]:+.5f}, K(3) = {Kav[3]:+.5f}, K(4) = {Kav[4]:+.5f}")
    for lam in (8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.03, 0.01, -0.01, -0.1, -1.0, -4.0):
        r = hfk.solve(lam, U)
        h = r["h"]; v = sum(i * h[i] for i in range(1, 200))
        print(f"   lam = {lam:+6.2f}: G1 = {r['G'][1]:+.4f} G3 = {r['G'][3]:+.4f}  h1 = {h[1]:+.4f} h3 = {h[3]:+.4f}  "
              f"h3/h1 = {h[3]/h[1]:+.4f}  v = {v:+.4f}  zeros = {r['zeros']}  c_MF = {r['c']}  (err {r['err']:.0e})")
print("pure RZFA term of the same sign (U_2 = -2 K(2) only), lam = t:")
K, _ = kernel_rows(L, np.ones(L))
for lam in (1.0, 0.1, 0.01, 0.001):
    r = hfk.solve(lam, {2: -2 * K[0][2]})
    h = r["h"]
    print(f"   lam = {lam}: h1 = {h[1]:+.4f} h3 = {h[3]:+.5f} h3/h1 = {h[3]/h[1]:+.4f} zeros = {r['zeros']}")
