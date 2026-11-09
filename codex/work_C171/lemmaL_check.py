"""C171: Lemma L (C163) with the composite boson's V_c at momentum pi (phi2 shift pi/2),
so that only conversion-type couplings can lock the two bosons."""
import math, itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from opcount import boson, host, spin0_invariant_count, ALPHA
b1 = math.pi + ALPHA
for n, delta in ((1, math.pi), (2, 0.0), (2, math.pi), (3, 0.0), (4, 0.0)):
    for K1, K2 in ((1.2, 0.8), (1.2, 1.5), (1.5, 1.0), (1.5, 1.5), (1.8, 1.5), (2.2, 1.5), (2.5, 2.0)):
        best = None
        for c2 in (0.0, math.pi / 2):
            b2 = boson(K2, +1, math.pi / 2, n * b1 + delta, +1, -1, c2, (n * math.pi) % (2 * math.pi), "2")
            k, labs = spin0_invariant_count([host(K1), b2])
            if best is None or k < best[0]:
                best = (k, labs)
        xc = n * n / (4 * K1) + 1 / (4 * K2)
        print(f"n={n} delta={delta:.2f} K1={K1} K2={K2} K1+K2={K1+K2:.2f} x_conv={xc:.2f}: k={best[0]} {best[1]}")
