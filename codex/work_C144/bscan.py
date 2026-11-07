"""ED of the effective B chain h_B(lam) (bmodel.py) on rings L, in the physical sector
(c, parity) = (+1, 0) for L = 0 mod 4 and (-1, 1) for L = 2 mod 4 (fixed by validate_bmodel.out).
Prints the lowest levels, the splitting of the lowest pair, the gap to the third level, and the
dimerization D = <psi+| b_1 - b_2 |psi+> of psi+ = (psi0 + psi1)/sqrt2."""
import sys, time
import numpy as np
from kernel import kernel_row
import bmodel

def phys(L):
    return (1, 0) if L % 4 == 0 else (-1, 1)

lams = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.0]
Ls = [int(x) for x in sys.argv[1].split(",")]
for lam in lams:
    for L in Ls:
        t0 = time.time()
        K, E0A = kernel_row(L, -1); K = K.real
        c, par = phys(L)
        w = bmodel.spectrum(L, lam, K, c, par, k=4)
        w2, V, st, index, _ = bmodel.ground_vec(L, lam, K, c, par)
        psi = (V[:, 0] + V[:, 1]) / np.sqrt(2)
        b = bmodel.expect_b(L, st, index, psi, c)
        Dm = np.mean(b[0::2]) - np.mean(b[1::2])
        print(f"lam={lam:+.3f} L={L:2d}: E/L={w[0]/L:+.5f} levels-E0={np.round(w-w[0],4).tolist()} "
              f"split={w[1]-w[0]:.3e} gap3={w[2]-w[0]:.4f} L*gap3={L*(w[2]-w[0]):.3f} D={Dm:+.4f} ({time.time()-t0:.1f}s)", flush=True)
