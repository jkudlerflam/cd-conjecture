"""Diagnostics of the effective B chain at given lam: more levels, <b_odd>, <b_even> in the
symmetry-broken combinations, product-state energies."""
import sys
import numpy as np
from kernel import kernel_row
import bmodel

lam = float(sys.argv[2])
for L in [int(x) for x in sys.argv[1].split(",")]:
    K, _ = kernel_row(L, -1); K = K.real
    c, par = (1, 0) if L % 4 == 0 else (-1, 1)
    w = bmodel.spectrum(L, lam, K, c, par, k=12)
    w2, V, st, index, _ = bmodel.ground_vec(L, lam, K, c, par)
    out = []
    for sgn in (1, -1):
        psi = (V[:, 0] + sgn * V[:, 1]) / np.sqrt(2)
        b = bmodel.expect_b(L, st, index, psi, c)
        out.append((b[0::2].mean(), b[1::2].mean()))
    # product-state energy (b_odd = 1) per site
    Eprod = -sum((0.5 if r == L // 2 else 1.0) * K[r] for r in range(2, L // 2 + 1, 2)) + lam * 0.5
    print(f"L={L} lam={lam}: E0/L={w[0]/L:.5f} (product b_odd=1: {Eprod:.5f}); "
          f"L*(E-E0)={np.round(L*(w-w[0]),2).tolist()}; <b_odd>,<b_even> in psi+-: {np.round(out,4).tolist()}", flush=True)
