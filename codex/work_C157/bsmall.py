"""ED of the second-order B chain (C144 bmodel.py, unmodified) at small lam = t/beta^2, where uniform Hartree-Fock
comes closest to a Lifshitz point (h3/h1 up to 0.94). Ring sector (c, P) = (+1, even), L = 0 mod 4 Majoranas.
Prints L (E - E0) of the 6 lowest levels and the staggered and uniform parts of <b_j> in the ground state."""
import numpy as np
from bmodel import spectrum, ground_vec, expect_b
from kernel import kernel_row
for lam in (0.1, 0.25):
    for L in (24, 32, 40):
        K, _ = kernel_row(L, -1); K = K.real
        w = spectrum(L, lam, K, 1, 0, k=6)
        e, V, st, index, _ = ground_vec(L, lam, K, 1, 0)
        b = expect_b(L, st, index, V[:, 0], 1)
        stg = np.mean(b * (-1.0) ** np.arange(L)); uni = np.mean(b)
        print(f"lam = {lam} L = {L}: L(E - E0) = {np.round((w - w[0]) * L, 3).tolist()}  <b> uniform {uni:+.4f} staggered {stg:+.4f}", flush=True)
