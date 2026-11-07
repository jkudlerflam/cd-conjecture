"""Validate the second-order effective B chain against full qubit ED on the tuned cone, small beta."""
import numpy as np
from xzlib import bond, ring_operator
from kernel import kernel_row
import bmodel

def full_spec(L, al, be, k):
    J = np.zeros((3, 3)); J[0, 2] = 1; J[2, 0] = al; J[1, 1] = be
    h = bond(J); N = 1 << L; op = ring_operator(h, L)
    H = np.zeros((N, N), complex)
    for x in range(N):
        e = np.zeros(N, complex); e[x] = 1; H[:, x] = op.matvec(e)
    return np.linalg.eigvalsh(H)[:k]

import sys
for L in [int(a) for a in sys.argv[1:]] or (8, 12):
    K, E0A = kernel_row(L, -1)
    K = K.real
    abar = E0A / L
    D = 2 ** (L // 2 - 1)
    for lam in (0.0, 1.0):
        be = 0.004
        al = -be * abar + lam * be**2
        w = full_spec(L, al, be, D)
        x = (w - E0A) / be**2 + L * K[0]
        print(f"L={L} lam={lam}: full (E-E0A)/beta^2 + L K0, lowest {min(D,12)}:", np.round(x[:12], 3).tolist())
        for c in (1, -1):
            for par in (0, 1):
                sp = bmodel.spectrum(L, lam, K, c, par, k=min(D, 12))
                print(f"     c={c:+d} parity={par}:", np.round(sp, 3).tolist())
