#!/usr/bin/env python3
"""At the boundary u = r: twist stiffness D_s (E(phi)-E0 = D_s phi^2/(2L)) of the Q = 0 kernel state,
and compressibility slope e''/(r-u) from E(N0+1)-E(N0) = e''/(2L) at u = r - eps.
Classical predictions: D_s = 2/(1+r), e''/(r-u) = 8 r."""
import numpy as np, sys
from ed import lowest
eps = 1e-5; phi = 0.02
for b in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]:
    r = np.sqrt(1 + b * b / 4)
    row = []
    for L in [8, 12, 16, 20]:
        N0 = L // 2
        E0 = lowest(L, N0, r, b)[0]; Et = lowest(L, N0, r, b, phi=phi)[0]
        Ds = 2 * L * (Et - E0) / phi ** 2
        u = r - eps
        e2 = 2 * L * (lowest(L, N0 + 1, u, b)[0] - lowest(L, N0, u, b)[0]) / eps
        row.append(f'L{L}: E0={E0:.2e} Ds={Ds:.5f} e2/(r-u)={e2:.4f}')
    print(f'b={b} r={r:.4f} Ds_cl={2/(1+r):.5f} 8r={8*r:.4f} | ' + ' | '.join(row), flush=True)
