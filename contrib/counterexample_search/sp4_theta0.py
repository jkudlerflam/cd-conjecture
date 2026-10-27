#!/usr/bin/env python3
"""Sp(4) chain h = cos(t) P5 + sin(t) P10 near t = 0.
(1) Kernel of H0 = sum P5 (dimension by momentum and Cartan sector).
(2) First-order degenerate perturbation theory: H_eff = Pi_K (sum P10) Pi_K.
Dense diagonalization, L <= 8.  Output: sp4_theta0.json"""
import json, sys, itertools
import numpy as np
from scipy import linalg
from ced import Sector
from models import sp4

m = sp4()
out = []
for L in range(2, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 8):
    kdim = 0
    total = 0
    heff_levels = []
    sym_dim = int(np.math.comb(L + 3, 3)) if hasattr(np, 'math') else None
    for q1 in range(-L, L + 1):
        for q2 in range(-L, L + 1):
            if abs(q1) + abs(q2) > L or (q1 + q2 + L) % 2:
                continue
            sec = Sector(4, L, m['hs'], m['charges'], (q1, q2))
            if len(sec.basis) == 0:
                continue
            total += len(sec.basis)
            for k in range(L):
                ids, (H5, H10) = sec.blocks(k)
                if len(ids) == 0:
                    continue
                w, v = linalg.eigh(H5.toarray())
                ker = v[:, w < 1e-9]
                kdim += ker.shape[1]
                if ker.shape[1]:
                    he = ker.conj().T @ H10.toarray() @ ker
                    ev = linalg.eigvalsh(he)
                    heff_levels += [(float(x), k, q1, q2) for x in ev]
    heff_levels.sort()
    from math import comb
    rec = dict(L=L, total=total, kernel_dim=kdim, sym_dim=comb(L + 3, 3),
               heff_lowest=heff_levels[:12])
    print(json.dumps(rec), flush=True)
    out.append(rec)
json.dump(out, open('data/sp4_theta0.json', 'w'), indent=1)
