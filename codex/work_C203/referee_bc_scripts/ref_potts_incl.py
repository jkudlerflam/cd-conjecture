"""Referee C203: check (C1) inclusion with the finite twist list, own code."""
import numpy as np
from ref_potts import potts_H, xxz_sector
for Q, Ls in [(2, [2, 3, 4]), (3, [2, 3, 4]), (4, [2, 3])]:
    g = np.arccos(np.sqrt(Q) / 2); a = np.cos(g)
    for sign in [+1, -1]:
        for L in Ls:
            ev = np.linalg.eigvalsh(potts_H(Q, L, sign).toarray())
            N = 2 * L
            pool = {}
            for phi in [2 * g, -2 * g, np.pi]:
                pool[(0, phi)] = np.linalg.eigvalsh(xxz_sector(N, sign * a, phi, 0).toarray())
            for j in range(1, L + 1):
                for k in range(j):
                    pool[(j, 2 * np.pi * k / j)] = np.linalg.eigvalsh(xxz_sector(N, sign * a, 2 * np.pi * k / j, 2 * j).toarray())
            allv = np.concatenate(list(pool.values()))
            nom1 = np.concatenate([v for (j, p), v in pool.items() if j != 1])
            d_all = max(np.min(np.abs(allv - e)) for e in ev)
            d_nom1 = max(np.min(np.abs(nom1 - e)) for e in ev)
            print(f"Q={Q} s={sign:+d} L={L}: {len(ev)} Potts levels, max dist to list={d_all:.1e}, without m=1: {d_nom1:.1e}")
