#!/usr/bin/env python3
"""Interval entanglement of the static (u = 0) selector in one marker pattern (all cells of
length 3), where only the payload chain is gapless: an in-model c = 1 calibration of c_half.
Writes data/ed/static_entropy.json.  Claude (Opus 5.5), 2026-10-27."""
import json, time, numpy as np
from ed import Ring
out = []
for L in (9, 12, 15, 18, 21, 24):
    M = L // 3
    t0 = time.time()
    r = Ring(L, 2 * M, M, t=1.0, u=0.0, mu=-1.0, J=0.5, cells=(3,) * M)
    w, V = r.lowest(k=3)
    # cuts at cell boundaries only (site 0 is a marker): ell = 3, 6, ..., L - 3
    ells = list(range(3, L, 3))
    S = [r.half_entropy(V[:, 0], ell) for ell in ells]
    out.append(dict(L=L, M=M, dim=r.dim, E=w.tolist(), ell=ells, S=S))
    print(L, r.dim, np.round(w, 6), dict(zip(ells, np.round(S, 5))), f'{time.time()-t0:.0f}s', flush=True)
json.dump(out, open('data/ed/static_entropy.json', 'w'), indent=1)
