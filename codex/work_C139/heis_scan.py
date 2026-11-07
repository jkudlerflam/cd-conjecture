"""Half-filled doubling defect Q=L(F_2L-2F_L)/M0 and global positivity of sigma_L as Delta->1,
including Delta=1 (XXX rapidities). Also sigma_L at roots. J=1, M0=1/2. Diagnostics only."""
import json, sys, numpy as np
from bethe_scan import solve
out = []
for delta in [0.5, 0.9, 0.99, 0.999, 1.0]:
    if delta < 1:
        g = np.arccos(delta)
        a1 = lambda x: np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
        K = lambda x: np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
        v = np.pi*np.sin(g)/(2*g)
    else:
        a1 = lambda x: 2/(np.pi*(1+4*x*x))
        K = lambda x: 1/(np.pi*(1+x*x))
        v = np.pi/2
    prev = {}
    for L in [2, 4, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
        for LL in [L, 2*L]:
            if LL not in prev:
                E, d, x = solve(LL, LL//2, delta, return_roots=True)
                X = np.max(np.abs(x))
                grid = np.concatenate([np.linspace(-X-8, X+8, 8001), x, np.linspace(X, 60 if delta<1 else 1e4, 2000)])
                sig = np.empty_like(grid)
                for i in range(0, len(grid), 2000):
                    gg = grid[i:i+2000]; sig[i:i+2000] = a1(gg)-K(gg[:, None]-x[None, :]).sum(1)/LL
                prev[LL] = (E, float(sig.min()), float(LL*sig[8001:8001+len(x)].min()))
        E, smin, sroot = prev[L]; E2 = prev[2*L][0]
        row = dict(delta=delta, L=L, Q=float(L*(E2-2*E)/0.5), casimir_limit=float(np.pi*v/2),
                   min_sigma=smin, min_L_sigma_roots=sroot)
        out.append(row); print(json.dumps(row), flush=True)
    json.dump(out, open('heis_scan.json', 'w'), indent=1)
