#!/usr/bin/env python3
"""One ED sector of ed.Ring (for the cluster): lowest k levels with momenta and the half-ring
entropy of the lowest state.  usage: ed_one.py OUT.json L W nup t u mu J [k]
Claude (Opus 5.5), 2026-10-27."""
import sys, json, time
import numpy as np
from ed import Ring
out, L, W, nup = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
t, u, mu, J = map(float, sys.argv[5:9])
k = int(sys.argv[9]) if len(sys.argv) > 9 else 6
t0 = time.time()
r = Ring(L, W, nup, t=t, u=u, mu=mu, J=J)
w, V = r.lowest(k=k)
P = r.momenta(w, V)
S = [r.half_entropy(V[:, 0], L // 2)] + ([r.half_entropy(V[:, 0], L // 2 + 1)] if L % 2 else [])
json.dump(dict(L=L, W=W, nup=nup, t=t, u=u, mu=mu, J=J, dim=r.dim, E=w.tolist(), P=P.tolist(), S=S,
               time_s=time.time() - t0), open(out, 'w'), indent=1)
print(L, W, nup, r.dim, np.round(w, 8), P, S, f'{time.time()-t0:.0f}s')
