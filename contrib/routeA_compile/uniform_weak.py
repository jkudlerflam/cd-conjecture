"""Uniform cell P D^n P in the weak-decay regime: c_max versus n."""
import json, time
import numpy as np
from onebody import logical_h, blocks, max_c
out = []
for delta in [0.3, 0.02]:
    for n in [2, 3, 4, 5]:
        w = (1 - delta) / 2
        p = dict(eT=1, eH=1, eD=1, wT=w, wH=w, wD=-w, tauA=0.1, tauB=0.1, kappa=1.0, phiT=0.0)
        E0, h = logical_h(n, p)
        T, M = blocks(h, n)
        c, mu = max_c(T, M, nk=4000)
        r = dict(delta=delta, n=n, c_max=c, T_abs_over_max=(np.abs(T) / np.abs(T).max()).tolist())
        print(json.dumps(dict(delta=delta, n=n, c_max=c)), flush=True); out.append(r)
json.dump(out, open('uniform_weak_results.json', 'w'), indent=1)
