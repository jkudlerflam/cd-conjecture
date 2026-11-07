"""C147 (a) Bethe vs ED at flux pi per bond, L<=12, both sides; (b) ED sector defect L(F_2L(2n)-2F_L(n))/M_0, L<=6;
(c) flux minimality spot check. Own code (ref_core.py). J=1, M_0=1/2."""
import json, sys
import numpy as np
from ref_core import ed_sector, bethe_rep, energy_rep, bethe_att, energy_att

out = dict(rep={}, att={}, defect={}, minhess_rep={})
for Delta in [0.3, 0.7, 0.9, 0.99, 0.999]:
    g = np.arccos(Delta); err = 0.; mh = np.inf
    for L in range(2, 13):
        for n in range(1, L//2+1):
            x, h = bethe_rep(L, n, g)
            E = energy_rep(L, n, g, x); Ee = ed_sector(L, n, Delta)
            err = max(err, abs(E-Ee)); mh = min(mh, h/L)
    out['rep'][str(Delta)] = err; out['minhess_rep'][str(Delta)] = mh
    print('rep', Delta, err, mh, flush=True)
for eps in [0.5, 0.1, 1e-2, 1e-4]:
    eta = np.arccos(1-eps); Delta = -np.cos(eta); err = 0.
    for L in range(2, 13):
        for n in range(1, L//2+1):
            x = bethe_att(L, n, eta)
            E = energy_att(L, n, eta, x); Ee = ed_sector(L, n, Delta)
            err = max(err, abs(E-Ee))
    out['att'][str(eps)] = err
    print('att', eps, err, flush=True)
for Delta in [0.9, 0.99, 0.999, 1.0, 0.0, -0.5, -0.99]:
    w = -np.inf; arg = None
    for L in range(2, 7):
        for n in range(1, L//2+1):
            D = L*(ed_sector(2*L, 2*n, Delta)-2*ed_sector(L, n, Delta))/0.5
            if D > w: w, arg = D, (L, n)
    out['defect'][str(Delta)] = dict(max=w, at=arg)
    print('defect', Delta, w, arg, flush=True)
json.dump(out, open('check_ed.json', 'w'), indent=1)
