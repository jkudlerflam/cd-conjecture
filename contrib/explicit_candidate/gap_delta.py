"""Cell gaps of the 23-letter cell (vacuum-walker gap and excited-token gap) for
n = 1..6, and the effect of a small excited on-site energy Delta on the
Hermitian part of the g^2 hopping coefficient.  Usage: python3 gap_delta.py
-> gap_delta_results.json"""
import json, os
import numpy as np
from bond23 import params, cell_from_bond, onebody_T, herm_eigs, HERE

out = dict(gaps=[], delta_scan=[])
for n in (1, 2, 3, 4, 5, 6):
    for (Dx, g, Delta) in ((4.0, 0.05, 0.0), (20.0, 0.05, 0.0), (20.0, 0.01, 0.0), (20.0, 0.05, 0.1)):
        P = params(Dx=Dx, DM=Dx, gp=g, gm=g, Delta=Delta)
        c = cell_from_bond(n, P)
        E0 = -2 * np.cos(np.pi / (n + 2))
        ev = np.linalg.eigvalsh(c['H'][np.ix_(~c['exc'], ~c['exc'])])
        ee = np.linalg.eigvalsh(c['H'][np.ix_(c['exc'], c['exc'])])
        rec = dict(n=n, Dx=Dx, g=g, Delta=Delta, E0_formula=E0, E0_vac=float(ev[0]),
                   vacuum_gap=float(ev[2 ** n] - ev[0]), excited_gap=float(ee[0] - ev[0]),
                   omega_formula=float(2 * (np.cos(np.pi / (n + 2)) - np.cos(2 * np.pi / (n + 2)))))
        print(json.dumps(rec), flush=True)
        out['gaps'].append(rec)
for n in (2, 3, 4, 5):
    for Delta in (0.0, 0.01, 0.1):
        P = params(Dx=20.0, DM=20.0, Delta=Delta)
        TC, _, _ = onebody_T(n, P, 'g2')
        ev, im, sc = herm_eigs(TC)
        rec = dict(n=n, Dx=20.0, Delta=Delta, T2_scale=sc, herm_eigs_over_T=list(ev))
        print(json.dumps(rec), flush=True)
        out['delta_scan'].append(rec)
json.dump(out, open(os.path.join(HERE, 'gap_delta_results.json'), 'w'), indent=1)
