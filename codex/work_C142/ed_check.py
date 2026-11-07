"""ED cross-check at L <= 12 (flux pi per bond, the sector minimizer). Floating point. J=1, M_0=1/2.
(a) Bethe (bethe_scan, repulsive; bethe_ferro, attractive) against ED sector energies.
(b) Sector doubling defect L(F_2L(2n)-2F_L(n))/M_0 by ED for L<=6, against the proved ceilings:
    4059 (Theorem I.9 assembled, 0<=Delta<=1) and (15pi/2) eta L (C138 Theorem A, attractive)."""
import json
import numpy as np
from check_c129 import energy
from bethe_scan import solve as solve_rep
from bethe_ferro import solve as solve_att

out = dict(rep=[], att=[], defect=[])
err_rep = err_att = 0.
for delta in [0.99, 0.999]:
    for L in range(2, 13):
        for n in range(1, L//2+1):
            E, d = solve_rep(L, n, delta)
            e = abs(float(E)-energy(L, n, delta)); err_rep = max(err_rep, e)
for eps in [1e-2, 1e-4]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); Delta = -np.cos(eta)
    for L in range(2, 13):
        for n in range(1, L//2+1):
            G, x, d = solve_att(L, n, eta)
            e = abs(Delta*L/4+G-energy(L, n, Delta)); err_att = max(err_att, e)
worst = {}
for delta in [0.99, 0.999, 1.0, -1+1e-2, -1+1e-4]:
    w = -np.inf
    for L in range(2, 7):
        for n in range(0, L//2+1):
            D = L*(energy(2*L, 2*n, delta)-2*energy(L, n, delta))/0.5 if n > 0 else 0.0
            w = max(w, D)
            out['defect'].append(dict(delta=delta, L=L, n=n, defect_over_M0=D))
    worst[str(delta)] = w
out['summary'] = dict(max_bethe_ed_error_repulsive=err_rep, max_bethe_ed_error_attractive=err_att,
                      max_defect_over_M0=worst)
json.dump(out, open('ed_check.json', 'w'), indent=1)
print(json.dumps(out['summary'], indent=1))
