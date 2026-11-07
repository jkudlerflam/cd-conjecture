"""ED cross-check of the Bethe branch (minimizing flux, sector ground energies), L<=12, all sectors,
0<Delta<1; Delta=1/2 supersymmetric identity; Theorem A ceilings. Output ed_check.json."""
import json, numpy as np
from bethe_scan import solve
from check_c129 import energy as ed
from continuum_rep import Model, e_bulk
from constants import K as Kconst
out = dict(rows=[], susy=[])
worst = 0.
for delta in [0.1, 0.5, 0.9, 0.99, 0.999]:
    M = Model(delta); Kc = Kconst(M.g)
    for L in range(2, 13):
        for n in range(1, L//2+1):
            E, d, x = solve(L, n, delta, return_roots=True)
            Eed = ed(L, n, delta)
            err = abs(float(E)-Eed); worst = max(worst, err)
            eb = e_bulk(M, n/L)
            out['rows'].append(dict(delta=delta, L=L, n=n, bethe=float(E), ed=Eed, err=err,
                                    L_F_minus_Le=float(L*(Eed-L*eb)), K=Kc))
    print(delta, 'max err so far', worst, flush=True)
# Delta=1/2: odd L exact -3L/8 at the minimizing flux (sector (L-1)/2), Bethe and ED
for L in [3, 5, 7, 9, 11, 33, 65, 129, 257]:
    E, d = solve(L, (L-1)//2, 0.5)
    row = dict(L=L, bethe_minus_exact=float(E+3*L/8))
    if L <= 11: row['ed_minus_exact'] = ed(L, (L-1)//2, 0.5)+3*L/8
    out['susy'].append(row); print(row)
# even L: twisted energy at flux pi+2pi/3 equals -3L/8 (ED)
for L in [4, 6, 8, 10, 12]:
    out['susy'].append(dict(L=L, twisted_ed_minus_exact=ed(L, L//2, 0.5, a=np.pi+2*np.pi/(3*L))+3*L/8))
out['max_bethe_ed_error'] = worst
out['max_ratio_to_K'] = max(abs(r['L_F_minus_Le'])/r['K'] for r in out['rows'])
print('max err', worst, 'max |L(F-Le)|/K', out['max_ratio_to_K'])
json.dump(out, open('ed_check.json', 'w'), indent=1)
