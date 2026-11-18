"""Bethe (circle) vs ED sector ground energy at flux pi per bond, every 1 <= n <= L/2, L = 2..14, both parities."""
import numpy as np, json, time
from ref148 import *
t0 = time.time()
Deltas = [1.0001, 1.001, 1.01, 1.1, 1.5, 2.0, cosh := float(np.cosh(1.36)), 3.0, 5.0, 10.0, 20.0, 100.0]
worst = 0; rows = []
for D in Deltas:
    lam = float(np.arccosh(D))
    for L in range(2, 15):
        for n in range(1, L // 2 + 1):
            Eb, x = bethe_energy(L, n, lam)
            Ee = sector_energy(L, n, D)
            err = abs(Eb - Ee)
            worst = max(worst, err)
            rows.append(dict(Delta=D, L=L, n=n, bethe=Eb, ed=Ee, err=err))
    print(D, 'max err so far', worst, flush=True)
bad = [r for r in rows if r['err'] > 1e-9]
json.dump(dict(n_cases=len(rows), max_err=worst, n_bad=len(bad), bad=bad[:20]), open('check_ed_bethe.json', 'w'), indent=1)
print('cases', len(rows), 'max err', worst, 'bad', len(bad), 'time', time.time() - t0)
