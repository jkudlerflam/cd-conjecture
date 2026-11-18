"""End-to-end ED of L(E_2L - 2E_L)/M for the full bond (field, DM), and sector defects L(F_2L(2n)-2F_L(n))/M_0 by ED."""
import numpy as np, json, time
from ref148 import *
t0 = time.time()
Deltas = [1.0001, 1.001, 1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 20.0]
# (a) sector defect by ED, L = 2..10 (2L up to 20), all 1 <= n <= L/2, stoquastic gauge
sec = []
for D in Deltas:
    for L in range(2, 11):
        for n in range(1, L // 2 + 1):
            d = L * (sector_energy(2 * L, 2 * n, D) - 2 * sector_energy(L, n, D)) / ((1 + D) / 4)
            sec.append(dict(Delta=D, L=L, n=n, defect=d))
    print('sector', D, max(r['defect'] for r in sec if r['Delta'] == D), time.time() - t0, flush=True)
# (b) full problem, J = 1, DM angle theta, field B, L = 2..6
full = []
for D in Deltas:
    hc2 = 1 + D
    lam = np.arccosh(D)
    for th in [0.0, 0.3, 1.0, np.pi / 2, 2.5, np.pi]:
        Jp, Dm = np.cos(th), np.sin(th)
        for B in np.r_[np.linspace(0, 1.1 * hc2, 23), -0.37 * hc2]:
            M = half_spread(1.0, D, B)
            for L in range(2, 7):
                E1 = full_energy(L, Jp, Dm, D, B); E2 = full_energy(2 * L, Jp, Dm, D, B)
                full.append(dict(Delta=D, theta=float(th), B=float(B), L=L, ratio=L * (E2 - 2 * E1) / M))
    print('full', D, max(r['ratio'] for r in full if r['Delta'] == D), time.time() - t0, flush=True)
ms = max(sec, key=lambda r: r['defect']); mf = max(full, key=lambda r: r['ratio'])
json.dump(dict(sector_cases=len(sec), sector_max=ms, full_cases=len(full), full_max=mf,
               full_min=min(full, key=lambda r: r['ratio']),
               sector_max_by_Delta={str(D): max(r['defect'] for r in sec if r['Delta'] == D) for D in Deltas},
               full_max_by_Delta={str(D): max(r['ratio'] for r in full if r['Delta'] == D) for D in Deltas}),
          open('check_full_ed.json', 'w'), indent=1)
print(ms, mf)
