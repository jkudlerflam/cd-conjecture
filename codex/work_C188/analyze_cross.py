"""C188: branch energies across the SU(2)_1/ordered boundary.
Branches per L (sigma-frame momenta m):  S = SU(2)_1 singlet ground state (m = 0 for L = 12, 24; m = L/2 for L = 18, 30),
N = A_4 triplet at {pi, pi +- 2pi/3} (read at m = L/6), F = A_4 triplet at {0, +-2pi/3} (read at m = L/3).
S is the lowest level in its sector that is not degenerate (1e-8) with the triplet member sharing that sector.
For each beta and L: the crossing t_x where S stops being the lowest branch, the per-site slopes of the two crossing
branches there (linear fits on each side, 3 points), and the slope jump.
Usage: python3 analyze_cross.py file.jsonl [...]"""
import sys, json
import numpy as np
from collections import defaultdict
D = defaultdict(dict)
for fn in sys.argv[1:]:
    for l in open(fn):
        r = json.loads(l)
        D[(r['b'], r['L'])].setdefault(round(-r['g'] / r['b'], 4), {})[r['m']] = r['levels']
def branches(sec, L):
    mS = 0 if L % 12 == 0 else L // 2
    mN, mF = L // 6, L // 3
    if any(m not in sec for m in (mS, mN, mF)):
        return None
    EN, EF = sec[mN][0], sec[mF][0]
    partner = EF if mS == 0 else EN
    cand = [x for x in sec[mS] if abs(x - partner) > 1e-8]
    ES = cand[0] if cand else np.nan
    return ES, EN, EF
for (b, L) in sorted(D):
    ts = sorted(D[(b, L)])
    rows = []
    for t in ts:
        br = branches(D[(b, L)][t], L)
        if br: rows.append((t,) + br)
    if len(rows) < 6: continue
    A = np.array(rows)
    lowest = np.argmin(A[:, 1:], axis=1)
    print(f"beta={b} L={L}: t from {A[0,0]} to {A[-1,0]}, lowest branch sequence: " + "".join("SNF"[k] for k in lowest))
    for i in range(1, len(A)):
        if lowest[i] != lowest[i - 1]:
            k1, k2 = lowest[i - 1], lowest[i]
            # crossing by linear interpolation of the difference
            d0 = A[i - 1, 1 + k1] - A[i - 1, 1 + k2]; d1 = A[i, 1 + k1] - A[i, 1 + k2]
            tx = A[i - 1, 0] + (A[i, 0] - A[i - 1, 0]) * d0 / (d0 - d1)
            lo = max(0, i - 3); hi = min(len(A), i + 3)
            s1 = np.polyfit(A[lo:hi, 0], A[lo:hi, 1 + k1], 1)[0] / L
            s2 = np.polyfit(A[lo:hi, 0], A[lo:hi, 1 + k2], 1)[0] / L
            print(f"   {'SNF'[k1]}->{'SNF'[k2]} at t_x = {tx:.4f}: slopes per site d(E/L)/dt = {s1:.5f} ({'SNF'[k1]}), {s2:.5f} ({'SNF'[k2]}), jump = {s1 - s2:.5f}")
