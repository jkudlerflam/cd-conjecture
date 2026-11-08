"""Open-chain Casimir analysis E0(L) = e L + e_b + f / L, f = -pi (sum_s c_s v_s) / 24, from DMRG energies.
Records: dmrg_run.py jsonl (Pauli units, field key 'label') and C156 kit110 jsonl (spin-1/2 units, E x 4 -> Pauli).
Fits all triples of consecutive sizes, and a least-squares fit over all sizes. Reference: XZ point (A critical,
c = 1/2, v = 4) has f = -pi/12 = -0.2618; a gapped A with only a slow sector would give |f| of order v_slow/24."""
import sys, json, glob
import numpy as np
from collections import defaultdict
data = defaultdict(dict)
for f in sys.argv[1:]:
    for line in open(f):
        d = json.loads(line)
        if "h" in d and "label" not in d:
            key = round(float(d["h"]), 4); E = 4 * d["E"]
        else:
            lab = d["label"]
            if not lab.startswith("bis_h"):
                continue
            key = round(float(lab[5:]), 4); E = d["E"]
        L = d["L"]
        data[key].setdefault(L, []).append(E)
for h in sorted(data):
    Ls = sorted(data[h]); Es = [np.mean(data[h][L]) for L in Ls]
    print(f"h = {h}: sizes {Ls}")
    for L in Ls:
        if len(data[h][L]) > 1:
            print(f"   L = {L}: duplicate runs differ by {np.ptp(data[h][L]):.2e}")
    for i in range(len(Ls) - 2):
        A = np.array([[L, 1, 1 / L] for L in Ls[i:i + 3]]); e, eb, ff = np.linalg.solve(A, Es[i:i + 3])
        print(f"   triple {Ls[i:i+3]}: e = {e:.8f}  e_b = {eb:+.5f}  f = {ff:+.4f}  -> sum c v = {-24 * ff / np.pi:+.3f}")
    if len(Ls) >= 4:
        A = np.array([[L, 1, 1 / L] for L in Ls]); sol, *_ = np.linalg.lstsq(A, Es, rcond=None)
        print(f"   lstsq all: e = {sol[0]:.8f}  f = {sol[2]:+.4f}  -> sum c v = {-24 * sol[2] / np.pi:+.3f}")
