"""Calibration: can our DMRG + entanglement pipeline see c = 3/2 when it is there?
RZFA spin chain (PRB 92, 235123, Eq. 7) with t = 1:  H = sum_i [ t Z_i - t X_i X_{i+1} - g Z_i Z_{i+1} - g X_i X_{i+2} ].
Spins blocked in pairs (d = 4), so the range-2 XX term is nearest-neighbour in blocks. Intra-block terms are
an on-site MPO term (dmrg_d.py), so the open chain is the plain RZFA open chain.
Usage: python3 rzfa_dmrg.py g Lspins chis out.jsonl.  Expected: c = 3/2 for -2.86 < g < -0.285, c = 1/2 for g > -0.285.
Entropy cuts are between blocks, i.e. at even spin cuts l = 2m."""
import sys, json, time
import numpy as np
from functools import reduce
from dmrg_d import dmrg
X = np.array([[0, 1], [1, 0]], float); Z = np.diag([1., -1]); I = np.eye(2)
def op4(d):  # operator on 4 spins (1,2 | 3,4)
    return reduce(np.kron, [d.get(k, I) for k in range(4)])
g = float(sys.argv[1]); Ls = int(sys.argv[2]); chis = tuple(int(c) for c in sys.argv[3].split(",")); out = sys.argv[4]
t = 1.0
intraL = t * op4({0: Z}) + t * op4({1: Z}) - t * op4({0: X, 1: X}) - g * op4({0: Z, 1: Z})
intraR = t * op4({2: Z}) + t * op4({3: Z}) - t * op4({2: X, 3: X}) - g * op4({2: Z, 3: Z})
inter = -t * op4({1: X, 2: X}) - g * op4({1: Z, 2: Z}) - g * op4({0: X, 2: X}) - g * op4({1: X, 3: X})
h = inter
Ib = np.eye(4)
onsite = t * np.kron(Z, I) + t * np.kron(I, Z) - t * np.kron(X, X) - g * np.kron(Z, Z)
n = Ls // 2
t0 = time.time()
r = dmrg(h, n, chis=chis, onsite=onsite)
S = np.asarray(r["S"]); m = np.arange(1, n); l = 2 * m
x = np.log((2 * Ls / np.pi) * np.sin(np.pi * l / Ls))
fits = {}
for w in (0.1, 0.2, 0.3):
    k = (l >= w * Ls) & (l <= (1 - w) * Ls)
    fits[w] = round(float(6 * np.polyfit(x[k], S[k], 1)[0]), 3)
rec = dict(g=g, L=Ls, chis=list(chis), E=r["E"], hist=r["hist"], S=S.tolist(), c_fit=fits, S_mid=float(S[n // 2 - 1]), sec=time.time() - t0)
open(out, "a").write(json.dumps(rec) + "\n")
print(f"RZFA g = {g} L = {Ls}: E/L = {r['E']/Ls:.8f} trunc = {r['hist'][-1][2]:.1e} dE = {r['hist'][-1][1]-r['hist'][-2][1]:.1e} "
      f"S_mid = {rec['S_mid']:.4f} c_fit(w) = {fits} ({rec['sec']:.0f}s)")
