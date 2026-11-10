"""C188 control: open-chain DMRG of the uniform XXZ chain XX + YY + Delta ZZ (exact c = 1 for -1 < Delta <= 1),
same estimators as dmrg_tau.py. Usage: python3 xxz_open.py L chis out.jsonl Delta [...]"""
import sys, json, time
import numpy as np
from dmrg3 import dmrg
from dmrg_tau import X, iY, Z, I2, cfit
L = int(sys.argv[1]); chis = tuple(int(c) for c in sys.argv[2].split(",")); fout = sys.argv[3]
for D in [float(x) for x in sys.argv[4:]]:
    ops = [X, iY, Z]; J = [1.0, -1.0, D]
    Ws = []
    for j in range(L):
        W = np.zeros((5, 5, 2, 2)); W[0, 0] = I2; W[4, 4] = I2
        for k in range(3):
            W[0, 1 + k] = J[k] * ops[k]; W[1 + k, 4] = ops[k]
        Ws.append(W[0:1] if j == 0 else (W[:, 4:5] if j == L - 1 else W))
    t0 = time.time(); r = dmrg(None, L, chis=chis, Ws=Ws); S = np.asarray(r["S"])
    rec = dict(Delta=D, L=L, E_per_site=r["E"] / L, trunc=r["hist"][-1][2], c_fit=cfit(S, L, 0.25, 0.75),
               S_mid4=float(np.mean(S[L // 2 - 3:L // 2 + 1])), sec=round(time.time() - t0, 1))
    open(fout, "a").write(json.dumps(rec) + "\n"); print(json.dumps(rec), flush=True)
