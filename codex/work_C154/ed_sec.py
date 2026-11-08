"""C154: lowest levels of h = a XZ + b YX + g ZY on a ring in chosen momentum sectors (Claude Opus agent, 2026-11-07).
Usage: python3 ed_sec.py L a,b,g m1,m2,... nev out.jsonl     (uses C155's RingTable, est_c155.py, unmodified)"""
import sys, json, time
import numpy as np
import est_c155 as E

L = int(sys.argv[1]); a, b, g = [float(x) for x in sys.argv[2].split(",")]
ms = [int(x) for x in sys.argv[3].split(",")]; nev = int(sys.argv[4]); fout = sys.argv[5]
J = np.zeros((3, 3)); J[0, 2] = a; J[1, 0] = b; J[2, 1] = g
h = E.bond_from_JB(J, np.zeros(3))
t0 = time.time()
T = E.RingTable(L, 2, keep=False)
for m in ms:
    Hk = T.sector_matrix(h, m)
    w, _ = E.lowest(Hk, nev)
    rec = dict(a=a, b=b, g=g, L=L, m=m, levels=[float(x) for x in w], dim=int(Hk.shape[0]), sec=round(time.time() - t0, 1))
    print(json.dumps(rec), flush=True)
    with open(fout, "a") as f:
        f.write(json.dumps(rec) + "\n")
