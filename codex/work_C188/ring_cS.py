"""C188: half-ring entanglement slope c_S = 3 [S(18, 9) - S(12, 6)] / ln(3/2) of the ring ground state, with C155's RingTable
(imported read-only from ../work_C154/est_c155.py, unmodified). Usage: python3 ring_cS.py out.jsonl a,b,g [...]"""
import sys, json, time
import numpy as np
sys.path.insert(0, "../work_C154")
import est_c155 as E
from ringk_run import bond_abg
out = sys.argv[1]
for spec in sys.argv[2:]:
    a, b, g = [float(x) for x in spec.split(",")]
    h = bond_abg(a, b, g); t0 = time.time()
    S = {}
    for L in (12, 18):
        d = E.ring_data(h, L, nev=2)
        S[L] = float(E.entropies(d['psi'], L)[L // 2]); S[f"m0_{L}"] = int(d['m0'])
        S[f"gap_{L}"] = float(L * (np.sort(d['lev'].ravel())[1] - d['E0']))
    cS = 3 * (S[18] - S[12]) / np.log(1.5)
    rec = dict(a=a, b=b, g=g, S12=S[12], S18=S[18], c_S=cS, m0_12=S["m0_12"], m0_18=S["m0_18"], xgap12=S["gap_12"], xgap18=S["gap_18"], sec=round(time.time() - t0, 1))
    open(out, "a").write(json.dumps(rec) + "\n"); print(json.dumps(rec), flush=True)
