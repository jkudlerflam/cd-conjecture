"""C188 control: half-ring c_S(12, 18) and Casimir c of the XXZ ring XX + YY + Delta ZZ near the ferromagnetic point
Delta -> -1 (exact c = 1 for -1 < Delta <= 1). Same estimator as ring_cS.py. Usage: python3 xxz_control.py out.jsonl D1 D2 ..."""
import sys, json
import numpy as np
sys.path.insert(0, "../work_C154")
import est_c155 as E
from ringk_run import SX, SY, SZ
out = sys.argv[1]
for D in [float(x) for x in sys.argv[2:]]:
    h = np.kron(SX, SX) + np.kron(SY, SY) + D * np.kron(SZ, SZ)
    S = {}; E0 = {}
    for L in (12, 18):
        d = E.ring_data(h, L, nev=2)
        S[L] = float(E.entropies(d['psi'], L)[L // 2]); E0[L] = float(d['E0'])
    gam = np.arccos(D); v = 2 * np.pi * np.sin(gam) / gam   # exact velocity for Pauli normalization (4x spin-1/2 units: v = 2 pi sin(gam)/gam)
    cS = 3 * (S[18] - S[12]) / np.log(1.5)
    ccas = 6 * (E0[18] / 18 - E0[12] / 12) / (np.pi * (1 / 144 - 1 / 324)) / v
    rec = dict(Delta=D, c_S=cS, c_cas_exact_v=ccas, v=v)
    open(out, "a").write(json.dumps(rec) + "\n"); print(json.dumps(rec), flush=True)
