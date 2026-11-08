"""Finite-size scaling of open-chain DMRG entropies (dmrg_pt.py records). For each label: L, chi-converged fits, central
average Sbar(L) = mean of S over the four central cuts (removes the parity oscillation to leading order), and
c_pair(L, L') = 6 [Sbar(L') - Sbar(L)] / ln(L'/L) for consecutive sizes. Usage: python3 fss.py files..."""
import json, sys
import numpy as np
from collections import defaultdict
R = defaultdict(dict)
for fn in sys.argv[1:]:
    for line in open(fn):
        try: r = json.loads(line)
        except Exception: continue
        R[r['label']][r['L']] = r
for lab in sorted(R):
    d = R[lab]; Ls = sorted(d)
    sb = {L: float(np.mean(np.array(d[L]['S'])[L // 2 - 3: L // 2 + 1])) for L in Ls}
    pairs = [f"c({Ls[i]},{Ls[i+1]})={6 * (sb[Ls[i+1]] - sb[Ls[i]]) / np.log(Ls[i+1] / Ls[i]):5.2f}" for i in range(len(Ls) - 1)]
    fits = " ".join(f"L{L}:{d[L]['c_25_75'][2]:.2f}/{d[L]['c_12_88'][2]:.2f}" for L in Ls)
    tr = max(d[L]['trunc'] for L in Ls); dE = max(abs(d[L]['dE_last']) for L in Ls)
    print(f"{lab:38s} Sbar " + " ".join(f"{sb[L]:.3f}" for L in Ls) + "  " + " ".join(pairs) + "  fits " + fits + f"  tr={tr:.0e} dE={dE:.0e}")
