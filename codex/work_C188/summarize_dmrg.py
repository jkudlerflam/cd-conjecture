"""C188: summary of dmrg_tau.py output. For each point: E/L, truncation, parity-resolved and averaged CC fits per L,
and the pair estimate c_pair(L1, L2) = 6 [Sbar(L2) - Sbar(L1)] / ln(L2/L1), Sbar = mean of the four central cuts.
Also the staggered and period-3 magnetization amplitudes in the bulk (tau frame). Usage: python3 summarize_dmrg.py file_or_dir [...]"""
import sys, json, glob
import numpy as np
from collections import defaultdict
D = defaultdict(dict)
files = []
for a in sys.argv[1:]:
    files += glob.glob(a + "/dmrg_L*.jsonl") if not a.endswith(".jsonl") else [a]
for fn in files:
    for l in open(fn):
        r = json.loads(l); D[(r['a'], r['b'], r['g'])][r['L']] = r
def key(p):
    return (p[1], -p[2] / p[1] if p[1] else 0) if p != (1.0, 1.0, 1.0) else (9, 9)
for p in sorted(D, key=key):
    t = -p[2] / p[1] if p[1] else 0
    s = f"({p[0]},{p[1]},{p[2]}) t={t:.2f}"
    Ls = sorted(D[p])
    for L in Ls:
        r = D[p][L]; mx = np.array(r['mx']); mz = np.array(r['mz']); mid = slice(L // 4, 3 * L // 4)
        m = max(np.max(np.abs(mx[mid])), np.max(np.abs(mz[mid])))
        s += f" | L{L}: E/L={r['E_per_site']:.7f} tr={r['trunc']:.0e} cfit={r['c_fit'][2]:.3f}({r['c_fit'][0]:.2f},{r['c_fit'][1]:.2f}) S4={r['S_mid4']:.4f} |m|={m:.3f}"
    for i in range(len(Ls) - 1):
        L1, L2 = Ls[i], Ls[i + 1]
        cp = 6 * (D[p][L2]['S_mid4'] - D[p][L1]['S_mid4']) / np.log(L2 / L1)
        s += f" | c_pair({L1},{L2})={cp:.3f}"
    print(s)
