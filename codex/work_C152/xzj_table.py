"""Table of the near-XZ family (job 17594701): averaged central-window fit c25 (c12) and S_mid versus h for each (phi, s, V)."""
import json, re, sys
from collections import defaultdict
T = defaultdict(dict)
for fn in sys.argv[1:]:
    for l in open(fn):
        r = json.loads(l); m = re.match(r"xzj_phi(\S+)_s(\S+)_V(\S+)_h(\S+)", r['label'])
        if m: T[(m[1], m[2], m[3])][float(m[4])] = r
hs = [round(0.05 * i, 2) for i in range(11)]
print("phi  s    V     " + "  ".join(f"h={h:<9}" for h in hs))
for k in sorted(T, key=lambda k: (k[0], k[1], float(k[2]))):
    row = []
    for h in hs:
        r = T[k].get(h)
        row.append("    .      " if r is None else f"{r['c_25_75'][2]:4.2f}({r['c_12_88'][2]:4.2f})")
    print(f"{k[0]:4s} {k[1]:3s} {k[2]:5s} " + " ".join(row))
