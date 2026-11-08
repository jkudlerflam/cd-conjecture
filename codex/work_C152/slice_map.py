"""Print a (Delta, mu) map of ring c_eff (L = 16 -> 20, and 12 -> 16 in brackets when they differ by > 0.3) per slice."""
import json, re, sys
from collections import defaultdict
R = [json.loads(l) for l in open(sys.argv[1])]
S = defaultdict(dict)
for r in R:
    if 'c' not in r: continue
    m = re.match(r"(S\d) V=(\S+) phi=(\S+) Delta=(\S+) mu=(\S+)", r['label'])
    S[m[1] + f" V={m[2]} phi={m[3]}"][(float(m[4]), float(m[5]))] = r
for name in sorted(S):
    d = S[name]; Ds = sorted({k[0] for k in d}); mus = sorted({k[1] for k in d})
    print(name, "  rows Delta, columns mu =", " ".join(f"{m:.2f}" for m in mus))
    for D in Ds:
        row = []
        for m in mus:
            r = d.get((D, m))
            if r is None: row.append("  .  "); continue
            c = r['c'][-1]; row.append(f"{c:5.2f}")
        print(f"{D:5.2f} " + " ".join(row))
    big = [(k, [round(x, 2) for x in r['c']]) for k, r in d.items() if max(r['c']) > 1.25]
    print("  points with c_eff > 1.25:", sorted(big)[:30])
