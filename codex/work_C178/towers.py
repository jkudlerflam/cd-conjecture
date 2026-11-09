"""Towers: for each pseudospin class (kappa = m mod N, M) of a twist-relaxed spectrum, the levels of the same
(m, theta_j) at which the class minimum is attained (lowest few, all M labels shown). Used to see whether the
kappa = 1 sector holds one boson level (next one at about 3x) or a dense pseudospin tower.
Usage: python3 towers.py file.jsonl"""
import sys, json, numpy as np
import analyze2 as A
for line in open(sys.argv[1]):
    r = json.loads(line)
    L, N = r['L'], r['N']
    ks, es, cs = A.valley_cos(r['phi'], r['wabs'])
    rel = A.relax(r['levels'], r['nth'], N, cs, r['model'], L)
    Eg = min(v[0] for v in rel.values())
    m0 = min(rel, key=lambda k: rel[k][0])[0]
    print(f"== L={L} N={N} |w|={r['wabs']:.4f}  E0={Eg:.5f}")
    for kap in range(N):
        key = ((m0 + kap) % N, 0)
        if key not in rel: continue
        (E, (m, j), d) = rel[key]
        lv = sorted((x[2], x[3] / cs) for x in r['levels'] if x[0] == j and x[1] == m)
        print(f"   kappa={kap} (m={m}, theta_j={j}): " + "  ".join(f"{e-Eg:.4f}[M{M:+.1f}]" for e, M in lv[:6]))
