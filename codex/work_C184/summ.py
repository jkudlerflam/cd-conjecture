"""Two-size central charge per label from DMRG jsonl files in a directory (dmrg_*.jsonl, ctrl_*.jsonl).
S-bar = mean of the four central cuts; c_pair(L1,L2) = 6 [Sbar(L2) - Sbar(L1)] / ln(L2/L1)."""
import json, glob, math, collections, sys, os
d = sys.argv[1] if len(sys.argv) > 1 else 'out'
rows = collections.defaultdict(dict)
for f in sorted(glob.glob(os.path.join(d, '*.jsonl'))):
    for line in open(f):
        try: r = json.loads(line)
        except Exception: continue
        rows[r['label']][r['L']] = r
def sbar(r):
    S = r['S']; m = r['L'] // 2
    return sum(S[m - 3:m + 1]) / 4.0
for lab in sorted(rows):
    dd = rows[lab]; Ls = sorted(dd)
    p = dd[Ls[0]].get('p')
    s = f"{lab:34s}" + (" p=[" + ", ".join('%.4f' % x for x in p) + "]" if p else "")
    s += "  " + "  ".join(f"L{L}:S={sbar(dd[L]):.4f},tr={dd[L]['trunc']:.0e},dE={dd[L].get('dE_last',0):.0e}" for L in Ls)
    s += "  | " + "  ".join(f"c({a},{b})={6*(sbar(dd[b])-sbar(dd[a]))/math.log(b/a):.3f}" for a, b in zip(Ls, Ls[1:]))
    print(s)
