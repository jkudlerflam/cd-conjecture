# Summarize C152 job 17597010 (copied to c152_j2_copy/): two-size central charge per label.
import json, glob, math, collections, sys
rows = collections.defaultdict(dict)
for f in sorted(glob.glob('c152_j2_copy/dmrg_*.jsonl')):
    for line in open(f):
        line=line.strip()
        if not line: continue
        try: r=json.loads(line)
        except Exception: continue
        rows[r['label']][r['L']]=r
def sbar(r):
    S=r['S']; L=r['L']; m=L//2
    # mean over 4 central cuts (bonds m-2..m+1), removes parity oscillation
    return sum(S[m-3:m+1])/4.0
for lab in sorted(rows):
    d=rows[lab]; Ls=sorted(d)
    p=d[Ls[0]]['p']
    out=[f"{lab:40s} p=[{', '.join('%.4g'%x for x in p)}]"]
    for L in Ls:
        out.append(f"  L={L:4d} Smid={d[L]['S_mid']:.4f} Sbar={sbar(d[L]):.4f} trunc={d[L]['trunc']:.1e} fit25={d[L]['c_25_75'][0]:.3f}")
    for a,b in zip(Ls,Ls[1:]):
        c=6*(sbar(d[b])-sbar(d[a]))/math.log(b/a)
        out.append(f"  c_pair({a},{b}) = {c:.3f}")
    print("\n".join(out))
