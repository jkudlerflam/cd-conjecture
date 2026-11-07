"""Classify search landing points."""
import json, sys, glob
import numpy as np
from local_ioms import comm3
from ed import params_to_JB

def xz_dist(x):
    """distance of normalized J to the set of rank-one a b^T with a orthogonal to b (B ignored), and |B|."""
    J, B = params_to_JB(x)
    u, s, vt = np.linalg.svd(J)
    a, b = u[:, 0], vt[0]
    return s[1] / s[0], abs(a @ b), np.linalg.norm(B) / s[0]

def summarize(d):
    x = np.array(d["x"]); J, B = params_to_JB(x)
    r = d["rings"]
    s = [f"g={d['g_final']:.3f}", f"Jsv={np.round(d['J_sv'],3).tolist()}", f"|DM|={np.linalg.norm(d['DM']):.3f}", f"|B|={np.linalg.norm(B):.3f}",
         f"comm3={comm3(x):.3f}", f"U1u={d['u1_uniform']:.3f}", f"U1s={d['u1_stag']:.3f}",
         "frustr=" + ",".join(f"{q['frustration']:.1e}" for q in r),
         "deg=" + ",".join(str(q['deg']) for q in r),
         "nlow(1M/L)=" + ",".join(str(q['n_low(0.5,1,2 M/L)'][1]) for q in r),
         "S=" + ",".join(f"{q['S']:.2f}" for q in r),
         "gap1*L=" + ",".join(f"{q['gaps'][0]*q['L']:.3f}" for q in r if q['gaps']),
         "xz(s2/s1, a.b, |B|/s1)=" + ",".join(f"{v:.3f}" for v in xz_dist(x))]
    if "S16_trunc" in d: s.append(f"S16={d['S16_trunc']}")
    return "  ".join(s)

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        for line in open(fn):
            d = json.loads(line)
            print(f"[{fn.split('/')[-1]} start {d.get('start', d.get('seed'))}]", summarize(d))
