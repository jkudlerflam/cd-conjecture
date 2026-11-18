"""Entry-by-entry comparison of two stored block files (author's scalar9.blocks vs referee's readme_scalar.blocks_mp)."""
import sys, json, os
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from item4_fermi import load
d1, T1, L1, R1 = load(sys.argv[1]); d2, T2, L2, R2 = load(sys.argv[2])
n = d1["n"]
sc = max(abs(T1[i, j]) for i in range(n) for j in range(n))
def md(A, B): return max(abs(A[i, j] - B[i, j]) for i in range(n) for j in range(n)) / sc
out = dict(files=[os.path.basename(sys.argv[1]), os.path.basename(sys.argv[2])], n=n, log10_maxT=mp.nstr(mp.log10(sc), 10),
           dT_over_T=mp.nstr(md(T1, T2), 5), dmuL_over_T=mp.nstr(md(L1, L2), 5), dmuR_over_T=mp.nstr(md(R1, R2), 5),
           muL_over_T=mp.nstr(max(abs(L1[i, j]) for i in range(n) for j in range(n)) / sc, 5),
           muR_over_T=mp.nstr(max(abs(R1[i, j]) for i in range(n) for j in range(n)) / sc, 5))
print(json.dumps(out))
json.dump(out, open(os.path.join(HERE, f"compare_blocks_n{n}.json"), "w"), indent=1)
