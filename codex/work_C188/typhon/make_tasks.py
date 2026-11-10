"""C188: task lists for the Typhon job. tasks_a.txt (single-core stream, 16 slots), tasks_b.txt (L = 30 rings, 8 threads, sequential)."""
import numpy as np
def pt(b, t):
    return f"1,{b},{round(-t*b, 6)}"
A = []
# long DMRG first: L = 144 at three points, then L = 96 along three cuts
for spec in ["1,1,1", pt(1.0, 0.45), pt(1.0, 0.5)]:
    A.append(f"dmrg 144 16,32,64,128,192,256,256 {spec}")
cut1 = [0.30, 0.40, 0.45, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.65, 0.70]
cut2 = [0.50, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
dm = [pt(1.0, t) for t in cut1] + [pt(0.7, t) for t in cut2] + [pt(0.4, t) for t in cut2] + ["1,1,1", "1,0.5,0"]
for spec in dm:
    A.append(f"dmrg 96 16,32,64,128,192,256,256 {spec}")
# L = 24 rings: grid and the fine beta = 1 cut, controls
ts = [0.30, 0.40] + [round(x, 3) for x in np.arange(0.45, 1.0001, 0.025)]
pts = ["1,1,1", "1,0.5,0"]
for b in [1.0, 0.7, 0.4, 0.9, 0.8, 0.6, 0.5, 0.3, 0.2, 0.1]:
    pts += [pt(b, t) for t in ts]
pts += [pt(1.0, t) for t in [0.51, 0.52, 0.53, 0.54, 0.56, 0.57, 0.58, 0.59, 0.61, 0.62]]
pts += [pt(0.7, t) for t in [0.62, 0.63, 0.64, 0.66, 0.67, 0.68, 0.69, 0.70]]
pts += [pt(0.4, t) for t in [0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.81, 0.82]]
seen = set()
for p in pts:
    if p not in seen:
        seen.add(p); A.append(f"ring24 {p}")
for spec in dm:
    A.append(f"dmrg 48 16,32,64,128,192,192 {spec}")
open("tasks_a.txt", "w").write("\n".join(A) + "\n")
B = [f"ring30 {pt(1.0, t)}" for t in [0.50, 0.55, 0.45, 0.53, 0.57, 0.40, 0.60]] + [f"ring30 {pt(0.7, t)}" for t in [0.64, 0.66]] + [f"ring30 {pt(0.4, t)}" for t in [0.77]]
open("tasks_b.txt", "w").write("\n".join(B) + "\n")
print(len(A), len(B))
