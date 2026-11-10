"""C188 job 2: task lists. tasks_a2.txt (16 single-core slots), ring30_queue.txt (L = 30, 8 threads, sequential)."""
def pt(b, t):
    return f"1,{b},{round(-t*b, 6)}"
dm = ["1,1,1"] + [pt(1.0, t) for t in [0.30, 0.45, 0.50, 0.54, 0.56, 0.58, 0.60, 0.70]] \
   + [pt(0.7, t) for t in [0.55, 0.62, 0.66, 0.70, 0.80]] + [pt(0.4, t) for t in [0.65, 0.75, 0.79, 0.85]] \
   + [pt(0.2, t) for t in [0.60, 0.85]]
A = [f"dmrg 96 16,32,64,128,160,160 {s}" for s in dm] + [f"dmrg 48 16,32,64,128,160,160 {s}" for s in dm]
ts = [0.30, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
for b in [1.0, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1]:
    A += [f"ring24 {pt(b, t)}" for t in ts]
open("tasks_a2.txt", "w").write("\n".join(A) + "\n")
open("ring30_queue.txt", "w").write("\n".join([pt(1.0, 0.5), pt(1.0, 0.45), pt(1.0, 0.55), pt(1.0, 0.40)]) + "\n")
print(len(A), len(dm))
