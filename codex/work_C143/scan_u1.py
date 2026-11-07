"""Grid scans of the singular-point indicator g (v2, L = 8 -> 10 -> 12, kappa = 3) on two exactly
characterized two-parameter families.
 uniform U(1):   h = cos(t)(XX+YY) + sin(t) ZZ + B (Z1+Z2)/2   (Delta = tan t)
   predicted singular set: saturation B = 2(cos t + sin t) for Delta > -1, lower critical field B_c1(Delta) for
   Delta > 1, FM point t = -pi/4, B = 0, and the Fibonacci points as t -> +pi/2 with B = 2 sin t.
 staggered U(1): h = cos(t)(XX-YY) + sin(t) ZZ + B (Z1+Z2)/2   (conserves sum_j (-1)^j Z_j)
   predicted singular set: the Bravyi-Gosset curve B = 2 sqrt(sin^2 t - cos^2 t) for sin t >= cos t > 0
   (from the FM Heisenberg point t = pi/4, B = 0 to the Fibonacci point t = pi/2, B = 2).
Output: rows of (family, t, B, g, S8, S10, S12)."""
import sys, json
import numpy as np
import search

def x_of(fam, t, B):
    J = np.zeros((3, 3)); J[2, 2] = np.sin(t)
    J[0, 0] = np.cos(t); J[1, 1] = np.cos(t) if fam == "uniform" else -np.cos(t)
    return np.concatenate([J.ravel(), [0, 0, B]])

fam = sys.argv[1]
ts = np.linspace(-np.pi / 2 + 0.05, np.pi / 2, 24)
Bs = np.linspace(0, 3.0, 25)
out = open(f"scan_{fam}.jsonl", "w")
for t in ts:
    for B in Bs:
        x = x_of(fam, t, B); x = x / np.linalg.norm(x)
        S = [search.softcount(x, L) for L in (8, 10, 12)]
        g = float(min(np.diff(np.log(S))))
        out.write(json.dumps({"t": float(t), "B": float(B), "g": g, "S": S}) + "\n"); out.flush()
