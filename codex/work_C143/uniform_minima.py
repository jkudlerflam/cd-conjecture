"""Distinct uniform product-state minima of a bond and classical two-wall cost on a ring of 2N sites."""
import json, sys
import numpy as np
from scipy.optimize import minimize
from ed import params_to_JB
from classical_landscape import energy
fn, idx = sys.argv[1], int(sys.argv[2])
d = [json.loads(l) for l in open(fn)][idx]
J, B = params_to_JB(np.array(d["x"]))
rng = np.random.default_rng(1)
sols = []
for s in range(200):
    r = minimize(energy, rng.uniform(0, 2 * np.pi, 2), args=(J, B, 1), method="BFGS")
    th, ph = r.x; n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    sols.append((r.fun, n))
e0 = min(s[0] for s in sols)
dist = []
for e, n in sols:
    if e < e0 + 1e-7 and all(np.linalg.norm(n - m) > 1e-3 for m in dist):
        dist.append(n)
print("uniform classical minimum e0 =", e0, "; distinct minimizing Bloch vectors:", [np.round(n, 4).tolist() for n in dist])
if len(dist) >= 2:
    a, b = dist[0], dist[1]
    print("overlap |<alpha|beta>|^2 = (1 + n_a.n_b)/2 =", (1 + a @ b) / 2)
    for N in [6, 10, 16]:
        p = 2 * N
        best = None
        for s in range(20):
            # start: half alpha, half beta
            ns = [a] * N + [b] * N
            th0 = [np.arccos(np.clip(v[2], -1, 1)) for v in ns]; ph0 = [np.arctan2(v[1], v[0]) for v in ns]
            a0 = np.array(th0 + ph0) + 0.05 * rng.normal(size=2 * p)
            r = minimize(energy, a0, args=(J, B, p), method="BFGS")
            best = r.fun if best is None else min(best, r.fun)
        print(f"  ring of {p}: classical energy with two walls minus bulk = {p * (best - e0):.5f}")
