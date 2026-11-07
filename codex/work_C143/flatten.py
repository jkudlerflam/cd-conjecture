"""From a landing point, minimize the relative width W = (E_D - E_0)/|x| of the lowest D = 2^(L/2 - 1) ring
states (the size of the XZ manifold), at L = 10, and report where it converges: does W -> 0 only on the XZ orbit
(rank-one J with orthogonal vectors and B = 0)?"""
import json, sys
import numpy as np
from scipy.optimize import minimize
from ed import bond, spectrum, all_levels, params_to_JB
fn, L = sys.argv[1], int(sys.argv[2])
D = 2 ** (L // 2 - 1)
x0 = np.array(json.loads(open(fn).readline())["x"])
def W(x):
    x = x / np.linalg.norm(x); J, B = params_to_JB(x)
    w, _ = all_levels(spectrum(bond(J, B), L)); w = np.sort(w)
    return (w[D - 1] - w[0]) / (w[D] - w[0] + 1e-12) if False else (w[D - 1] - w[0])
def report(x, tag):
    x = x / np.linalg.norm(x); J, B = params_to_JB(x)
    u, s, vt = np.linalg.svd(J)
    print(f"{tag}: W={W(x):.3e}  J singular values {np.round(s,4).tolist()}  a.b={abs(u[:,0]@vt[0]):.4f}  |B|/s1={np.linalg.norm(B)/s[0]:.4f}", flush=True)
report(x0, "start")
x = x0.copy()
for rnd in range(3):
    r = minimize(W, x, method="Nelder-Mead", options=dict(maxfev=3000, xatol=1e-7, fatol=1e-10, adaptive=True))
    x = r.x / np.linalg.norm(r.x)
    report(x, f"round {rnd}")
