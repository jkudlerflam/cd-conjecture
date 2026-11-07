"""Continue the maximization of g (v2, L = 8,10,12) from a landing point, with a smaller kappa to sharpen
near-degeneracies, and print the trajectory's diagnostics."""
import json, sys
import numpy as np
from scipy.optimize import minimize
import search
from analyze import summarize
fn, idx, kappa, maxfev = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
d = [json.loads(l) for l in open(fn)][idx]
x0 = np.array(d["x"])
search.KAPPA = kappa
f = lambda x: -search.objective(x / np.linalg.norm(x), 8, 12)
r = minimize(f, x0, method="Nelder-Mead", options=dict(maxfev=maxfev, xatol=1e-6, fatol=1e-8, adaptive=True))
x = r.x / np.linalg.norm(r.x)
search.KAPPA = 3.0
out = {"g_final": -r.fun, "kappa_follow": kappa}
out.update(search.diagnostics(x))
print(json.dumps(out))
print(summarize(out))
J, B = search.params_to_JB(x)
print("J=\n", np.round(J, 4), "\nB=", np.round(B, 4))
print("dist from start", np.linalg.norm(x - x0 / np.linalg.norm(x0)))
