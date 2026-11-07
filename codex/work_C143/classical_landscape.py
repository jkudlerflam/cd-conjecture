"""Product-state (classical) energy landscape of a bond: minimize the energy per site of periodic product states
with period p = 1..6 (unit Bloch vectors), and report the minimum per period and the spread of near-optimal
solutions (a degenerate classical manifold shows up as many distinct optima at equal energy)."""
import json, sys
import numpy as np
from scipy.optimize import minimize
from ed import params_to_JB

def energy(angles, J, B, p):
    th, ph = angles[:p], angles[p:]
    n = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    e = 0.0
    for i in range(p):
        e += n[i] @ J @ n[(i + 1) % p] + B @ n[i]
    return e / p

def landscape(x, rng, pmax=6, nst=40):
    J, B = params_to_JB(np.asarray(x))
    res = {}
    for p in range(1, pmax + 1):
        vals = []
        for s in range(nst):
            a0 = rng.uniform(0, 2 * np.pi, 2 * p)
            r = minimize(energy, a0, args=(J, B, p), method="BFGS")
            vals.append(r.fun)
        vals = np.sort(vals)
        res[p] = (vals[0], int(np.sum(vals < vals[0] + 1e-6)))
    return res

if __name__ == "__main__":
    fn, idx = sys.argv[1], int(sys.argv[2])
    d = [json.loads(l) for l in open(fn)][idx]
    rng = np.random.default_rng(0)
    for p, (e, n) in landscape(d["x"], rng).items():
        print(f"  period {p}: min classical energy per site {e:.6f} (hits within 1e-6: {n}/40)")
    print("  quantum E0/L at L = 8, 10, 12:", [round(r["E0_per_L"], 6) for r in d["rings"]])
