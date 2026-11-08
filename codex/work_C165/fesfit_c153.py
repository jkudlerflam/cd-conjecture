#!/usr/bin/env python3
"""Central charge from fes.py output (S per bond and xi per chi): pairwise c = 6 dS/dln xi, least squares over
the last 3 and 4 chi, and the Pollmann-Moore kappa = dln xi/dln chi inverted (as in C158's fesfit.py).
Also prints g = b^3 (r - u). Written for codex/work_C153 (Claude Opus agents, 2026-11-07).
usage: fesfit.py files.json ..."""
import sys, json, numpy as np
from scipy.optimize import brentq
def k2c(k):
    try: return brentq(lambda c: 6 / (c * (np.sqrt(12 / c) + 1)) - k, 1e-3, 50)
    except ValueError: return float('nan')
for f in sys.argv[1:]:
    try:
        d = json.load(open(f)); R = d['results']
        chi = np.array([r['chi'] for r in R]); S = np.array([np.mean(r['S']) for r in R]); xi = np.array([r['xi'] for r in R])
        lx = np.log(xi); b, u = d['b'], d['u']; r = np.sqrt(1 + b * b / 4)
        pairs = [6 * (S[i + 1] - S[i]) / (lx[i + 1] - lx[i]) if abs(lx[i + 1] - lx[i]) > 1e-3 else np.nan for i in range(len(S) - 1)]
        ls = {k: 6 * np.polyfit(lx[-k:], S[-k:], 1)[0] if len(S) >= k else np.nan for k in (3, 4)}
        kap = np.polyfit(np.log(chi[-3:]), lx[-3:], 1)[0] if len(S) >= 3 else np.nan
        print(f"{f}: b={b} r-u={r-u:.5f} g={b**3*(r-u):.3f} chi={chi.max()} S={S[-1]:.4f} xi={xi[-1]:.1f} "
              f"c_pair={','.join(f'{p:.2f}' for p in pairs)} c_ls3={ls[3]:.3f} c_ls4={ls[4]:.3f} kappa={kap:.3f} c_kappa={k2c(kap):.2f}")
    except Exception as e:
        print('fail', f, e)
