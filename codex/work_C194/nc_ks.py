"""C194: n_c from the gapless side (registered definition 5). K_s(n) at the smallest measured L >= 120 of each
density (L = 120 to 128, so all densities are compared at the same scale), weighted linear and quadratic fits,
crossings of K_s = 1 and of K_s = 1 + |y1(n)|/2 with C191's bare y1 (predict194.json, interpolated).
Errors by Monte Carlo over the K_s errors (2000 draws). Also the unregistered (K_s - 1)^2 linear fit."""
import json, numpy as np
ks = json.load(open('ks_table.json')); pr = json.load(open('predict194.json'))
pn = np.array([p['n'] for p in pr]); py = np.array([abs(p['y1']) for p in pr])
y1 = lambda n: np.interp(n, pn, py)
n = []; K = []; E = []
for v in ks.values():
    r = v['rows'][0]
    n.append(v['n']); K.append(r['Ks']); E.append(max(r['Ks_err'], 0.01))
n, K, E = map(np.array, (n, K, E)); o = np.argsort(n); n, K, E = n[o], K[o], E[o]
print('input:', ' '.join(f'n={a:.4f}:{b:.3f}+-{c:.3f}' for a, b, c in zip(n, K, E)))
grid = np.linspace(0.0, 0.2, 4001)
def cross(f, target):
    d = f(grid) - target(grid); s = np.where(np.diff(np.sign(d)) != 0)[0]
    return grid[s[-1]] if len(s) else np.nan
rng = np.random.default_rng(1)
for deg in (1, 2):
    res = {'1': [], 'sep': [], 'sq': []}
    for k in range(2001):
        Kk = K if k == 0 else K + E * rng.standard_normal(len(K))
        p = np.polyfit(n, Kk, deg, w=1 / E)
        f = lambda x: np.polyval(p, x)
        res['1'].append(cross(f, lambda x: np.ones_like(x)))
        res['sep'].append(cross(f, lambda x: 1 + y1(x) / 2))
        if deg == 1:
            q = np.polyfit(n, (Kk - 1) ** 2 * np.sign(Kk - 1), 1, w=1 / (2 * np.abs(K - 1) * E + 1e-3))
            res['sq'].append(-q[1] / q[0])
    for key, v in res.items():
        if not v: continue
        v = np.array(v); c = v[0]; s = np.nanstd(v[1:])
        print(f"deg {deg} crossing {'K_s = 1' if key == '1' else ('K_s = 1 + |y1|/2' if key == 'sep' else '(K_s-1)^2 -> 0 (unregistered)')}: n = {c:.4f} +- {s:.4f}")
