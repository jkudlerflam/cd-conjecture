"""C194: the boundary n_c from the gapped side, KT form (PROOF section 3), registered definition 5.
f(n) = A exp(-B / sqrt(n_c - n)) for n < n_c, 0 for n >= n_c.
Data: (n, value, error). Weighted chi^2 over all densities (gapless ones enter as value ~ 0 +- error).
n_c is profiled on a grid (A, B minimized at each n_c); the 1-sigma interval is chi^2 <= min + 1. A leave-one-out
spread over the densities with value > 3 err is added in quadrature.
Usage: python3 -E -s nc_fit.py   (reads gap_table.json and inf_table.json)."""
import json, numpy as np
from scipy.optimize import minimize


def model(p, n, nc):
    A, B = np.exp(p[0]), p[1]
    x = np.clip(nc - n, 1e-12, None)
    return np.where(n < nc, A * np.exp(-B / np.sqrt(x)), 0.0)


def chi2_at(nc, n, y, e):
    best = None
    for lA in (np.log(max(y.max(), 1e-6)) + s for s in (-2, 0, 2, 4)):
        for B in (0.05, 0.2, 0.5, 1.0):
            r = minimize(lambda p: np.sum(((model(p, n, nc) - y) / e) ** 2), [lA, B], method='Nelder-Mead',
                         options=dict(xatol=1e-8, fatol=1e-10, maxiter=4000))
            if best is None or r.fun < best.fun:
                best = r
    return best.fun, best.x


def fit(n, y, e, grid=np.linspace(0.035, 0.2, 331)):
    n, y, e = map(np.asarray, (n, y, e))
    prof = np.array([chi2_at(nc, n, y, e)[0] for nc in grid])
    i = int(np.argmin(prof)); c0 = prof[i]
    ok = grid[prof <= c0 + 1]
    return dict(nc=float(grid[i]), lo=float(ok.min()), hi=float(ok.max()), chi2=float(c0), dof=len(n) - 3,
                params=chi2_at(grid[i], n, y, e)[1].tolist(), prof=list(zip(grid.tolist(), prof.tolist())))


def run(label, n, y, e):
    n, y, e = map(np.asarray, (n, y, e))
    print(f"--- {label}")
    for a, b, c in zip(n, y, e):
        print(f"   n={a:.4f} value={b:+.5f} +- {c:.5f} ({b / c:+.1f} err)")
    if np.sum(y > 3 * e) < 2:
        print("   fewer than two densities with a value above three errors: the KT fit is not constrained")
        return None
    f = fit(n, y, e)
    loo = []
    for k in range(len(n)):
        if y[k] > 3 * e[k]:
            m = np.arange(len(n)) != k
            if np.sum(y[m] > 3 * e[m]) >= 2:
                loo.append(fit(n[m], y[m], e[m])['nc'])
    spread = max([abs(x - f['nc']) for x in loo], default=0.0)
    err = np.hypot(max(f['hi'] - f['nc'], f['nc'] - f['lo']), spread)
    A, B = np.exp(f['params'][0]), f['params'][1]
    print(f"   n_c = {f['nc']:.4f}  (chi2+1 interval {f['lo']:.4f} to {f['hi']:.4f}; leave-one-out {['%.4f' % x for x in loo]}) "
          f"-> n_c = {f['nc']:.3f} +- {err:.3f};  A = {A:.3g}, B = {B:.3f}, chi2 = {f['chi2']:.2f} for {f['dof']} dof")
    f.update(err=float(err), loo=loo)
    return f


if __name__ == '__main__':
    out = {}
    g = json.load(open('gap_table.json'))
    pts = [(v['n'], v['fit']['D'], v['fit']['err']) for v in g.values() if 'fit' in v]
    if pts:
        out['gap'] = run('Delta_pb(L -> infinity), registered fit', *zip(*sorted(pts)))
    try:
        t = json.load(open('inf_table.json'))
        pts = []
        for tag, v in t.items():
            if tag == 'inf_c1' or 'growth' not in v:
                continue
            xi1 = v['rows'][-1]['xi1']
            # saturating xi_1 (G1-type): 1/xi_1 with the last-doubling change as its error; otherwise 0 +- 1/xi_1
            if v['growth']['xi1'] < 1.15 and v['growth']['xi2'] > 1.5:
                pts.append((v['n'], 1 / xi1, abs(1 / xi1 - 1 / (xi1 / v['growth']['xi1'])) + 1e-4))
            else:
                pts.append((v['n'], 0.0, 1 / xi1))
        if pts:
            out['xi1'] = run('1/xi_1 at the largest chi (saturating densities), 0 +- 1/xi_1 otherwise', *zip(*sorted(pts)))
    except FileNotFoundError:
        pass
    json.dump(out, open('nc_fit.json', 'w'), indent=1, default=float)
