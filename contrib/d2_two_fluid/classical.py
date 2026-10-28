#!/usr/bin/env python3
"""Variational product states for a qubit bond: the lowest energy per site over product
states of period p <= pmax (two angles per site), an upper bound on the true ground-state
energy per site.  Used to locate where a dense, ordered state undercuts the vacuum (energy 0
for the exact-vacuum bonds of twobody.py) and so preempts the dilute regime.

usage: classical.py IN.json OUT.json      (entries with 'p' and 'mu_star', e.g. tda_survivors.json)
"""
import sys, json
import numpy as np
from scipy.optimize import minimize, brentq
from twobody import bond


def energy(x, h, p):
    th, ph = x[:p], x[p:]
    s = np.stack([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)], 1)      # |0> = down, |1> = up
    e = 0.0
    for j in range(p):
        v = np.kron(s[j], s[(j + 1) % p])
        e += np.vdot(v, h @ v).real
    return e / p


def min_classical(h, pmax=6, restarts=6, seed=0):
    rng = np.random.default_rng(seed)
    best = (0.0, None)          # the vacuum (all theta = 0) has energy 0 for exact-vacuum bonds
    best = (energy(np.zeros(2), h, 1), np.zeros(2))
    for p in range(1, pmax + 1):
        for _ in range(restarts):
            x0 = np.concatenate([rng.uniform(0, np.pi, p), rng.uniform(0, 2 * np.pi, p)])
            r = minimize(energy, x0, args=(h, p), method='L-BFGS-B')
            if r.fun < best[0] - 1e-12:
                best = (float(r.fun), r.x)
    return best


def mu_dense(p, lo, hi, tol=1e-4):
    """Smallest mu in [lo, hi] at which some product state goes below the vacuum (energy 0)."""
    f = lambda mu: min_classical(bond(1.0, p['phi'], p['V'], mu, p['gL'], p['gR']))[0] + 1e-7
    if f(lo) < 0:
        return lo
    if f(hi) >= 0:
        return hi
    return brentq(f, lo, hi, xtol=tol)


if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    R = json.load(open(src))
    if len(sys.argv) > 4:                      # chunk i of n
        i, n = int(sys.argv[3]), int(sys.argv[4])
        R = R[i::n]
    res = []
    for r in R:
        p = dict(r['p']); p['gR'] = complex(*p['gR'])
        ms = r['mu_star']
        md = mu_dense(p, ms - 2.0, ms + 3.0)
        res.append(dict(p=r['p'], mu_star=ms, mu_dense_cl=md, window=md - ms,
                        dmu_two_sea=r.get('two_sea', {}).get('dmu')))
        print(json.dumps(res[-1]), flush=True)
        json.dump(res, open(out, 'w'))
