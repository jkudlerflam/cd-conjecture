#!/usr/bin/env python3
"""Ring-ED estimators for data/ring_NAME_L*.json.

  c_S1, c_S2   fits S_alpha(l) = (c/6)(1 + 1/alpha) ln[(L/pi) sin(pi l/L)] + b over l = 2..L/2
  c_half       3 Delta S1(L/2) / ln(L2/L1) between consecutive sizes
  gaps         lowest levels at momentum k0 +- 1, k0 +- 2 relative to the ground state (k0),
               with L * gap / (2 pi) (constant for z = 1) and L^2 * gap
  Casimir      fit of E0/L = e + a / L^2 over the sizes given, a = -pi sum_a c_a v_a / 6
usage: analyze_ring.py NAME [NAME ...]
"""
import sys, glob, json, re
import numpy as np


def fit_c(S, L, alpha):
    l = np.arange(1, L)
    m = (l >= 2) & (l <= L // 2)
    x = np.log(L / np.pi * np.sin(np.pi * l[m] / L))
    A = np.vstack([x, np.ones_like(x)]).T
    slope = np.linalg.lstsq(A, np.asarray(S)[m], rcond=None)[0][0]
    return slope * 6 / (1 + 1 / alpha)


def load(name):
    out = {}
    for f in glob.glob(f'data/ring_{name}_L*.json'):
        r = json.load(open(f)); out[r['L']] = r
    return dict(sorted(out.items()))


def report(name):
    R = load(name)
    print(f"== {name}  sizes {list(R)}")
    Ls, E0s, Sh = [], [], []
    for L, r in R.items():
        lev = {int(k): np.array(v) for k, v in r['levels'].items()}
        E0 = r['E0']; k0 = r['k0']
        g = {}
        for dk in (1, 2):
            for s in (+1, -1):
                k = (k0 + s * dk) % L
                g[s * dk] = lev[k][0] - E0 if k != k0 else lev[k][1] - E0
        N0 = r['N'][str(k0)][0] if r.get('N') else float('nan')
        c1 = fit_c(r['S1'], L, 1); c2 = fit_c(r['S2'], L, 2)
        print(f"  L={L:2d} k0={k0:2d} E0/L={E0/L:.8f} <N>/L={N0/L:.4f} c_S1={c1:.3f} c_S2={c2:.3f} "
              f"L*gap(+-1)/2pi={L*g[1]/2/np.pi:.4f},{L*g[-1]/2/np.pi:.4f} L*gap(+-2)/2pi={L*g[2]/2/np.pi:.4f},{L*g[-2]/2/np.pi:.4f} "
              f"L^2*gap(1)={L*L*min(g[1], g[-1]):.3f}")
        # second level at k0 +- 1 (a second branch)
        Ls.append(L); E0s.append(E0); Sh.append(r['S1'][L // 2 - 1])
    for i in range(1, len(Ls)):
        print(f"  c_half({Ls[i-1]}->{Ls[i]}) = {3*(Sh[i]-Sh[i-1])/np.log(Ls[i]/Ls[i-1]):.3f}")
    if len(Ls) >= 3:
        Ls_ = np.array(Ls, float); e = np.array(E0s) / Ls_
        A = np.vstack([np.ones_like(Ls_), 1 / Ls_ ** 2]).T
        co = np.linalg.lstsq(A, e, rcond=None)[0]
        print(f"  Casimir: e_inf={co[0]:.8f}, sum c_a v_a = {-6*co[1]/np.pi:.4f}")


if __name__ == '__main__':
    for n in sys.argv[1:]:
        report(n)
