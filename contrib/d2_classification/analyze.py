#!/usr/bin/env python3
"""Finite-size estimators for the D2 classification.

For every chain (results/NAME.json from run_point.py) and every L:
  gap      lowest level above the ground multiplet (any momentum)
  D1, D2   lowest level at relative momentum dm = +-1, +-2 (units 2 pi/L,
           measured from the ground momentum m0)
  v2       L*D2/(4 pi).  In a unitary CFT whose finite-L ground state is the
           vacuum, the lowest spin-2 state has x = 2 exactly (x >= |s| by
           unitarity, and L_{-2}|0> saturates it), so v2 -> v.  With several
           velocities it gives the smallest one.
  r21      D2/D1.  CFT: 2/x_min(spin 1) <= 2.  Quadratic (z = 2) modes: 4.
  c_S      periodic-chain entanglement fit S = (c/3) ln[(L/pi) sin(pi l/L)] + b
           over l in [lmin, L/2], for lmin = 2 and lmin = ceil(L/4).
Across L:
  z        from gaps at (L, L+2):  z = -ln(D(L+2)/D(L)) / ln((L+2)/L)
  c_E      Casimir fit E0/L = e - pi c v /(6 L^2) + B/L^4 over the chosen L set,
           with v = v2 at the largest L (single-velocity assumption).
  c_half   slope of S(L/2) against (1/3) ln L over the last even sizes.
"""
import json, os, sys, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    return json.load(open(os.path.join(HERE, 'results', f'{name}.json')))


def levels(rec, L):
    out = []
    for m, es in rec['spectrum'].items():
        for e in es:
            out.append((e, int(m)))
    out.sort()
    return out


def rel(m, m0, L):
    d = (m - m0) % L
    return d - L if d > L // 2 else d


def per_size(rec, L):
    E0 = rec['E0']; m0 = rec['ground_momenta'][0]
    lv = levels(rec, L)
    tol = 1e-7 * max(1, abs(E0))
    ground = [x for x in lv if x[0] < E0 + tol]
    exc = [x for x in lv if x[0] >= E0 + tol]
    gap = exc[0][0] - E0 if exc else np.nan
    def lowest_at(d):
        c = [e for e, m in exc if abs(rel(m, m0, L)) == d]
        # also ground-multiplet partners at that momentum count as excitations of 0 energy
        c += [e for e, m in ground if abs(rel(m, m0, L)) == d]
        return min(c) - E0 if c else np.nan
    D1, D2, D3 = lowest_at(1), lowest_at(2), lowest_at(3)
    out = dict(L=L, E0=E0, m0=m0, gms=rec['ground_momenta'], gdeg=len(ground), gap=gap,
               D1=D1, D2=D2, D3=D3, v1=L * D1 / (2 * np.pi), v2=L * D2 / (4 * np.pi),
               r21=D2 / D1 if D1 > 0 else np.nan, gapL=gap * L, gapL2=gap * L * L)
    S = rec.get('S')
    if S is not None:
        S = np.array(S)
        for tag, lmin in (('cS2', 2), ('cSq', max(2, math.ceil(L / 4)))):
            ls = np.arange(lmin, L // 2 + 1)
            if len(ls) >= 2:
                x = np.log(L / np.pi * np.sin(np.pi * ls / L)) / 3
                A = np.vstack([x, np.ones_like(x)]).T
                c, b = np.linalg.lstsq(A, S[ls - 1], rcond=None)[0]
                out[tag] = c
            else:
                out[tag] = np.nan
        out['Shalf'] = S[L // 2 - 1]
    return out


def zfit(rows, key='gap', step=2):
    zs = {}
    byL = {r['L']: r for r in rows}
    for L in sorted(byL):
        if L + step in byL:
            a, b = byL[L][key], byL[L + step][key]
            if a > 0 and b > 0:
                zs[L] = -math.log(b / a) / math.log((L + step) / L)
    return zs


def casimir(rows, Ls, v):
    Ls = [L for L in Ls if any(r['L'] == L for r in rows)]
    if len(Ls) < 3:
        return np.nan, np.nan, np.nan
    E = np.array([next(r['E0'] for r in rows if r['L'] == L) for L in Ls])
    L = np.array(Ls, float)
    A = np.vstack([np.ones_like(L), L ** -2, L ** -4]).T
    (e, a, b), *_ = np.linalg.lstsq(A, E / L, rcond=None)
    cv = -6 * a / np.pi
    return e, cv, cv / v


def chalf(rows, n=4):
    ev = sorted([r for r in rows if r['L'] % 2 == 0 and 'Shalf' in r], key=lambda r: r['L'])[-n:]
    if len(ev) < 2:
        return np.nan
    x = np.log([r['L'] for r in ev]) / 3
    y = np.array([r['Shalf'] for r in ev])
    return np.polyfit(x, y, 1)[0]


def analyse(name, Lmin_fit=12, parity=None):
    d = load(name)
    rows = [per_size(d['sizes'][str(L)], int(L)) for L in sorted(map(int, d['sizes']))]
    Lmax = max(r['L'] for r in rows)
    last = rows[-1]
    sel = [r['L'] for r in rows if r['L'] >= Lmin_fit and (parity is None or r['L'] % 2 == parity)]
    e, cv, cE = casimir(rows, sel, last['v2'])
    zg = zfit(rows, 'gap'); z2 = zfit(rows, 'D2')
    return dict(name=name, rows=rows, Lmax=Lmax, e=e, cv=cv, cE=cE, z_gap=zg, z_D2=z2,
                c_half=chalf(rows), meta=d.get('meta', {}))


def show(name, parity=None):
    a = analyse(name, parity=parity)
    print(f"=== {name}  {a['meta'].get('note', a['meta'].get('kind', ''))}")
    print(f"{'L':>3} {'m0/L':>6} {'deg':>3} {'gap*L':>8} {'gap*L^2':>8} {'v1':>6} {'v2':>6} {'r21':>5} {'cS2':>6} {'cSq':>6}")
    for r in a['rows']:
        print(f"{r['L']:3d} {r['m0']/r['L']:6.3f} {r['gdeg']:3d} {r['gapL']:8.4f} {r['gapL2']:8.3f} {r['v1']:6.3f} {r['v2']:6.3f} {r['r21']:5.2f} {r.get('cS2', np.nan):6.3f} {r.get('cSq', np.nan):6.3f}")
    zg = a['z_gap']; z2 = a['z_D2']
    print('z(gap):', {k: round(v, 2) for k, v in zg.items() if k >= 14})
    print('z(D2): ', {k: round(v, 2) for k, v in z2.items() if k >= 14})
    print(f"Casimir: e={a['e']:.8f} cv={a['cv']:.4f} c_E(v2)={a['cE']:.3f}   c_half={a['c_half']:.3f}")
    return a


if __name__ == '__main__':
    for n in sys.argv[1:]:
        show(n)


def level_spacing(name, ns=(4, 10, 20), Lmin=12):
    """eps_n(L): n-th excitation energy above E0 among all momentum sectors (k=4 per sector
    stored, so n <= 4L).  Fit ln eps_n = -z ln L + const over L >= Lmin.  Robust to the
    ground-state level crossings that spoil single-gap ratios in incommensurate chains."""
    d = load(name)
    Ls = sorted(int(L) for L in d['sizes'] if int(L) >= Lmin)
    eps = {n: [] for n in ns}
    for L in Ls:
        rec = d['sizes'][str(L)]
        lv = sorted(e for es in rec['spectrum'].values() for e in es)
        E0 = lv[0]
        for n in ns:
            eps[n].append(lv[n] - E0 if n < len(lv) else np.nan)
    out = {}
    for n in ns:
        y = np.log(np.array(eps[n])); x = np.log(np.array(Ls, float))
        ok = np.isfinite(y)
        p, cov = np.polyfit(x[ok], y[ok], 1, cov=True)
        # local exponents from the last 4 sizes
        p2 = np.polyfit(x[ok][-4:], y[ok][-4:], 1)
        out[n] = dict(z=-p[0], z_err=float(np.sqrt(cov[0, 0])), z_last4=-p2[0],
                      epsL=[e * L for e, L in zip(eps[n], Ls)], Ls=Ls)
    return out
