"""Central-charge estimators on the JSON records written by run.py.

c_ent     periodic ground-state entanglement, S(l) = (c/3) log chord(l) + s0, l in [lmin, L/2]
c_pair    entropy difference between two sizes at fixed l/L (l = L/2 and l = L/4 when integer)
v_adj     L * (lowest level at k0 +/- 2pi/L  - E0) / (2 pi)  (assumes an x = 1 current)
c_cas     from E0(L)/L = e - pi c v / (6 L^2): two-size and three-size (with 1/L^4) fits
z_eff     -log(gap(L2)/gap(L1)) / log(L2/L1) for the lowest gap at any momentum
"""
import json
import numpy as np
from ced import chord


def load(path, **filt):
    out = []
    for line in open(path):
        r = json.loads(line)
        if all(r.get(k) == v for k, v in filt.items()):
            out.append(r)
    return out


def by_L(recs):
    return {r['L']: r for r in recs}


def c_ent(r, lmin=2, renyi=False):
    S = np.array(r['S2'] if renyi else r['S'])
    L = r['L']
    ell = np.arange(1, len(S) + 1)
    m = ell >= lmin
    if m.sum() < 2:
        return np.nan
    x = np.log(chord(ell[m], L))
    A = np.stack([x, np.ones_like(x)], 1)
    coef = np.linalg.lstsq(A, S[m], rcond=None)[0]
    f = 3.0 if not renyi else 4.0   # Renyi-2: S2 = (c/4) log chord
    return float(f * coef[0])


def c_half(r1, r2):
    """3 [S(L2/2) - S(L1/2)] / log(L2/L1)."""
    return float(3 * (r2['S'][-1] - r1['S'][-1]) / np.log(r2['L'] / r1['L']))


def levels(r):
    """List of (E - E0, k) over all momenta, k folded into 0..L/2."""
    out = []
    for k, es in r['spec'].items():
        for e in es:
            out.append((e - r['E0'], int(k)))
    out.sort()
    return out


def fold(k, L):
    k %= L
    return min(k, L - k)


def v_adj(r):
    L, k0 = r['L'], r['k0']
    ks = {fold(k0 + 1, L), fold(k0 - 1, L)}
    g = min(r['spec'][str(k)][0] - r['E0'] for k in ks if str(k) in r['spec'])
    return float(L * g / (2 * np.pi))


def lowest_gap(r):
    lv = [x for x in levels(r) if x[0] > 1e-9]
    return lv[0][0] if lv else np.nan


def z_eff(r1, r2):
    g1, g2 = lowest_gap(r1), lowest_gap(r2)
    return float(-np.log(g2 / g1) / np.log(r2['L'] / r1['L']))


def c_cas(rs, v):
    """rs: records sorted by L (2 or 3). Returns c using velocity v."""
    L = np.array([r['L'] for r in rs], float)
    e = np.array([r['E0'] for r in rs]) / L
    if len(rs) == 2:
        a = (e[0] - e[1]) / (1 / L[0] ** 2 - 1 / L[1] ** 2)
    else:
        A = np.stack([np.ones_like(L), 1 / L ** 2, 1 / L ** 4], 1)
        a = np.linalg.solve(A, e)[1]
    return float(-6 * a / (np.pi * v))


def v_extrap(rs):
    """Velocity extrapolated linearly in 1/L^2 from the last two records."""
    if len(rs) == 1:
        return v_adj(rs[0])
    L1, L2 = rs[-2]['L'], rs[-1]['L']
    v1, v2 = v_adj(rs[-2]), v_adj(rs[-1])
    return float(v2 + (v2 - v1) * (1 / L2 ** 2) / (1 / L1 ** 2 - 1 / L2 ** 2))


def casimir_states(r, k, value, tol=1e-6):
    """Energies (minus E0) at momentum k of states with <C2> = value and small variance."""
    out = []
    cas = r.get('casimir', {}).get(str(k))
    if not cas:
        return out
    for e, (c, var) in zip(r['spec'][str(k)], cas):
        if abs(c - value) < tol and var < 1e-6:
            out.append(e - r['E0'])
    return out


def summary(recs, lmin=2):
    recs = sorted(recs, key=lambda r: r['L'])
    rows = []
    for i, r in enumerate(recs):
        row = dict(L=r['L'], k0=r['k0'], E0L=r['E0'] / r['L'], v=v_adj(r),
                   cent=c_ent(r, lmin) if 'S' in r else np.nan,
                   cent3=c_ent(r, 3) if 'S' in r else np.nan,
                   cren=c_ent(r, lmin, renyi=True) if 'S2' in r else np.nan,
                   gapL=lowest_gap(r) * r['L'] / (2 * np.pi))
        if i >= 1:
            row['z'] = z_eff(recs[i - 1], r)
            row['ccas2'] = c_cas(recs[i - 1:i + 1], v_adj(r))
            if 'S' in r:
                row['chalf'] = c_half(recs[i - 1], r)
        if i >= 2:
            row['ccas3'] = c_cas(recs[i - 2:i + 1], v_extrap(recs[:i + 1]))
        rows.append(row)
    return rows


def fmt(rows):
    keys = ['L', 'k0', 'E0L', 'v', 'gapL', 'z', 'cent', 'cent3', 'cren', 'chalf', 'ccas2', 'ccas3']
    s = ' '.join('%8s' % k for k in keys) + '\n'
    for r in rows:
        s += ' '.join(('%8.4f' % r[k]) if isinstance(r.get(k), float) else '%8s' % r.get(k, '') for k in keys) + '\n'
    return s
