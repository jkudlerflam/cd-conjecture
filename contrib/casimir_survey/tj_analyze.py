#!/usr/bin/env python3
"""Supersymmetric t-J chain as a bosonic chain on C^3: h_mu = -P_g + (mu/2)(n_e(x)1 + 1(x)n_e),
P_g the graded swap with the hole bosonic and both electron colours fermionic. Canonical ED
sectors (N_e electrons, S_z = 0 or 1/2) from data/tj.jsonl. For a density n, mu(n) is minus
the centred difference of E_0 in N_e at the largest ring, ||h_mu|| is computed with the bond's
chemical-potential term included, and A comes from E_0(L) at fixed density (1/L^2 fit), with
v_c, v_s read off as the two lowest levels at k0 +- 2 pi/L in the ground-state sector."""
import json
from collections import defaultdict
import numpy as np
import core, chains as ch

E = defaultdict(dict); R = {}
for line in open('data/tj.jsonl'):
    r = json.loads(line)
    L = r['L']; Ne = L - r['target'][0]
    if Ne not in E[L] or r['E0'] < E[L][Ne]:
        E[L][Ne] = r['E0']; R[(L, Ne)] = r


def lev(r, k, n):
    for tg, kk, e in r['spec']:
        if kk % r['L'] == k % r['L'] and len(e) > n:
            return e[n]
    return np.nan


Lbig = max(L for L in E if len(E[L]) == L + 1)
Pg = ch.graded_swap([0, 1, 1])
ne = np.diag([0.0, 1, 1])
for n_num, n_den in ((1, 4), (1, 3), (1, 2), (2, 3), (3, 4)):
    n = n_num / n_den
    Lm = max(L for L in E if len(E[L]) == L + 1 and (n * L) % 1 == 0)
    N0 = round(n * Lm)
    mu = -(E[Lm][N0 + 1] - E[Lm][N0 - 1]) / 2
    h = -Pg + mu * ch.onsite(ne, 3)
    hn = core.hnorm(h, 3)
    rows = []
    for L in sorted(E):
        if (n * L) % 1 or round(n * L) not in E[L]:
            continue
        Ne = round(n * L)
        Eg = {N: e + mu * N for N, e in E[L].items()}
        gs = min(Eg, key=Eg.get)
        r = R[(L, Ne)]
        k0 = r['k0']
        g = sorted([lev(r, k0 + 1, 0), lev(r, k0 + 1, 1), lev(r, k0 - 1, 0), lev(r, k0 - 1, 1)])
        rows.append((L, Ne, E[L][Ne], gs, k0, [(x - r['E0']) * L / (2 * np.pi) for x in g[:2]]))
    print('n=%.3f mu=%.4f ||h_mu||=%.4f (Lmu=%d)' % (n, mu, hn, Lm))
    for L, Ne, e, gs, k0, vv in rows:
        print('   L=%2d Ne=%2d E0=%.8f gc-ground-sector Ne=%d k0=%d  L*gaps/2pi=%s' % (L, Ne, e, gs, k0, np.round(vv, 3)))
    for par in (0, 2):
        sel = [(L, e) for L, Ne, e, gs, k0, vv in rows if Ne % 4 == par and L >= 6]
        if len(sel) >= 3:
            A, e0, b = core.fit_cas([s[0] for s in sel], [s[1] for s in sel])
            print('   Ne = %d mod 4, L=%s: A_cas=%.3f  R=A/(2pi||h_mu||)=%.3f' % (par, [s[0] for s in sel], A, A / (2 * np.pi * hn)))
        elif len(sel) == 2:
            A = core.fit_cas([s[0] for s in sel], [s[1] for s in sel], quartic=False)[0]
            print('   Ne = %d mod 4, L=%s: A_2pt=%.3f  R=%.3f' % (par, [s[0] for s in sel], A, A / (2 * np.pi * hn)))
