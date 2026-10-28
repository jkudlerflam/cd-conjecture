#!/usr/bin/env python3
"""Exact ring check of two-body (TDA) candidates, numpy/scipy only.

For a particle bond (twobody.bond) on a ring of L sites, scan mu and record, per momentum
sector, the lowest levels and their particle number <N> (N = number of up spins; not
conserved when gL, gR != 0).  From this:
  * mu_ex: first mu on the grid at which the ring ground state is not the exact vacuum,
    and <N>/L of the new ground state there (a jump to large <N>/L is a first-order,
    phase-separated onset; a dilute onset has <N> of order 1 or 2);
  * at the last mu below mu_ex: the lowest excitation in every sector with its <N>
    (the exact few-body band structure above the vacuum).

usage: screen_ed.py IN_TDA.json OUT.json [L]      (all two-sea entries of a TDA file)
       screen_ed.py --point phi V gL |gR| argR mu_lo mu_hi [L]
"""
import sys, json
import numpy as np
from twobody import bond
from ring_ed import sector_matrix

POP = {}


def popcount(x):
    return np.array([bin(int(v)).count('1') for v in x])


def ring_levels(h, L, nev=3):
    out = []
    for k in range(L):
        H, reps, per = sector_matrix(h, L, k)
        Hd = H.toarray(); Hd = (Hd + Hd.conj().T) / 2
        w, v = np.linalg.eigh(Hd)
        if (L, k) not in POP:
            POP[(L, k)] = popcount(reps)
        n = (np.abs(v[:, :nev]) ** 2 * POP[(L, k)][:, None]).sum(0)
        par = (np.abs(v[:, :nev]) ** 2 * (1 - 2 * (POP[(L, k)][:, None] % 2))).sum(0)
        out.append((w[:nev], n, par))
    return out


def scan(p, mus, L=12):
    rec = []
    for mu in mus:
        h = bond(1.0, p['phi'], p['V'], mu, p['gL'], p['gR'])
        lev = ring_levels(h, L)
        E = np.array([l[0] for l in lev]); N = np.array([l[1] for l in lev]); P = np.array([l[2] for l in lev])
        rec.append(dict(mu=float(mu), E=E.tolist(), N=N.tolist(), P=P.tolist()))
    # vacuum has E = 0 exactly at k = 0
    res = dict(L=L, scan=rec)
    prev = None
    for r in rec:
        E = np.array(r['E']); N = np.array(r['N'])
        k, i = np.unravel_index(np.argmin(E), E.shape)
        if E[k, i] < -1e-9:
            res['mu_ex'] = r['mu']; res['N_gs'] = float(N[k, i]); res['k_gs'] = int(k)
            res['E_gs'] = float(E[k, i])
            break
        prev = r
    if prev is not None:
        E = np.array(prev['E']); N = np.array(prev['N'])
        # lowest non-vacuum level per sector
        exc = []
        for k in range(L):
            j = 1 if k == 0 else 0
            exc.append((float(E[k, j]), float(N[k, j])))
        res['below'] = dict(mu=prev['mu'], exc=exc)
    return res


if __name__ == '__main__':
    if sys.argv[1] == '--point':
        phi, V, gL, gRa, gRp, lo, hi = map(float, sys.argv[2:9])
        L = int(sys.argv[9]) if len(sys.argv) > 9 else 12
        p = dict(phi=phi, V=V, gL=gL, gR=gRa * np.exp(1j * gRp))
        r = scan(p, np.arange(lo, hi + 1e-9, 0.02), L)
        b = r.get('below')
        print('mu_ex', r.get('mu_ex'), 'N_gs', r.get('N_gs'), 'k_gs', r.get('k_gs'))
        if b:
            print('below mu =', b['mu'])
            for k, (e, n) in enumerate(b['exc']):
                print(f"  k={k:2d} K={2*np.pi*k/L:+.3f} E={e:.4f} <N>={n:.3f}")
        sys.exit()
    src, out = sys.argv[1], sys.argv[2]
    L = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    R = json.load(open(src))
    res = []
    for r in R:
        if 'two_sea' not in r:
            continue
        p = dict(r['p']); p['gR'] = complex(*p['gR'])
        ms = r['mu_star']
        s = scan(p, np.arange(ms - 0.6, ms + 0.6 + 1e-9, 0.02), L)
        s.pop('scan')
        res.append(dict(p=r['p'], mu_star_tda=ms, two_sea=r['two_sea'], ed=s))
        json.dump(res, open(out, 'w'))
    json.dump(res, open(out, 'w'))
