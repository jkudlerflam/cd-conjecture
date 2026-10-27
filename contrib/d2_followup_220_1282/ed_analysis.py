#!/usr/bin/env python3
"""Ring ED estimators from data/ed/NAME_k6.json (output of d2_classification/run_point.py, six
levels in every momentum sector m, k = 2 pi m / L).

  c_S(L)   S(l) = (c/3) ln[(L/pi) sin(pi l/L)] + b over l in [ceil(L/4), L/2]
  z_n      eps_n(L) = n-th excitation energy above E0 among all sectors; fit ln eps_n = -z ln L
           over L >= Lmin (all L, and the last six sizes)
  soft     the lowest excitations with their momentum relative to the ground state, in units of
           2 pi/L, and x_L = L (E - E0) / (2 pi) (velocity not divided out)
  q_ED     N0/L, where +-N0 is the momentum offset of the lowest "rotor" pair
usage: ed_analysis.py NAME [NAME ...] [--json out.json]
"""
import sys, json, os, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def rel(m, m0, L):
    d = (m - m0) % L
    return d - L if d > L // 2 else d


def analyse(name, Lmin=12, ns=(10, 20, 40)):
    d = json.load(open(os.path.join(HERE, 'data', 'ed', f'{name}_k6.json')))
    Ls = sorted(int(L) for L in d['sizes'])
    rows = []
    for L in Ls:
        r = d['sizes'][str(L)]
        E0, m0 = r['E0'], r['ground_momenta'][0]
        lv = sorted((e, int(m)) for m, es in r['spectrum'].items() for e in es)
        row = dict(L=L, E0=E0, e0=E0 / L, m0=m0, gms=r['ground_momenta'])
        row['soft'] = [(rel(m, m0, L), L * (e - E0) / (2 * np.pi)) for e, m in lv[1:9]]
        row['eps'] = {n: lv[n][0] - E0 for n in ns if n < len(lv)}
        # rotor pair: lowest excitation with |dm| >= 2, and the lowest level at -dm
        exc = [(rel(m, m0, L), L * (e - E0) / (2 * np.pi)) for e, m in lv[1:]]
        cand = [x for x in exc if abs(x[0]) >= 2]
        if cand:
            dm1 = cand[0][0]
            partner = min((x[1] for x in exc if x[0] == -dm1), default=np.nan)
            row['rotor'] = (abs(dm1), cand[0][1], partner)
        S = r.get('S')
        if S is not None:
            S = np.array(S)
            ls = np.arange(max(2, math.ceil(L / 4)), L // 2 + 1)
            x = np.log(L / np.pi * np.sin(np.pi * ls / L)) / 3
            row['c_S'] = float(np.polyfit(x, S[ls - 1], 1)[0])
            row['S_half'] = float(S[L // 2 - 1])
        rows.append(row)
    z = {}
    for n in ns:
        sel = [(r['L'], r['eps'][n]) for r in rows if r['L'] >= Lmin and n in r['eps']]
        if len(sel) < 3:
            continue
        x = np.log([s[0] for s in sel]); y = np.log([s[1] for s in sel])
        p, cov = np.polyfit(x, y, 1, cov=True)
        p6 = np.polyfit(x[-6:], y[-6:], 1)
        z[n] = dict(z=-p[0], err=float(np.sqrt(cov[0, 0])), z_last6=-p6[0])
    return dict(name=name, rows=rows, z=z, note=d.get('meta', {}).get('note', ''))


def show(a):
    print(f"=== {a['name']} {a['note']}")
    for r in a['rows']:
        if r['L'] < 14:
            continue
        rot = r.get('rotor')
        rs = f"rotor dm=+-{rot[0]:2d} x={rot[1]:.3f},{rot[2]:.3f} N0/L={rot[0]/r['L']:.3f}" if rot else ''
        cs = f"c_S={r['c_S']:.3f}" if 'c_S' in r else ' ' * 9
        print(f"L={r['L']:2d} E0/L={r['e0']:.8f} m0={r['m0']:2d} {cs} {rs} | " +
              ' '.join(f"{dm:+d}:{x:.2f}" for dm, x in r['soft'][:6]))
    for n, v in a['z'].items():
        print(f"   z(eps_{n}) = {v['z']:.3f} +- {v['err']:.3f}   last six sizes: {v['z_last6']:.3f}")


if __name__ == '__main__':
    args = sys.argv[1:]
    jout = None
    if '--json' in args:
        i = args.index('--json'); jout = args[i + 1]; del args[i:i + 2]
    res = []
    for n in args:
        a = analyse(n)
        show(a)
        res.append(a)
    if jout:
        json.dump(res, open(jout, 'w'), indent=1, default=float)
