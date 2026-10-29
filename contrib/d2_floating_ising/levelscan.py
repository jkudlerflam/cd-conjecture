#!/usr/bin/env python3
"""Neutral q = 0 levels for every liquid-like or Ising-like point of data/ring_summary.json.
For the ground-state momentum sector at L, the lowest --k levels with their parity <P> and magnetization
<M> = sum_j <Z_j>.  x0n = s/sJ of the lowest excited q = 0 level with |M - M0| < 1 and, in the Z2 cuts,
the ground-state parity: the neutral channel where an extra Ising energy field (x = v_I/v) would sit.
One boson gives x0n ~ 2 (J Jbar).  Writes data/levelscan.json.
usage: levelscan.py --L 16 18 --procs 10
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'd2_classification'))
import numpy as np
import ed, models_fi as MF


def one(args):
    t, fam, p, L, m0, sJ, k = args
    J, f = MF.family(fam, p); h = ed.bond(J, f)
    e, u = ed.lowest(ed.sector_hamiltonian(L, h, m0), k=k)
    pc = np.array([bin(i).count('1') for i in range(1 << L)]); P = (-1.0) ** pc; M = L - 2 * pc
    lv = []
    for n in range(len(e)):
        w = np.abs(ed.expand(L, m0, u[:, n])) ** 2
        lv.append(dict(s=float(L * (e[n] - e[0]) / (2 * np.pi)), P=float(w @ P), M=float(w @ M)))
    z2 = fam in ('Z2',) or fam == 'J67'
    neut = [x for x in lv[1:] if abs(x['M'] - lv[0]['M']) < 1.0 and (not z2 or x['P'] * lv[0]['P'] > 0.5)]
    x0n = neut[0]['s'] / sJ if neut else None
    return t, L, dict(levels=lv, x0n=x0n)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--L', type=int, nargs='+', default=[16, 18])
    ap.add_argument('--k', type=int, default=4); ap.add_argument('--procs', type=int, default=4)
    a = ap.parse_args()
    S = json.load(open(os.path.join(HERE, 'data', 'ring_summary.json')))
    tasks = [(t, s['fam'], s['p'], L, s[f'L{L}']['m0'], s[f'L{L}']['sJ'], a.k)
             for t, s in S.items() if s['gapless'] or s['ising'] for L in a.L]
    out = {}
    from multiprocessing import Pool
    with Pool(a.procs) as pool:
        for t, L, r in pool.imap_unordered(one, tasks, chunksize=2):
            out.setdefault(t, {})[f'L{L}'] = r
    json.dump(out, open(os.path.join(HERE, 'data', 'levelscan.json'), 'w'))
    xs = [(t, [out[t].get(f'L{L}', {}).get('x0n') for L in a.L]) for t in out]
    xs = [(t, v) for t, v in xs if all(x is not None for x in v)]
    xs.sort(key=lambda z: min(z[1]))
    print(len(out), 'points; lowest neutral q = 0 levels x0n:')
    for t, v in xs[:25]:
        print(' ', t, S[t]['p'], [round(x, 3) for x in v], 'cS20', round(S[t]['L20']['cS'], 3), 'gapless' if S[t]['gapless'] else 'ising-like')
