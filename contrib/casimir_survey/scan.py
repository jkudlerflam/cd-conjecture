#!/usr/bin/env python3
"""Scan a linear family h = sum_i x_i h_i at fixed ring length, reusing momentum blocks.
Writes one JSON line per point in the format of run.py (so analyze.py reads it).

  python3 scan.py --family sp4 --points "[[c,s],...]" --L 8 12 --out data/sp4_scan.jsonl
  python3 scan.py --family sp4 --theta 1.6 3.9 31 --L 8 12 --out ...     (circle families)
"""
import argparse, json, time
import numpy as np
import core, chains as ch
from ced import Sector
import models as cm
from run import sz_charges, colors


def family(name):
    if name == 'sp4':
        m = cm.sp4()
        return dict(d=4, hs=m['hs'], charges=m['charges'], tg=lambda L: [[L % 2, 0]])
    if name == 'blbq':
        x = ch.ss(1.0)
        return dict(d=3, hs=[x, x @ x], charges=sz_charges(1.0), tg=lambda L: [[0]])
    if name == 'spin32':          # P0, P1, P2, P3
        return dict(d=4, hs=ch.spin_projectors(1.5), charges=sz_charges(1.5), tg=lambda L: [[L % 2]])
    if name == 'orbital':         # AB, A, B
        m = cm.spin_orbital()
        return dict(d=4, hs=[x.astype(complex) for x in m['hs']], charges=m['charges'], tg=lambda L: [[L % 2, L % 2]])
    if name.startswith('so'):     # SO(N) vector, (P, K)
        N = int(name[2:]); m = cm.son_vector(N); r = N // 2
        tg = (lambda L: [[0] * r]) if N % 2 else (lambda L: [[L % 2] + [0] * (r - 1)])
        return dict(d=N, hs=m['hs'], charges=m['charges'], tg=tg)
    raise ValueError(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True)
    ap.add_argument('--points', default=None, help='JSON list of coefficient vectors')
    ap.add_argument('--theta', type=float, nargs=3, default=None, help='start stop n on the circle (cos, sin)')
    ap.add_argument('--L', type=int, nargs='+', required=True)
    ap.add_argument('--nlev', type=int, default=2)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    f = family(a.family)
    if a.points:
        pts = json.loads(a.points)
    else:
        th = np.linspace(a.theta[0], a.theta[1], int(a.theta[2]))
        pts = [[float(np.cos(t)), float(np.sin(t))] for t in th]
    norms = [core.hnorm(sum(c * h for c, h in zip(p, f['hs'])), f['d']) for p in pts]
    for L in a.L:
        t0 = time.time()
        secs = [Sector(f['d'], L, f['hs'], charges=f['charges'], target=t) for t in f['tg'](L)]
        for p, hn in zip(pts, norms):
            spec, best = [], (np.inf,)
            for t, sec in zip(f['tg'](L), secs):
                for k in range(L // 2 + 1):
                    e, v, res = sec.solve(p, k, a.nlev)
                    if len(e) == 0:
                        continue
                    spec.append([t, k, [float(x) for x in e]])
                    if e[0] < best[0] - 1e-12:
                        best = (float(e[0]), k, t, v[:, 0], sec)
            psi = best[4].expand(best[1], best[3])
            S = best[4].entropies(psi, cuts=[L // 2])[0]
            r = dict(L=L, k0=best[1], E0=best[0], target=best[2], spec=[[s[0], s[1], s[2]] for s in spec])
            rr = {'L': L, 'k0': best[1], 'E0': best[0], 'target': best[2],
                  'spec': {(tuple(s[0]), s[1]): s[2] for s in spec}}
            rec = dict(model=a.family, par=list(p), L=L, hnorm=hn, E0=best[0], k0=best[1], target=best[2],
                       gap1=core.gap_at(rr, 1), S_half=S, spec=r['spec'])
            with open(a.out, 'a') as fo:
                fo.write(json.dumps(rec) + '\n')
            print(a.family, ['%.4f' % x for x in p], L, 'E0/L=%.6f k0=%d S=%.4f |h|=%.4f' % (best[0] / L, best[1], S, hn), flush=True)
        print('L', L, 'done in %.1fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
