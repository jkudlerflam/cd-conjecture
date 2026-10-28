#!/usr/bin/env python3
"""Ising crossings from a scan.py map, column by column.

For each value of XNAME, scans YNAME and reports every y where
L1*gap_s(L1) - L2*gap_s(L2) changes sign, for the sector s (default (0,-),
the band-edge Ising sigma), together with the interpolated L2*gap_s, the
ratio sigma/eps at L2, and c_S(L2) at the nearest grid point.  Points whose
ground state is not in (0,+) at both sizes are skipped (incommensurate or
ordered regions), so roots are reported only where the crossing is clean.

usage: maplines.py data/maps/X.jsonl XNAME YNAME [--sig 0,-1]
"""
import sys, json, argparse, gzip
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fn'); ap.add_argument('xn'); ap.add_argument('yn')
    ap.add_argument('--sig', default='0,-1')
    a = ap.parse_args()
    ms, Ps = map(int, a.sig.split(','))
    recs = {}
    for line in (gzip.open(a.fn, 'rt') if a.fn.endswith('.gz') else open(a.fn)):
        r = json.loads(line)
        recs[(r['p'][a.xn], r['p'][a.yn], r['L'])] = r
    Ls = sorted({k[2] for k in recs}); L1, L2 = Ls[0], Ls[-1]
    xs = sorted({k[0] for k in recs}); ys = sorted({k[1] for k in recs})

    def g(r, which):
        L = r['L']; m = 0 if ms == 0 else L // 2
        key = f'{m},{Ps}' if f'{m},{Ps}' in r['lev'] else f'{m},0'
        e = r['lev'][key]
        e0 = r['E0']
        if which == 'sig':
            return L * ((e[0] if (m, Ps) != (0, 1) else e[1]) - e0)
        return L * (r['lev']['0,1' if '0,1' in r['lev'] else '0,0'][1] - e0)
    for x in xs:
        out = []
        prev = None
        for y in ys:
            if (x, y, L1) not in recs or (x, y, L2) not in recs:
                prev = None; continue
            r1, r2 = recs[(x, y, L1)], recs[(x, y, L2)]
            f = g(r1, 'sig') - g(r2, 'sig')
            if prev is not None and np.sign(f) != np.sign(prev[1]) and r1['gs'][0] == 0 and r2['gs'][0] == 0:
                t = prev[1] / (prev[1] - f)
                yc = prev[0] + t * (y - prev[0])
                out.append((round(yc, 3), round(g(r2, 'sig'), 3), round(g(r2, 'sig') / g(r2, 'eps'), 3), round(r2['cS'], 3)))
            prev = (y, f)
        print(f'{a.xn}={x:+.3f}', out)


if __name__ == '__main__':
    main()
