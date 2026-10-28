#!/usr/bin/env python3
"""Single-quasiparticle dispersion from the parity-resolved ring spectra of a
Z2-symmetric scan (scan.py output).

In a gapped P-even vacuum at momentum 0, the lowest P-odd level at momentum
2 pi m / L is one Bogoliubov quasiparticle, E_qp(k) = E(m, -) - E(0, +), up to
corrections exponentially small in L.  Near the trivial side of the band-edge
Ising line its minimum over k is the Ising mass (at k = 0) or the onset of a
chiral Fermi sea (at k* != 0).  A second local minimum at k* != 0 that reaches
zero before the k = 0 gap closes is what an Ising line continuing inside a
chiral c = 1 region would need.

usage: qp.py data/maps/X.jsonl --L 16 --at D=0.3,h=1.2 [--at ...]
       qp.py data/maps/X.jsonl --L 16 --line XNAME     (minimum structure along the
                                                       lowest-gap point of each column)
"""
import sys, json, argparse, gzip
import numpy as np


def disp(r):
    L = r['L']
    e0 = r['lev']['0,1'][0]
    return np.array([r['lev'][f'{m},-1'][0] - e0 for m in range(L)])


def minima(E):
    """Local minima of a periodic array: (index, value), sorted by value."""
    n = len(E)
    out = [(i, E[i]) for i in range(n) if E[i] <= E[(i - 1) % n] and E[i] <= E[(i + 1) % n]]
    return sorted(out, key=lambda x: x[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fn')
    ap.add_argument('--L', type=int, default=16)
    ap.add_argument('--at', action='append', default=[])
    ap.add_argument('--line', default=None)
    a = ap.parse_args()
    recs = [json.loads(l) for l in (gzip.open(a.fn, 'rt') if a.fn.endswith('.gz') else open(a.fn))]
    recs = [r for r in recs if r['L'] == a.L]
    for spec in a.at:
        want = {s.split('=')[0]: float(s.split('=')[1]) for s in spec.split(',')}
        best = min(recs, key=lambda r: sum((r['p'].get(k, 0) - v) ** 2 for k, v in want.items()))
        E = disp(best)
        ks = (np.arange(a.L) + a.L // 2) % a.L
        print(best['p'], 'GS', best['gs'])
        print('  m :', ' '.join(f'{(m if m <= a.L // 2 else m - a.L):6d}' for m in ks))
        print('  E :', ' '.join(f'{E[m]:6.3f}' for m in ks))
        print('  minima (m, E):', [((i if i <= a.L // 2 else i - a.L), round(v, 4)) for i, v in minima(E)[:3]])
    if a.line:
        cols = sorted({r['p'][a.line] for r in recs})
        for c in cols:
            col = [r for r in recs if r['p'][a.line] == c and r['gs'] == [0, 1]]
            if not col:
                continue
            # the column point with the smallest quasiparticle gap on the vacuum side
            best = min(col, key=lambda r: disp(r).min())
            E = disp(best)
            mins = minima(E)
            print(f'{a.line}={c:.3f}', {k: v for k, v in best['p'].items() if k != a.line},
                  'Emin=%.4f' % E.min(), 'minima', [((i if i <= a.L // 2 else i - a.L), round(v, 3)) for i, v in mins[:3]])


if __name__ == '__main__':
    main()
