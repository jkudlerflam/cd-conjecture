#!/usr/bin/env python3
"""Maps from a scan.py output: c_S at the two sizes, the ground-state sector,
and two gap indicators.

usage: mapplot.py data/maps/X.jsonl XNAME YNAME [--png out.png]

Indicators (L1 < L2 the two sizes):
  cS(L)      chord-fit entanglement slope of the ring ground state.
  rgap       ratio L2*g(L2) / (L1*g(L1)) of the lowest excitation above the
             quasi-degenerate ground manifold (levels within 1e-6 * L of E0 are
             counted in the manifold).  About 1 for a critical point, about
             L2/L1 or larger in a gapped phase.
  nsoft(L)   number of levels (over all sectors) with L*(E-E0) < 2.
matplotlib is used only for the picture.
"""
import sys, json, argparse, gzip
import numpy as np


def load(fn):
    recs = {}
    for line in (gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn)):
        r = json.loads(line)
        recs[(json.dumps(r['p'], sort_keys=True), r['L'])] = r
    return recs


def levels(r):
    out = []
    for key, v in r['lev'].items():
        m, P = map(int, key.split(','))
        for e in v:
            if e == e:
                out.append((e - r['E0'], m, P))
    out.sort()
    return out


def indicators(r):
    L = r['L']
    lv = levels(r)
    gaps = np.array([x[0] for x in lv])
    manifold = gaps < 1e-6 * L
    exc = gaps[~manifold]
    g = exc[0] if len(exc) else np.nan
    return dict(cS=r['cS'], g=L * g, ndeg=int(manifold.sum()), nsoft=int((L * gaps < 2.0).sum()),
                gs=tuple(r['gs']))


def grid(recs, xn, yn):
    Ls = sorted({k[1] for k in recs})
    pts = sorted({k[0] for k in recs})
    xs = sorted({json.loads(p)[xn] for p in pts}); ys = sorted({json.loads(p)[yn] for p in pts})
    return Ls, pts, np.array(xs), np.array(ys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fn'); ap.add_argument('xn'); ap.add_argument('yn')
    ap.add_argument('--png', default=None)
    a = ap.parse_args()
    recs = load(a.fn)
    Ls, pts, xs, ys = grid(recs, a.xn, a.yn)
    L1, L2 = Ls[0], Ls[-1]
    shape = (len(ys), len(xs))
    M = {k: np.full(shape, np.nan) for k in ('cS1', 'cS2', 'rgap', 'gsm', 'nsoft')}
    for p in pts:
        d = json.loads(p)
        i, j = np.searchsorted(ys, d[a.yn]), np.searchsorted(xs, d[a.xn])
        if (p, L1) in recs and (p, L2) in recs:
            A, B = indicators(recs[(p, L1)]), indicators(recs[(p, L2)])
            M['cS1'][i, j] = A['cS']; M['cS2'][i, j] = B['cS']
            M['rgap'][i, j] = B['g'] / A['g']
            M['gsm'][i, j] = B['gs'][0] / L2 * (1 if B['gs'][1] >= 0 else -1)
            M['nsoft'][i, j] = B['nsoft']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    ext = [xs[0], xs[-1], ys[0], ys[-1]]
    for axx, key, title, kw in ((ax[0, 0], 'cS2', f'c_S L={L2}', dict(vmin=0, vmax=2, cmap='viridis')),
                                (ax[0, 1], 'cS1', f'c_S L={L1}', dict(vmin=0, vmax=2, cmap='viridis')),
                                (ax[1, 0], 'rgap', f'gap ratio {L2}/{L1}', dict(vmin=0.5, vmax=2, cmap='coolwarm')),
                                (ax[1, 1], 'gsm', f'GS momentum/2pi (sign = parity) L={L2}', dict(vmin=-0.5, vmax=0.5, cmap='twilight'))):
        im = axx.imshow(M[key], origin='lower', extent=ext, aspect='auto', interpolation='nearest', **kw)
        axx.set_title(title); axx.set_xlabel(a.xn); axx.set_ylabel(a.yn); plt.colorbar(im, ax=axx)
    fig.suptitle(a.fn)
    fig.tight_layout()
    fig.savefig(a.png or a.fn.replace('.gz', '').replace('.jsonl', '.png'), dpi=80)


if __name__ == '__main__':
    main()
