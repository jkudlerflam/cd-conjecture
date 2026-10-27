#!/usr/bin/env python3
"""Order parameters of the finite-chi iMPS states (data/fes_*.json, field 'mag' = <sigma^a_k> on
the n sites of the cell).

The cell magnetization is split into a uniform part m_0 = |mean_k <sigma_k>| and modulated
parts |M_p| = |(1/n) sum_k <sigma_k> e^{-2 pi i p k/n}| (vector norm).  The dominant p gives the
spiral wave vector q = p/n (mod 1, folded to [0, 1/2]).  At a z = 1 critical point a
symmetry-breaking order parameter of scaling dimension x vanishes as xi^{-x}; we fit
ln|M_q| = -x ln xi + const over the last `nfit` bond dimensions.
usage: analyze_order.py [files...] [--nfit 4]
"""
import sys, json, glob, os
import numpy as np


def order(path, nfit=4, errmax=1e-3):
    d = json.load(open(path))
    n = d['n']
    rows = []
    byc = {}
    for r in d['results']:
        byc[r['chi']] = r
    for r in [byc[c] for c in sorted(byc)]:
        if 'mag' not in r or r['err'] > errmax or not np.isfinite(r['xi']) or r['xi'] > 1e7:
            continue
        m = np.array(r['mag'])                      # (3, n)
        F = np.fft.fft(m, axis=1) / n               # (3, n)
        amp = np.linalg.norm(F, axis=0)             # vector norm of each Fourier component
        m0 = amp[0]
        if n > 1:
            p = 1 + int(np.argmax(amp[1:]))
            q = min(p / n, 1 - p / n)
            Mq = amp[p] * (2 if p != n - p else 1)  # +-p pair (real field)
        else:
            q, Mq = 0.0, 0.0
        rows.append(dict(chi=r['chi'], xi=r['xi'], m0=float(m0), q=q, Mq=float(Mq), m0vec=np.real(F[:, 0]).tolist()))
    out = dict(file=os.path.basename(path), model=d['model'], n=n, rows=rows)
    sel = [x for x in rows if x['Mq'] > 1e-8][-nfit:]
    if len(sel) >= 3:
        out['x_q'] = -np.polyfit(np.log([x['xi'] for x in sel]), np.log([x['Mq'] for x in sel]), 1)[0]
    sel0 = rows[-nfit:]
    if len(sel0) >= 3 and all(x['m0'] > 1e-8 for x in sel0):
        out['x_0'] = -np.polyfit(np.log([x['xi'] for x in sel0]), np.log([x['m0'] for x in sel0]), 1)[0]
    return out


if __name__ == '__main__':
    args = sys.argv[1:]
    nfit = 4
    if '--nfit' in args:
        i = args.index('--nfit'); nfit = int(args[i + 1]); del args[i:i + 2]
    files = args or sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'fes_*.json')))
    for f in files:
        o = order(f, nfit)
        print(f"{o['model']:14s} n={o['n']:2d}  x_q={o.get('x_q', float('nan')):.3f}  x_0={o.get('x_0', float('nan')):.3f}")
        for x in o['rows']:
            print(f"    chi={x['chi']:4d} xi={x['xi']:9.2f} q={x['q']:.4f} |M_q|={x['Mq']:.4f} m_0={x['m0']:.4f}")
