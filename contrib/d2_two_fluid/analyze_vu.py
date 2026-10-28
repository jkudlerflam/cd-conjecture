#!/usr/bin/env python3
"""VUMPS finite-entanglement scaling summary for data/vu_NAME.json (run_vumps.py).
c(chi) = 6 Delta S / Delta ln xi between consecutive chi (S averaged over the cell bonds),
and a least-squares c over the last three chi.  For d = 2 also the mean density
n = (1 - <sigma^z>)/2 and the in-plane moment |(<sigma^x>, <sigma^y>)| over the cell.
usage: analyze_vu.py NAME [NAME ...]   (NAME may contain a glob)"""
import sys, glob, json
import numpy as np


def report(path, quiet=False):
    r = json.load(open(path)); res = r['results']
    out = []
    head = f"== {r['spec']}  n={r['n']}"
    prev = None
    for x in res:
        S = np.mean(x['S']); xi = x['xi']
        line = f"  chi={x['chi']:3d} e={x['e']:.10f} err={x['err']:.0e} S={S:.4f} xi={xi:9.2f}"
        if 'mag' in x:
            m = np.array(x['mag'])          # (3, n)
            line += f" n={np.mean((1 - m[2]) / 2):.4f} |m_xy|={np.mean(np.hypot(m[0], m[1])):.4f}"
        if prev is not None and xi > prev[1] * 1.0001:
            line += f"  c(chi)={6 * (S - prev[0]) / np.log(xi / prev[1]):.3f}"
        out.append(line); prev = (S, xi)
    cfit = None
    if len(res) >= 3:
        S = np.array([np.mean(x['S']) for x in res[-3:]]); X = np.array([x['xi'] for x in res[-3:]])
        A = np.vstack([np.log(X), np.ones(3)]).T
        cfit = 6 * np.linalg.lstsq(A, S, rcond=None)[0][0]
        out.append(f"  c (fit, last three chi) = {cfit:.3f}")
    if not quiet:
        print(head); print('\n'.join(out))
    return cfit


if __name__ == '__main__':
    for pat in sys.argv[1:]:
        for f in sorted(glob.glob(f'data/vu_{pat}.json')):
            report(f)
