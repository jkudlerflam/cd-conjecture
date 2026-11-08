#!/usr/bin/env python3
"""C165: summary of fesx.py output for a two-site cell: staggered moment |m0 - m1|/2, uniform moment,
mean S, xi, leading transfer-eigenvalue phase, pairwise c = 6 dS/dln xi."""
import sys, json, numpy as np
for f in sys.argv[1:]:
    d = json.load(open(f)); R = d['results']
    print('==', f)
    Sp = None
    for r in R:
        m = np.array(r['mag']).real.T  # rows = sites of the cell, columns = (x, y, z)
        N = np.linalg.norm(m[0] - m[1]) / 2; M = np.linalg.norm(m[0] + m[1]) / 2
        S = np.mean(r['S']); lx = np.log(r['xi'])
        lam = np.array([complex(*z) for z in r['lam']]); ang = np.abs(np.angle(lam[1:6])) / np.pi
        cp = '' if Sp is None else f' c_pair={6*(S-Sp[0])/(lx-Sp[1]):.3f}'
        print(f"  chi={r['chi']:4d} |N|={N:.4f} |M|={M:.4f} Nvec={np.round((m[0]-m[1])/2,3)} S={S:.4f} dS={np.ptp(r['S']):.1e} xi={r['xi']:.1f} args/pi={np.round(ang,3)}{cp}")
        Sp = (S, lx)
