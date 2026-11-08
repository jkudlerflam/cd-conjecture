"""C163 (Claude Opus agent, 2026-11-08): scaled levels of C159's ring scans (L = 16), relative to the ground state.
x_pi = lowest level at relative momentum pi, x_0 = second level at the ground-state momentum, both divided by v
(v = L dE/2pi from the lowest level at relative momentum +-2pi/L, est.py).  For a half-floating liquid in T or K_I,
x_pi <= K (the commensurate vertex e^{2i phi} sits at momentum pi), so x_pi > 15/8 certifies the Ising x boson window."""
import json, glob, sys, numpy as np
from tk import luttinger_K
rows = []
for fn in sorted(glob.glob('c159_scan/*.jsonl')):
    for l in open(fn):
        r = json.loads(l)
        if 'lev16' not in r or 'crit' not in r: continue
        lev = np.array(r['lev16'], float); L = 16; m0 = r['m0'][-1]
        v = r['v']
        if not (v > 0.05): continue
        xpi = lev[(m0 + 8) % L, 0] / v
        x0 = lev[m0, 1] / v
        x1 = min(lev[(m0 + 1) % L, 0], lev[(m0 - 1) % L, 0]) / v
        xmin_other = min(lev[(m0 + j) % L, 0] for j in range(2, 15) if j != 8) / v
        rows.append((r['mode'], r['Delta'], r['harm'], r['crit'], r['c_S'], r['c_cas'], xpi, x0, xmin_other, r['vec']))
modes = sorted({x[0] for x in rows})
for mo in modes:
    R = [x for x in rows if x[0] == mo and x[3] < 0.15 and 0.8 < x[4] < 1.25 and 0.8 < x[5] < 1.3]
    print(f'== {mo}: liquid-like points (crit<0.15, 0.8<c_S<1.25, 0.8<c_cas<1.3): {len(R)}')
    if not R: continue
    xp = np.array([x[6] for x in R]); D = np.array([x[1] for x in R]); hm = np.array([x[2] for x in R])
    for lo, hi in [(-1.0, -0.67), (-0.67, -0.3), (-0.3, 0.3), (0.3, 1.0), (1.0, 99)]:
        s = (D >= lo) & (D < hi)
        if s.sum() == 0: continue
        print(f'  parent Delta [{lo},{hi}): n={s.sum():5d}  x_pi pct 10/50/90 = {np.percentile(xp[s],10):.2f} {np.percentile(xp[s],50):.2f} {np.percentile(xp[s],90):.2f}   n(x_pi>15/8)={np.sum(xp[s]>15/8)}')
    W = [x for x in R if x[6] > 15 / 8]
    print(f'  in window x_pi > 15/8: {len(W)}; their c_S pct 50/95/max = ' + (f'{np.median([x[4] for x in W]):.3f} {np.percentile([x[4] for x in W],95):.3f} {max(x[4] for x in W):.3f}' if W else '-'))
    print(f'     x_0 min over window = ' + (f'{min(x[7] for x in W):.3f}; x_other min = {min(x[8] for x in W):.3f}' if W else '-'))
    # small-harmonic calibration: x_pi against exact K of the parent
    C = [x for x in R if x[2] < 0.15 and -0.95 < x[1] < 0.95]
    if C:
        err = [x[6] - luttinger_K(x[1]) for x in C]
        print(f'  calibration harm<0.15: n={len(C)}, x_pi - K_parent median {np.median(err):+.3f}, pct 10/90 {np.percentile(err,10):+.3f} {np.percentile(err,90):+.3f}')
