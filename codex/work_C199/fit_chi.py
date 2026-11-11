"""C199b: finite-bond-dimension extrapolation of S̄₆(L) on the ground branch, and c_two from the extrapolated entropies.
Data (S̄₆, max discarded weight of the last sweep of each χ stage) from codex/work_C190/gpu typhon_out (outG, outX1024,
outX1536) and the C195/C195b runs at L = 512 as of 2026-11-11 07:00 EST (sym40_n316 at L = 512 is interim, inside its
χ = 1536 stage). Fit forms: (a) S = S∞ − a χ^(−b); (b) S = S∞ − A ε^α; (c) linear in ε through the last two stages
(the lower estimate of R8). Stages that were not run to a stop rule at L = 512 (χ ≤ 512) are not used there."""
import numpy as np
from scipy.optimize import brentq
D = {
 'p153c': {128: [(256, 1.8385, 3.6e-9), (384, 1.8545, 7.9e-10), (512, 1.8591, 2.6e-10)],
           256: [(512, 2.0559, 1.2e-9), (768, 2.0828, 3.5e-10), (1024, 2.0899, 1.3e-10), (1536, 2.0938, 2.6e-11)],
           512: [(768, 2.2327, 8.9e-10), (1024, 2.2729, 4.4e-10), (1536, 2.3077, 1.5e-10)]},
 'sym40_n316': {128: [(256, 1.8138, 2.9e-9), (384, 1.8274, 6.3e-10), (512, 1.8310, 2.0e-10)],
           256: [(512, 2.0322, 1.0e-9), (768, 2.0559, 2.9e-10), (1024, 2.0619, 1.0e-10), (1536, 2.0650, 2.0e-11)],
           512: [(768, 2.2140, 7.9e-10), (1024, 2.2524, 3.8e-10), (1536, 2.2837, 1.3e-10)]}}
def three_point(x, y):
    """S = S∞ − a x^(−b) through three points (x increasing): solve for b from the ratio of differences."""
    (x1, x2, x3), (y1, y2, y3) = x, y
    r = (y3 - y2) / (y2 - y1)
    f = lambda b: (x2 ** -b - x3 ** -b) / (x1 ** -b - x2 ** -b) - r
    try: b = brentq(f, 1e-3, 20)
    except ValueError: return None
    a = (y2 - y1) / (x1 ** -b - x2 ** -b)
    return y3 + a * x3 ** -b, b
def fits(rows):
    chi = np.array([r[0] for r in rows], float); S = np.array([r[1] for r in rows]); e = np.array([r[2] for r in rows])
    out = {}
    t = three_point(chi[-3:], S[-3:]); out['chi^-b'] = t
    t = three_point(1 / e[-3:], S[-3:]); out['eps^a'] = t  # x = 1/eps, so S = S∞ − A eps^a
    out['lin eps'] = (S[-1] + (S[-1] - S[-2]) * e[-1] / (e[-2] - e[-1]), 1.0)
    return out, S[-1]
k = 6 / np.log(2)
for key, d in D.items():
    print(f'== {key}')
    ext = {}
    for L, rows in d.items():
        f, raw = fits(rows); ext[L] = (f, raw)
        print(f'  L = {L:3d}: raw {raw:.4f}  ' + '  '.join(f'{n}: {v[0]:.4f} (exp {v[1]:.2f})' if v else f'{n}: no fit' for n, v in f.items()))
    for a, b in ((128, 256), (256, 512)):
        s = f'  c_two({a}->{b}): raw {k * (ext[b][1] - ext[a][1]):.3f}'
        for n in ('lin eps', 'eps^a', 'chi^-b'):
            va, vb = ext[a][0][n], ext[b][0][n]
            s += f', {n} {k * (vb[0] - va[0]):.3f}' if va and vb else f', {n} n/a'
        print(s)
