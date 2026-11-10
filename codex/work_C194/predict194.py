"""C194 step 1: C191's KT flow (copied ktflow.py conventions) at the C194 density grid.
For each n: bare (K_s, y1) from C191's hybrid-ABA rows (four g4 variants), flow
dK/dl = -K^2 y^2/2, dy/dl = (2-2K) y from l = 0 at a0 = 1/k_F, k_F = pi n/2.
Runaway scale l* where |y| = 1 gives xi_s ~ a0 e^{l*}; gap Delta_s ~ u_s / xi_s with u_s from C191's velocities."""
import numpy as np, json, sys
from predict import run
def flow(K0, y0, lmax=60, dl=1e-3):
    K, y, l = K0, y0, 0.0
    while l < lmax:
        if abs(y) > 1.0: return l, K
        K, y, l = K - 0.5*K*K*y*y*dl, y + (2-2*K)*y*dl, l + dl
    return None, K
out = []
for n in [float(eval(s)) for s in sys.argv[1:]]:
    o = run(n); kF = np.pi*n/2; a0 = 1/kF
    rows = [r for r in o['rows'] if r['method'] == 'hybrid-ABA']
    rec = dict(n=n, y1=o['y1'], Ks=[r['Ks'] for r in rows], Kc=[r['Kc'] for r in rows], xi=[], us=[min(r['vel']) for r in rows])
    for r in rows:
        ls, Kf = flow(r['Ks'], o['y1'])
        rec['xi'].append(None if ls is None else a0*np.exp(ls))
    margin = [2*(k-1)-abs(o['y1']) for k in rec['Ks']]
    rec['margin'] = margin
    out.append(rec)
    xs = [x for x in rec['xi'] if x is not None]
    print(f"n={n:.4f} y1={o['y1']:.3f} K_s0={min(rec['Ks']):.3f}..{max(rec['Ks']):.3f} margin={min(margin):+.3f}..{max(margin):+.3f} "
          f"xi_s={'gapless' if not xs else f'{min(xs):.3g}..{max(xs):.3g}'}{' (some variants gapless)' if xs and len(xs)<len(rec['xi']) else ''} "
          f"u_s~{np.mean(rec['us']):.2f}", flush=True)
json.dump(out, open('predict194.json', 'w'), indent=1, default=float)
