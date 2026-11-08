"""C163 entanglement-peak-first path scan of the attractive window (Claude Opus agent, 2026-11-08).
A path starts at an attractive U(1) parent (XXZ about z with Delta in (-1, -0.3), spiral angle theta, zero field for
T and K_I, a field along z for N) and switches on non-axial harmonics with amplitude lam (c163lib.path_bond).
Stage 1: open-chain DMRG at L = L1 on a lam grid; Smid(lam) and the profile fit c_prof; ring K_S(12) (kest.py).
Stage 2: an Ising sector that becomes critical inside the liquid adds (1/12) ln min(xi_I, L) to Smid, a bump in
Smid(lam) of relative width ~ 1/L; the two most prominent interior local maxima of Smid (and the maximum of c_prof)
are zoomed twice (lam +- d/3, +- 2d/3, then d/9).
Stage 3: the best zoomed point (largest Smid excess over the linear interpolation of its neighbours, or largest c_prof)
gets DMRG at L = 64, 128 (chi to 128) and K_S(16).
usage: path_scan.py TAG [--L1 48] [--n 16] [--lmax 1.2] [--out paths]"""
import sys, os, json, time, argparse
import numpy as np
import c163lib as C, dmrg_c as DC, tk
from kest import k_est

ap = argparse.ArgumentParser()
ap.add_argument('tag'); ap.add_argument('--L1', type=int, default=48); ap.add_argument('--n', type=int, default=16)
ap.add_argument('--lmax', type=float, default=1.2); ap.add_argument('--out', default='paths')
ap.add_argument('--final', type=int, nargs='+', default=[64, 128]); ap.add_argument('--chi', type=int, default=128)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
spec = [s for s in C.make_specs() if s['tag'] == a.tag][0]
fn = os.path.join(a.out, a.tag + '.json')
res = dict(spec=spec, pts={})
t0 = time.time()


def point(lam, L=None, chis=None, ks=12):
    key = f'{lam:.6f}'
    if key in res['pts'] and L is None:
        return res['pts'][key]
    J, f = C.path_bond(spec, lam)
    v = C.J_to_v(J, f)
    par = tk.spiral_parent(J, f)
    L = L or a.L1
    chis = chis or (16, 32, 48, 48)
    r = DC.run(v, (L,), chis=chis)['L'][L]
    try:
        K = k_est(v, ks)[0]
    except Exception:
        K = float('nan')
    rec = dict(lam=lam, v=v.tolist(), Delta=float(par['Delta']), theta=float(par['theta']), harm=float(par['harm']),
               Smid=r['Smid'], c_prof=r['c_prof'], c_sm=r['c_sm'], trunc=r['trunc'], K_S12=K)
    res['pts'][key] = rec
    json.dump(res, open(fn, 'w'), default=float)
    return rec


grid = np.linspace(0.0, a.lmax, a.n)
for lam in grid:
    point(float(lam))


def excess(lams):
    """Smid minus the linear interpolation of the two neighbours, for interior points of the sorted list."""
    P = sorted(res['pts'].values(), key=lambda r: r['lam'])
    out = []
    for i in range(1, len(P) - 1):
        l0, l1, l2 = P[i - 1]['lam'], P[i]['lam'], P[i + 1]['lam']
        s0, s1, s2 = P[i - 1]['Smid'], P[i]['Smid'], P[i + 1]['Smid']
        lin = s0 + (s2 - s0) * (l1 - l0) / (l2 - l0)
        out.append((s1 - lin, l1, i))
    return out


d = grid[1] - grid[0]
for rnd in range(2):
    ex = sorted(excess(None), reverse=True)
    peaks = [e for e in ex if e[0] > 0.002][:2]
    P = sorted(res['pts'].values(), key=lambda r: r['lam'])
    cmax = max(P, key=lambda r: r['c_prof'])
    centres = [p[1] for p in peaks] + ([cmax['lam']] if cmax['c_prof'] > 1.1 else [])
    dd = d / 3 ** (rnd + 1)
    for c in centres:
        for s in (-2, -1, 1, 2):
            lam = c + s * dd
            if 0 <= lam <= a.lmax:
                point(float(lam))
ex = sorted(excess(None), reverse=True)
P = sorted(res['pts'].values(), key=lambda r: r['lam'])
best = ex[0][1] if ex and ex[0][0] > 0.002 else max(P, key=lambda r: r['c_prof'])['lam']
res['best_lam'] = best
res['excess_top'] = ex[:4]
json.dump(res, open(fn, 'w'), default=float)
J, f = C.path_bond(spec, best)
v = C.J_to_v(J, f)
chis = tuple(c for c in (16, 32, 64, 96, a.chi, a.chi) if c <= a.chi)
res['final'] = DC.run(v, tuple(a.final), chis=chis)
try:
    res['final']['K_S16'] = k_est(v, 16)[0]
except Exception as ex:
    res['final']['K_S16'] = str(ex)
res['time_s'] = time.time() - t0
json.dump(res, open(fn, 'w'), default=float)
print('done', a.tag, f'{time.time() - t0:.0f}s', flush=True)
