#!/usr/bin/env python3
"""C158: central charge from finite-entanglement scaling data written by fesx.py.
Estimators, for each file:
  c_pair   6 dS/dlog xi between consecutive chi (xi from the second transfer eigenvalue, in sites)
  c_ls(k)  least-squares slope 6 dS/dlog xi over the last k chi values
  c_d      the same with xi replaced by 1/delta, delta = eps_2 - eps_1 among distinct moduli (Rams et al.)
  c_kappa  from the Pollmann-Moore exponent kappa = dlog xi/dlog chi over the last 3 chi, inverting
           kappa = 6/(c (sqrt(12/c) + 1))
S is averaged over the bonds of the cell (and the bond spread is reported).
usage: fesfit.py files.json ... [--k 3]"""
import sys, json
import numpy as np
from scipy.optimize import brentq


def kappa_to_c(k):
    f = lambda c: 6 / (c * (np.sqrt(12 / c) + 1)) - k
    try:
        return brentq(f, 1e-3, 50)
    except ValueError:
        return float('nan')


def eps_levels(lam, n):
    m = np.abs(np.array([complex(*z) for z in lam]))
    m = m / m[0]
    e = -np.log(np.maximum(m, 1e-300)) / n
    # distinct moduli
    out = [e[0]]
    for x in e[1:]:
        if x - out[-1] > 1e-7:
            out.append(x)
    return np.array(out)


def analyse(fn, k=3, quiet=False):
    d = json.load(open(fn))
    n = d['args']['n']
    R = d['results']
    chi = np.array([r['chi'] for r in R], float)
    S = np.array([np.mean(r['S']) for r in R])
    Ssp = np.array([np.ptp(r['S']) for r in R])
    xi = np.array([r['xi'] for r in R])
    dl = []
    for r in R:
        e = eps_levels(r['lam'], n)
        dl.append(e[2] - e[1] if len(e) > 2 else np.nan)
    dl = np.array(dl)
    e = np.array([r['e'] for r in R])
    # momentum-resolved correlation lengths: k = 0, k = pi (per cell) and incommensurate
    xik = {'0': [], 'pi': [], 'ic': []}
    for r in R:
        lam = np.array([complex(*z) for z in r['lam']])
        m = np.abs(lam) / abs(lam[0]); ar = np.abs(np.angle(lam[1:])) / np.pi
        xs = -n / np.log(np.maximum(m[1:], 1e-300))
        for key, msk in (('0', ar < 0.05), ('pi', ar > 0.95), ('ic', (ar >= 0.05) & (ar <= 0.95))):
            xik[key].append(float(xs[msk].max()) if msk.any() else np.nan)
    xik = {kk: np.array(v) for kk, v in xik.items()}
    lx = np.log(xi)
    cp = 6 * np.diff(S) / np.diff(lx)
    res = dict(file=fn, chi=chi.tolist(), S=S.tolist(), xi=xi.tolist(), c_pair=cp.tolist())
    for kk in (k, k + 1):
        if len(chi) >= kk:
            res[f'c_ls{kk}'] = float(6 * np.polyfit(lx[-kk:], S[-kk:], 1)[0])
            if np.all(np.isfinite(dl[-kk:])):
                res[f'c_d{kk}'] = float(6 * np.polyfit(-np.log(dl[-kk:]), S[-kk:], 1)[0])
    for key in ('ic',):
        v = xik[key]
        if len(chi) >= k and np.all(np.isfinite(v[-k:])):
            res[f'c_{key}{k}'] = float(6 * np.polyfit(np.log(v[-k:]), S[-k:], 1)[0])
    res['xi_k'] = {kk: v.tolist() for kk, v in xik.items()}
    if len(chi) >= 3:
        kap = np.polyfit(np.log(chi[-3:]), lx[-3:], 1)[0]
        res['kappa'] = float(kap); res['c_kappa'] = kappa_to_c(kap)
    if not quiet:
        ex = d.get('exact_e')
        print(f"== {fn}  n={n}  exact_e={ex}")
        for i, r in enumerate(R):
            ang = np.angle(complex(*r['lam'][1])) / np.pi
            print(f"  chi={int(chi[i]):4d} e={e[i]:.10f}{'' if ex is None else f' (de={e[i]-ex:.1e})'} S={S[i]:.5f}"
                  f"(spread {Ssp[i]:.1e}) xi={xi[i]:9.2f} 1/delta={1/dl[i]:9.2f} arg1/pi={ang:+.4f} err={r['err']:.0e} it={r['iters']}"
                  + f" xi0={xik['0'][i]:.1f} xipi={xik['pi'][i]:.1f} xiic={xik['ic'][i]:.1f}"
                  + (f"  c_pair={cp[i-1]:.3f}" if i > 0 else ''))
        print('  ' + '  '.join(f"{kk}={v:.3f}" for kk, v in res.items() if kk.startswith('c_') and kk != 'c_pair' and np.isscalar(v))
              + (f"  kappa={res['kappa']:.3f}" if 'kappa' in res else ''))
    return res


if __name__ == '__main__':
    k = 3
    args = sys.argv[1:]
    if '--k' in args:
        i = args.index('--k'); k = int(args[i + 1]); del args[i:i + 2]
    for f in args:
        try:
            analyse(f, k)
        except Exception as ex:
            print('fail', f, ex)
