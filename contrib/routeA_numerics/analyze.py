#!/usr/bin/env python3
"""Central-charge estimators from the data/ files written by effchain.py.

  FES  : S against ln xi over the VUMPS chi sequence, c = 6 dS/dln(xi) (last four chi, and the
         local slope between the last two), plus the bond-energy modulation inside the cell.
  OBC  : open-chain profile S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + s_{l mod p}, fitted over
         L/8 <= l <= 7L/8 with a separate constant for each residue l mod p (p = 3 for d = 3, else 2),
         which absorbs the parity/trimer oscillation.
  ED   : Casimir c = 6 L^2 (e_inf - E0/L) / (pi v), with v from the k = 2 pi/L level (x = 1) and
         e_inf the chi-extrapolated VUMPS energy (exact for Heisenberg and ULS), and the two-point
         version c_2pt from consecutive sizes, which needs no e_inf; also the ratio of the k = 2 pi/3 level to it
         (2/3 for SU(3)_1) as an operator-content check and gap * L as a z = 1 check.
usage: analyze.py [names ...]   (default: everything in data/)
Written by Claude (Opus 5.5) for contrib/routeA_numerics, 2026-10-27.
"""
import os, sys, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')


def fes(path, last=4):
    d = json.load(open(path))
    R = [r for r in d['results'] if np.isfinite(r['xi']) and r['xi'] > 0]
    if len(R) < 2:
        return None
    S = np.array([np.mean(r['S']) for r in R])
    xi = np.array([r['xi'] for r in R])
    e = np.array([r['e'] for r in R])
    lx = np.log(xi)
    k = min(last, len(R))
    c_fit = 6 * np.polyfit(lx[-k:], S[-k:], 1)[0] if k >= 2 else np.nan
    c_loc = 6 * (S[-1] - S[-2]) / (lx[-1] - lx[-2]) if lx[-1] != lx[-2] else np.nan
    mod = [float(np.ptp(r['ebond'])) for r in R] if 'ebond' in R[0] else []
    # e_inf from e = e_inf + a / xi^2 over the last k points
    A = np.vstack([np.ones(k), xi[-k:] ** -2.0]).T
    einf = np.linalg.lstsq(A, e[-k:], rcond=None)[0][0]
    return dict(chi=[r['chi'] for r in R], S=S.tolist(), xi=xi.tolist(), e=e.tolist(), c_fit=float(c_fit),
                c_last=float(c_loc), mod=mod, e_inf=float(einf), err=[r['err'] for r in R])


def obc(path, p=None):
    d = json.load(open(path))
    r = d['results'][-1]
    S = np.array(r['S'])
    L = len(S) + 1
    if p is None:
        p = 3 if any(t in path for t in ('uls', 'eff', 'blbq', 'amp')) else 2
    l = np.arange(1, L)
    sel = (l >= L / 8) & (l <= 7 * L / 8)
    X = np.log(2 * L / np.pi * np.sin(np.pi * l / L)) / 6
    cols = [X[sel]] + [((l[sel] % p) == q).astype(float) for q in range(p)]
    coef = np.linalg.lstsq(np.vstack(cols).T, S[sel], rcond=None)[0]
    # also the plain fit on the residue class with the largest S (no oscillation term)
    return dict(L=L, chi=r['chi'], trunc=r['trunc'], c=float(coef[0]), chis=[x['chi'] for x in d['results']],
                Smid=[x['S'][L // 2 - 1] for x in d['results']])


def ed(path, einf, kx=1):
    d = json.load(open(path))
    out = []
    for r in d['results']:
        L = r['L']
        E = np.array(r['E']); k = np.array(r['k'])
        g = (E - E[0]) * L / (2 * np.pi)
        k0 = k[0] if np.isfinite(k[0]) else 0.0
        dk = lambda q: abs((q - k0 + L / 2) % L - L / 2)          # momentum relative to the ground state
        cur = [i for i in range(1, len(E)) if np.isfinite(k[i]) and abs(dk(k[i]) - kx) < 0.15]
        pri = [i for i in range(1, len(E)) if np.isfinite(k[i]) and abs(dk(k[i]) - L / 3) < 0.15]
        v = g[cur[0]] if cur else np.nan
        x23 = g[pri[0]] / v if (pri and cur) else np.nan
        c = 6 * L ** 2 * (einf - E[0] / L) / (np.pi * v) if einf is not None and np.isfinite(v) else np.nan
        out.append(dict(L=L, E0L=E[0] / L, v=float(v), c=float(c), x_primary=float(x23),
                        g1=float(g[1]), k1=float(k[1]) if np.isfinite(k[1]) else None))
    # two-point Casimir fit (no e_inf needed): E0/L = e_inf - pi v c / (6 L^2) on consecutive sizes
    for a, b in zip(out[:-1], out[1:]):
        if np.isfinite(b['v']):
            b['c2pt'] = float(-6 / (np.pi * b['v']) * (a['E0L'] - b['E0L']) / (1 / a['L'] ** 2 - 1 / b['L'] ** 2))
    return out


def main():
    names = sys.argv[1:]
    res = {}
    for f in sorted(glob.glob(os.path.join(D, 'fes_*.json'))):
        n = os.path.basename(f)[4:-5]
        if names and n not in names:
            continue
        r = fes(f)
        if r is None:
            continue
        res.setdefault(n, {})['fes'] = r
        print(f"FES {n:22s} chi<= {r['chi'][-1]:4d} xi={r['xi'][-1]:8.2f} S={r['S'][-1]:.4f} "
              f"c_fit(last4)={r['c_fit']:.3f} c_local={r['c_last']:.3f} e_inf={r['e_inf']:.9f} "
              f"bond-mod={r['mod'][-1] if r['mod'] else float('nan'):.2e} err={r['err'][-1]:.0e}")
    for f in sorted(glob.glob(os.path.join(D, 'obc_*.json'))):
        n = os.path.basename(f)[4:-5]
        r = obc(f)
        res.setdefault(n, {})['obc'] = r
        print(f"OBC {n:22s} L={r['L']} chi={r['chi']} trunc={r['trunc']:.1e} c={r['c']:.3f}  S_mid(chi)={np.round(r['Smid'], 4).tolist()}")
    for f in sorted(glob.glob(os.path.join(D, 'ed*_*.json'))):
        n = os.path.basename(f).split('_', 1)[1][:-5]
        key = 'fes_' + n.replace('eff3_', 'eff3_').replace('eff_', 'eff_')
        einf = res.get(n, {}).get('fes', {}).get('e_inf')
        if n == 'uls':
            einf = 1 - np.log(3) - np.pi / (3 * np.sqrt(3))      # Sutherland, h = P_3
        if n == 'heis':
            einf = 0.25 - np.log(2)
        r = ed(f, einf)
        res.setdefault(n, {})['ed_' + os.path.basename(f).split('_')[0]] = r
        for x in r:
            print(f"ED  {n:22s} L={x['L']:3d} E0/L={x['E0L']:.8f} v={x['v']:.4f} x_prim/x_cur={x['x_primary']:.4f} "
                  f"c_Casimir={x['c']:.3f} c_2pt={x.get('c2pt', float('nan')):.3f} first gap*L/2pi={x['g1']:.4f} (k={x['k1']})")
    json.dump(res, open(os.path.join(D, 'summary.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
