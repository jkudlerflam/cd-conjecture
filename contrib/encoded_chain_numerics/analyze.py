#!/usr/bin/env python3
"""Analysis of the VUMPS (fes_*.json, eff_*.json) and open-chain (obc_*.json) runs.

FES: S = cell-boundary entropy averaged over the cells of the unit cell, xi in cells.
  c_pair = 6 dS / d ln xi between consecutive chi;  c_fit = 6 x slope over the last `k` chi.
  kappa  = d ln xi / d ln chi (CFT: 6 / (c (sqrt(12/c) + 1)), Pollmann et al. 2009).
  p      = exponent of e(chi) - e_inf ~ xi^-p from a 3-parameter fit to the last 4 chi (z = p - 1).
  dS     = max over cells of |S_j - mean| (one-cell translation breaking), d_eb same for bond energies.
OBC: Calabrese-Cardy fit S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + s0 on cell-boundary bonds with
  L/4 <= l <= 3L/4, plain and with a staggered term (-1)^l A/[...]; gaps L (E1 - E0).

usage: analyze.py fes data/fes_*.json      analyze.py obc data/obc_*.json
Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import sys, json, glob, os
import numpy as np
from scipy.optimize import curve_fit


def fes_table(path, k=3):
    d = json.load(open(path))
    a = d['args']
    nc = a['ncells']
    sub = len(d['dims']) // nc
    rows = []
    for r in d['results']:
        S0 = [r['S'][sub * j] for j in range(nc)]
        row = dict(chi=r['chi'], S=float(np.mean(S0)), dS=float(np.max(np.abs(np.array(S0) - np.mean(S0)))),
                   xi=r['xi_cells'], e=r['e_cell'], err=r['err'], it=r['iters'])
        if 'e_bond' in r:
            eb = np.array(r['e_bond']); row['d_eb'] = float(np.max(np.abs(eb - eb.mean())))
            row['eb'] = float(eb.mean())
        if 'charge' in r:
            row['Q'] = float(np.mean(r['charge']))
        if 'leakL' in r:
            row['leak'] = float(np.mean(r['leakL']))
        lam1 = r['lam'][1]
        row['phase'] = float(np.arctan2(lam1[1], lam1[0]) / np.pi)
        rows.append(row)
    out = dict(path=os.path.basename(path), args=a, rows=rows)
    if len(rows) >= 2:
        S = np.array([x['S'] for x in rows]); lx = np.log([x['xi'] for x in rows])
        out['c_pair'] = (6 * np.diff(S) / np.diff(lx)).tolist()
        kk = min(k, len(rows))
        out['c_fit'] = float(6 * np.polyfit(lx[-kk:], S[-kk:], 1)[0])
        lc = np.log([x['chi'] for x in rows])
        out['kappa'] = float(np.polyfit(lc[-kk:], lx[-kk:], 1)[0])
        if len(rows) >= 4:
            xi = np.array([x['xi'] for x in rows[-4:]]); e = np.array([x['e'] for x in rows[-4:]])
            try:
                f = lambda x, einf, A, p: einf + A * x ** (-p)
                popt, _ = curve_fit(f, xi, e, p0=(e[-1], (e[0] - e[-1]) * xi[0] ** 2, 2.0), maxfev=20000)
                out['p'] = float(popt[2]); out['e_inf'] = float(popt[0])
            except Exception:
                pass
    return out


def fmt_fes(o):
    a = o['args']
    lab = f"{a['model']} q={a.get('q')} g={a.get('g')} lam={a.get('lam')}" + (f" kappa={a.get('kappa')}" if a['model'] == 'effective' else '')
    s = [f"== {o['path']}: {lab}, cell={a['ncells']}"]
    for r, x in enumerate(o['rows']):
        cp = '' if r == 0 else f" c_pair={o['c_pair'][r-1]:.3f}"
        extra = ''.join(f" {k}={x[k]:.2e}" if k in ('dS', 'd_eb', 'err') else f" {k}={x[k]:.5f}"
                        for k in ('Q', 'leak', 'dS', 'd_eb', 'err') if k in x)
        s.append(f"  chi={x['chi']:4d} S={x['S']:.5f} xi={x['xi']:9.3f} e={x['e']:.10f} ph={x['phase']:+.4f}{cp}{extra}")
    if 'c_fit' in o:
        s.append(f"  c_fit(last 3)={o['c_fit']:.3f} kappa={o['kappa']:.3f}" + (f" p={o['p']:.2f}" if 'p' in o else ''))
    return '\n'.join(s)


def obc_fit(path):
    d = json.load(open(path))
    a = d['args']
    sub = d['sub']
    L = a['L']
    st0 = [s for s in d['states'] if s['state'] == 0][0]['recs']
    r = st0[-1]
    S = np.array(r['S'])
    Sb = S[sub - 1::sub] if sub > 1 else S          # bond after cell l, l = 1..L-1
    l = np.arange(1, len(Sb) + 1)
    w = (l >= L / 4) & (l <= 3 * L / 4)
    X = np.log(2 * L / np.pi * np.sin(np.pi * l / L))
    c_plain = 6 * np.polyfit(X[w], Sb[w], 1)[0]
    A = np.vstack([X[w] / 6, np.ones(w.sum()), (-1.0) ** l[w] / np.exp(X[w])]).T
    c_stag = np.linalg.lstsq(A, Sb[w], rcond=None)[0][0]
    ce = [6 * np.polyfit(X[w & (l % 2 == p)], Sb[w & (l % 2 == p)], 1)[0] for p in (0, 1)]
    out = dict(path=os.path.basename(path), args=a, L=L, chi=r['chi'], E0=r['E'], c_plain=float(c_plain),
               c_stag=float(c_stag), c_even=float(ce[0]), c_odd=float(ce[1]), trunc=r['trunc'], Smid=float(Sb[len(Sb) // 2]))
    if 'e_bond' in r:
        eb = np.array(r['e_bond'])
        mid = eb[len(eb) // 4: 3 * len(eb) // 4]
        out['dimer'] = float(np.abs(np.diff(mid)).max())
    ex = [s for s in d['states'] if s['state'] == 1]
    if ex:
        r1 = ex[0]['recs'][-1]
        out['gap'] = r1['E'] - r['E']
        out['Lgap'] = L * out['gap']
        out['ov'] = r1['ov']
    return out


if __name__ == '__main__':
    mode = sys.argv[1]
    files = sorted(sum([glob.glob(p) for p in sys.argv[2:]], []))
    if mode == 'fes':
        for f in files:
            try:
                print(fmt_fes(fes_table(f)))
            except Exception as exc:
                print('!!', f, exc)
    else:
        for f in files:
            try:
                o = obc_fit(f)
                a = o['args']
                print(f"{o['path']}: L={o['L']} chi={o['chi']} E0={o['E0']:.10f} trunc={o['trunc']:.1e} "
                      f"c_plain={o['c_plain']:.3f} c_stag={o['c_stag']:.3f} c_even={o['c_even']:.3f} c_odd={o['c_odd']:.3f} Smid={o['Smid']:.4f}"
                      + (f" gap={o['gap']:.6f} L*gap={o['Lgap']:.4f}" if 'gap' in o else '')
                      + (f" dimer={o['dimer']:.2e}" if 'dimer' in o else ''))
            except Exception as exc:
                print('!!', f, exc)


def md_rows(files, k=3):
    """Markdown rows: g, lam, lam/(alpha^2 Delta), chi range, xi range, c_fit, spread, p, leak, Q."""
    rows = []
    for f in files:
        o = fes_table(f, k)
        a = o['args']
        R = o['rows']
        if len(R) < 3:
            continue
        cp = o['c_pair'][-k + 1:] if len(o['c_pair']) >= k - 1 else o['c_pair']
        spread = max(abs(x - o['c_fit']) for x in cp)
        al = json.load(open(f))['info'].get('alpha', float('nan'))
        ratio = a.get('lam', 0) / al ** 2 if al == al else float('nan')
        rows.append(f"| {a.get('g')} | {a.get('lam')} | {ratio:.2f} | {R[-k]['chi']}–{R[-1]['chi']} | "
                    f"{R[-k]['xi']:.0f}–{R[-1]['xi']:.0f} | {o['c_fit']:.3f} ± {spread:.3f} | "
                    f"{o.get('p', float('nan')):.2f} | {R[-1].get('leak', float('nan')):.2f} | {R[-1].get('Q', 0):+.3f} |")
    return rows
