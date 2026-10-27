#!/usr/bin/env python3
"""Collect all estimators into results/summary.json and print a Markdown table.

ED (periodic):   c_S = entanglement slope with blocks l >= L/4, mean and spread over the
                 three largest even L; drift = c_S(Lmax_even) - c_S(Lmax_even - 4).
                 c_half = slope of S(L/2) vs (1/3) ln L over the last four even L.
                 z10, z20 = level-spacing exponents (analyze.level_spacing), all L >= 12.
DMRG (open):     c_mid = 6 [S_win(96) - S_win(48)] / ln 2, with S_win the mean of S(l) over
                 l in [L/4, 3L/4] (removes the non-universal constant and averages
                 incommensurate oscillations); c_fit = OBC fit with the bond-energy term over l in [L/4, 3L/4].
"""
import json, os, glob, math
import numpy as np
import analyze
from dmrg import fit_obc

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')


def dmrg_rec(name, L):
    best = None
    for f in glob.glob(os.path.join(R, f'dmrg_{name}_L{L}_s*.json')) + glob.glob(os.path.join(HERE, 'results_dmrg_local', f'dmrg_{name}_L{L}_s*.json')):
        r = json.load(open(f))
        if best is None or r['E'] < best['E'] - 1e-9:
            best = r
    return best


def smid(r):
    """Mean of S(l) over l in [L/4, 3L/4].  For a CFT the L-dependence of this window
    average is (c/6) ln L exactly (the rest depends only on l/L), and averaging over L/2
    bonds suppresses incommensurate oscillations."""
    S = np.array(r['S']); L = r['L']
    return float(S[L // 4 - 1:3 * L // 4].mean())


def row(name):
    a = analyze.analyse(name)
    rows = [r for r in a['rows'] if r['L'] % 2 == 0]
    cs = [r['cSq'] for r in rows[-3:]]
    drift = rows[-1]['cSq'] - rows[-3]['cSq'] if len(rows) >= 3 else np.nan
    ls = analyze.level_spacing(name, ns=(10, 20))
    out = dict(name=name, Lmax=a['Lmax'], cS=float(np.mean(cs)), cS_spread=float(np.ptp(cs)),
               cS_last=rows[-1]['cSq'], cS_drift=float(drift), c_half=float(a['c_half']),
               z10=ls[10]['z'], z10_err=ls[10]['z_err'], z20=ls[20]['z'], z20_err=ls[20]['z_err'],
               v2_last=rows[-1]['v2'])
    d48, d96 = dmrg_rec(name, 48), dmrg_rec(name, 96)
    if d48 is not None and d96 is not None:
        out['c_mid'] = 6 * (smid(d96) - smid(d48)) / math.log(2)
        out['c_fit96'] = fit_obc(d96['S'], 96, 0.25, d96['Eb'])
        out['c_fit48'] = fit_obc(d48['S'], 48, 0.25, d48['Eb'])
        out['trunc96'] = max(h['maxtrunc'] for h in d96['hist'][-2:])
        out['Swin48'] = smid(d48); out['Swin96'] = smid(d96); out['chi96'] = d96['chis'][-1]
    elif d48 is not None:
        out['c_fit48'] = fit_obc(d48['S'], 48, 0.25, d48['Eb'])
    return out


if __name__ == '__main__':
    import sys
    names = sys.argv[1:] or (list(json.load(open(os.path.join(HERE, 'candidates_params.json')))) +
                             ['ising', 'xx', 'heis', 'xxz05', 'xxzm08', 'xxzm095', 'xxz_dm', 'xxzm07_dm_h', 'cluster_ising'])
    allrows = []
    for n in names:
        try:
            allrows.append(row(n))
        except Exception as e:
            print('skip', n, e)
    json.dump(allrows, open(os.path.join(R, 'summary.json'), 'w'), indent=1)
    f = lambda x, p=2: ('%.*f' % (p, x)) if (x is not None and np.isfinite(x)) else '-'
    print('| point | ED Lmax | c_S (l>=L/4) | drift | c_half | z10 | z20 | DMRG c_mid(48->96) | c_fit(96) |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for r in allrows:
        print(f"| {r['name']} | {r['Lmax']} | {f(r['cS'])} ± {f(r['cS_spread'])} | {f(r['cS_drift'])} | {f(r['c_half'])} | "
              f"{f(r['z10'],1)} ± {f(r['z10_err'],1)} | {f(r['z20'],1)} ± {f(r['z20_err'],1)} | {f(r.get('c_mid'))} | {f(r.get('c_fit96'))} |")
