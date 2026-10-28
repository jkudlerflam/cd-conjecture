#!/usr/bin/env python3
"""Tables for the two-flavour README: FES and block-entropy central charges, fillings, and the M = 3 ring
spectra against the effective chain.  Writes data/summary2.json.  Claude (Opus 5.5), 2026-10-28.
"""
import os, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}


def fes(name, maxerr=1e-4, chimin=16, bond=0):
    f = os.path.join(HERE, 'data', 'fes', name + '.json')
    if not os.path.exists(f):
        return None
    r = [x for x in json.load(open(f))['results'] if x['err'] < maxerr and x['chi'] >= chimin]
    if len(r) < 2:
        return None
    S = np.array([x['S'][bond] for x in r]); xi = np.array([x['xi_cells'] for x in r])
    x = np.log(xi); A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, S, rcond=None)
    n = len(x)
    se = float(np.sqrt(res[0] / (n - 2) / np.sum((x - x.mean()) ** 2))) if n > 2 and len(res) else float('nan')
    loc = r[-1]['local']
    fill = [max(v) for v in loc] if isinstance(loc[0], list) else loc
    return dict(chis=[x['chi'] for x in r], xi=xi.tolist(), S=S.tolist(), c=float(6 * coef[0]), se=6 * se,
                pairwise=(6 * np.diff(S) / np.diff(x)).tolist(), err=[x['err'] for x in r], fill=fill)


def blockc(fname):
    """Block estimator used in the README: c_b = 3 [S(2) - S(1)] / ln 2 for blocks of whole cells."""
    d = json.load(open(fname))
    S = dict((int(l), s) for l, s in d['S'])
    if 1 not in S or 2 not in S:
        return None
    return dict(xi=d['xi_cells'], chi=d['chi'], c=float(3 * (S[2] - S[1]) / np.log(2)))


if __name__ == '__main__':
    print('== FES ==')
    out['fes'] = {}
    for nm in ['xx1s', 'xx2', 'effg_l0.15', 'effg_l0.20', 'micg_l0.15', 'micg_l0.20']:
        t = fes(nm)
        if t:
            out['fes'][nm] = t
            print(f"{nm:12s} chi {t['chis']} xi {np.round(t['xi'], 1).tolist()} c = {t['c']:.3f} +- {t['se']:.3f} "
                  f"pairwise {np.round(t['pairwise'], 2).tolist()} fill {np.round(t['fill'], 4).tolist()}")
    # the microscopic fits quoted in the README: chi >= 32 and chi >= 48, all converged points (err < 1e-2)
    for nm in ['micg_l0.15', 'micg_l0.20']:
        for cm in (32, 48):
            t = fes(nm, maxerr=1e-2, chimin=cm)
            if t:
                out['fes'][f'{nm}_chi>={cm}'] = t
                print(f"{nm} chi >= {cm}: chi {t['chis']} xi {np.round(t['xi'], 2).tolist()} c = {t['c']:.3f} +- {t['se']:.3f}")
    print('== block entropies ==')
    out['block'] = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'block', '*.json'))):
        b = blockc(f)
        if b:
            out['block'][os.path.basename(f)] = b
            print(os.path.basename(f), b)
    print('== M = 3 rings: micro vs effective (exact two-cell bond) ==')
    out['rings'] = {}
    lines = open(os.path.join(HERE, 'tasks_ring2.txt')).read().strip().split('\n')
    import effective2 as E, ring2 as R
    Heffs = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'ring', 'ring2_*.json'))):
        i = int(os.path.basename(f)[6:-5]); a = lines[i].split()
        lam = float(a[a.index('--lam') + 1]); N = (int(a[a.index('--N') + 1]), int(a[a.index('--N') + 2]))
        if lam not in Heffs:
            Heffs[lam], _, _ = E.exact_eff(dict(E.DEFAULT), lam)
        lv = json.load(open(f))['levels']
        for key, w in lv.items():
            m = int(key.split('_')[2])
            we = R.eff_ring(3, N, m, Heffs[lam], len(w))
            out['rings'][f'{lam}_{N[0]}_{N[1]}_{m}'] = dict(micro=w, eff=we.tolist())
    for lam in sorted(Heffs):
        rows = {k: v for k, v in out['rings'].items() if k.startswith(f'{lam}_')}
        e0m = min(min(v['micro']) for v in rows.values()); e0e = min(min(v['eff']) for v in rows.values())
        dev = []
        print(f'lam = {lam}: E0 micro {e0m:.8f} eff {e0e:.8f}')
        for k, v in sorted(rows.items(), key=lambda kv: min(kv[1]['micro'])):
            gm = np.array(v['micro']) - e0m; ge = np.array(v['eff'][:len(v['micro'])]) - e0e
            print(f'   {k}: gaps micro {np.round(gm, 6).tolist()} eff {np.round(ge, 6).tolist()}')
            dev += [abs(x - y) for x, y in zip(gm, ge) if y > 1e-9]
        out['rings'][f'summary_{lam}'] = dict(E0_micro=e0m, E0_eff=e0e, max_gap_dev=max(dev) if dev else None)
    json.dump(out, open(os.path.join(HERE, 'data', 'summary2.json'), 'w'), indent=1, default=float)
