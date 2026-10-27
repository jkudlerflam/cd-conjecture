#!/usr/bin/env python3
"""Open-chain entropy profiles S(l): fits with oscillation corrections, and window-averaged
doubling estimates.  Uses data/obc_NAME_L*_s*.json (fdmrg.py, largest chi) and, for L = 48 and 96,
the runs of contrib/d2_classification (results_cluster/dmrg_NAME_L*_s*.json, dmrg.py).

For each profile, with x(l) = (1/6) ln[(2L/pi) sin(pi l/L)] and the window l in [L/4, 3L/4]:
  c_plain   S = c x + b
  c_Eb      S = c x + b + beta (E_b(l) - mean)            (Laflorencie et al. 2006)
  c_osc     S = c x + b + A cos(2 pi k l + phi) / [(2L/pi) sin(pi l/L)]^p, scanned over k in
            (0, 1/2] and p in {0.5, 1, 1.5, 2}, best least-squares fit
  c_smooth  S averaged over a sliding window of one oscillation period 1/k*, then S = c x + b
Across L: c_dbl(L -> 2L) = 6 [Sbar(2L) - Sbar(L)] / ln 2 with Sbar = mean of S over [L/4, 3L/4].
usage: analyze_obc.py NAME [NAME ...]
"""
import sys, os, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(HERE, '..', 'd2_classification', 'results_cluster')


def profiles(name, conv):
    """conv = 'full': full field on the end sites (d2_classification/dmrg.py, fdmrg.py --edgefield);
    conv = 'half': half field on the end sites (sum of bond terms, fdmrg.py default)."""
    out = {}

    def take(L, rec):
        if L not in out or rec['E'] < out[L]['E'] - 1e-9:
            out[L] = rec
    if conv == 'full':
        old_files = glob.glob(os.path.join(OLD, f'dmrg_{name}_L*_s*.json')) + \
            glob.glob(os.path.join(OLD, '..', 'results_dmrg_local', f'dmrg_{name}_L*_s*.json'))
        for f in old_files:
            d = json.load(open(f))
            take(d['L'], dict(S=np.array(d['S']), Eb=np.array(d['Eb']), E=d['E'], chi=d['chis'][-1],
                              src=os.path.basename(f)))
    pre = 'obcedge' if conv == 'full' else 'obc'
    for f in glob.glob(os.path.join(HERE, 'data', f'{pre}_{name}_L*_s*.json')):
        d = json.load(open(f))
        if not d['results']:
            continue
        r = d['results'][-1]
        L = d['L']
        take(L, dict(S=np.array(r['S']), Eb=np.array(r.get('Eb', [np.nan] * (L - 1))), E=r['E'], chi=r['chi'],
                     src=os.path.basename(f), trunc=r['trunc']))
    return dict(sorted(out.items()))


def fits(S, Eb, L):
    l = np.arange(1, L)
    sel = (l >= L / 4) & (l <= 3 * L / 4)
    g = (2 * L / np.pi) * np.sin(np.pi * l / L)
    x = np.log(g) / 6
    res = {}
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    res['c_plain'] = np.linalg.lstsq(A, S[sel], rcond=None)[0][0]
    if np.all(np.isfinite(Eb)):
        e = Eb[sel] - Eb[sel].mean()
        A2 = np.vstack([x[sel], np.ones(sel.sum()), e]).T
        res['c_Eb'] = np.linalg.lstsq(A2, S[sel], rcond=None)[0][0]
    best = (np.inf, None)
    for p in (0.5, 1.0, 1.5, 2.0):
        for k in np.linspace(0.01, 0.5, 491):
            A3 = np.vstack([x[sel], np.ones(sel.sum()), np.cos(2 * np.pi * k * l[sel]) / g[sel] ** p,
                            np.sin(2 * np.pi * k * l[sel]) / g[sel] ** p]).T
            coef, r, *_ = np.linalg.lstsq(A3, S[sel], rcond=None)
            rr = float(r[0]) if len(r) else np.inf
            if rr < best[0]:
                best = (rr, (coef[0], k, p))
    res['c_osc'], res['k_osc'], res['p_osc'] = best[1]
    per = min(max(1, int(round(1 / best[1][1]))), L // 8)
    Ssm = np.convolve(S, np.ones(per) / per, mode='same')
    ok = sel & (l > per) & (l < L - per)
    A4 = np.vstack([x[ok], np.ones(ok.sum())]).T
    res['c_smooth'] = np.linalg.lstsq(A4, Ssm[ok], rcond=None)[0][0]
    res['Sbar'] = float(S[sel].mean())
    res['Smid'] = float(S[L // 2 - 1])
    return res


if __name__ == '__main__':
    for name in sys.argv[1:]:
        for conv in ('full', 'half'):
            P = profiles(name, conv)
            if not P:
                continue
            print(f"=== {name}, {conv} field on the end sites")
            prev = None
            for L, rec in P.items():
                r = fits(rec['S'], rec['Eb'], L)
                dbl = ''
                if prev is not None and L == 2 * prev[0]:
                    dbl = (f"c_dbl={6 * (r['Sbar'] - prev[1]['Sbar']) / np.log(2):.3f} "
                           f"c_mid_cut={6 * (r['Smid'] - prev[1]['Smid']) / np.log(2):.3f}")
                print(f"L={L:4d} chi={rec['chi']:4d} E/L={rec['E']/L:.10f} c_plain={r['c_plain']:.3f} "
                      f"c_Eb={r.get('c_Eb', np.nan):.3f} c_osc={r['c_osc']:.3f} (k={r['k_osc']:.3f}, p={r['p_osc']}) "
                      f"c_smooth={r['c_smooth']:.3f} {dbl}  [{rec['src']}]")
                prev = (L, r)
