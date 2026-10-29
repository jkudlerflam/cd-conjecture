#!/usr/bin/env python3
"""Calibrated c from the open-chain DMRG runs (data/obc/*.json) and the FES runs (data/fes/*.json).

Profile estimator: per run, fit S(l) to (c/6) ln[(2(L+1)/pi) sin(pi l/(L+1))] + s0 + Friedel terms
(fitlib.profile_fit) at the largest converged chi; the identical fit is applied to the exact
reference (independent open free-fermion chains with the same L and particle numbers, one per
gapless mode), and c_cal = m c / c_ref with m the number of reference modes (3 for LV_2, 2 for the
Neel control).  Variants: lmin in (6, 8, 12), with and without sum/difference wavevectors; the
quoted error combines the variant spread with the change of c between the last two chi.
Centre estimator: c = 6 [S(L2/2) - S(L1/2)] / ln[(L2+1)/(L1+1)] for consecutive L at fixed fillings,
calibrated the same way.
FES: S = (c/6) ln xi over chi windows (xi from the leading transfer-matrix eigenvalue), and the
multi-scale variant S = (c/6) (1/m) sum_{j<=m} ln xi_j with xi_j from the m leading distinct
transfer-matrix magnitudes, each calibrated by the matched free reference at the same chi.
usage: analyze.py  (writes data/summary.json).  Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import models as MD
import fitlib as F

VARIANTS = [(lmin, cb) for lmin in (6, 8, 12) for cb in (False, True)]
MAIN = (8, True)


def load(pattern):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'obc', pattern))):
        d = json.load(open(f))
        if d['results']:
            out.append(d)
    return out


def fits(S, L, Ns):
    qs = {cb: F.wavevectors(L, Ns, cb) for cb in (False, True)}
    return {v: F.profile_fit(S, L, qs[v[1]], v[0])[0] for v in VARIANTS}


def dfit(dS, L, Ns):
    """slope of a profile difference: c_model - c_reference, over the fit variants."""
    qs = {cb: F.wavevectors(L, Ns, cb) for cb in (False, True)}
    v = np.array([F.profile_fit(dS, L, qs[cb], lmin)[0] for (lmin, cb) in VARIANTS])
    return float(np.median(v)), float(v.min()), float(v.max())


def run_table(runs, refs, nref):
    """Per run and chi: (i) ratio calibration of the raw profile c by the exact reference,
    (ii) difference estimator c = m + slope(S_model - S_exact ref),
    (iii) the same with the matched-velocity reference computed by DMRG at the same chi and L."""
    rows = []
    for d in runs:
        a = d['args']; L = a['L']; Ns = [a['M']] + a['N']
        Sref = MD.ref_S(L, Ns[:nref])
        cref = fits(Sref, L, Ns[:nref])
        rd = refs.get(L)
        rchi = {r['chi']: np.array(r['S']) for r in rd['results']} if rd else {}
        for r in d['results']:
            S = np.array(r['S'])
            cm = fits(S, L, Ns)
            ratio = nref * cm[MAIN] / cref[MAIN]
            dx = dfit(S - Sref, L, Ns)
            row = dict(name=d['name'], L=L, Ns=Ns, chi=r['chi'], trunc=r['trunc'], E=r['E'], sweeps=r['sweeps'],
                       c_raw=cm[MAIN], c_ref_exact=cref[MAIN], c_ratio=ratio,
                       c_diff_exact=nref + dx[0], c_diff_exact_range=[nref + dx[1], nref + dx[2]],
                       S_mid=float(S[L // 2]), S_mid_exact_ref=float(Sref[L // 2]))
            if r['chi'] in rchi:
                dc = dfit(S - rchi[r['chi']], L, Ns)
                row['c_diff_dmrgref'] = nref + dc[0]
                row['c_diff_dmrgref_range'] = [nref + dc[1], nref + dc[2]]
                row['S_mid_dmrg_ref'] = float(rchi[r['chi']][L // 2])
            rows.append(row)
    return rows


def centre_pairs(runs, nref):
    rows = []
    runs = sorted(runs, key=lambda d: d['args']['L'])
    for d1, d2 in zip(runs[:-1], runs[1:]):
        L1, L2 = d1['args']['L'], d2['args']['L']
        N1 = [d1['args']['M']] + d1['args']['N']; N2 = [d2['args']['M']] + d2['args']['N']
        lg = np.log((L2 + 1) / (L1 + 1))
        c = 6 * (d2['results'][-1]['S'][L2 // 2] - d1['results'][-1]['S'][L1 // 2]) / lg
        cr = 6 * (MD.ref_S(L2, N2[:nref])[L2 // 2] - MD.ref_S(L1, N1[:nref])[L1 // 2]) / lg
        rows.append(dict(L1=L1, L2=L2, chi1=d1['results'][-1]['chi'], chi2=d2['results'][-1]['chi'],
                         c_raw=float(c), c_ref=float(cr), c_cal=float(nref * c / cr)))
    return rows


def fes_fit(res, chimin, m=1):
    """c from S vs mean ln xi over the m leading distinct transfer-matrix scales."""
    X, Y = [], []
    for r in res:
        if r['chi'] < chimin or 'lamN' not in r:
            continue
        lam = np.abs(np.array([complex(*z) for z in r['lamN']]))
        lam = lam / lam[0]
        mags = []
        for v in lam[1:]:
            if v < 1 - 1e-12 and all(abs(v - u) > 1e-6 * u for u in mags):
                mags.append(v)
        if len(mags) < m:
            continue
        X.append(np.mean([np.log(-1 / np.log(v)) for v in mags[:m]]))
        Y.append(r['S'][0])
    if len(X) < 3:
        return None, None, len(X)
    A = np.column_stack([np.array(X) / 6, np.ones(len(X))])
    coef, res_, *_ = np.linalg.lstsq(A, np.array(Y), rcond=None)
    r = np.array(Y) - A @ coef
    se = np.sqrt(np.sum(r ** 2) / max(1, len(X) - 2) / np.sum((np.array(X) - np.mean(X)) ** 2)) * 6
    return float(coef[0]), float(se), len(X)


def fes_table():
    out = {}
    fs = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob(os.path.join(HERE, 'data', 'fes', '*.json'))}
    for name, d in fs.items():
        res = d['results']
        out[name] = dict(chis=[r['chi'] for r in res], S=[r['S'][0] for r in res], xi=[r['xi'] for r in res],
                         err=[r['err'] for r in res], local=[r['local'] for r in res])
        for chimin in (32, 64, 96, 128):
            for m in (1, 2, 3):
                c, se, n = fes_fit(res, chimin, m)
                out[name][f'c_chi>={chimin}_m{m}'] = (c, se, n)
    return out


def main():
    model = load('lv2_L*.json'); ctrl = load('ctrl*_L*.json')
    ref3 = {d['args']['L']: d for d in load('ref3_L*.json')}
    ref2 = {d['args']['L']: d for d in load('ref2_L*.json')}
    tm = run_table(model, ref3, 3); tc = run_table(ctrl, ref2, 2)
    summ = dict(model=tm, control=tc, centre_model=centre_pairs(model, 3), centre_control=centre_pairs(
        [d for d in ctrl if d['name'].startswith('ctrl_')], 2))
    for name, rows in (('LV_2 (target 3)', tm), ('control (target 2)', tc)):
        print(name)
        for r in rows:
            extra = ''
            if 'c_diff_dmrgref' in r:
                extra = f" | vs DMRG ref same chi {r['c_diff_dmrgref']:.3f} [{r['c_diff_dmrgref_range'][0]:.3f},{r['c_diff_dmrgref_range'][1]:.3f}] Sref_mid={r['S_mid_dmrg_ref']:.4f}"
            print(f"  {r['name']:11s} chi={r['chi']:5d} trunc={r['trunc']:.1e} ratio={r['c_ratio']:.3f} "
                  f"diff_exact={r['c_diff_exact']:.3f} [{r['c_diff_exact_range'][0]:.3f},{r['c_diff_exact_range'][1]:.3f}] "
                  f"Smid={r['S_mid']:.4f} (exact ref {r['S_mid_exact_ref']:.4f}){extra}")
    for name in ('centre_model', 'centre_control'):
        print(name)
        for r in summ[name]:
            print(f"  L={r['L1']}(chi {r['chi1']})->{r['L2']}(chi {r['chi2']}): c_raw={r['c_raw']:.4f} c_ref={r['c_ref']:.4f} c_cal={r['c_cal']:.4f}")
    # matched-velocity references by DMRG against their exact profiles
    for d in load('ref*_L*.json'):
        a = d['args']; L = a['L']; Ns = [a['M']] + a['N']; Sx = MD.ref_S(L, Ns)
        chk = []
        for r in d['results']:
            dc = dfit(np.array(r['S']) - Sx, L, Ns)
            chk.append(dict(chi=r['chi'], max_dS=float(np.max(np.abs(np.array(r['S']) - Sx))),
                            dS_mid=float(r['S'][L // 2] - Sx[L // 2]), dc=dc[0], trunc=r['trunc'], dE=None))
        summ['refcheck_' + d['name']] = chk
        print('DMRG reference vs exact', d['name'], [(c['chi'], round(c['dS_mid'], 4), round(c['dc'], 3)) for c in chk])
    # final table: value at the largest chi, variant range, and the linear-in-truncation extrapolation from the
    # last two chi (the same extrapolation applied to the DMRG references shows its residual bias)
    final = []
    for d in load('lv2_L*.json') + load('ctrl_L*.json') + load('ref*_L*.json'):
        a = d['args']; L = a['L']; Ns = [a['M']] + a['N']
        nref = 2 if (a.get('vz') or d['name'].startswith('ref2')) else 3
        Sx = MD.ref_S(L, Ns[:nref])
        pts = [(r['chi'], r['trunc'], dfit(np.array(r['S']) - Sx, L, Ns)) for r in d['results']]
        (c1, t1, v1), (c2, t2, v2) = pts[-2], pts[-1]
        ext = v2[0] - t2 * (v2[0] - v1[0]) / (t2 - t1)
        final.append(dict(name=d['name'], L=L, chi=c2, trunc=t2, c=nref + v2[0], lo=nref + v2[1], hi=nref + v2[2],
                          c_prev=nref + v1[0], chi_prev=c1, c_extrap=nref + ext))
        print(f"FINAL {d['name']:10s} L={L} chi={c2} trunc={t2:.1e} c={nref + v2[0]:.3f} [{nref + v2[1]:.3f},{nref + v2[2]:.3f}] "
              f"(chi {c1}: {nref + v1[0]:.3f}) extrap={nref + ext:.3f}")
    summ['final'] = final
    if os.path.isdir(os.path.join(HERE, 'data', 'fes')):
        ft = fes_table()
        summ['fes'] = ft
        for name, v in sorted(ft.items()):
            print('FES', name, 'chi_max', v['chis'][-1], {k: (round(x[0], 3), round(x[1], 3), x[2])
                                                          for k, x in v.items() if k.startswith('c_') and x[0] is not None})
    json.dump(summ, open(os.path.join(HERE, 'data', 'summary.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
