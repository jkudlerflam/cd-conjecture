"""Central-charge estimators for the projected n = 2 chain, with calibration.

(a) FES: one-site VUMPS, S = (c/6) ln xi + b over chi (xi from the leading
    transfer-matrix gap); FES3 uses instead the mean of -ln over the three
    leading distinct transfer-matrix scales (one per species).
(b) Profile: open-chain S(l) = (c/6) ln x + s0 + sum_k [A_k cos(k l) + B_k sin(k l)] / x,
    x = (2L/pi) sin(pi l / L), with Friedel wavevectors k from the Fermi points
    of the quadratic payload part (all pairwise differences) and 2 pi N_l / L.
Calibration: c_cal = 3 c(model)/c(M0) (c = 2 for C1 against C0) at identical
chi (FES) or identical L and chi (profile).  python3 -B fit2.py -> c_results.json
"""
import glob
import json
import sys
import numpy as np

ERRMAX = 1e-6
FOLDERS = ('results', 'results_typhon')


def wavevectors(fermi, nl_frac=None):
    ks = [k for _, k, _ in fermi]
    out = []
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            q = abs(np.angle(np.exp(1j * (ks[i] - ks[j]))))
            if q > 1e-6 and all(abs(q - x) > 1e-3 for x in out):
                out.append(q)
    if nl_frac is not None and 0 < nl_frac < 1:
        q = abs(np.angle(np.exp(2j * np.pi * nl_frac)))
        if all(abs(q - x) > 1e-3 for x in out):
            out.append(q)
    return out


def fit_profile(S, ks, l0=4):
    S = np.asarray(S)
    L = len(S) + 1
    l = np.arange(1, L)
    sel = (l >= l0) & (l <= L - l0)
    x = (2 * L / np.pi) * np.sin(np.pi * l / L)
    cols = [np.log(x) / 6, np.ones_like(x)]
    for k in ks:
        cols += [np.cos(k * l) / x, np.sin(k * l) / x]
        if abs(k - np.pi) < 1e-9:
            cols.pop()
    A = np.array(cols).T[sel]
    coef, *_ = np.linalg.lstsq(A, S[sel], rcond=None)
    return float(coef[0]), float(np.sqrt(np.mean((S[sel] - A @ coef) ** 2)))


def slope_c(lx, S):
    A = np.vstack([lx, np.ones_like(lx)]).T
    (m, b), *_ = np.linalg.lstsq(A, np.asarray(S), rcond=None)
    return float(6 * m)


def fes_points(m):
    pts = {}
    for folder in FOLDERS:
        for fn in glob.glob('%s/%s_vumps_chi*.json' % (folder, m)):
            for p in json.load(open(fn))['runs']:
                if p['err'] < ERRMAX and (p['chi'] not in pts or p['e_site'] < pts[p['chi']]['e_site']):
                    pts[p['chi']] = p
    return [pts[k] for k in sorted(pts)]


def tm_scales(p):
    ev = [e for e in p['tm_eigs'][1:]]
    distinct = []
    for e in ev:
        if all(abs(e - x) > 1e-6 * x for x in distinct):
            distinct.append(e)
    return distinct[:3]


if __name__ == '__main__':
    out = {}
    models = ('M0', 'M1', 'M2', 'C0', 'C1')
    allP = {m: fes_points(m) for m in models}
    groups = {'three': ('M0', 'M1', 'M2'), 'two': ('C0', 'C1')}
    for g, ms in groups.items():
        sets = [set(p['chi'] for p in allP[m]) for m in ms]
        common = sorted(set.intersection(*sets)) if all(sets) else []
        for m in ms:
            P = [p for p in allP[m] if p['chi'] in common]
            rec = dict(points=[dict(chi=p['chi'], xi=p['xi'], S=p['S'], e_site=p['e_site'], err=p['err'],
                                    nl=p['nl'], nz=p['nz'], tm=p['tm_eigs'][:6], iters=p['iters']) for p in allP[m]],
                       common_chi=common)
            for lo in (16, 32, 48):
                Q = [p for p in P if p['chi'] >= lo]
                if len(Q) >= 3:
                    rec['c_fes_chi>=%d' % lo] = slope_c(np.log([p['xi'] for p in Q]), [p['S'] for p in Q])
                    sc = [tm_scales(p) for p in Q]
                    if all(len(s) == 3 for s in sc):
                        lx = [-np.mean([np.log(-np.log(e)) for e in s]) for s in sc]
                        rec['c_fes3_chi>=%d' % lo] = slope_c(np.array(lx), [p['S'] for p in Q])
            out[m] = rec
    # spread of the calibrated FES value over all subsets of at least three common chi values
    import itertools
    for ref, ms, target in (('M0', ('M1', 'M2'), 3.0), ('C0', ('C1',), 2.0)):
        P = {m: {p['chi']: p for p in allP[m]} for m in (ref,) + ms}
        for m in ms:
            common = sorted(set(P[ref]) & set(P[m]))
            vals = []
            for r in range(3, len(common) + 1):
                for sub in itertools.combinations(common, r):
                    cm = slope_c(np.log([P[m][c]['xi'] for c in sub]), [P[m][c]['S'] for c in sub])
                    cr = slope_c(np.log([P[ref][c]['xi'] for c in sub]), [P[ref][c]['S'] for c in sub])
                    vals.append(target * cm / cr)
            if vals:
                full = target * slope_c(np.log([P[m][c]['xi'] for c in common]), [P[m][c]['S'] for c in common]) / \
                    slope_c(np.log([P[ref][c]['xi'] for c in common]), [P[ref][c]['S'] for c in common])
                raw = slope_c(np.log([P[m][c]['xi'] for c in common]), [P[m][c]['S'] for c in common])
                rawref = slope_c(np.log([P[ref][c]['xi'] for c in common]), [P[ref][c]['S'] for c in common])
                out[m]['c_fes_calibrated_subsets'] = dict(chi=common, all_points=full, raw=raw, raw_ref=rawref, n=len(vals),
                                                          mean=float(np.mean(vals)), min=float(np.min(vals)), max=float(np.max(vals)))
    for ref, ms, target in (('M0', ('M1', 'M2'), 3.0), ('C0', ('C1',), 2.0)):
        for key in list(out.get(ref, {})):
            if key.startswith('c_fes'):
                for m in ms:
                    if key in out[m]:
                        out[m][key + '_calibrated'] = target * out[m][key] / out[ref][key]
    # profiles
    for m in models:
        for folder in FOLDERS:
            for fn in sorted(glob.glob('%s/%s_L*_chi*.json' % (folder, m))):
                d = json.load(open(fn))
                L = len(d['S']) + 1
                ks = wavevectors(d['info']['fermi'], d['Nl'] / L)
                fits = {l0: fit_profile(d['S'], ks, l0)[0] for l0 in (3, 4, 6)}
                rec = dict(file=fn, L=L, Nl=d['Nl'], ks=ks, fits=fits, max_trunc=max(h['max_trunc'] for h in d['hist'][-2:]),
                           dE_last_sweep=d['hist'][-1]['E'] - d['hist'][-2]['E'], E=d['E'],
                           nz_mid=float(np.mean(d['meas']['nz'][L // 4: 3 * L // 4])),
                           nl_mid=float(np.mean(d['meas']['nl'][L // 4: 3 * L // 4])))
                if 'S_exact' in d:
                    rec['max_dev_from_exact'] = d['max_dev_from_exact']
                    rec['dE_from_exact'] = d['E'] - d['E_exact']
                    rec['exact_fits'] = {l0: fit_profile(d['S_exact'], ks, l0)[0] for l0 in (3, 4, 6)}
                chi = int(fn.split('chi')[1].split('.')[0])
                out[m]['profile_L%d_chi%d' % (L, chi)] = rec
    for m, ref, target in (('M1', 'M0', 3.0), ('M2', 'M0', 3.0), ('C1', 'C0', 2.0)):
        for k, rec in out[m].items():
            if k.startswith('profile') and k in out[ref]:
                rec['calibrated'] = {l0: target * rec['fits'][l0] / out[ref][k]['fits'][l0] for l0 in rec['fits']}
                # difference form: c = c_ref [1 + (f_model - f_ref^DMRG) / f_ref^exact], i.e. the
                # truncation deficit of the reference is removed additively before the bias correction
                rf = out[ref][k]
                rec['calibrated_diff'] = {l0: target * (1 + (rec['fits'][l0] - rf['fits'][l0]) / rf['exact_fits'][l0]) for l0 in rec['fits']}
    json.dump(out, open('c_results.json', 'w'), indent=1, default=float)
    for m, v in out.items():
        print(m, 'chi', v['common_chi'], {k: (round(x, 4) if isinstance(x, float) else x) for k, x in v.items() if k.startswith('c_fes')})
        for k, x in v.items():
            if k.startswith('profile'):
                print('   ', k, 'nz %.3f nl %.3f' % (x['nz_mid'], x['nl_mid']), 'fits', {a: round(b, 3) for a, b in x['fits'].items()},
                      'exact', {a: round(b, 3) for a, b in x.get('exact_fits', {}).items()},
                      'dev %.1e' % x.get('max_dev_from_exact', float('nan')), 'trunc %.1e' % x['max_trunc'],
                      'cal', {a: round(b, 3) for a, b in x.get('calibrated', {}).items()},
                      'diff', {a: round(b, 3) for a, b in x.get('calibrated_diff', {}).items()})
