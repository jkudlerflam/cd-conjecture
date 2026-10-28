"""Central-charge estimators for the projected n=1 chain, with calibration.

(a) finite-entanglement scaling (one-site VUMPS): S = (c/6) ln xi + b over chi.
    Two variants: xi from the leading transfer-matrix gap, and, for the
    two-species models, the geometric mean of the two leading scales.
(b) open-chain interval entanglement (finite DMRG):
    S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + s0
           + sum_a [A_a cos(2 k_a l) + B_a sin(2 k_a l)] / [(2L/pi) sin(pi l/L)],
    k_a = pi n_a from the measured bulk densities.
Calibration: the same estimators applied to the decoupled free reference M0
(c = 2 exactly, profile also computed exactly) and to the payload-only C0 (c = 1).
Usage: python3 -B fitc.py   -> c_results.json
"""
import glob
import json
import numpy as np

FOLDERS = ('results', 'results_typhon')


def fit_profile(S, dens, l0=8, osc=True):
    S = np.asarray(S)
    L = len(S) + 1
    l = np.arange(1, L)
    sel = (l >= l0) & (l <= L - l0)
    x = (2 * L / np.pi) * np.sin(np.pi * l / L)
    cols = [np.log(x) / 6, np.ones_like(x)]
    if osc:
        for n in dens:
            k = np.pi * n
            cols += [np.cos(2 * k * l) / x, np.sin(2 * k * l) / x]
    A = np.array(cols).T[sel]
    coef, *_ = np.linalg.lstsq(A, S[sel], rcond=None)
    resid = S[sel] - A @ coef
    return float(coef[0]), float(np.sqrt(np.mean(resid ** 2)))


def slope_c(lx, S):
    A = np.vstack([lx, np.ones_like(lx)]).T
    (m, b), *_ = np.linalg.lstsq(A, np.asarray(S), rcond=None)
    return float(6 * m)


def free_profile(L, hop, onsite):
    H = np.zeros((L, L), complex)
    for i in range(L - 1):
        H[i + 1, i] = hop
        H[i, i + 1] = np.conj(hop)
    H += onsite * np.eye(L)
    e, v = np.linalg.eigh(H)
    occ = v[:, e < 0]
    C = occ @ occ.conj().T
    prof = []
    for l in range(1, L):
        ev = np.linalg.eigvalsh(C[:l, :l])
        ev = ev[(ev > 1e-15) & (ev < 1 - 1e-15)]
        prof.append(float(-(ev * np.log(ev) + (1 - ev) * np.log(1 - ev)).sum()))
    return np.array(prof), float((e < 0).sum() / L)


def fes_points(m):
    pts = {}
    for folder in FOLDERS:
        for fn in glob.glob('%s/%s_vumps*.json' % (folder, m)):
            for p in json.load(open(fn))['points']:
                if p['err'] < 1e-6 and (p['chi'] not in pts or p['e_site'] < pts[p['chi']]['e_site']):
                    pts[p['chi']] = p
    return [pts[k] for k in sorted(pts)]


if __name__ == '__main__':
    out = {}
    two = {'M0', 'M1', 'M2'}
    allP = {m: fes_points(m) for m in ('M0', 'M1', 'M2', 'C0', 'C1')}
    # identical chi sets inside each calibration group
    groups = {'two': ('M0', 'M1', 'M2'), 'one': ('C0', 'C1')}
    common = {g: sorted(set.intersection(*[set(p['chi'] for p in allP[m]) for m in ms])) for g, ms in groups.items()}
    for m in ('M0', 'M1', 'M2', 'C0', 'C1'):
        g = 'two' if m in two else 'one'
        P = [p for p in allP[m] if p['chi'] in common[g]]
        rec = dict(points=[dict(chi=p['chi'], xi=p['xi'], S=p['S'], e_site=p['e_site'], err=p['err'],
                                nl=p['nl'], nz=p['nz'], tm=p['tm_eigs'][:4]) for p in allP[m]],
                   common_chi=common[g])
        for lo in (16, 32):
            Q = [p for p in P if p['chi'] >= lo]
            if len(Q) >= 3:
                rec['c_fes_chi>=%d' % lo] = slope_c(np.log([p['xi'] for p in Q]), [p['S'] for p in Q])
                if m in two:
                    l2 = [-np.log(p['tm_eigs'][1]) for p in Q]
                    l3 = [-np.log(p['tm_eigs'][2]) for p in Q]
                    lxe = [-(np.log(a) + np.log(b)) / 2 for a, b in zip(l2, l3)]
                    rec['c_fes2_chi>=%d' % lo] = slope_c(np.array(lxe), [p['S'] for p in Q])
        out[m] = rec
    # calibration of FES on identical chi sets
    for key in list(out['M0']):
        if key.startswith('c_fes'):
            for m in ('M1', 'M2'):
                if key in out[m]:
                    out[m][key + '_calibrated'] = 2 * out[m][key] / out['M0'][key]
    for key in list(out['C0']):
        if key.startswith('c_fes') and key in out['C1']:
            out['C1'][key + '_calibrated'] = out['C1'][key] / out['C0'][key]
    # profiles
    for m in ('M0', 'M1', 'M2', 'C0', 'C1'):
        for folder in FOLDERS:
            for fn in sorted(glob.glob('%s/%s_finite_L*_chi*.json' % (folder, m))):
                d = json.load(open(fn))
                L = len(d['S']) + 1
                nl = float(np.mean(d['meas']['nl'][L // 4: 3 * L // 4]))
                nz = float(np.mean(d['meas']['nz'][L // 4: 3 * L // 4]))
                dens = [x for x in (nl, nz) if 0.02 < x < 0.98]
                fits = {l0: fit_profile(d['S'], dens, l0, True)[0] for l0 in (4, 8, 12)}
                rec = dict(file=fn, nl=nl, nz=nz, fits=fits, max_trunc=max(h['max_trunc'] for h in d['hist'][-2:]),
                           dE_last_sweep=d['hist'][-1]['E'] - d['hist'][-2]['E'])
                if m == 'M0':
                    info = d['info']
                    Sp, _ = free_profile(L, -1.0, 0.0)
                    Sl, _ = free_profile(L, -info['jl'], -info['kappa'])
                    Sex = Sp + Sl
                    rec['max_dev_from_exact'] = float(np.abs(np.array(d['S']) - Sex).max())
                    rec['exact_fits'] = {l0: fit_profile(Sex, dens, l0, True)[0] for l0 in (4, 8, 12)}
                if m == 'C0':
                    Sp, _ = free_profile(L, -1.0, 0.0)
                    rec['max_dev_from_exact'] = float(np.abs(np.array(d['S']) - Sp).max())
                    rec['exact_fits'] = {l0: fit_profile(Sp, dens, l0, True)[0] for l0 in (4, 8, 12)}
                chi = int(fn.split('chi')[1].split('.')[0])
                out[m]['profile_L%d_chi%d' % (L, chi)] = rec
    # equal-numerics calibration of the profile estimator (same L and chi as M0, or C0)
    for m, ref in (('M1', 'M0'), ('M2', 'M0'), ('C1', 'C0')):
        target = 2.0 if ref == 'M0' else 1.0
        for k, rec in out[m].items():
            if k.startswith('profile') and k in out[ref]:
                rec['calibrated'] = {l0: target * rec['fits'][l0] / out[ref][k]['fits'][l0] for l0 in rec['fits']}
    json.dump(out, open('c_results.json', 'w'), indent=1, default=float)
    for m, v in out.items():
        print(m, 'chi', v['common_chi'], {k: round(x, 4) for k, x in v.items() if k.startswith('c_fes')})
        for k, x in v.items():
            if k.startswith('profile'):
                print('   ', k, 'nl %.4f nz %.4f' % (x['nl'], x['nz']), {a: round(b, 4) for a, b in x['fits'].items()},
                      'exact' if 'exact_fits' in x else '', {a: round(b, 4) for a, b in x.get('exact_fits', {}).items()},
                      'dev %.1e' % x.get('max_dev_from_exact', float('nan')), 'trunc %.1e' % x['max_trunc'],
                      'cal', {a: round(b, 4) for a, b in x.get('calibrated', {}).items()})
