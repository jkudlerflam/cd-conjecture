#!/usr/bin/env python3
"""Calibrated c for the logical vacancy model LV_n (logical_vac.py): ring interval entanglement and FES.
Writes data/lv_summary.json.

Ring estimator.  For each ring, S(l) for l = 2..L-2 is fitted to a + (c/3) ln x + sum_k b_k cos(q_k l)/x,
x = (L/pi) sin(pi l/L), with q_k the Friedel wavevectors 2 pi M/L (charge) and 2 pi N_f/L (flavours).
The same fit is applied to the reference, a sum of free-fermion XX rings with the same L and particle
numbers (one per gapless mode: M on L sites for the charge, N_f on L sites per occupied flavour).
Calibrated c = m c_raw / c_ref (m = number of modes in the reference); the additive version
c_raw - c_ref + m is also reported.  A second ring estimator uses half-ring entropies at fixed density,
c = 3 [S(L2/2) - S(L1/2)] / ln(L2/L1) (mean of the two middle cuts when L is odd), calibrated the same way.

FES estimator.  c from S = (c/6) ln xi + s0 fitted over the last four chi, calibrated by decoupled XX chains
with the first-order velocity ratios (run_fes.py kxx), c_cal = m c_raw / c_ref.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, glob, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import logical_vac as LV  # noqa: E402
D = os.path.join(HERE, 'data')


def fit_c(L, S, qs):
    l = np.arange(2, L - 1)
    x = (L / np.pi) * np.sin(np.pi * l / L)
    cols = [np.ones_like(x), np.log(x) / 3]
    for q in qs:
        cols.append(np.cos(q * l) / x)
    A = np.array(cols).T
    y = np.array(S)[l - 1]
    if A.shape[0] < A.shape[1] + 3:
        return None
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(coef[1])


def half(S, L):
    S = np.asarray(S)
    return float(S[L // 2 - 1]) if L % 2 == 0 else float(0.5 * (S[L // 2 - 1] + S[L // 2]))


def qlist(L, M, N):
    qs = []
    for p in [M] + [x for x in N if x > 0]:
        q = (2 * np.pi * p / L) % (2 * np.pi)
        q = min(q, 2 * np.pi - q)
        if q > 1e-9 and all(abs(q - r) > 1e-9 for r in qs):
            qs.append(q)
    return qs


rows = {}
for f in sorted(glob.glob(os.path.join(D, 'lv', '*.json'))):
    d = json.load(open(f))
    a = d['args']
    name = os.path.basename(f)
    series = re.sub(r'_L\d+\.json$', '', name)
    L, M, N = a['L'], a['M'], a['N']
    Sref = np.array(LV.ff_ring_entropies(L, M))
    m = 1
    vz = a.get('vz') or [0.0] * len(N)
    gapless = [x for x, v, tau in zip(N, vz, a['taus']) if v <= 2 * tau]     # XXZ Delta = vz/(2 tau) > 1: gapped
    N = gapless
    for x in N:
        if 0 < x < L:
            Sref = Sref + np.array(LV.ff_ring_entropies(L, x)); m += 1
    qs = qlist(L, M, N)
    rec = dict(L=L, M=M, N=N, J=a['J'], dim=d['dim'], gap=d['gap'], modes=m,
               c_raw=fit_c(L, d['S'], qs), c_ref=fit_c(L, Sref, qs), S_half=half(d['S'], L), Sref_half=half(Sref, L))
    if rec['c_raw'] is not None:
        rec['c_cal'] = m * rec['c_raw'] / rec['c_ref']
        rec['c_cal_add'] = rec['c_raw'] - rec['c_ref'] + m
    rows.setdefault(series, []).append(rec)

out = dict(ring={}, fes={})
print('== ring interval entanglement (fit per L, calibrated by free-fermion rings with the same L and fillings)')
for s, rr in rows.items():
    rr.sort(key=lambda r: r['L'])
    pairs = []
    from fractions import Fraction as Fr
    groups = {}
    for r in rr:
        groups.setdefault((Fr(r['M'], r['L']),) + tuple(Fr(x, r['M']) for x in r['N']), []).append(r)
    for g in groups.values():
        for r1, r2 in zip(g[:-1], g[1:]):
            lr = np.log(r2['L'] / r1['L'])
            cp = 3 * (r2['S_half'] - r1['S_half']) / lr
            cr = 3 * (r2['Sref_half'] - r1['Sref_half']) / lr
            pairs.append(dict(L=(r1['L'], r2['L']), c_raw=cp, c_ref=cr, c_cal=r1['modes'] * cp / cr))
    out['ring'][s] = dict(per_L=rr, half_pairs=pairs)
    print(f'  {s} (modes {rr[0]["modes"]})')
    for r in rr:
        if r['c_raw'] is None:
            continue
        print(f"    L={r['L']:2d} M={r['M']} N={r['N']} dim={r['dim']} gap={r['gap']:.4f} c_raw={r['c_raw']:.4f} "
              f"c_ref={r['c_ref']:.4f} c_cal={r['c_cal']:.4f} (additive {r['c_cal_add']:.4f})")
    for p in pairs:
        print(f"    half-ring pair L={p['L']}: c_raw={p['c_raw']:.4f} c_ref={p['c_ref']:.4f} c_cal={p['c_cal']:.4f}")


def fes_fit(name):
    # the same run may exist twice (compute node, and a login-node copy with suffix L): take the longer one
    R = []
    for nm in (name, name + 'L'):
        p = os.path.join(D, 'fes', nm + '.json')
        if os.path.exists(p):
            r = json.load(open(p))['results']
            if len(r) > len(R):
                R = r
    if len(R) < 4:
        return None
    S = np.array([np.mean(r['S']) for r in R]); xi = np.array([r['xi'] for r in R])
    pw = [6 * (S[i + 1] - S[i]) / np.log(xi[i + 1] / xi[i]) for i in range(len(S) - 1)]
    fits = [float(6 * np.polyfit(np.log(xi[k:k + 4]), S[k:k + 4], 1)[0]) for k in range(max(0, len(S) - 6), len(S) - 3)]
    return dict(chis=[r['chi'] for r in R], xi=xi.tolist(), S=S.tolist(), pairwise=pw, fit_last4=fits[-1],
                fits_sliding=fits, err=[r['err'] for r in R], local_last=R[-1].get('local'))


print('\n== FES (fit of the last four chi; sliding windows in brackets), calibrated by decoupled XX chains')
for model, ref, m in [('lv_n1_J0.5', 'kxx_1_0.48', 2), ('lv_n2_J0.5', 'kxx_1_0.48_0.24', 3),
                      ('lv_ctrl1_J0.5', 'xx', 1), ('lv_ctrl2_J0.5', 'kxx_1_0.48', 2), ('lv_n1_J0.25', 'kxx_1_0.24', 2),
                      ('kxx_1_0.48', None, 2), ('kxx_1_0.48_0.24', None, 3), ('kxx_1_0.24', None, 2), ('xx', None, 1)]:
    fm = fes_fit(model)
    if fm is None:
        print(f'  {model}: not available'); continue
    rec = dict(model=fm, ref=ref, modes=m)
    fr = fes_fit(ref) if ref else None
    if fr:
        # calibrate on the same chi window as the model's last four chi
        win = fm['chis'][-4:]
        idx = [fr['chis'].index(c) for c in win if c in fr['chis']]
        if len(idx) < 4:      # reference stopped earlier: use its last four chi
            idx = list(range(len(fr['chis']) - 4, len(fr['chis'])))
            rec['ref_window_mismatch'] = True
        if len(idx) == 4:
            S = np.array(fr['S'])[idx]; xi = np.array(fr['xi'])[idx]
            rec['ref_fit_same_chi'] = float(6 * np.polyfit(np.log(xi), S, 1)[0])
            rec['c_cal'] = m * fm['fit_last4'] / rec['ref_fit_same_chi']
        rec['ref_fit'] = fr
    out['fes'][model] = rec
    print(f"  {model}: chi={fm['chis']} xi={np.round(fm['xi'], 1).tolist()}\n     pairwise={np.round(fm['pairwise'], 3).tolist()} "
          f"fit={fm['fit_last4']:.4f} {np.round(fm['fits_sliding'], 3).tolist()} local={np.round(fm['local_last'], 4).tolist() if fm['local_last'] else None}"
          + (f"\n     ref {ref}: same-chi fit={rec.get('ref_fit_same_chi')}  c_cal={rec.get('c_cal')}" if fr else ''))



def fit_se(S, xi):
    x = np.log(xi); S = np.asarray(S); p = np.polyfit(x, S, 1); res = S - np.polyval(p, x)
    se = np.sqrt(res @ res / max(1, len(S) - 2) / ((x - x.mean()) @ (x - x.mean())))
    return 6 * p[0], 6 * se


print('\n== FES, common-window fit over all chi >= 16 shared by model and reference (the README table)')
for model, ref, m in [('lv_n1_J0.5', 'kxx_1_0.48', 2), ('lv_n1_J0.25', 'kxx_1_0.24', 2),
                      ('lv_n2_J0.5', 'kxx_1_0.48_0.24', 3), ('lv_n2g_J0.5', 'kxx_1_0.48_0.22', 3),
                      ('lv_ctrl1_J0.5', 'xx', 1), ('lv_ctrl2_J0.5', 'kxx_1_0.48', 2)]:
    fm, fr = fes_fit(model), fes_fit(ref)
    if not fm or not fr:
        continue
    ch = [c for c in fm['chis'] if c in fr['chis'] and c >= 16]
    if len(ch) < 3:
        continue
    im = [fm['chis'].index(c) for c in ch]; ir = [fr['chis'].index(c) for c in ch]
    cm, sm = fit_se(np.array(fm['S'])[im], np.array(fm['xi'])[im])
    cr, sr = fit_se(np.array(fr['S'])[ir], np.array(fr['xi'])[ir])
    cal = m * cm / cr; dcal = cal * np.hypot(sm / cm, sr / cr)
    out['fes'].setdefault(model, {})['common_window'] = dict(chis=ch, c_raw=cm, se_raw=sm, c_ref=cr, se_ref=sr, c_cal=cal, se_cal=dcal)
    print(f"  {model} vs {ref}: chi={ch} raw {cm:.3f}+-{sm:.3f}, ref {cr:.3f}+-{sr:.3f}, calibrated {cal:.3f}+-{dcal:.3f}")

json.dump(out, open(os.path.join(D, 'lv_summary.json'), 'w'), indent=1, default=float)
