#!/usr/bin/env python3
"""Tables for README.md from data/obc/*.json and data/fes/*.json; writes data/summary.json.

Estimators (all at down-spin density n = 1/3, lengths in cells):
  K_gap(M)  = E1 / Dc, E1 the lowest neutral gap, Dc = E(Q0+1) + E(Q0-1) - 2 E(Q0)   (open chain of M cells)
  K_inf     : K_gap(M) = K + a/M + b/M^2 through M = 12, 24, 48
  v         : 1/E1 = (M + delta)/(pi v), linear fit in M through M = 24, 48 (and M = 12 as a check)
  c_int, K_fl : left-block entropy and charge variance at inter-cell cuts, fitted to
              (c/6 or K/2pi^2) ln[(2(M+1)/pi) sin(pi l/(M+1))] + a + (b1 cos 2kF l + b2 sin 2kF l) / chord^K,
              kF = pi N/(M+1), exponent K fixed to K_gap(M), cuts 4 <= l <= M-4
  c_FES     : 6 dS/d ln xi, least squares over the listed chi window (VUMPS, one-cell unit cell)
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chain as CH
import bethe as BA

LAMS = [0.0125, 0.025, 0.05, 0.10, 0.15]
LSTR = {0.0125: '0.0125', 0.025: '0.025', 0.05: '0.05', 0.10: '0.10', 0.15: '0.15'}
RHOS = [2.0, 1.0, 0.5]


def load(name):
    f = os.path.join(HERE, 'data', 'obc', name + '.json')
    if not os.path.exists(f):
        return None
    d = json.load(open(f))
    return d if 'summary' in d else None


def interval_fits(d, lmin=4):
    s = d['sectors']['0']; cl = d['cell_len']; M = d['args']['M']
    S = np.array(s['S']); QA = np.array(s['QA']); QA2 = np.array(s['QA2'])
    idx = np.arange(1, M) * cl - 1
    l = np.arange(1, M).astype(float)
    F = QA2[idx] - QA[idx] ** 2
    Sl = S[idx]
    N = M // 3
    kF = np.pi * N / (M + 1)
    K = d['summary']['K_gap']
    ch = 2 * (M + 1) / np.pi * np.sin(np.pi * l / (M + 1))
    sel = (l >= lmin) & (l <= M - lmin)
    A = np.vstack([np.log(ch), np.ones_like(l), np.cos(2 * kF * l) / ch ** K, np.sin(2 * kF * l) / ch ** K]).T[sel]
    c = 6 * np.linalg.lstsq(A, Sl[sel], rcond=None)[0][0]
    Kf = 2 * np.pi ** 2 * np.linalg.lstsq(A, F[sel], rcond=None)[0][0]
    return float(c), float(Kf), float(Sl[M // 2 - 1])


def kinf(Ms, Ks):
    Ms = np.asarray(Ms, float); Ks = np.asarray(Ks)
    if len(Ms) >= 3:
        A = np.vstack([np.ones_like(Ms), 1 / Ms, 1 / Ms ** 2]).T
    else:
        A = np.vstack([np.ones_like(Ms), 1 / Ms]).T
    return float(np.linalg.lstsq(A, Ks, rcond=None)[0][0])


def vfit(Ms, E1s):
    Ms = np.asarray(Ms, float)
    p = np.polyfit(Ms, 1 / np.asarray(E1s), 1)       # 1/E1 = M/(pi v) + delta/(pi v)
    return float(1 / (np.pi * p[0])), float(p[1] / p[0])


def fes_fit(f, xi_max=None, chi_min=12):
    if not os.path.exists(f):
        return None
    d = json.load(open(f))
    r = [x for x in d['results'] if x['chi'] >= chi_min and x['err'] < 1e-4 and (xi_max is None or x['xi'] <= xi_max)]
    if len(r) < 3:
        return None
    S = np.array([x['S'][0] for x in r]); lx = np.log([x['xi'] for x in r])
    A = np.vstack([lx, np.ones_like(lx)]).T
    p, res, *_ = np.linalg.lstsq(A, S, rcond=None)
    resid = S - A @ p
    err = np.sqrt(np.sum(resid ** 2) / max(len(S) - 2, 1) / np.sum((lx - lx.mean()) ** 2)) if len(S) > 2 else np.nan
    return dict(c=float(6 * p[0]), err=float(6 * err), chis=[x['chi'] for x in r], xi=[x['xi'] for x in r],
                n=[x['n'] for x in r], maxerr=max(x['err'] for x in r))


def fes_points(files, errmax=1e-4):
    pts = []
    for f in files:
        if os.path.exists(f):
            pts += [(x['xi'], x['S'][0], x['chi'], x['err'], x['n']) for x in json.load(open(f))['results'] if x['err'] < errmax]
    pts.sort()
    return pts


def fes_matched(micro_f, log_fs, xi_min=5.0):
    """c from 6 dS/d ln xi for the micro chain (points with xi >= xi_min) and for the logical_x chain over
    the same xi window; offset = S_micro - S_logical(xi) by interpolation in ln xi."""
    mp = [p for p in fes_points([micro_f]) if p[0] >= xi_min]
    lp = fes_points(log_fs)
    if len(mp) < 2 or len(lp) < 3:
        return None
    lo, hi = mp[0][0], mp[-1][0]
    lw = [p for p in lp if 0.8 * lo <= p[0] <= 1.25 * hi]

    def fit(pts):
        x = np.log([p[0] for p in pts]); y = np.array([p[1] for p in pts])
        A = np.vstack([x, np.ones_like(x)]).T
        coef = np.linalg.lstsq(A, y, rcond=None)[0]
        if len(pts) > 2:
            r = y - A @ coef
            err = np.sqrt(np.sum(r ** 2) / (len(pts) - 2) / np.sum((x - x.mean()) ** 2))
        else:
            err = np.nan
        return 6 * coef[0], 6 * err
    cm, em = fit(mp)
    cl, el = fit(lw)
    lx = np.log([p[0] for p in lp]); ly = np.array([p[1] for p in lp])
    off = [p[1] - np.interp(np.log(p[0]), lx, ly) for p in mp]
    # consecutive-pair slopes of the micro chain
    pairs = [6 * (b[1] - a[1]) / np.log(b[0] / a[0]) for a, b in zip(mp[:-1], mp[1:])]
    return dict(c_micro=float(cm), err_micro=float(em), c_log=float(cl), err_log=float(el),
                xi_micro=[round(p[0], 2) for p in mp], chi_micro=[p[2] for p in mp], xi_log=[round(p[0], 2) for p in lw],
                offset=[float(o) for o in off], pair_slopes=[float(x) for x in pairs], n_micro=[p[4] for p in mp],
                err_grad=max(p[3] for p in mp))


def combined(pt):
    """Thermodynamic K and v of the micro chain = logical_x limit corrected by the micro/logical_x ratio at the
    largest size where both exist (the finite-size corrections are common to both chains to the shown digits)."""
    Ml = [M for M in (12, 24, 48, 96) if 'logx' in pt['M'].get(str(M), {})]
    Mc = [M for M in (12, 24, 48) if 'micro' in pt['M'].get(str(M), {}) and 'logx' in pt['M'][str(M)]]
    if len(Ml) < 3 or not Mc:
        return
    K = lambda M, k: pt['M'][str(M)][k]['K_gap']
    E = lambda M, k: pt['M'][str(M)][k]['E1']
    Kl_a = kinf([12, 24, 48], [K(M, 'logx') for M in (12, 24, 48)])
    Kl_b = kinf(Ml[-3:], [K(M, 'logx') for M in Ml[-3:]])
    dK = [K(M, 'micro') - K(M, 'logx') for M in Mc]
    rE = [E(M, 'micro') / E(M, 'logx') for M in Mc]
    vl_a = vfit([24, 48], [E(M, 'logx') for M in (24, 48)])[0]
    vl_b = vfit(Ml[-2:], [E(M, 'logx') for M in Ml[-2:]])[0]
    pt['dK'] = dict(zip(map(str, Mc), dK)); pt['rE'] = dict(zip(map(str, Mc), rE))
    # primary size for the micro/logical_x comparison: M = 24 (energies converged to 1e-12); M = 48 runs are
    # truncation-limited at 1e-9 (ground) to 1e-8 (excited), so they enter only through the error
    i0 = Mc.index(24) if 24 in Mc else len(Mc) - 1
    spread_K = max(abs(x - dK[i0]) for x in dK) if len(dK) > 1 else abs(dK[i0])
    spread_E = max(abs(x - rE[i0]) for x in rE) if len(rE) > 1 else abs(1 - rE[i0])
    pt['K_logx_inf'] = Kl_a
    pt['Kc'] = Kl_a + dK[i0]
    pt['Kc_err'] = abs(Kl_a - Kl_b) + spread_K + 5e-5
    pt['v_logx_inf'] = vl_b
    pt['vc'] = vl_b * rE[i0]   # logical_x velocity from (48, 96) when M = 96 exists
    pt['vc_err'] = abs(vl_a - vl_b) * rE[i0] + vl_b * (spread_E + 1e-3)
    pt['M_micro_max'] = Mc[-1]
    if 'K_inf_micro' in pt:                       # direct micro extrapolation available: include its spread
        pt['Kc_err'] += abs(pt['K_inf_micro'] - pt['Kc'])
    if 'v_micro' in pt:
        pt['vc_err'] += abs(pt['v_micro'] - pt['vc'])


def md_tables(out):
    P = out['points']
    print('\n### Table 1. Bare logical XXZ chain (calibration of the open-chain estimators), n = 1/3\n')
    print('| rho | g | K_BA | K_gap M=24 | 48 | 96 | 192 | K_inf (12,24,48) | v_BA/J | v(48,96)/J | c_int M=48 | 96 | 192 | K_fl M=192 |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for l in out['logical']:
        m = l['M']
        g = lambda M, k: (f"{m[M][k]:.5f}" if k == 'K_gap' else f"{m[M][k]:.3f}") if M in m else '-'
        vv = vfit([48, 96], [m[M]['E1'] for M in (48, 96)])[0] / (m[48]['E1'] * 48 / (np.pi * m[48]['v_M'])) if 96 in m else np.nan
        Ks = [M for M in (12, 24, 48) if M in m]
        Kinf = kinf(Ks, [m[M]['K_gap'] for M in Ks])
        print(f"| {l['rho']:g} | {l['g']:.4f} | {l['K_BA']:.5f} | {g(24,'K_gap')} | {g(48,'K_gap')} | {g(96,'K_gap')} | {g(192,'K_gap')} | "
              f"{Kinf:.5f} | {l['v_BA']:.4f} | {vv:.4f} | {g(48,'c_int')} | {g(96,'c_int')} | {g(192,'c_int')} | {g(192,'K_fl')} |")
    print('\n### Table 2. Micro chain against the logical chain with exact two-cell couplings, equal M\n')
    print('| lam | rho | M | K_gap micro | K_gap logical_x | dK | E1 micro / E1 logical_x | c_int micro | c_int logical_x | K_fl micro | K_fl logical_x |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for p in P:
        for M, x in sorted(p['M'].items(), key=lambda kv: int(kv[0])):
            if 'micro' not in x:
                continue
            a, b = x['micro'], x['logx']
            ci = lambda y: f"{y:.3f}" if int(M) >= 24 else '-'
            print(f"| {p['lam']:g} | {p['rho']:g} | {M} | {a['K_gap']:.5f} | {b['K_gap']:.5f} | {a['K_gap'] - b['K_gap']:+.1e} | "
                  f"{a['E1'] / b['E1']:.5f} | {ci(a['c_int'])} | {ci(b['c_int'])} | {ci(a['K_fl'])} | {ci(b['K_fl'])} |")
    print('\n### Table 3. Thermodynamic K and v of the micro chain against the XXZ predictions\n')
    print('| lam | rho | g | g_x(lam) | K micro | K_BA(g) | K_BA(g_x) | K micro - K_BA(g) | v micro / J2 | v_BA(g)/J2 | v_BA(g_x)/J2 | M max |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for p in P:
        if 'Kc' not in p:
            continue
        J2 = p['J2']
        print(f"| {p['lam']:g} | {p['rho']:g} | {p['g']:.4f} | {p['g_x']:+.4f} | {p['Kc']:.4f}({int(round(p['Kc_err'] * 1e4))}) | "
              f"{p['K_BA']:.4f} | {p['K_BA_x']:.4f} | {p['Kc'] - p['K_BA']:+.4f} | "
              f"{p['vc'] / J2:.4f}({int(round(p['vc_err'] / J2 * 1e4))}) | {p['v_BA'] / J2:.4f} | {p['v_BA_x'] / J2:.4f} | {p['M_micro_max']} |")
    print('\n### Table 4. Direct micro extrapolations (where M = 12, 24, 48 all exist)\n')
    print('| lam | rho | K_inf micro (12,24,48) | K_inf logical_x (12,24,48) | v micro (24,48)/J2 | v logical_x (24,48)/J2 |')
    print('|---|---|---|---|---|---|')
    for p in P:
        if 'K_inf_micro' in p:
            print(f"| {p['lam']:g} | {p['rho']:g} | {p['K_inf_micro']:.5f} | {p['K_inf_logx']:.5f} | {p.get('v_micro', np.nan) / p['J2']:.4f} | "
                  f"{p['v_logx'] / p['J2']:.4f} |")
    print('\n### Table 5. Finite-entanglement scaling (VUMPS), matched xi window\n')
    print('| lam | rho | chi (micro) | xi micro (cells) | c micro | c logical_x, same xi window | S micro - S logical_x at equal xi | n reached |')
    print('|---|---|---|---|---|---|---|---|')
    for p in P:
        f = p.get('fes_matched')
        if not f:
            continue
        print(f"| {p['lam']:g} | {p['rho']:g} | {f['chi_micro'][0]} to {f['chi_micro'][-1]} | {f['xi_micro'][0]:.1f} to {f['xi_micro'][-1]:.1f} | "
              f"{f['c_micro']:.3f} +- {f['err_micro']:.3f} | {f['c_log']:.3f} +- {f['err_log']:.3f} | "
              f"{min(f['offset']):.3f} to {max(f['offset']):.3f} | {min(f['n_micro']):.4f} to {max(f['n_micro']):.4f} |")


def fill_readme(out):
    """Replace the blocks between <!-- Tk --> and <!-- /Tk --> in README.md with the generated tables."""
    import io, contextlib, re
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        md_tables(out)
    parts = re.split(r'\n### Table (\d)\. [^\n]*\n\n', buf.getvalue())
    tabs = {parts[i]: parts[i + 1].strip() for i in range(1, len(parts) - 1, 2)}
    f = os.path.join(HERE, 'README.md')
    txt = open(f).read()
    for k, t in tabs.items():
        txt = re.sub(rf'<!-- T{k} -->.*?<!-- /T{k} -->', lambda m: f'<!-- T{k} -->\n{t}\n<!-- /T{k} -->', txt, flags=re.S)
    open(f, 'w').write(txt)


def main():
    out = dict(points=[], logical=[])
    # bare logical chains: finite-size calibration of the estimators against the Bethe ansatz
    print('== bare logical chain (second order), estimator calibration vs Bethe ansatz, n = 1/3')
    for rho in RHOS:
        Jp, Jz = CH.couplings(rho)
        Jxy, Jzz, g = CH.second_order(Jp, Jz)
        ba = BA.at_density(g, 1 / 3)
        rows = {}
        for M in (12, 24, 48, 96, 192):
            d = load(f'obc_logical_M{M}_r{rho:g}')
            if d:
                rows[M] = d
        if not rows:
            continue
        J = rows[min(rows)]['summary']['J']
        line = dict(rho=rho, g=g, K_BA=ba['K'], v_BA=ba['v'], M={})
        for M, d in rows.items():
            c, Kf, Smid = interval_fits(d)
            line['M'][M] = dict(K_gap=d['summary']['K_gap'], E1=d['summary']['E1'], v_M=d['summary']['E1'] * M / (np.pi * J),
                                c_int=c, K_fl=Kf)
        Ms = [M for M in (12, 24, 48) if M in rows]
        if len(Ms) == 3:
            line['K_inf_12_48'] = kinf(Ms, [rows[M]['summary']['K_gap'] for M in Ms])
            line['v_fit_24_48'] = vfit([24, 48], [rows[M]['summary']['E1'] for M in (24, 48)])[0] / J
        print(f"rho={rho} g={g:.5f} K_BA={ba['K']:.5f} v_BA/J={ba['v']:.5f}  K_inf(12,24,48)={line.get('K_inf_12_48', np.nan):.5f} "
              f"v(24,48)/J={line.get('v_fit_24_48', np.nan):.5f}")
        for M, x in sorted(line['M'].items()):
            print(f"    M={M:4d} K_gap={x['K_gap']:.5f} c_int={x['c_int']:.4f} K_fl={x['K_fl']:.4f}")
        fl = fes_fit(os.path.join(HERE, 'data', 'fes', f'fes_logical_r{rho:g}.json'), xi_max=60)
        line['fes'] = fl
        if fl:
            print(f"    FES (xi <= 60): c = {fl['c']:.3f} +- {fl['err']:.3f}, chis {fl['chis']}, xi {np.round(fl['xi'], 1).tolist()}")
        out['logical'].append(line)

    print('\n== micro chain vs logical chains')
    for lam in LAMS:
        for rho in RHOS:
            Jp, Jz = CH.couplings(rho)
            Jxy, Jzz, g = CH.second_order(Jp, Jz)
            Jx, Jzx, _ = CH.exact_two_cell(lam, Jp, Jz)
            gx = 2 * Jzx / abs(Jx)
            J2, Jxd = 2 * abs(Jxy) * lam ** 2, 2 * abs(Jx)
            ba, bax = BA.at_density(g, 1 / 3), BA.at_density(gx, 1 / 3)
            pt = dict(lam=lam, rho=rho, Jp=Jp, Jz=Jz, g=g, g_x=gx, J2=J2, J_x=Jxd,
                      K_BA=ba['K'], v_BA=ba['v'] * J2, K_BA_x=bax['K'], v_BA_x=bax['v'] * Jxd, M={})
            for kind, pre in (('micro', 'obc_micro'), ('logx', 'obc_logx')):
                for M in (12, 24, 48, 96):
                    d = load(f'{pre}_M{M}_l{LSTR[lam]}_r{rho:g}')
                    if d is None:
                        continue
                    c, Kf, Smid = interval_fits(d)
                    s = d['summary']
                    pt['M'].setdefault(str(M), {})[kind] = dict(K_gap=s['K_gap'], E1=s['E1'], Dc=s['charge_gap'], c_int=c,
                                                               K_fl=Kf, S_mid=Smid,
                                                               trunc=d['sectors']['0']['hist'][-1][2])
            for kind in ('micro', 'logx'):
                Ms = [M for M in (12, 24, 48) if kind in pt['M'].get(str(M), {})]
                if len(Ms) == 3:
                    pt[f'K_inf_{kind}'] = kinf(Ms, [pt['M'][str(M)][kind]['K_gap'] for M in Ms])
                if 24 in Ms and 48 in Ms:
                    pt[f'v_{kind}'], pt[f'delta_{kind}'] = vfit([24, 48], [pt['M'][str(M)][kind]['E1'] for M in (24, 48)])
            for kind in ('micro', 'logx'):
                pre = 'fes_micro' if kind == 'micro' else 'fes_logx'
                pt[f'fes_{kind}'] = fes_fit(os.path.join(HERE, 'data', 'fes', f'{pre}_l{LSTR[lam]}_r{rho:g}.json'))
            pt['fes_matched'] = fes_matched(os.path.join(HERE, 'data', 'fes', f'fes_micro_l{LSTR[lam]}_r{rho:g}.json'),
                                            [os.path.join(HERE, 'data', 'fes', f'{p}_l{LSTR[lam]}_r{rho:g}.json')
                                             for p in ('fes_logx', 'fesS_logx')])
            combined(pt)
            out['points'].append(pt)
            print(f"lam={lam:g} rho={rho:g} g={g:.5f} g_x={gx:.5f} K_BA={ba['K']:.5f} K_BA_x={bax['K']:.5f} "
                  f"v_BA={pt['v_BA']:.4e} v_BA_x={pt['v_BA_x']:.4e}")
            for M, x in sorted(pt['M'].items(), key=lambda kv: int(kv[0])):
                for kind, y in x.items():
                    print(f"    M={M:>3} {kind:5s} K_gap={y['K_gap']:.5f} E1={y['E1']:.4e} c_int={y['c_int']:.4f} "
                          f"K_fl={y['K_fl']:.4f} trunc={y['trunc']:.0e}")
            for kind in ('micro', 'logx'):
                if f'K_inf_{kind}' in pt:
                    print(f"    {kind}: K_inf={pt[f'K_inf_{kind}']:.5f}", end='')
                if f'v_{kind}' in pt:
                    print(f"  v={pt[f'v_{kind}']:.4e} (delta={pt[f'delta_{kind}']:.2f})", end='')
                f = pt[f'fes_{kind}']
                if f:
                    print(f"  FES c={f['c']:.3f}+-{f['err']:.3f} xi={np.round(f['xi'], 1).tolist()} n={np.round(f['n'], 4).tolist()}", end='')
                print()
            fm = pt['fes_matched']
            if fm:
                print(f"    FES matched: c_micro={fm['c_micro']:.3f}+-{fm['err_micro']:.3f} (xi {fm['xi_micro']}), "
                      f"c_logx={fm['c_log']:.3f}+-{fm['err_log']:.3f} (xi {fm['xi_log']}), "
                      f"offset S_micro-S_logx={np.round(fm['offset'], 4).tolist()}")
    json.dump(out, open(os.path.join(HERE, 'data', 'summary.json'), 'w'), indent=1)
    if '--md' in sys.argv:
        md_tables(out)
    if '--readme' in sys.argv:
        fill_readme(out)


if __name__ == '__main__':
    main()
