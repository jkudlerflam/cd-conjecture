#!/usr/bin/env python3
"""Summaries for contrib/mobile_cells: central-charge estimators and ring degeneracies.

  FES:    c_pair = 6 dS/d ln xi between consecutive chi; c_fit = 6 x slope of S against ln xi over
          the last nfit points (base run and its *_hi continuation merged by name).
  ED:     c_half(L1 -> L2) = 3 [S(L2/2) - S(L1/2)] / ln(L2/L1) on the fixed wall density W/L = 2/3
          (L = 9, 12, 15, and 18 if the cluster files exist); static calibration from cuts at
          cell boundaries, S(L, L/2) for L = 12, 18, 24.
Writes data/summary.json.  Claude (Opus 5.5), 2026-10-27.
"""
import os, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')


def fes(nfit=4, errmax=1e-4):
    runs = {}
    for f in sorted(glob.glob(os.path.join(D, 'fes', '*.json'))):
        d = json.load(open(f))
        base = d['name'].replace('_hi', '')
        for r in d['results']:
            runs.setdefault(base, {})[r['chi']] = r            # later (hi) runs override equal chi
        runs[base]['_args'] = d['args']
    out = {}
    for name, byc in runs.items():
        args = byc.pop('_args')
        rs = [byc[c] for c in sorted(byc) if byc[c]['err'] < errmax and np.isfinite(byc[c]['xi']) and byc[c]['xi'] < 1e6]
        chi = np.array([r['chi'] for r in rs], float)
        S = np.array([np.mean(r['S']) for r in rs]); xi = np.array([r['xi'] for r in rs])
        o = dict(chi=chi.tolist(), S=S.tolist(), xi=xi.tolist(), err=[r['err'] for r in rs],
                 illegal=[r.get('illegal') for r in rs], local_last=rs[-1]['local'] if rs else None,
                 args={k: args.get(k) for k in ('model', 't', 'u', 'mu', 'J', 'Jz', 'dJ', 'n')})
        if len(rs) >= 2:
            o['c_pair'] = (6 * np.diff(S) / np.diff(np.log(xi))).tolist()
        if len(rs) >= 3:
            k = min(nfit, len(rs))
            o['c_fit'] = float(6 * np.polyfit(np.log(xi[-k:]), S[-k:], 1)[0])
            o['c_fit_all'] = float(6 * np.polyfit(np.log(xi), S, 1)[0])
            o['nfit'] = k
        out[name] = o
    return out


def ed():
    out = {}
    for name in ('mobile_A_u1', 'mobile_B_u03', 'mobile_C_u01'):
        f = os.path.join(D, 'ed', name + '.json')
        if not os.path.exists(f):
            continue
        res = json.load(open(f))['results']
        tag = name.split('_')[1]
        S = {r['L']: (r['W0'], r['S_half'][0]) for r in res}
        # cluster L = 18 sectors (W = 12 is the 2/3 density)
        for W in (10, 12, 14):
            g = os.path.join(D, 'ed', f'L18_{tag}_W{W}.json')
            if os.path.exists(g):
                x = json.load(open(g))
                S.setdefault(('L18', W), (W, x['S'][0], x['E'][0]))
        seq = [L for L in (9, 12, 15) if L in S and S[L][0] == 2 * L // 3]
        if ('L18', 12) in S:
            E18 = {W: S[('L18', W)][2] for W in (10, 12, 14) if ('L18', W) in S}
            if min(E18, key=E18.get) == 12:
                S[18] = (12, S[('L18', 12)][1]); seq.append(18)
        pairs = []
        for a, b in zip(seq[:-1], seq[1:]):
            pairs.append(dict(L1=a, L2=b, c_half=3 * (S[b][1] - S[a][1]) / np.log(b / a)))
        if len(seq) >= 2:
            a, b = seq[0], seq[-1]
            pairs.append(dict(L1=a, L2=b, c_half=3 * (S[b][1] - S[a][1]) / np.log(b / a)))
        degs = {r['L']: r['degeneracy'] for r in res}
        gaps = {}
        for r in res:
            lv = sorted(l['E'] for l in r['levels'])
            gaps[r['L']] = r['L'] * (lv[1] - lv[0]) / (2 * np.pi) if len(lv) > 1 else None
        out[name] = dict(c_half=pairs, degeneracy=degs, L_gap_over_2pi=gaps)
    f = os.path.join(D, 'ed', 'static_entropy.json')
    if os.path.exists(f):
        st = json.load(open(f))
        Sh = {r['L']: dict(zip(r['ell'], r['S'])).get(r['L'] // 2) for r in st}
        Ls = [L for L in (12, 18, 24) if Sh.get(L) is not None]
        out['static_payload_only'] = [dict(L1=a, L2=b, c_half=3 * (Sh[b] - Sh[a]) / np.log(b / a)) for a, b in zip(Ls[:-1], Ls[1:])]
    f = os.path.join(D, 'ed', 'static_mu-1_J05.json')
    if os.path.exists(f):
        st = json.load(open(f))['results']
        out['static_degeneracy'] = {r['L']: dict(deg=r['degeneracy'], ground=[g['cells'] for g in r['ground']],
                                                 next_pattern_gap=r['next_pattern_gap']) for r in st}
    return out


def modes():
    """Velocities from run_modes.py and the Casimir test on the fixed density W/L = 2/3."""
    out = {}
    for tag in ('A_u1', 'B_u03', 'C_u01'):
        recs = []
        for f in sorted(glob.glob(os.path.join(D, 'ed', f'modes_{tag}*.json'))):
            recs += json.load(open(f))
        if not recs:
            continue
        byLW = {(r['L'], r['W']): r for r in recs}
        rows = [dict(L=r['L'], W=r['W'], E0=r['E0'], v_rho=r['rho']['v'], w_rho=r['rho']['peak_weight'],
                     v_S=r['Sz']['v'], w_S=r['Sz']['peak_weight']) for r in sorted(recs, key=lambda r: (r['L'], r['W']))]
        o = dict(rows=rows)
        seq = [(L, 2 * L // 3) for L in (9, 12, 15, 18) if (L, 2 * L // 3) in byLW]
        if len(seq) >= 2:
            Ls = np.array([L for L, W in seq], float)
            e = np.array([byLW[k]['E0'] for k in seq]) / Ls
            pairs = []
            for i in range(len(seq) - 1):
                L1, L2 = Ls[i], Ls[i + 1]
                C = 6 / np.pi * (e[i + 1] - e[i]) / (1 / L1 ** 2 - 1 / L2 ** 2)
                vs = byLW[seq[i + 1]]['rho']['v'] + byLW[seq[i + 1]]['Sz']['v']
                pairs.append(dict(L1=int(L1), L2=int(L2), C=float(C), v_sum_at_L2=float(vs), ratio=float(C / vs)))
            o['casimir_pairs'] = pairs
            if len(seq) >= 3:
                X = np.column_stack([np.ones_like(Ls), -np.pi / (6 * Ls ** 2)])
                coef = np.linalg.lstsq(X, e, rcond=None)[0]
                vs = byLW[seq[-1]]['rho']['v'] + byLW[seq[-1]]['Sz']['v']
                o['casimir_fit'] = dict(L=Ls.tolist(), e_inf=float(coef[0]), C=float(coef[1]), v_sum_at_Lmax=float(vs),
                                        ratio=float(coef[1] / vs))
        out[tag] = o
    return out


if __name__ == '__main__':
    F = fes(); E = ed(); Mo = modes()
    for k, o in Mo.items():
        print('MODES', k)
        for r in o['rows']:
            print(f"   L={r['L']} W={r['W']} v_rho={r['v_rho']:.4f} (w {r['w_rho']:.2f}) v_S={r['v_S']:.4f} (w {r['w_S']:.2f})")
        for p in o.get('casimir_pairs', []):
            print(f"   Casimir ({p['L1']},{p['L2']}): C={p['C']:.4f} v_rho+v_S={p['v_sum_at_L2']:.4f} ratio={p['ratio']:.3f}")
        if 'casimir_fit' in o:
            cf = o['casimir_fit']; print(f"   Casimir fit L={cf['L']}: C={cf['C']:.4f} ratio={cf['ratio']:.3f}")
    for k, o in sorted(F.items()):
        cp = np.round(o.get('c_pair', []), 2).tolist()
        print(f"FES {k:12s} chi={o['chi']} xi={np.round(o['xi'], 1).tolist()}\n    c_pair={cp} c_fit(last {o.get('nfit')})={o.get('c_fit', float('nan')):.2f} "
              f"c_fit(all)={o.get('c_fit_all', float('nan')):.2f} illegal_max={max([x or 0 for x in o['illegal']]):.1e}")
    for k, o in E.items():
        print('ED', k, json.dumps(o)[:600])
    json.dump(dict(fes=F, ed=E, modes=Mo), open(os.path.join(D, 'summary.json'), 'w'), indent=1)
