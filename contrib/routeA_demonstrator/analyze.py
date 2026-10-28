#!/usr/bin/env python3
"""Tables for the README: FES central charges, ring spectra (Casimir c, velocities, z), open-chain fits,
and the comparison with the effective logical chain.  Writes data/summary.json and prints the tables.
Claude (Opus 5.5), 2026-10-27.
"""
import os, json, glob
import numpy as np
import effective as E
import eff_ring as ER

HERE = os.path.dirname(os.path.abspath(__file__))
P = dict(D0=1.0, D1=1.0, kap=0.6, Jp=2.0, Jz=-1.0)
LAMS = [0.05, 0.1, 0.15]
out = {}


def fes_table(name, bond=None):
    f = os.path.join(HERE, 'data', 'fes', name + '.json')
    if not os.path.exists(f):
        return None
    r = json.load(open(f))['results']
    if bond is None:
        # inter-cell bond: for the d = 10 runs, the bond between two P sites
        nB = np.array(r[-1]['local'][0])
        n = len(nB)
        if n == 1:
            bond = 0
        else:
            bsite = int(np.argmax(nB))
            bond = (bsite + 2) % n
    S = np.array([x['S'][bond] for x in r]); xi = np.array([x['xi'] for x in r])
    chis = [x['chi'] for x in r]
    pw = 6 * np.diff(S) / np.diff(np.log(xi))
    k = min(4, len(S))
    fit = 6 * np.polyfit(np.log(xi[-k:]), S[-k:], 1)[0] if len(S) >= 2 else np.nan
    return dict(chis=chis, xi=xi.tolist(), S=S.tolist(), pairwise=pw.tolist(), c_fit_last4=float(fit),
                e=[x['e'] for x in r], err=[x['err'] for x in r], bond=bond,
                local=r[-1]['local'], BB=r[-1].get('BB'))


def two_cell_couplings(lam):
    He, w, sm = E.exact_eff(lam, P['D0'], P['D1'], P['kap'], P['Jp'], P['Jz'])
    s = E.summary(He)
    return s['Jxy'] / 2 * 2, s['Jzz'], s['a0']     # Jxy in H = Jxy (s+s- + h.c.) convention


def ring_micro(M, lam):
    f = os.path.join(HERE, 'data', 'ring', f'micro_M{M}_l{lam:.2f}.json')
    alt = glob.glob(os.path.join(HERE, 'data', 'ring', f'*M{M}*l{lam:.2f}*.json'))
    lv = {}
    for g in ([f] if os.path.exists(f) else []) + alt:
        lv.update(json.load(open(g))['levels'])
    return lv


def ring_summary(lv, M):
    """E0, current gap (Q = 0 or lowest-|Q| sector, m = 1), charge gap."""
    if not lv:
        return None
    q0 = 0 if M % 2 == 0 else 1
    q1 = q0 + 2
    sec = lambda Q2: {int(k.split('_')[1]): v for k, v in lv.items() if int(k.split('_')[0]) == Q2 and len(v)}
    s0, s1 = sec(q0), sec(q1)
    if 0 not in s0:
        return None
    E0 = min(min(v) for v in s0.values())
    res = dict(E0=E0)
    if 1 in s0:
        res['gapJ'] = s0[1][0] - E0
    if s1:
        res['gapQ'] = min(min(v) for v in s1.values()) - E0
    return res


if __name__ == '__main__':
    print('== FES (S at the inter-cell bond vs log xi) ==')
    out['fes'] = {}
    for nm in ['xx', 'eff2', 'effx_0.10', 'cell_l0.05', 'cell_l0.10', 'cell_l0.15', 'micro_l0.05', 'micro_l0.10',
               'micro_l0.15', 'sel_n3_s1', 'sel_n6']:
        t = fes_table(nm)
        if t is None:
            continue
        out['fes'][nm] = t
        print(f"{nm:12s} chi {t['chis']}\n   xi {np.round(t['xi'], 1).tolist()}\n   pairwise c {np.round(t['pairwise'], 3).tolist()}"
              f"  fit(last 4) {t['c_fit_last4']:.3f}  local(last) {np.round(t['local'], 3).tolist()}")
    print('\n== exact two-cell couplings (H = Jxy (s+s- + h.c.) + Jzz ZZ + a0) ==')
    out['couplings'] = {}
    for lam in LAMS:
        J, Jz, a0 = two_cell_couplings(lam)
        g = np.arccos(-2 * Jz / J)     # AFM-equivalent anisotropy Delta' = -2 Jzz / Jxy  (Jxy < 0)
        Dp = -2 * Jz / J
        gam = np.arccos(Dp)
        v = 2 * abs(J) * (np.pi / 2) * np.sin(gam) / gam
        out['couplings'][lam] = dict(Jxy=J, Jzz=Jz, a0=a0, Delta=Dp, v_bethe=v)
        print(f'lam={lam}: Jxy={J:.6e} Jzz={Jz:.6e} a0={a0:.6e} Delta\'={Dp:.4f} v(Bethe ansatz)={v:.6e}')
    print('\n== rings: microscopic vs effective (exact two-cell bond) vs XX ==')
    out['rings'] = {}
    for lam in LAMS:
        J, Jz, a0 = out['couplings'][lam]['Jxy'], out['couplings'][lam]['Jzz'], out['couplings'][lam]['a0']
        rows = {}
        for M in (3, 4, 5, 6):
            mic = ring_summary(ring_micro(M, lam), M)
            q0 = 0 if M % 2 == 0 else 1
            effl = {}
            for Q2 in (q0, q0 + 2):
                for m in range(M // 2 + 1):
                    effl[f'{Q2}_{m}'] = ER.levels(M, Q2, m, J, Jz, a0=a0, nev=2).tolist()
            eff = ring_summary(effl, M)
            rows[M] = dict(micro=mic, eff=eff)
            fm = lambda d, k: (f"{d[k]:.6e}" if d and k in d else '   -   ')
            print(f"lam={lam} M={M}: E0 micro {fm(mic, 'E0')} eff {fm(eff, 'E0')} | gapJ micro {fm(mic, 'gapJ')} eff {fm(eff, 'gapJ')}"
                  f" | gapQ micro {fm(mic, 'gapQ')} eff {fm(eff, 'gapQ')}")
        out['rings'][lam] = rows
    json.dump(out, open(os.path.join(HERE, 'data', 'summary.json'), 'w'), indent=1, default=float)


def casimir(E0a, Ma, E0b, Mb, gapJb):
    """c from E0/M = e_inf - pi v c / (6 M^2) on two sizes, v = Mb gapJ(Mb) / (2 pi) (x = 1 current level)."""
    v = Mb * gapJb / (2 * np.pi)
    return 6 * (E0a / Ma - E0b / Mb) / (np.pi * v * (1 / Mb ** 2 - 1 / Ma ** 2)), v


if __name__ == '__main__':
    print('\n== Casimir c from rings (M1, M2), v from the m = 1 current level at M2 ==')
    out['casimir'] = {}
    for name, J, Jz, a0 in [('XX', -1.0, 0.0, 0.0)] + [(f'eff {lam}', out['couplings'][lam]['Jxy'], out['couplings'][lam]['Jzz'],
                                                      out['couplings'][lam]['a0']) for lam in LAMS]:
        row = {}
        for (Ma, Mb) in [(4, 6), (6, 8), (8, 10), (10, 12), (12, 14)]:
            E = {M: ER.levels(M, 0, 0, J, Jz, a0=a0, nev=1)[0] for M in (Ma, Mb)}
            g = ER.levels(Mb, 0, 1, J, Jz, a0=a0, nev=1)[0] - E[Mb]
            c, v = casimir(E[Ma], Ma, E[Mb], Mb, g)
            row[f'{Ma},{Mb}'] = dict(c=c, v=v)
        out['casimir'][name] = row
        print(name, {k: (round(v['c'], 4), float('%.4g' % v['v'])) for k, v in row.items()})
    for lam in LAMS:
        r = out['rings'][lam]
        if r[4]['micro'] and r[6]['micro'] and 'gapJ' in r[6]['micro']:
            c, v = casimir(r[4]['micro']['E0'], 4, r[6]['micro']['E0'], 6, r[6]['micro']['gapJ'])
            out['casimir'][f'micro {lam}'] = {'4,6': dict(c=c, v=v)}
            print(f'micro {lam}: (4,6) c = {c:.4f}, v = {v:.4g}')
    print('\n== z from the current gap, M = 4, 5, 6 (micro) ==')
    out['z'] = {}
    for lam in LAMS:
        r = out['rings'][lam]
        Ms = [M for M in (3, 4, 5, 6) if r[M]['micro'] and 'gapJ' in r[M]['micro']]
        if len(Ms) >= 2:
            g = np.array([r[M]['micro']['gapJ'] for M in Ms])
            ge = np.array([r[M]['eff']['gapJ'] for M in Ms])
            Ml = np.array(Ms, float)
            sel = Ml >= 4
            z = -np.polyfit(np.log(Ml[sel]), np.log(g[sel]), 1)[0] if sel.sum() >= 2 else np.nan
            ze = -np.polyfit(np.log(Ml[sel]), np.log(ge[sel]), 1)[0] if sel.sum() >= 2 else np.nan
            out['z'][lam] = dict(M=Ms, gapJ=g.tolist(), gapJ_eff=ge.tolist(), Mgap_over_2pi=(Ml * g / 2 / np.pi).tolist(), z=z, z_eff=ze)
            print(f'lam={lam}: M={Ms} M*gapJ/2pi micro={np.round(Ml * g / 2 / np.pi, 6).tolist()} eff={np.round(Ml * ge / 2 / np.pi, 6).tolist()} '
                  f'z(M>=4) micro={z:.3f} eff={ze:.3f}')
    json.dump(out, open(os.path.join(HERE, 'data', 'summary.json'), 'w'), indent=1, default=float)
