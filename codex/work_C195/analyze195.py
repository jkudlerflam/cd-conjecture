#!/usr/bin/env python3
"""C195 analysis: T1 to T4 of RESULT.md (Claude Opus 5.5 agent, 2026-11-10).  numpy only.

Usage (from codex/work_C195, after the collection commands of RESULT.md):
    python3 -E -s analyze195.py [--partial]
Defaults: L = 512 outputs in typhon_out/out512, L = 256 outputs in ../work_C190/gpu/typhon_out/outG and
../work_C190/gpu/typhon_out/outX1024 (the chi = 1024 continuation, used where it exists), key p153c.
General form (used for the test on other sizes):
    python3 -E -s analyze195.py --big DIR --big-L 512 --small DIR [DIR ...] --small-L 256 --key p153c [--partial]

Definitions are those of ../work_C190/RESULT.md (S_W mean over cuts |l - L/2| <= W, c_two = 6 dS / ln(L2/L1),
c_prof of C176 with a 12-cut moving average and fit window [L/6, 5L/6], n_mid over [L/4, 3L/4), ground branch =
the lower-energy parity sector at each size).  For a finished run the final exact entropies and the measured
density are used, as in analyze190.py; for a running one (--partial) the last mixer-off sweep.

T1 is read at the final chi of the L = 512 run: residual <= 1e-8 in the last sweep, between the last two sweeps
of that chi |dE| <= 2e-8 (1e-8 L/256), |dS_W| <= 1e-3 for W = 1, 6, 12, |dn_mid| <= 5e-4, discarded weight of
the last sweep <= 1e-7; C4 (step of S6 between the last sweeps of the last two chi stages, against 0.01) is
reported and does not enter T1.  T2 needs T1 for the ground-branch run at the large size.
S(256) is taken, per parity sector, from the run with the highest chi that has status done (the chi = 1024
continuation where it is done, otherwise the chi = 768 run of outG).  The ground sector at 256 is the lower
energy of these two; if they are at different chi, the sector order at the common chi is printed as well.
"""
import os
import sys
import glob
import json
import math
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--big', default=os.path.join(HERE, 'typhon_out', 'out512'))
ap.add_argument('--big-L', type=int, default=512)
ap.add_argument('--small', nargs='+', default=[os.path.join(HERE, '..', 'work_C190', 'gpu', 'typhon_out', 'outG'),
                                               os.path.join(HERE, '..', 'work_C190', 'gpu', 'typhon_out', 'outX1024')])
ap.add_argument('--small-L', type=int, default=256)
ap.add_argument('--key', default='p153c')
ap.add_argument('--partial', action='store_true', help='include runs that are not done (their last mixer-off sweep)')
a = ap.parse_args()
WS = (1, 6, 12)


def sbar(S, L, W):
    S = np.asarray(S)
    l = np.arange(1, L)
    sel = (l >= L // 2 - W) & (l <= L // 2 + W)
    return float(S[sel].mean())


def c_prof(S, L, sm=12, lo=1 / 6, hi=5 / 6):
    """analyze190.c_prof (C176), unchanged."""
    S = np.asarray(S)
    l = np.arange(1, L)
    k = np.ones(sm) / sm
    Ss = np.convolve(S, k, mode='valid')
    ls = np.convolve(l.astype(float), k, mode='valid')
    x = np.log((2 * L / np.pi) * np.sin(np.pi * ls / L))
    sel = (ls >= lo * L) & (ls <= hi * L)
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    return 6 * np.linalg.lstsq(A, Ss[sel], rcond=None)[0][0]


def stages(h, cont=None):
    """Last mixer-off sweep of each chi, in order: chi -> record (with E_phys, S_W, n_mid, trunc).  cont: the
    'continues' entry of a continuation config (values of the chi it started from)."""
    out = {}
    if cont:
        out[cont['chi']] = dict(chi=cont['chi'], E_phys=cont['E_phys'], Sbar1=cont['Sbar1'], Sbar6=cont['Sbar6'],
                                Sbar12=cont['Sbar12'], trunc=cont['trunc'], n_mid=None, src='continued from')
    for r in h:
        if 'Sbar6' in r:
            out[r['chi']] = r
    return out


def load_run(f, L):
    try:
        r = json.load(open(f))
    except Exception as e:
        print(f'  unreadable {f}: {e}')
        return None
    if r.get('status') != 'done' and not a.partial:
        return None
    h = r['hist']
    hs = [x for x in h if 'Sbar6' in x]
    if not hs:
        print(f'  {os.path.basename(f)}: {len(h)} sweeps, no mixer-off sweep yet (last chi {h[-1]["chi"] if h else None})')
        return None
    if r.get('status') == 'done':
        S = np.array(r['S'])
        dens = np.array(r['meas']['dens'])
        E = r.get('E_phys', h[-1]['E_phys'])
        src = 'final measurement'
    else:
        S = np.array(hs[-1]['Sprof'])
        dens = (1 - np.array(hs[-1]['Zprof'])) / 2
        E = hs[-1]['E_phys']
        src = f'sweep record (running, chi {hs[-1]["chi"]})'
    cfg = r['cfg']
    run = dict(file=f, tag=os.path.basename(f)[:-5], P=cfg['P'], status=r.get('status'), E=E, S=S,
               SW={W: sbar(S, L, W) for W in WS}, cprof=c_prof(S, L), n_mid=float(np.mean(dens[L // 4: 3 * L // 4])),
               chi=h[-1]['chi'], chi_final=cfg['final']['chi'], hist=h, src=src,
               st=stages(h, cfg.get('continues')), cfg=cfg)
    return run


def t1(run, L):
    """T1 at the run's last chi.  Returns (ok, dict of numbers)."""
    h = run['hist']
    chi = h[-1]['chi']
    fin = [x for x in h if x['chi'] == chi]
    x = fin[-1]
    d = dict(chi=chi, nsweeps_at_chi=len(fin), res=x['res_max'], trunc=x['trunc'])
    tolE = 1e-8 * L / 256
    if len(fin) >= 2:
        y = fin[-2]
        d['dE'] = abs(x['E'] - y['E'])
        d['dS'] = max(abs(x['Sbar%d' % W] - y['Sbar%d' % W]) for W in WS)
        d['dn'] = abs(x['n_mid'] - y['n_mid'])
        ok = (d['res'] <= 1e-8 and d['dE'] <= tolE and d['dS'] <= 1e-3 and d['dn'] <= 5e-4 and d['trunc'] <= 1e-7)
    else:
        d['dE'] = d['dS'] = d['dn'] = float('nan')
        ok = False
    prev = [c for c in run['st'] if c < chi]
    d['C4_step'] = abs(x['Sbar6'] - run['st'][max(prev)]['Sbar6']) if prev else float('nan')
    d['C4_from_chi'] = max(prev) if prev else None
    at_final = chi == run['chi_final'] and run['status'] == 'done'
    d['tolE'] = tolE
    return ok and at_final, ok, at_final, d


def ctwo(Sb, Ss, Lb, Ls):
    return 6 * (Sb - Ss) / math.log(Lb / Ls)


def extrap(st, chis):
    """R8 extrapolation in the discarded weight through the last two converged stages chis = (lo, hi)."""
    if len(chis) != 2:
        return None
    lo, hi = chis
    if lo not in st or hi not in st:
        return None
    el, eh = st[lo]['trunc'], st[hi]['trunc']
    if not (eh < el):
        return None
    return {W: st[hi]['Sbar%d' % W] + (st[hi]['Sbar%d' % W] - st[lo]['Sbar%d' % W]) * eh / (el - eh) for W in WS}


def main():
    Lb, Ls, key = a.big_L, a.small_L, a.key
    print(f'C195 analysis: key {key}, L = {Ls} -> {Lb}{"  (--partial: running runs at their last mixer-off sweep)" if a.partial else ""}')
    print(f'large size from {os.path.relpath(a.big, HERE)}, small size from ' + ', '.join(os.path.relpath(s, HERE) for s in a.small))
    big = {}
    for P in (0, 1):
        fs = sorted(glob.glob(os.path.join(a.big, f'{key}_L{Lb}_P{P}_*.json')))
        rs = [x for x in (load_run(f, Lb) for f in fs) if x is not None]
        if rs:
            big[P] = min(rs, key=lambda r: r['E'])
    small = {}
    for P in (0, 1):
        cands = []
        for d in a.small:
            for f in sorted(glob.glob(os.path.join(d, f'{key}_L{Ls}_P{P}_*.json'))):
                r = load_run(f, Ls)
                if r is not None and (r['status'] == 'done' or a.partial):
                    cands.append(r)
        if cands:
            done = [r for r in cands if r['status'] == 'done'] or cands
            small[P] = max(done, key=lambda r: (r['chi'], -r['E']))
            small[P]['all'] = cands
    print('\n-- runs')
    for lab, D, L in (('large', big, Lb), ('small', small, Ls)):
        for P, r in sorted(D.items()):
            print(f"  L={L} P={P} {r['tag']} [{os.path.basename(os.path.dirname(r['file']))}] status={r['status']} chi={r['chi']} "
                  f"(final {r['chi_final']}) E={r['E']:.10f} S1,S6,S12 = {r['SW'][1]:.4f}, {r['SW'][6]:.4f}, {r['SW'][12]:.4f} "
                  f"n_mid={r['n_mid']:.4f} c_prof={r['cprof']:.3f}  [{r['src']}]")
            for c, x in sorted(r['st'].items()):
                print(f"      stage chi={c:5d}: E={x['E_phys']:.10f} S6={x['Sbar6']:.4f} trunc={x['trunc']:.1e}"
                      + (f" n_mid={x['n_mid']:.4f}" if x.get('n_mid') is not None else '') + (f" ({x['src']})" if x.get('src') else ''))
    if not big:
        print('\nno usable large-size run yet (T1 to T4 not evaluated)')
        return
    # T1
    print('\n-- T1 (convergence at the final chi), per sector')
    T1 = {}
    for P, r in sorted(big.items()):
        ok, ok_num, at_final, d = t1(r, Lb)
        T1[P] = ok
        print(f"  P={P}: chi={d['chi']} sweeps at chi={d['nsweeps_at_chi']} res={d['res']:.1e} (<=1e-8) dE={d['dE']:.1e} "
              f"(<={d['tolE']:.0e}) dS_W={d['dS']:.1e} (<=1e-3) dn_mid={d['dn']:.1e} (<=5e-4) trunc={d['trunc']:.1e} (<=1e-7) "
              f"-> numbers {'meet' if ok_num else 'miss'} T1, {'at the final chi and done' if at_final else 'NOT at the final chi / not done'}"
              f"; T1 {'holds' if ok else 'not established'}")
        print(f"        C4 step of S6 from chi={d['C4_from_chi']}: {d['C4_step']:.4f} against 0.01 "
              f"-> {'ok' if d['C4_step'] <= 0.01 else 'fails (reported, does not void T2/T3 alone)'}")
    # ground branches
    def ground(D):
        if 0 in D and 1 in D:
            return min((0, 1), key=lambda p: D[p]['E'])
        return next(iter(D)) if D else None
    gb, gs = ground(big), ground(small)
    print('\n-- ground branch')
    if len(big) == 2:
        print(f"  L={Lb}: E(P0)-E(P1) = {big[0]['E'] - big[1]['E']:+.3e} -> P={gb}"
              + ('' if big[0]['chi'] == big[1]['chi'] else f"  (sectors at different chi {big[0]['chi']}, {big[1]['chi']})"))
    else:
        print(f'  L={Lb}: only sector P={gb} available, ground branch provisional')
    if len(small) == 2:
        print(f"  L={Ls}: E(P0)-E(P1) = {small[0]['E'] - small[1]['E']:+.3e} -> P={gs}, chi {small[0]['chi']}/{small[1]['chi']}")
        if small[0]['chi'] != small[1]['chi']:
            cc = min(small[0]['chi'], small[1]['chi'])
            e = {}
            for P in (0, 1):
                for r in small[P]['all']:
                    if cc in r['st']:
                        e[P] = r['st'][cc]['E_phys']
            if len(e) == 2:
                print(f"        at the common chi {cc}: E(P0)-E(P1) = {e[0] - e[1]:+.3e}")
    elif small:
        print(f'  L={Ls}: only sector P={gs} available')
    if gs is None:
        print('no small-size run: c_two not evaluated')
        return
    # T2 / T3
    print(f'\n-- c_two({Ls} -> {Lb}), W = 1, 6, 12')
    rows = [('ground branch', big[gb], small[gs])]
    for P in (0, 1):
        if P in big and P in small:
            rows.append((f'sector P={P}', big[P], small[P]))
    for P in (0, 1):
        for Q in (0, 1):
            if P in big and Q in small:
                rows.append((f'pairing {Ls}:P={Q} -> {Lb}:P={P}', big[P], small[Q]))
    ct = None
    for lab, rb, rs in rows:
        c = {W: ctwo(rb['SW'][W], rs['SW'][W], Lb, Ls) for W in WS}
        if ct is None:
            ct = c
        print(f"  {lab:32s}: {c[1]:.3f} {c[6]:.3f} {c[12]:.3f}   (chi {rs['chi']} -> {rb['chi']})")
    t1g = T1.get(gb, False)
    if t1g and ct[6] >= 1.6 and ct[12] >= 1.6:
        verdict = 'T2 PASS'
    elif t1g and ct[6] <= 1.3 and ct[12] <= 1.3:
        verdict = 'T3 FAIL (falls toward 1)'
    elif not t1g:
        verdict = 'not evaluable (T1 not established for the ground-branch run)'
    else:
        verdict = 'undecided (between T2 and T3)'
    print(f"  -> ground-branch c_two at W = 6, 12: {ct[6]:.3f}, {ct[12]:.3f}; T1 for the ground-branch run: "
          f"{'holds' if t1g else 'not established'}; verdict: {verdict}")
    print(f"  c_prof({Lb}) ground branch: {big[gb]['cprof']:.3f} (sectors: " +
          ', '.join(f"P={P} {r['cprof']:.3f}" for P, r in sorted(big.items())) + ')')
    # matched chi
    print('\n-- c_two at matched chi (last sweep of each chi stage at both sizes), W = 1, 6, 12 [reported, not deciding]')
    allst_s = {P: {} for P in (0, 1)}
    for P in small:
        for r in small[P]['all']:
            for c, x in r['st'].items():
                if x.get('src') is None or c not in allst_s[P]:
                    allst_s[P][c] = x
    for c in sorted(set().union(*[set(big[P]['st']) for P in big])):
        parts = []
        for P in (0, 1):
            if P in big and c in big[P]['st'] and c in allst_s.get(P, {}):
                xb, xs = big[P]['st'][c], allst_s[P][c]
                parts.append(f"P={P}: " + ' '.join(f"{ctwo(xb['Sbar%d' % W], xs['Sbar%d' % W], Lb, Ls):.3f}" for W in WS))
        eb = {P: big[P]['st'][c]['E_phys'] for P in big if c in big[P]['st']}
        es = {P: allst_s[P][c]['E_phys'] for P in allst_s if c in allst_s[P]}
        if len(eb) == 2 and len(es) == 2:
            pb, ps = min(eb, key=eb.get), min(es, key=es.get)
            xb, xs = big[pb]['st'][c], allst_s[ps][c]
            parts.append(f"ground ({ps}->{pb}): " + ' '.join(f"{ctwo(xb['Sbar%d' % W], xs['Sbar%d' % W], Lb, Ls):.3f}" for W in WS))
        if parts:
            print(f'  chi={c}: ' + '; '.join(parts))
    # extrapolation in the discarded weight (R8 formula), last two converged stages
    print('\n-- c_two from S extrapolated in the discarded weight (R8 formula, last two chi stages run to the stop rule) [reported]')
    conv_b = [s['chi'] for s in big[gb]['cfg'].get('conv', [])] + [big[gb]['cfg']['final']['chi']]
    pair_b = tuple(conv_b[-2:])
    for P in sorted(big):
        if P not in small:
            continue
        eb = extrap(big[P]['st'], pair_b)
        cs = sorted(allst_s[P])
        es = extrap(allst_s[P], tuple(cs[-2:])) if len(cs) >= 2 and cs[-2:] != [512, 768] else None
        if eb is None or es is None:
            print(f"  P={P}: not available (large size stages {pair_b} {'present' if eb else 'missing or no decrease of the discarded weight'}, "
                  f"small size needs two converged stages, have {cs[-2:]})")
            continue
        c = {W: ctwo(eb[W], es[W], Lb, Ls) for W in WS}
        print(f"  P={P}: S6(eps->0) {Ls}: {es[6]:.4f}, {Lb}: {eb[6]:.4f}; c_two = {c[1]:.3f} {c[6]:.3f} {c[12]:.3f} "
              f"(an upper estimate, see R8)")
    # T4
    print('\n-- T4 (density)')
    nb, ns = big[gb]['n_mid'], small[gs]['n_mid']
    print(f"  n_mid({Lb}) ground = {nb:.4f}, n_mid({Ls}) ground = {ns:.4f}, difference {nb - ns:+.4f} against 0.006 -> "
          f"{'T4 holds' if abs(nb - ns) <= 0.006 else 'T4 fails'}")
    for P in (0, 1):
        if P in big and P in small:
            print(f"  sector P={P}: {big[P]['n_mid']:.4f} against {small[P]['n_mid']:.4f} ({big[P]['n_mid'] - small[P]['n_mid']:+.4f})")
    st = 'all runs done' if all(r['status'] == 'done' for r in big.values()) and len(big) == 2 else 'PROVISIONAL: large-size runs not all done'
    print(f'\n{st}')


if __name__ == '__main__':
    main()
