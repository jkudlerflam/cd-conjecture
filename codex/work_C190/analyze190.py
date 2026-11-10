#!/usr/bin/env python3
"""C190 analysis (Claude Opus 5.5 agent, 2026-11-09).  Usage: python3 analyze190.py <outdir> [--partial]

Evaluates the registered estimators, convergence targets C1 to C6 and verdict criteria V1 to V6 of RESULT.md.
With --partial, runs still in progress are included at their last completed sweep (flagged 'running')."""
import sys
import os
import re
import json
import glob
import numpy as np

T2 = {'p153': 0.0125, 'rat': 0.0125, 'sym40': 0.0125, 'sym80': 0.00625}
TARGET = {'n316': 0.1957, 'n18': 0.1326}
REF_FAMS = ('sym18', 'sym532', 'sym16', 'sym316', 'sym14')
HERE = os.path.dirname(os.path.abspath(__file__))


def sbar(S, L, W):
    S = np.asarray(S)
    l = np.arange(1, L)
    sel = (l >= L // 2 - W) & (l <= L // 2 + W)
    return float(S[sel].mean())


def c_prof(S, L, sm=12, lo=1 / 6, hi=5 / 6):
    S = np.asarray(S)
    l = np.arange(1, L)
    k = np.ones(sm) / sm
    Ss = np.convolve(S, k, mode='valid')
    ls = np.convolve(l.astype(float), k, mode='valid')
    x = np.log((2 * L / np.pi) * np.sin(np.pi * ls / L))
    sel = (ls >= lo * L) & (ls <= hi * L)
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    return 6 * np.linalg.lstsq(A, Ss[sel], rcond=None)[0][0]


def valley_weights(meas):
    G = np.array(meas['Gf_re']) + 1j * np.array(meas['Gf_im'])
    W = G.shape[0]
    ks = np.linspace(-np.pi, np.pi, 256, endpoint=False)
    x = np.arange(W)
    ph = np.exp(1j * np.outer(ks, x))
    nk = np.real(np.einsum('ki,ij,kj->k', ph, G, ph.conj())) / W
    return float(nk[np.cos(ks) > 0].sum()), float(nk[np.cos(ks) < 0].sum())


def ff_entropies(L, mu, legs):
    """Exact entanglement at every cut of the free-fermion control (open chains, hopping -1, chemical pot. mu)."""
    def chain(n):
        h = -np.eye(n, k=1) - np.eye(n, k=-1) - mu * np.eye(n)
        e, U = np.linalg.eigh(h)
        occ = U[:, e < 0]
        return occ @ occ.conj().T

    def S_block(C, l):
        if l == 0 or l == C.shape[0]:
            return 0.0
        nu = np.linalg.eigvalsh(C[:l, :l])
        nu = nu[(nu > 1e-14) & (nu < 1 - 1e-14)]
        return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))
    if legs == 1:
        C = chain(L)
        return np.array([S_block(C, l) for l in range(1, L)])
    C = chain(L // 2)
    return np.array([S_block(C, (l + 1) // 2) + S_block(C, l // 2) for l in range(1, L)])


def load(outdir, partial):
    runs = []
    for f in sorted(glob.glob(os.path.join(outdir, '*.json'))):
        tag = os.path.basename(f)[:-5]
        m = re.match(r'^(.*)_L(\d+)_P(\d)_(.*)$', tag)
        if not m or tag.startswith('tune_'):
            continue
        try:
            r = json.load(open(f))
        except Exception:
            continue
        if r.get('status') != 'done' and not partial:
            continue
        key, L, P, init = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        h = r['hist']
        hs = [x for x in h if 'Sbar6' in x]
        if not hs:
            continue
        if r.get('status') == 'done':
            S = np.array(r['S'])
            dens = np.array(r['meas']['dens'])
            E = r.get('E_phys', h[-1]['E_phys'])
        else:
            S = np.array(hs[-1]['Sprof'])
            dens = (1 - np.array(hs[-1]['Zprof'])) / 2
            E = hs[-1]['E_phys']
        fin = [x for x in h if x.get('final')]
        chi_f = r['cfg']['final']['chi']
        run = dict(tag=tag, key=key, L=L, P=P, init=init, status=r.get('status'), E=E, S=S,
                   S1=sbar(S, L, 1), S6=sbar(S, L, 6), S12=sbar(S, L, 12), cprof=c_prof(S, L),
                   n_mid=float(np.mean(dens[L // 4: 3 * L // 4])), chi_last=h[-1]['chi'], chi_f=chi_f,
                   nsweeps=len(h), var=r.get('variance'), cfg=r['cfg'], last=h[-1])
        last = h[-1]
        run['C1'] = last['res_max'] <= 1e-8
        run['res_max'] = last['res_max']
        run['C3'] = last['trunc'] <= 1e-7
        run['trunc'] = last['trunc']
        if len(fin) >= 2:
            a, b = fin[-2], fin[-1]
            dE = abs(b['E'] - a['E'])
            dS = max(abs(b[k] - a[k]) for k in ('Sbar1', 'Sbar6', 'Sbar12'))
            dn = abs(b['n_mid'] - a['n_mid'])
            run['C2'] = dE <= 1e-8 and dS <= 1e-3 and dn <= 5e-4
            run['dE'], run['dS'], run['dn'] = dE, dS, dn
        else:
            run['C2'] = False
            run['dE'] = run['dS'] = run['dn'] = None
        prev = [x for x in hs if x['chi'] < chi_f]
        if fin and prev and 'Sbar6' in fin[-1]:
            run['dS_chi'] = abs(fin[-1]['Sbar6'] - prev[-1]['Sbar6'])
            run['chi_prev'] = prev[-1]['chi']
            run['C4'] = run['dS_chi'] <= 0.01
        else:
            run['dS_chi'] = None
            run['C4'] = False
        run['conv'] = run['C1'] and run['C2'] and run['C3'] and run['C4'] and r.get('status') == 'done'
        if r.get('status') == 'done':
            try:
                run['Wa'], run['Wb'] = valley_weights(r['meas'])
            except Exception:
                pass
            run['parity_meas'] = r.get('parity')
        # chi history: last sweep of each chi stage with alpha = 0
        stages = {}
        for x in hs:
            stages[x['chi']] = (x['E_phys'], x['Sbar6'], x['n_mid'])
        run['stages'] = stages
        runs.append(run)
    return runs


def ref_interp(L, n, W='Sbar6'):
    ref = json.load(open(os.path.join(HERE, 'c176_ref.json')))
    pts = []
    for f in REF_FAMS:
        e = ref.get(f, {}).get(str(L))
        if e and 'n_mid' in e:
            pts.append((e['n_mid'], e[W]))
    pts.sort()
    if len(pts) < 2:
        return None
    xs = np.log([p[0] for p in pts])
    ys = [p[1] for p in pts]
    if not xs[0] <= np.log(n) <= xs[-1]:
        return None
    return float(np.interp(np.log(n), xs, ys))


def main():
    outdir = sys.argv[1]
    partial = '--partial' in sys.argv
    runs = load(outdir, partial)
    keys = sorted(set(r['key'] for r in runs))
    print('C190 analysis' + (' (partial: running jobs at their last sweep)' if partial else ''))
    print()
    print('== per run ==')
    for r in sorted(runs, key=lambda r: (r['key'], r['L'], r['P'], r['init'])):
        fl = ''.join(c if r.get(c) else '-' for c in ('C1', 'C2', 'C3', 'C4'))
        print(f"{r['tag']:28s} {r['status']:7s} sw={r['nsweeps']:2d} chi={r['chi_last']:3d} E={r['E']:.10f} "
              f"S1={r['S1']:.4f} S6={r['S6']:.4f} S12={r['S12']:.4f} cprof={r['cprof']:.3f} n_mid={r['n_mid']:.4f} "
              f"res={r['res_max']:.1e} tr={r['trunc']:.1e} dE={r['dE'] if r['dE'] is None else format(r['dE'], '.1e')} "
              f"dS={r['dS'] if r['dS'] is None else format(r['dS'], '.1e')} dSchi={r['dS_chi'] if r['dS_chi'] is None else format(r['dS_chi'], '.4f')} "
              f"[{fl}]" + (f" Wa/Wb={min(r['Wa'], r['Wb']) / max(r['Wa'], r['Wb']):.3f}" if 'Wa' in r else '') +
              (f" var={r['var']:.2e}" if r.get('var') is not None else ''))
        print('      stages (chi: E, S6, n_mid): ' + '; '.join(f"{c}: {v[0]:.9f}, {v[1]:.4f}, {v[2]:.4f}"
                                                            for c, v in sorted(r['stages'].items())))
    print()
    for key in keys:
        rs = [r for r in runs if r['key'] == key]
        print(f'== {key} ==')
        best = {}
        for L in sorted(set(r['L'] for r in rs)):
            for P in (0, 1):
                cand = [r for r in rs if r['L'] == L and r['P'] == P]
                if not cand:
                    continue
                cand.sort(key=lambda r: r['E'])
                best[(L, P)] = cand[0]
                if len(cand) > 1:
                    a, b = cand[0], cand[1]
                    agree = abs(a['E'] - b['E']) <= 1e-7
                    print(f"  C5 L={L} P={P}: inits {a['init']} (E={a['E']:.10f}, S6={a['S6']:.4f}) vs {b['init']} "
                          f"(E={b['E']:.10f}, S6={b['S6']:.4f}): dE={b['E'] - a['E']:.2e}, dS6={b['S6'] - a['S6']:+.4f} "
                          + ('agree in E; S6 within 0.01: ' + str(abs(a['S6'] - b['S6']) <= 0.01) if agree else
                             'higher one metastable'))
        Ls = sorted(set(L for (L, P) in best))
        gb = {}
        for L in Ls:
            if (L, 0) in best and (L, 1) in best:
                e0, e1 = best[(L, 0)]['E'], best[(L, 1)]['E']
                gb[L] = best[(L, 0)] if e0 <= e1 else best[(L, 1)]
                bond = best[(L, 0)]['cfg'].get('model')
                t2 = T2.get(bond)
                txt = f"L*dGC/|t2| = {L * (e1 - e0) / t2:+.3f}" if t2 else f"E_odd - E_even = {e1 - e0:+.3e}"
                print(f"  L={L}: E_even={e0:.10f} E_odd={e1:.10f} {txt}  ground P={gb[L]['P']} "
                      f"n_mid={gb[L]['n_mid']:.4f}")
            else:
                gb[L] = best.get((L, 0)) or best.get((L, 1))
        for (a, b) in ((64, 128), (128, 256), (64, 256)):
            if a in gb and b in gb:
                cs = [6 * (gb[b][f'S{W}'] - gb[a][f'S{W}']) / np.log(b / a) for W in (1, 6, 12)]
                conv = gb[a]['conv'] and gb[b]['conv']
                print(f"  c_two {a}->{b} (ground branch P={gb[a]['P']},{gb[b]['P']}) W=1,6,12: "
                      + ', '.join(f'{c:.3f}' for c in cs) + ('' if conv else '   [not all C1-C4]'))
                allp = []
                for Pa in (0, 1):
                    for Pb in (0, 1):
                        if (a, Pa) in best and (b, Pb) in best:
                            allp.append(6 * (best[(b, Pb)]['S1'] - best[(a, Pa)]['S1']) / np.log(b / a))
                if allp:
                    print(f"      all parity pairings W=1: {min(allp):.3f} to {max(allp):.3f}")
        # density-corrected secondary estimator (RESULT.md): interpolate S_W between the two parity sectors
        if Ls and all((L, 0) in best and (L, 1) in best for L in Ls):
            nstar = gb[max(Ls)]['n_mid']
            sdc = {}
            for L in Ls:
                a0, a1 = best[(L, 0)], best[(L, 1)]
                if abs(a1['n_mid'] - a0['n_mid']) < 1e-6:
                    continue
                sdc[L] = {W: a0[f'S{W}'] + (a1[f'S{W}'] - a0[f'S{W}']) * (nstar - a0['n_mid']) /
                          (a1['n_mid'] - a0['n_mid']) for W in (1, 6, 12)}
            print(f"  density-corrected S_W at n*={nstar:.4f}: " + '; '.join(
                f"L={L}: " + ', '.join(f'{sdc[L][W]:.4f}' for W in (1, 6, 12)) for L in sorted(sdc)))
            for (a, b) in ((64, 128), (128, 256), (64, 256)):
                if a in sdc and b in sdc:
                    print(f"  c_dc {a}->{b} W=1,6,12: " + ', '.join(
                        f"{6 * (sdc[b][W] - sdc[a][W]) / np.log(b / a):.3f}" for W in (1, 6, 12)))
        for L in Ls:
            r = gb[L]
            ref = ref_interp(L, r['n_mid'])
            extra = f" S6_eff(n_mid)={ref:.4f} dS6={r['S6'] - ref:+.4f}" if ref is not None else ''
            print(f"  L={L}: cprof={r['cprof']:.3f} S6={r['S6']:.4f}{extra}"
                  + (f" Wa/Wb={min(r['Wa'], r['Wb']) / max(r['Wa'], r['Wb']):.3f}" if 'Wa' in r else ''))
        if key.startswith('xx'):
            legs = int(key[2])
            for L in Ls:
                r = gb[L]
                Se = ff_entropies(L, r['cfg']['mu'], legs)
                print(f"  exact L={L}: S1={sbar(Se, L, 1):.4f} S6={sbar(Se, L, 6):.4f} S12={sbar(Se, L, 12):.4f}; "
                      f"DMRG - exact: {r['S1'] - sbar(Se, L, 1):+.2e}, {r['S6'] - sbar(Se, L, 6):+.2e}, "
                      f"{r['S12'] - sbar(Se, L, 12):+.2e}")
            for (a, b) in ((64, 128), (128, 256)):
                if a in gb and b in gb:
                    Sa = ff_entropies(a, gb[a]['cfg']['mu'], legs)
                    Sb = ff_entropies(b, gb[b]['cfg']['mu'], legs)
                    print(f"  exact c_two {a}->{b} W=1,6,12: " + ', '.join(
                        f"{6 * (sbar(Sb, b, W) - sbar(Sa, a, W)) / np.log(b / a):.3f}" for W in (1, 6, 12)))
        for tk, tv in TARGET.items():
            if key.endswith(tk) and 256 in gb:
                print(f"  C6 density at L=256: n_mid={gb[256]['n_mid']:.4f} target {tv} -> "
                      f"{abs(gb[256]['n_mid'] - tv) <= 0.006}")
        if 64 in gb and 256 in gb:
            r64, r256 = gb[64], gb[256]
            a, b = ref_interp(64, r64['n_mid']), ref_interp(256, r256['n_mid'])
            if a is not None and b is not None:
                d = (r256['S6'] - b) - (r64['S6'] - a)
                print(f"  V4 dS6(256) - dS6(64) = {d:+.4f} (faithful if within +-0.05)")
        print()
    tf = sorted(glob.glob(os.path.join(outdir, 'tune', '*.json')))
    if tf:
        print('== field tuning ==')
        for f in tf:
            t = json.load(open(f))
            print(f"  {t['key']}: points (mu_eff, n_mid, hz) = " +
                  '; '.join(f"({p[0]:.3f}, {p[1]:.4f}, {p[2]:.9f})" for p in t['points']) +
                  f" -> mu_eff={t['mu_eff']:.4f}, hz={t['hz']:.9f}")


if __name__ == '__main__':
    main()
