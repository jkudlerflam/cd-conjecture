#!/usr/bin/env python3
"""Tables for README.md from data/: ground-state degeneracies, spin-charge separation errors, velocities,
twist levels, Casimir estimates, FES fits, d = 24 rings.  Writes data/summary.json.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, glob, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sqed  # noqa: E402

D = os.path.join(HERE, 'data')
out = {}


def load(p):
    try:
        return json.load(open(p))
    except Exception:  # noqa: BLE001
        return None


# ---------------------------------------------------------------- ground-state scans
def scans(pattern):
    rows = []
    for f in glob.glob(os.path.join(D, 'scan', pattern)):
        d = load(f)
        if d and 'summary' not in d and d.get('sectors'):
            # partial run: report the lowest computed level and its multiplicity, flagged
            secs = d['sectors']
            E0 = min(min(x['levels']) for x in secs)
            mult = sum((2 if x['Q2'] else 1) for x in secs for e in x['levels'] if e < E0 + 1e-8 * max(1, abs(E0)))
            g = min(secs, key=lambda x: min(x['levels']))
            rows.append(dict(L=d['args']['L'], deg=mult, M=g['M'], word=g['word'], Q2=g['Q2'], m=g['m'], gap=None,
                             E0=E0, sat=0, n_sectors=len(secs), partial=sorted(set((x['M'], x['word'], x['Q2']) for x in secs))))
            continue
        if not d or 'summary' not in d:
            continue
        s = d['summary']
        L = d['args']['L']
        g = s['gs_sectors'][0]
        rows.append(dict(L=L, deg=s['degeneracy'], M=g['M'], word=g['word'], Q2=g['Q2'], m=g['m'],
                         gap=s['gaps'][0] if s['gaps'] else None, E0=s['E0'],
                         sat=len(s.get('saturated_sectors', [])), n_sectors=len(d['sectors'])))
    return sorted(rows, key=lambda r: r['L'])


for name, pat in [('vac_J0.25', 'vac_J0.25_mu0.6_L*.json'), ('frozen', 'frozen_L*.json'),
                  ('vac_J1.0', 'vac_J1.0_mu0.34_L*.json'), ('vac_J2.0', 'vac_J2.0_mu0.07_L*.json'),
                  ('vac_J4.0', 'vac_J4.0_mu-0.48_L*.json')]:
    rows = scans(pat)
    out['scan_' + name] = rows
    print(f'\n== {name}: L, degeneracy, GS (M, word, 2Q, m), first gap')
    for r in rows:
        print(f"  L={r['L']:2d} deg={r['deg']:3d}  M={r['M']:2d} {r['word']:<14s} 2Q={r['Q2']} m={r['m']} "
              f"gap={r['gap'] if r['gap'] is None else round(r['gap'], 6)} sat={r['sat']}"
              + (f" PARTIAL, sectors done: {r['partial']}" if 'partial' in r else ''))

# ---------------------------------------------------------------- spin-charge separation
print('\n== first-order separation: max |E_ED - E_pred| and max gap error over the lowest levels per sector')
sep = []
for f in sorted(glob.glob(os.path.join(D, 'sep', '*.json'))):
    d = load(f)
    if not d:
        continue
    L = d['args']['L']; word = d['args']['word']
    for run in d['runs']:
        J = run['J']
        ed0 = min(min(s['ed']) for s in run['sectors'].values())
        pr0 = min(np.nanmin(s['pred']) for s in run['sectors'].values())
        errs, gerrs = [], []
        for m, s in run['sectors'].items():
            if 2 * int(m) == L:      # P = pi: degenerate charge states mix at first order (not in the prediction)
                continue
            ed = np.array(s['ed']); pr = np.array(s['pred'])
            ok = np.isfinite(pr)
            errs.append(np.abs(ed[ok] - pr[ok]).max())
            gerrs.append(np.abs((ed[ok] - ed0) - (pr[ok] - pr0)).max())
        # charge velocity: lowest level at m = 1 labelled as a charge particle-hole state (r = 0, i = 0)
        s1 = run['sectors'].get('1')
        vc = vc_pred = None
        if s1:
            for e, p, lab in zip(s1['ed'], s1['pred'], s1['labels']):
                if lab['r'] == 0 and lab['i'] == 0 and lab['charge_exc'] > 1e-9:
                    vc = L * (e - ed0) / (2 * np.pi); vc_pred = L * (p - pr0) / (2 * np.pi); break
        # pattern-twist level for a period-3 word: internal sector r = M/3, charge flux 2 pi/3
        M = len(word)
        tw = None
        if sqed.period(word) == 3 and M % 3 == 0:
            for m, s in run['sectors'].items():
                for e, lab in zip(s['ed'], s['labels']):
                    if lab['r'] == M // 3 and lab['i'] == 0 and lab['charge_exc'] < 1e-9:
                        tw = (int(m), e - ed0); break
                if tw:
                    break
        rec = dict(file=os.path.basename(f), L=L, word=word, J=J, max_abs_err=float(max(errs)),
                   max_gap_err=float(max(gerrs)), v_c=vc, v_c_pred=vc_pred,
                   twist=None if tw is None else dict(m=tw[0], gap=tw[1],
                                                      x=(tw[1] * L / (2 * np.pi * vc)) if vc else None))
        sep.append(rec)
        print(f"  L={L} {word} J={J}: |dE|max={rec['max_abs_err']:.2e} gap err={rec['max_gap_err']:.2e} "
              f"v_c={vc} (first order {vc_pred}) twist={rec['twist']}")
out['separation'] = sep

# ---------------------------------------------------------------- velocities from structure factors
print('\n== modes: charge velocity (density probe), q-probe peak, twist level x = L gap / (2 pi v_c)')
modes = []
for f in glob.glob(os.path.join(D, 'modes', 'J*_L*.json')):
    d = load(f)
    if not d or 'n' not in d.get('probes', {}):
        continue
    a = d['args']; L = a['L']; M = len(a['word'])
    vc = d['probes']['n']['dominant']['v']
    wfrac = d['probes']['n']['dominant']['weight'] / sum(p['weight'] for p in d['probes']['n']['peaks'])
    # free hard-core bosons at the same (L, M): particle-hole velocity and twist gap
    n = np.arange(L); k = (2 * np.pi * n + np.pi * (M - 1)) / L; e = -2 * np.cos(k); occ = set(np.argsort(e)[:M])
    vfree = L * min(e[(i + dd) % L] - e[i] for i in occ for dd in (1, -1) if (i + dd) % L not in occ) / (2 * np.pi)
    tw = d.get('twist')
    twf = sqed.lowest_free(L, M, 2 * np.pi / 3) - sqed.lowest_free(L, M, 0.0)
    rec = dict(L=L, M=M, rho=M / L, J=a['J'], E0=d['E0'], v_c=vc, v_free=vfree, weight=wfrac,
               q_peak=d['probes'].get('q', {}).get('dominant'),
               twist_gap=None if tw is None else tw['gap'], twist_gap_free=twf,
               x_twist=None if tw is None else tw['gap'] * L / (2 * np.pi * vc))
    modes.append(rec)
modes.sort(key=lambda r: (r['J'], r['M'], r['L']))
for r in modes:
    print(f"  J={r['J']} L={r['L']:2d} M={r['M']} rho={r['rho']:.3f} v_c={r['v_c']:.5f} (free {r['v_free']:.5f}, weight {r['weight']:.3f}) "
          f"twist gap={r['twist_gap']} (free {r['twist_gap_free']:.6f}) x={r['x_twist']}")
out['modes'] = modes

# ---------------------------------------------------------------- FES
print('\n== FES: c from S = (c/6) ln xi + s0, fits over the last 4 chi; pairwise slopes')
fes = {}
for f in sorted(glob.glob(os.path.join(D, 'fes', '*.json'))):
    d = load(f)
    if not d or not d['results']:
        continue
    R = d['results']
    chis = [r['chi'] for r in R]
    S = np.array([np.mean(r['S']) for r in R]); xi = np.array([r['xi'] for r in R])
    pw = [6 * (S[i + 1] - S[i]) / np.log(xi[i + 1] / xi[i]) for i in range(len(S) - 1)]
    fit = None
    if len(S) >= 3:
        k = min(4, len(S))
        fit = float(6 * np.polyfit(np.log(xi[-k:]), S[-k:], 1)[0])
    loc = R[-1].get('local')
    fes[d['name']] = dict(chis=chis, S=S.tolist(), xi=xi.tolist(), pairwise=pw, c_fit_last4=fit,
                          e=[r['e'] for r in R], err=[r['err'] for r in R], local_last=loc,
                          adjacent_last=R[-1].get('adjacent_occupied'))
    print(f"  {d['name']}: chi={chis} xi={np.round(xi, 2).tolist()} S={np.round(S, 4).tolist()}\n"
          f"      pairwise={np.round(pw, 3).tolist()} fit={fit} local={loc}")
out['fes'] = fes

# ---------------------------------------------------------------- d = 24
print('\n== d = 24 (23-letter bond, n = 1, plus vacancy)')
d24 = []
for f in sorted(glob.glob(os.path.join(D, 'd24', '*_L*.json'))):
    d = load(f)
    if not d:
        continue
    r = d['result']; s = r['summary']; a = r['args']
    g = s['gs'][0]
    rec = dict(file=os.path.basename(f), L=a['L'], mode=a['mode'], deg=s['degeneracy'], M=g['M'], word=g['word'],
               Q=g['Q'], ms=sorted(set(x['m'] for x in s['gs'])), gaps=s['gaps'][:4])
    d24.append(rec)
    print(f"  {rec['file']}: deg={rec['deg']} GS M={rec['M']} {rec['word']} Q={rec['Q']} m={rec['ms']} gaps={np.round(rec['gaps'], 7).tolist()}")
out['d24'] = d24

json.dump(out, open(os.path.join(D, 'summary.json'), 'w'), indent=1, default=float)
