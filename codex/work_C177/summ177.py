#!/usr/bin/env python3
"""C177 summary (Claude Opus 5.5 agent, 2026-11-09).  usage: summ177.py DIR  (DIR contains oc/ and vc/)
(1) open chains: per tag and L, Sbar (three middle cuts), S at the cuts L/2 - 3, L/2, L/2 + 3 (the block cuts of C173 for
    L divisible by 3 are reported separately), c_prof, c_sm, last-sweep truncation, energy history, n_mid, emid,
    the valley-polarization profile max |P_j| in the middle half and the correlator <P_i0 P_j> at |j - i0| = L/4, L/2.
(2) two-size slopes 6 dSbar / dln L for every pair of lengths of a tag.
(3) VUMPS FES (fesfit_c158.analyse) and the order parameters.
(4) C177 pass criterion P1 to P3 (RESULT.md) for every tag with L = 128 and 256 and a VUMPS file."""
import sys, os, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fesfit_c158 as FF

D = sys.argv[1] if len(sys.argv) > 1 else '.'
oc = {}
for f in sorted(glob.glob(os.path.join(D, 'oc', '*.json'))):
    d = json.load(open(f))
    oc.setdefault(d['tag'], {})[d['L']] = d
print('== (1) open chains')
for tag in sorted(oc):
    for L in sorted(oc[tag]):
        d = oc[tag][L]
        S = np.array(d['S']); c = L // 2
        P = np.array(d.get('P', [np.nan])); PP = np.array(d.get('PP', [np.nan]))
        mid = slice(L // 4, 3 * L // 4)
        pmax = float(np.max(np.abs(P[mid]))) if P.size > 1 else np.nan
        def ppat(r):
            return float(PP[r - 2]) if PP.size > r - 2 >= 0 else np.nan
        hist = d['hist']
        dE = hist[-1][1] - hist[-2][1] if len(hist) > 1 else np.nan
        print(f"{tag:8s} L={L:4d} Sbar={d['Sbar']:.4f} S(L/2)={S[c - 1]:.4f} c_prof={d['c_prof']:.3f} c_sm={d['c_sm']:.3f} "
              f"trunc={d['trunc']:.1e} E={d['E']:.10f} dE_last={dE:.1e} n_mid={d['n_mid']:.4f} emid={d['emid']:.10f} "
              f"|P|max={pmax:.1e} PP(L/4)={ppat(L // 4):+.2e} PP(L/2-2)={ppat(L // 2 - 2):+.2e} chi={d['chis'][-1]}")
print('== (2) two-size slopes 6 dSbar/dlnL')
for tag in sorted(oc):
    Ls = sorted(oc[tag])
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            a, b = oc[tag][Ls[i]], oc[tag][Ls[j]]
            print(f"{tag:8s} {Ls[i]:4d} -> {Ls[j]:4d}: c_two = {6 * (b['Sbar'] - a['Sbar']) / np.log(Ls[j] / Ls[i]):.3f}")
print('== (3) VUMPS')
vc = {}
for f in sorted(glob.glob(os.path.join(D, 'vc', '*.json'))):
    try:
        res = FF.analyse(f, 3)
    except Exception as ex:
        print('fail', f, ex); continue
    d = json.load(open(f))
    tag = d['args']['tag']
    R = d['results']
    last = R[-1]
    print(f"   order parameters at chi={last['chi']}: n={np.round(last['n'], 5).tolist()} P={np.round(last['xxyy'], 6).tolist()} "
          f"e={last['e']:.12f}")
    vc.setdefault(tag, []).append(dict(file=f, n=d['args']['n'], e=last['e'], chi=last['chi'], res=res,
                                       echi=[r['e'] for r in R], chis=[r['chi'] for r in R]))
print('== (4) pass criterion P1-P3')
for tag in sorted(oc):
    if 128 not in oc[tag] or 256 not in oc[tag]:
        continue
    a, b = oc[tag][128], oc[tag][256]
    ctwo = 6 * (b['Sbar'] - a['Sbar']) / np.log(2)
    p1 = (ctwo > 1.2 and a['c_prof'] > 1.2 and b['c_prof'] > 1.2 and b['c_prof'] >= a['c_prof'] - 0.03
          and a['trunc'] < 1e-8 and b['trunc'] < 1e-8)
    line = (f"{tag}: P1 {'pass' if p1 else 'FAIL'} (c_two 128->256 = {ctwo:.3f}, c_prof {a['c_prof']:.3f} -> {b['c_prof']:.3f}, "
            f"trunc {a['trunc']:.1e}, {b['trunc']:.1e})")
    if tag in vc:
        best = min(vc[tag], key=lambda x: x['e'])
        r = best['res']
        fes = r.get('c_ls3', np.nan)
        p2 = best['chi'] >= 256 and fes > 1.2
        e0 = best['e']
        echi = best['echi']
        err_e0 = abs(echi[-1] - echi[len(echi) // 2]) if len(echi) > 2 else np.nan
        allow = 2 * abs(b['emid'] - a['emid']) + 1e-6
        de = b['emid'] - e0
        p3 = abs(de) <= allow
        line += (f"; P2 {'pass' if p2 else 'FAIL'} (n={best['n']}, chi_max={best['chi']}, FES c_ls3={fes:.3f})"
                 f"; P3 {'pass' if p3 else 'FAIL'} (emid(256) - e0 = {de:+.2e}, allowance {allow:.1e}, chi-error of e0 ~ {err_e0:.1e})")
    else:
        line += '; P2, P3 n/a (no VUMPS file)'
    print(line)
