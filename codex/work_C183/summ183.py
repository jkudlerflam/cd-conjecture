#!/usr/bin/env python3
"""C183 summary of the Typhon output (Claude Opus 5.5 agent, 2026-11-09), written before any output existed.
usage: summ183.py DIR   (DIR holds screen_summary.json, conf/, vc/, ctrl/ as written by pipeline183.py)
(1) the stage-S table (screen_summary.txt is reprinted, with the controls);
(2) stage C per selected tag: S-bar, c_prof, truncation at L = 128 and 256, c2(128 -> 256), e_mid;
    VUMPS FES slope (fesfit_c158.analyse, c_ls3) and e0; the parity gap L * (E_odd - E_even) at L = 64, 128;
(3) the registered criterion P1 to P4 and the control calibration (RESULT.md, stage C)."""
import sys, os, json, glob, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fesfit_c158 as FF

D = sys.argv[1] if len(sys.argv) > 1 else '.'


def load(fn):
    try:
        return json.load(open(fn))
    except Exception:
        return None


print('== (1) screen')
t = os.path.join(D, 'screen_summary.txt')
if os.path.exists(t):
    print(open(t).read())
SS = load(os.path.join(D, 'screen_summary.json')) or {}
sel = SS.get('selected', [])
rows = {r['tag']: r for r in SS.get('rows', [])}
print('== (2, 3) confirmation')
ctl = {}
for kind in ('xx', 'leg2'):
    a, b = load(os.path.join(D, 'ctrl', f'{kind}_128.json')), load(os.path.join(D, 'ctrl', f'{kind}_256.json'))
    if a and b:
        ctl[kind] = 6 * (b['Sbar'] - a['Sbar']) / math.log(2)
        print(f"control {kind}: c2(128->256) = {ctl[kind]:.3f}, c_prof {a['c_prof']:.3f} -> {b['c_prof']:.3f}, trunc {a['trunc']:.1e}, {b['trunc']:.1e}")
calib = ('xx' in ctl and abs(ctl['xx'] - 1) <= 0.2 and 'leg2' in ctl and abs(ctl['leg2'] - 2) <= 0.3)
print('controls calibrated' if calib else 'controls NOT calibrated (or missing)')
for tag in sel:
    a, b = load(os.path.join(D, 'conf', f'{tag}_128.json')), load(os.path.join(D, 'conf', f'{tag}_256.json'))
    line = f"{tag} [{rows.get(tag, {}).get('group', '?')}]: "
    p1 = None
    if a and b:
        c2 = 6 * (b['Sbar'] - a['Sbar']) / math.log(2)
        p1 = (c2 > 1.2 and a['c_prof'] > 1.2 and b['c_prof'] > 1.2 and b['c_prof'] >= a['c_prof'] - 0.03
              and a['trunc'] < 1e-7 and b['trunc'] < 1e-7)
        line += (f"Sbar {a['Sbar']:.4f} -> {b['Sbar']:.4f}, c2 = {c2:.3f}, c_prof {a['c_prof']:.3f} -> {b['c_prof']:.3f}, "
                 f"trunc {a['trunc']:.1e}, {b['trunc']:.1e}, n_mid {b['n_mid']:.4f}; P1 {'pass' if p1 else 'FAIL'}")
    else:
        line += 'open chains missing; P1 n/a'
    vf = os.path.join(D, 'vc', f'{tag}_n2.json')
    v = load(vf)
    if v and v.get('results'):
        try:
            res = FF.analyse(vf, 3)
            fes = res.get('c_ls3', np.nan)
        except Exception as ex:
            fes = np.nan
        R = v['results']; e0 = R[-1]['e']; chim = R[-1]['chi']
        p2 = chim >= 256 and fes > 1.2
        line += f"; P2 {'pass' if p2 else 'FAIL'} (chi_max {chim}, FES {fes:.3f})"
        if a and b:
            allow = 2 * abs(b['emid'] - a['emid']) + 1e-6
            de = b['emid'] - e0
            line += f"; P3 {'pass' if abs(de) <= allow else 'FAIL'} (emid - e0 = {de:+.2e}, allowance {allow:.1e})"
    else:
        line += '; VUMPS missing'
    gaps = {}
    for L in (64, 128):
        o, e = load(os.path.join(D, 'conf', f'{tag}_odd_{L}.json')), load(os.path.join(D, 'conf', f'{tag}_even_{L}.json'))
        if o and e:
            gaps[L] = L * (o['E'] - e['E'])
    if gaps:
        line += '; L*Dp ' + ', '.join(f'{L}: {g:+.4f}' for L, g in sorted(gaps.items()))
        if 64 in gaps and 128 in gaps:
            p4 = abs(gaps[128]) < 1.5 * abs(gaps[64])
            line += f"; P4 {'pass' if p4 else 'FAIL'}"
    print(line)
