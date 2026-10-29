#!/usr/bin/env python3
"""Ring-spectrum diagnostics for the floating-Ising search (reads data/ring/*.jsonl).

Per point and L, with s = L (E - E0) / (2 pi) and q the momentum relative to the ground state:
  sJ   lowest level at q = +-2 pi/L. In a one-boson liquid this is the U(1) current (x = 1), so sJ ~ v.
  s0   first excited level at q = 0.  One boson: x0 = s0/sJ -> 2 (J Jbar).  An extra neutral Ising sector
       puts its energy field (x = 1, always at q = 0) at x0 = v_I/v.  Ising (x) XX gives x0 = 1/2.
  spi  lowest level at q = pi;  xpi = spi/sJ (one boson: K at the staggered vertex, Ising sigma adds v_I/(8v)).
  smin lowest gap over all momenta;  gap ratio rg = smin(20)/smin(16) (> 1.2: gapped at these sizes).
Gapless (liquid-like): cS(20) > 0.7 and sJ(20)/sJ(16) < 1.15.  Ising-like: 0.3 < cS(20) < 0.7, same sJ test.
Flagged (extra neutral level): gapless, x0 < 1.3 at L = 16, 18, 20, and x0(20) <= x0(16) + 0.05.
Stable flag: flagged and max - min of x0 over L = 16..20 below 0.15 (a moving zero mode is not stable).
usage: analyze_ring.py [--json data/ring_summary.json] [--maps]
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, glob, argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models_fi as MF


def diag(r, L):
    E = np.array(r[f'L{L}']['E']); m0 = r[f'L{L}']['m0']; e0 = E[m0, 0]
    s = L * (E - e0) / (2 * np.pi)
    sJ = min(s[(m0 + 1) % L, 0], s[(m0 - 1) % L, 0])
    s0 = s[m0, 1]
    spi = s[(m0 + L // 2) % L, 0]
    allv = np.sort(s.ravel())
    smin = allv[1]
    return dict(m0=m0, sJ=float(sJ), s0=float(s0), spi=float(spi), smin=float(smin),
                x0=float(s0 / sJ) if sJ > 1e-9 else np.inf, xpi=float(spi / sJ) if sJ > 1e-9 else np.inf,
                cS=r[f'L{L}'].get('cS'), low=[float(x) for x in allv[1:9]])


def load():
    recs = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'ring', '*.jsonl'))):
        for line in open(f):
            try:
                r = json.loads(line)
            except Exception:
                continue
            recs[r['tag']] = r
    return recs


def summarize(recs, sizes=(16, 18, 20)):
    out = {}
    for t, r in recs.items():
        if not all(f'L{L}' in r for L in sizes):
            continue
        d = {L: diag(r, L) for L in sizes}
        rg = d[20]['smin'] / d[16]['smin'] if d[16]['smin'] > 1e-12 else np.inf
        rJ = d[20]['sJ'] / d[16]['sJ'] if d[16]['sJ'] > 1e-12 else np.inf
        cS = d[20]['cS']
        gapless = bool(cS > 0.7 and rJ < 1.15)
        ising = bool(0.3 < cS < 0.7 and rJ < 1.15)
        x0s = [d[L]['x0'] for L in sizes]
        flag = bool(gapless and max(x0s) < 1.3 and x0s[-1] <= x0s[0] + 0.05)
        stable = bool(flag and max(x0s) - min(x0s) < 0.15)
        out[t] = dict(fam=r['fam'], p=r['p'], rg=float(rg), rJ=float(rJ), gapless=gapless, ising=ising, flag=flag, stable=stable,
                      **{f'L{L}': d[L] for L in sizes})
    return out


def maps(S):
    for cut in MF.CUTS:
        fam, fixed, (xn, x0, x1), (yn, y0, y1) = MF.CUTS[cut]
        keys = sorted(k for k in S if k.startswith(cut + '_'))
        if not keys:
            continue
        n = 1 + max(int(k.split('_')[1]) for k in keys)
        print(f'\ncut {cut} ({fam}, {fixed}); rows {yn} from {y1} (top) to {y0}, columns {xn} {x0}..{x1}')
        print('  symbols: . neither, 1 liquid (cS20<1.5), h liquid with cS20>=1.5, i Ising-like, + flagged, * stable flag')
        for j in range(n - 1, -1, -1):
            row = ''
            for i in range(n):
                k = f'{cut}_{i:02d}_{j:02d}'
                if k not in S:
                    row += ' '; continue
                s = S[k]
                if s['stable']:
                    row += '*'
                elif s['flag']:
                    row += '+'
                elif s['ising']:
                    row += 'i'
                elif not s['gapless']:
                    row += '.'
                else:
                    row += '1' if s['L20']['cS'] < 1.5 else 'h'
            print('  ' + row)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', default=os.path.join(HERE, 'data', 'ring_summary.json'))
    ap.add_argument('--maps', action='store_true')
    a = ap.parse_args()
    S = summarize(load())
    json.dump(S, open(a.json, 'w'), indent=0)
    fl = [k for k in S if S[k]['flag']]
    print('stable flags:', sum(S[k]['stable'] for k in S), ' Ising-like:', sum(S[k]['ising'] for k in S))
    print(f'{len(S)} points, {sum(S[k]["gapless"] for k in S)} gapless, {len(fl)} flagged')
    for k in sorted(fl, key=lambda k: S[k]['L20']['x0']):
        s = S[k]
        print(k, s['p'], 'x0 =', [round(s[f'L{L}']['x0'], 3) for L in (16, 18, 20)],
              'sJ =', [round(s[f'L{L}']['sJ'], 3) for L in (16, 18, 20)], 'cS20 =', round(s['L20']['cS'], 3))
    if a.maps:
        maps(S)
