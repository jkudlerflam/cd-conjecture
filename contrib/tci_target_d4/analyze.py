#!/usr/bin/env python3
"""Spectroscopy of chart.py output.  For each point and size, levels in units of 2 pi / L:
  Sgs    ground-state spin
  cur    lowest spin-1 level at k0 +- 2pi/L     (SU(2) current, x = 1, velocity of the current sector)
  pi0    lowest spin-0 level at k0 + pi         (Majorana masses: x = 1 in the c = 3 phase;
                                                  x = 1/5 (TCI energy) at the glued point; -> 0 if dimerized)
  pi0b   second spin-0 level at k0 + pi
  pis    lowest level at k0 + pi of any spin, and its spin
  k00    lowest excited spin-0 level at k0       (x = 2 in the c = 3 phase; 6/5 twice at the glued point)
  q      lowest level at k0 +- pi/2 and its spin (x = 3/4 in SU(4)_1, 7/10 at the glued point)
  gap    lowest excitation, its momentum offset and spin
c_half(8 -> 12) = 3 [S(6) - S(4)] / ln(3/2)."""
import json, glob, sys
import numpy as np
from collections import defaultdict


def spin_of(c):
    return (-1 + np.sqrt(1 + 4 * max(c, 0))) / 2


def fold(k, L):
    k %= L
    return min(k, L - k)


def load(files):
    recs = defaultdict(dict)
    for f in files:
        for line in open(f):
            r = json.loads(line)
            recs[(round(r['theta'], 6), round(r['eta'], 6))][r['L']] = r
    return recs


def levels(r):
    """list of (E - E0 in units 2pi/L, k offset from k0 folded, spin) for all labelled states"""
    L, k0, e0 = r['L'], r['k0'], r['E0']
    f = L / (2 * np.pi)
    out = []
    for k, es in r['spec'].items():
        cas = r['casimir'][k]
        k = int(k)
        dk = min(fold(k - k0, L), fold(k + k0, L))   # k and -k are degenerate
        for e, (c, var) in zip(es, cas):
            out.append(((e - e0) * f, dk, round(2 * spin_of(c)) / 2 if var < 1e-5 else -1.0, var))
    out.sort()
    return out


def summary(r):
    L = r['L']
    lv = levels(r)
    exc = [x for x in lv if x[0] > 1e-7]
    pick = lambda dk, s=None: [x[0] for x in exc if x[1] == dk and (s is None or x[2] == s)]
    first = lambda l: l[0] if l else np.nan
    out = dict(L=L, k0=r['k0'], Sgs=lv[0][2], e0L=r['E0'] / L)
    out['cur'] = first(pick(1, 1))
    p0 = pick(L // 2, 0)
    out['pi0'] = first(p0)
    out['pi0b'] = p0[1] if len(p0) > 1 else np.nan
    ps = [x for x in exc if x[1] == L // 2]
    out['pis'], out['pis_S'] = (ps[0][0], ps[0][2]) if ps else (np.nan, np.nan)
    out['k00'] = first(pick(0, 0))
    q = [x for x in exc if x[1] == L // 4] if L % 4 == 0 else []
    out['q'], out['q_S'] = (q[0][0], q[0][2]) if q else (np.nan, np.nan)
    out['gap'], out['gap_dk'], out['gap_S'] = exc[0][:3] if exc else (np.nan, np.nan, np.nan)
    out['S_half'] = r['S'][-1]
    return out


def table(recs, sizes=(8, 12)):
    rows = []
    for key in sorted(recs):
        b = recs[key]
        if not all(L in b for L in sizes):
            continue
        s = {L: summary(b[L]) for L in sizes}
        L1, L2 = sizes
        row = dict(theta=key[0] / np.pi, eta=key[1], chalf=3 * (s[L2]['S_half'] - s[L1]['S_half']) / np.log(L2 / L1))
        for L in sizes:
            for k, v in s[L].items():
                row['%s_%d' % (k, L)] = v
        rows.append(row)
    return rows


if __name__ == '__main__':
    files = sys.argv[1:] or glob.glob('data/grid/*.jsonl')
    rows = table(load(files))
    hdr = '%6s %6s | %5s | %3s %3s | %5s %5s %5s %5s | %5s %5s %5s | %5s %5s %3s'
    print(hdr % ('th/pi', 'eta', 'chalf', 'S8', 'S12', 'cur', 'pi0', 'pi0b', 'k00', 'q', 'pis', 'pS', 'gap', 'gdk', 'gS'))
    for r in rows:
        print('%6.3f %6.2f | %5.2f | %3.1f %3.1f | %5.3f %5.3f %5.3f %5.3f | %5.3f %5.3f %5.1f | %5.3f %5.0f %3.1f' % (
            r['theta'], r['eta'], r['chalf'], r['Sgs_8'], r['Sgs_12'], r['cur_12'], r['pi0_12'], r['pi0b_12'], r['k00_12'],
            r['q_12'], r['pis_12'], r['pis_S_12'], r['gap_12'], r['gap_dk_12'], r['gap_S_12']))
