#!/usr/bin/env python3
"""Tables for the spin-3/2 chart (tci_scan.py output).
Per (theta, eta): ground spin, c_half(8->12), c_ent(12), and spin-resolved levels
in units 2 pi / L at L = 12:
  e1   lowest level at k0 +/- 1 (any spin)          ~ v * 1 (current)
  s0pi lowest spin-0 level at k0 + pi                (Ising energy x = 1, or dimer partner -> 0)
  s0k0 lowest spin-0 excited level at k0             (x = 2 in the c = 3 phase)
A tricritical point would show s0pi/e1 well below its c = 3 value without a dimer doublet,
and c_half above the SU(4) calibration (3.05 for 8 -> 12)."""
import json, glob, sys
import numpy as np
from collections import defaultdict
from est import c_ent, c_half

recs = defaultdict(dict)
for f in sys.argv[1:] or glob.glob('data/tci/tci_chart_c*.jsonl'):
    for line in open(f):
        r = json.loads(line)
        recs[(round(r['theta'], 6), round(r['eta'], 6))][r['L']] = r
S = lambda c: (-1 + np.sqrt(1 + 4 * max(c, 0))) / 2


def lev(r):
    L, k0, e0 = r['L'], r['k0'], r['E0']
    f = L / (2 * np.pi)
    out = dict(Sgs=S(r['casimir'][str(k0)][0][0]))
    kpi = (k0 + L // 2) % L
    kpi = min(kpi, L - kpi)
    ka = min((k0 + 1) % L, L - (k0 + 1) % L)
    out['e1'] = (r['spec'][str(ka)][0] - e0) * f
    s0 = [e for e, c in zip(r['spec'][str(kpi)], r['casimir'][str(kpi)]) if abs(c[0]) < 1e-6 and abs(e - e0) > 1e-9]
    out['s0pi'] = (s0[0] - e0) * f if s0 else np.nan
    s0k = [e for e, c in zip(r['spec'][str(k0)], r['casimir'][str(k0)]) if abs(c[0]) < 1e-6 and abs(e - e0) > 1e-9]
    out['s0k0'] = (s0k[0] - e0) * f if s0k else np.nan
    allv = sorted(e for es in r['spec'].values() for e in es)
    out['gap'] = (allv[1] - e0) * f
    return out


rows = []
print('%6s %6s | %4s %4s | %6s %6s | %6s %6s %6s %6s' % ('th/pi', 'eta', 'S8', 'S12', 'chalf', 'cent12', 'e1', 's0pi', 's0k0', 'gap'))
for (t, e) in sorted(recs):
    b = recs[(t, e)]
    if 8 not in b or 12 not in b:
        continue
    l8, l12 = lev(b[8]), lev(b[12])
    ch = c_half(b[8], b[12])
    row = dict(theta=t, eta=e, S8=l8['Sgs'], S12=l12['Sgs'], chalf=ch, cent12=c_ent(b[12], 3), **{k: l12[k] for k in ('e1', 's0pi', 's0k0', 'gap')})
    rows.append(row)
    print('%6.3f %6.2f | %4.1f %4.1f | %6.2f %6.2f | %6.3f %6.3f %6.3f %6.3f' % (t / np.pi, e, row['S8'], row['S12'], ch, row['cent12'],
                                                                            row['e1'], row['s0pi'], row['s0k0'], row['gap']))
json.dump(rows, open('data/tci_chart_table.json', 'w'), indent=0)
