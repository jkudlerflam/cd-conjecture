#!/usr/bin/env python3
"""Sp(4) small-angle analysis: dimer doublet versus CFT tower.
For each (theta, L): ground momentum, Casimir of the ground state, the splitting
D0pi between the lowest singlets at k=0 and k=pi, the lowest gap above that doublet
(Dgap, any momentum and Casimir), and entanglement slope. In a dimerized phase
D0pi/Dgap -> 0 exponentially in L; in a CFT the ratio tends to a constant."""
import json, sys
import numpy as np
from est import c_ent

recs = []
for f in sys.argv[1:]:
    for line in open(f):
        r = json.loads(line)
        if 'coeff' in r:
            r['theta'] = float(np.arctan2(r['coeff'][1], r['coeff'][0]))
        recs.append(r)


def singlets(r, k):
    es = r['spec'].get(str(k), [])
    cas = r.get('casimir', {}).get(str(k))
    if cas is None:
        return list(es)
    return [e for e, (c, v) in zip(es, cas) if abs(c) < 1e-6]


rows = []
for r in recs:
    L = r['L']
    s0 = singlets(r, 0)
    sp = singlets(r, L // 2)
    e0 = r['E0']
    allv = sorted([(e, int(k)) for k, es in r['spec'].items() for e in es])
    if not s0 or not sp:
        continue
    pair = sorted([s0[0], sp[0]])
    D0pi = pair[1] - pair[0]
    # lowest level that is neither of the two doublet states
    rest = [e for e, k in allv if not (abs(e - s0[0]) < 1e-10 and k == 0) and not (abs(e - sp[0]) < 1e-10 and k == L // 2)]
    Dgap = rest[0] - pair[1] if rest else np.nan
    cas0 = r.get('casimir', {}).get(str(r['k0']), [[np.nan, 0]])[0][0]
    rows.append(dict(theta=r['theta'], L=L, Lth=L * r['theta'], k0=r['k0'], C2=cas0, D0pi=D0pi * L, Dgap=Dgap * L,
                     ratio=D0pi / Dgap if Dgap > 0 else np.nan, cent=c_ent(r, 2), e0=e0 / L))
rows.sort(key=lambda x: (x['theta'], x['L']))
print('%7s %3s %6s %3s %5s %9s %9s %8s %7s %9s' % ('theta', 'L', 'L*th', 'k0', 'C2', 'L*D0pi', 'L*Dgap', 'ratio', 'c_ent', 'E0/L'))
for x in rows:
    print('%7.4f %3d %6.2f %3d %5.2f %9.4f %9.4f %8.4f %7.3f %9.5f' % (x['theta'], x['L'], x['Lth'], x['k0'], x['C2'],
                                                                       x['D0pi'], x['Dgap'], x['ratio'], x['cent'], x['e0']))
json.dump(rows, open('data/sp4_small_table.json', 'w'), indent=0)
