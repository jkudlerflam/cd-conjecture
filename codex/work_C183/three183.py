#!/usr/bin/env python3
"""C183: three-body check (C175 nbody_box.py, bare band, |t2| = 1) with the full-chain fitted V2, V3 (twob_map.json,
model A) at unbound stage-T points.  Threshold: 3 eps_min (no dimer below threshold at these points).
A trimer is declared if E0 - 3 eps_min < -1e-3 at R and stays within 10 percent at the larger R.
K: the aab and abb momenta (2 k_a + k_b, k_a + 2 k_b) and a few others."""
import json, sys, numpy as np
import nbody_box175 as NB
from scat175 import valleys
R = json.load(open('twob_map.json'))
pts = [(float(a), float(b)) for a, b in (x.split(':') for x in sys.argv[1].split(','))]
Rs = [int(x) for x in sys.argv[2].split(',')]
for r, wa in pts:
    rec = [q for q in R if q['r'] == r and q['wa'] == wa][0]
    V2, V3 = rec['A']['V']['2'], rec['A']['V']['3']
    w = 1j * wa
    vs = valleys(w); em = vs[0][1]
    ka, kb = vs[0][0], vs[1][0]
    Ks = sorted({round(float(np.angle(np.exp(1j * K))), 4) for K in (2 * ka + kb, ka + 2 * kb, 3 * ka, 3 * kb, 0.0, np.pi)})
    for K in Ks:
        out = []
        for Rb in Rs:
            e = NB.lowest(3, K, w, Rb, V2=V2, V3=V3)
            e0 = float(np.real(np.atleast_1d(e)[0]))
            out.append(e0 - 3 * em)
        print(f"r={r:g} |w|={wa:g} V2={V2:+.4f} V3={V3:+.4f} K={K:+.4f}: E0 - 3 eps_min = " + ', '.join(f'R={Rb}: {v:+.5f}' for Rb, v in zip(Rs, out)), flush=True)
