#!/usr/bin/env python3
"""C177: open-chain DMRG (fc.py path: dmrg2.py, real gauge where it applies + shift + edge field) against sparse ED (c166lib.H_sparse) at L = 10, 12, 14,
for the tuned F_8 and F_16 bonds and the generic point G_8."""
import json, subprocess, sys, numpy as np
import c166lib as X, band as B
T = {d['r']: d['couplings'] for d in json.load(open('tune.json'))}
cases = [('F8', 8.0, {}), ('F16', 16.0, {}), ('G8', 8.0, dict(jxpy=0.02*0.0579, jxz=0.2*0.0579, hx=0.2*0.0579))]
for tag, r, g in cases:
    cp = dict(T[r]); cp['Jx'] = (1 + g.get('jxpy', 0)) / 2; cp['Jy'] = (g.get('jxpy', 0) - 1) / 2
    Jxz = g.get('jxz', 0.0); hx = g.get('hx', 0.0)
    J = np.array([[cp['Jx'], cp['D'], Jxz], [-cp['D'], cp['Jy'], 0], [Jxz, 0, cp['Jz']]]); f = np.array([hx, 0, cp['hz']])
    for L in (10, 12, 14):
        e_ed = X.ed_e0(J, f, L)[0]
        args = ['python3', 'fc.py', f'chk_{tag}', str(L), '--r', str(r), '--chis', '16', '32', '64', '64', '64', '--outdir', 'edchk2']
        for k, v in g.items(): args += [f'--{k}', str(v)]
        subprocess.run(args, capture_output=True)
        d = json.load(open(f'edchk2/chk_{tag}_{L}.json'))
        print(tag, L, f'ED {e_ed:.12f} DMRG {d["E"]:.12f} diff {d["E"]-e_ed:.1e}', flush=True)
