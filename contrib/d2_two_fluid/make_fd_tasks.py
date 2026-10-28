#!/usr/bin/env python3
"""Write slurm/fd_tasks.txt: calibrations, candidate C, and the three TDA survivors with the
widest product-state windows (data/cl_survivors_*.json), each at two chemical potentials
inside its window."""
import json, glob
import numpy as np
R = []
for f in glob.glob('data/cl_survivors_*.json'):
    R += json.load(open(f))
R = sorted(R, key=lambda r: -r['window'])[:3]
lines = ["xx xx 500", "xxh xxh 500", "xx2 xx2 300", "u1ctl p:1,0,-1.03,-1.8,0,0,0 500",
         "C_m2.40 p:1,1.40,-1.03,-2.40,1.43,0.38,1.40 500", "C_m2.30 p:1,1.40,-1.03,-2.30,1.43,0.38,1.40 500",
         "C_m2.22 p:1,1.40,-1.03,-2.22,1.43,0.38,1.40 500"]
for i, r in enumerate(R):
    p = r['p']; g = complex(*p['gR']); ms = r['mu_star']; md = r['mu_dense_cl']
    for j, f in enumerate([0.35, 0.7]):
        mu = ms + f * (md - ms)
        lines.append(f"S{i+1}_{j} p:1,{p['phi']:.6f},{p['V']:.6f},{mu:.4f},{p['gL']:.6f},{abs(g):.6f},{np.angle(g):.6f} 500")
open('slurm/fd_tasks.txt', 'w').write('\n'.join(lines) + '\n')
print(open('slurm/fd_tasks.txt').read())
