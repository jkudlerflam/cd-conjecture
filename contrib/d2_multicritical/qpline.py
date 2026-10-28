#!/usr/bin/env python3
"""Quasiparticle dispersion E_qp(k) = E0(m, P=-1) - E0(0, P=+1) on a ring of L sites
at given Kitaev-family parameters (trivial side of the band-edge Ising line).
usage: qpline.py L 'tau=1,delta=0.5,V=-1,h=0.6' D1 D2 ...
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring
from scan import bond_of

L = int(sys.argv[1])
base = {s.split('=')[0]: float(s.split('=')[1]) for s in sys.argv[2].split(',')}
out = []
for D in map(float, sys.argv[3:]):
    kw = dict(base); kw['D'] = D
    h = bond_of('kit', **kw)
    e0 = ring.lowest(ring.sector_hamiltonian(L, h, 0, 1), k=1)[0][0]
    E = np.array([ring.lowest(ring.sector_hamiltonian(L, h, m, -1), k=1)[0][0] - e0 for m in range(L)])
    ms = np.array([m if m <= L // 2 else m - L for m in range(L)])
    o = np.argsort(ms)
    print(f'D={D:.3f}  E_qp(k=2pi m/L), m={ms[o].tolist()}')
    print('   ', ' '.join(f'{x:6.3f}' for x in E[o]))
    out.append(dict(D=D, m=ms[o].tolist(), E=E[o].tolist()))
json.dump(dict(L=L, base=base, res=out), open(f"data/qp_L{L}_{sys.argv[2].replace(',', '_')}.json", 'w'))
