#!/usr/bin/env python3
"""Dynamical exponent at a candidate endpoint: the lowest few excitation energies
over all (momentum, parity) sectors on rings of length L, and z from
ln eps_n = -z ln L + const (n-th excitation, n = 1..4).  A relativistic point gives
z = 1 (L * eps_n constant); a Lifshitz point gives z > 1.

usage: endpoint.py FAMILY 'k=v,...' --L 10 12 14 16 18 [--out data/endpoint_X.json]
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring
from scan import bond_of

ap = argparse.ArgumentParser()
ap.add_argument('fam'); ap.add_argument('params')
ap.add_argument('--L', type=int, nargs='+', default=[10, 12, 14, 16, 18])
ap.add_argument('--out', default=None)
a = ap.parse_args()
kw = {s.split('=')[0]: float(s.split('=')[1]) for s in a.params.split(',')}
h = bond_of(a.fam, **kw)
rows = []
for L in a.L:
    sp, key, psi = ring.spectrum(L, h, k=3, want_gs=True)
    e = np.sort(np.concatenate([np.asarray(v)[np.isfinite(v)] for v in sp.values()]))
    eps = e[1:6] - e[0]
    S = ring.entropies(psi, L)
    rows.append(dict(L=L, gs=list(key), eps=eps.tolist(), cS=ring.c_chord(S, L)))
    print(L, key, 'L*eps =', np.round(L * eps, 4).tolist(), 'cS=%.3f' % rows[-1]['cS'], flush=True)
Ls = np.array(a.L, float)
for n in range(4):
    y = np.log([r['eps'][n] for r in rows])
    z = -np.polyfit(np.log(Ls[-4:]), y[-4:], 1)[0]
    print(f'n={n + 1}: z (last four L) = {z:.3f}')
if a.out:
    json.dump(dict(fam=a.fam, params=kw, rows=rows), open(a.out, 'w'), indent=1)
