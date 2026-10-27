#!/usr/bin/env python3
"""Spin-3/2 chart around the Sp(4) line, h = cos(t) P2 + (sin(t)+eta) P1 + (sin(t)-eta) P3
(P_J projects on total pair spin J; eta = 0 is the Sp(4) line, eta != 0 keeps only SU(2)).

Records, per point and size, the low spectrum at every momentum with total-spin labels,
for the tricritical-Ising test: the lowest spin-0 level at the ground momentum, in
units of the lowest adjacent-momentum level, drops from about 2 (marginal operators of the
c = 3 phase) towards x = 1/5 times a velocity ratio near a TCI point.
Usage: python3 tci_scan.py --sizes 8 12 --thetas ... --etas ... --chunk i n --out f.jsonl"""
import argparse, json, time
import numpy as np
from ced import Sector, measure
import models

p = argparse.ArgumentParser()
p.add_argument('--sizes', type=int, nargs='+', default=[8, 12])
p.add_argument('--thetas', type=float, nargs='+', required=True)
p.add_argument('--etas', type=float, nargs='+', required=True)
p.add_argument('--chunk', type=int, nargs=2, default=[0, 1])
p.add_argument('--nlev', type=int, default=4)
p.add_argument('--out', required=True)
a = p.parse_args()

m = models.spin_s(1.5)
P0, P1, P2, P3 = m['P']
hs = [P2, P1, P3]
gens = m['gens']
pts = [(t, e) for t in a.thetas for e in a.etas]
pts = pts[a.chunk[0]::a.chunk[1]]
for L in a.sizes:
    sec = Sector(4, L, hs, m['charges'], [L % 2], cache=True)
    t0 = time.time()
    for (t, e) in pts:
        coeff = [np.cos(t), np.sin(t) + e, np.sin(t) - e]
        r = measure(sec, coeff, nlev=a.nlev, entropy=True, gens=gens, ncas=a.nlev)
        r.update(theta=t, eta=e, coeff=coeff, family='spin32chart')
        with open(a.out, 'a') as f:
            f.write(json.dumps(r) + '\n')
    print('L=%d %d points %.0fs' % (L, len(pts), time.time() - t0), flush=True)
