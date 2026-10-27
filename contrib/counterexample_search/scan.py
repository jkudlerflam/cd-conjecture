#!/usr/bin/env python3
"""Dense scans of symmetric NN families. Bond components are built once per size
(momentum blocks cached) and every parameter point reuses them.

Families (identity direction removed, overall positive scale fixed):
  orbital  (AB, A, B) on S^2, SU(2)xSU(2) spin-orbital; --half keeps c_A >= c_B
  spin32   SU(2) spin 3/2, projectors P0..P3 in the D3 frame orthogonal to (1,1,1,1), on S^2
  g2       G2 on the 7, (P7, P14, P27) on S^2
  sp4      Sp(4) on C^4, (P5, P10) on S^1
  so5      SO(5) vector, (P, K) on S^1
  so4      SO(4) vector = spin-orbital line c_A = c_B, (P, K) on S^1
Usage: python3 scan.py --family orbital --npts 600 --half --sizes 8 10 --chunk 0 1 --out f.jsonl
"""
import argparse, json, time
import numpy as np
from ced import Sector, measure
import models


def fib_sphere(n):
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    r = np.sqrt(1 - z * z)
    ph = i * np.pi * (3 - np.sqrt(5))
    return np.stack([r * np.cos(ph), r * np.sin(ph), z], 1)


def family(name):
    if name == 'orbital':
        m = models.spin_orbital()
        return dict(d=4, hs=m['hs'], charges=m['charges'], gens=m['gens'], tgt=lambda L: [L % 2, L % 2], dim=3)
    if name == 'spin32':
        m = models.spin_s(1.5)
        B = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float)
        B /= np.linalg.norm(B, axis=1)[:, None]
        hs = [sum(B[i, J] * m['P'][J] for J in range(4)) for i in range(3)]
        return dict(d=4, hs=hs, charges=m['charges'], gens=m['gens'], tgt=lambda L: [L % 2], dim=3, frame=B.tolist())
    if name == 'g2':
        m = models.g2()
        return dict(d=7, hs=m['hs'], charges=m['charges'], gens=m['gens'], tgt=lambda L: [0, 0], dim=3)
    if name == 'sp4':
        m = models.sp4()
        return dict(d=4, hs=m['hs'], charges=m['charges'], gens=m['gens'], tgt=lambda L: [L % 2, 0], dim=2)
    if name in ('so5', 'so4', 'so6', 'so7'):
        N = int(name[2])
        m = models.son_vector(N)
        r = N // 2
        tgt = (lambda L: [0] * r) if N % 2 else (lambda L: [L % 2] + [0] * (r - 1))
        return dict(d=N, hs=m['hs'], charges=m['charges'], gens=m['gens'], tgt=tgt, dim=2)
    raise ValueError(name)


def points(f, a):
    if f['dim'] == 2:
        th = 2 * np.pi * np.arange(a.npts) / a.npts
        return [[np.cos(t), np.sin(t)] for t in th]
    P = fib_sphere(a.npts)
    if a.half:
        P = P[P[:, 1] >= P[:, 2] - 1e-12]
    return P.tolist()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--family', required=True)
    p.add_argument('--npts', type=int, default=200)
    p.add_argument('--half', action='store_true')
    p.add_argument('--points', type=str, default=None, help='JSON file with explicit coefficient list')
    p.add_argument('--sizes', type=int, nargs='+', required=True)
    p.add_argument('--chunk', type=int, nargs=2, default=[0, 1])
    p.add_argument('--nlev', type=int, default=3)
    p.add_argument('--cas', type=int, default=0)
    p.add_argument('--noent', action='store_true')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    f = family(a.family)
    pts = json.load(open(a.points)) if a.points else points(f, a)
    idx = list(range(len(pts)))[a.chunk[0]::a.chunk[1]]
    for L in a.sizes:
        t0 = time.time()
        sec = Sector(f['d'], L, f['hs'], f['charges'], f['tgt'](L), cache=True)
        for i in idx:
            t = time.time()
            r = measure(sec, pts[i], nlev=a.nlev, entropy=not a.noent,
                        gens=f['gens'] if a.cas else None, ncas=a.cas)
            r.update(family=a.family, point=i, coeff=[float(c) for c in pts[i]], seconds=time.time() - t)
            with open(a.out, 'a') as fh:
                fh.write(json.dumps(r) + '\n')
        print('L=%d done %d points in %.0fs (dim %d)' % (L, len(idx), time.time() - t0, len(sec.basis)), flush=True)


if __name__ == '__main__':
    main()
