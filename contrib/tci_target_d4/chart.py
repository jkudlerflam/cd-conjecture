#!/usr/bin/env python3
"""Spin-3/2 SU(2) chain, h = a1 P1 + a2 P2 + a3 P3 (a0 = 0), in the chart
   (a1, a2, a3) = (sin t + eta, cos t, sin t - eta),
so eta = 0 is the Sp(4) circle (cos t P5 + sin t P10) and t = pi/2, eta = 0 is SU(4).
Any direction (a1, a2, a3) except +-(1, 0, -1) has a unique (t, eta).

ED with Sz = 0 (or L mod 2) and momentum resolution, spin labels from the total Casimir
of the lowest nlev states in each momentum block, half-chain entanglement of the ground state.
Uses ced.py and models.py from ../counterexample_search.

usage: chart.py --sizes 8 12 --points pts.json --chunk i n --nlev 8 --out f.jsonl
       pts.json = [[t, eta], ...]   (t in radians)
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'counterexample_search'))
import numpy as np
from ced import Sector, measure
import ced
import models
from scipy.sparse.linalg import ArpackNoConvergence

_solve = ced.Sector.solve


def _robust_solve(self, coeff, k, n=3, tol=1e-10):
    """ced.Sector.solve, retried with a larger Krylov space when ARPACK fails to converge
    (this happens next to the ferromagnetic boundary, where many levels nearly cross)."""
    try:
        return _solve(self, coeff, k, n, tol)
    except ArpackNoConvergence:
        from scipy.sparse.linalg import eigsh
        ids, hs = self.blocks(k)
        H = sum(float(c) * h for c, h in zip(coeff, hs) if c != 0)
        rng = np.random.default_rng(99 + k)
        v0 = rng.normal(size=H.shape[0]).astype(H.dtype)
        e, v = eigsh(H, k=n, which='SA', tol=1e-9, maxiter=100000, v0=v0, ncv=min(H.shape[0] - 1, 120))
        ix = np.argsort(e)
        e, v = e[ix], v[:, ix]
        res = float(np.max(np.linalg.norm(H @ v - v * e, axis=0)))
        return e, v, res


ced.Sector.solve = _robust_solve


def bond_coeff(t, eta):
    """coefficients of (P0, P1, P2, P3)"""
    return [0.0, np.sin(t) + eta, np.cos(t), np.sin(t) - eta]


def to_chart(a):
    """(a0..a3) -> (t, eta)"""
    a = np.asarray(a, float) - a[0]
    s = np.hypot(a[2], (a[1] + a[3]) / 2)
    return float(np.arctan2((a[1] + a[3]) / 2, a[2]) % (2 * np.pi)), float((a[1] - a[3]) / 2 / s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sizes', type=int, nargs='+', default=[8, 12])
    p.add_argument('--points', required=True)
    p.add_argument('--chunk', type=int, nargs=2, default=[0, 1])
    p.add_argument('--nlev', type=int, default=8)
    p.add_argument('--momenta', type=int, nargs='*', default=None)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    m = models.spin_s(1.5)
    P0, P1, P2, P3 = m['P']
    hs = [P1, P2, P3]
    pts = json.load(open(a.points))[a.chunk[0]::a.chunk[1]]
    for L in a.sizes:
        t0 = time.time()
        sec = Sector(4, L, hs, m['charges'], [L % 2], cache=True)
        for (t, e) in pts:
            c = bond_coeff(t, e)
            ts = time.time()
            r = measure(sec, c[1:], nlev=a.nlev, momenta=a.momenta, entropy=True, gens=m['gens'], ncas=a.nlev)
            r.update(theta=t, eta=e, a=c, family='spin32chart', seconds=time.time() - ts)
            with open(a.out, 'a') as f:
                f.write(json.dumps(r) + '\n')
        print('L=%d %d points %.0fs' % (L, len(pts), time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
