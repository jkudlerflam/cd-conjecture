#!/usr/bin/env python3
"""Run one model at one parameter point for several sizes; append JSON lines.

Examples
  python3 run.py --model uls --d 3 --sizes 6 9 12 --out cal_uls3.jsonl
  python3 run.py --model sp4 --theta 0.0654 --sizes 8 10 12 --cas 4 --out sp4.jsonl
  python3 run.py --model orbital --coeff 1 0.25 0.25 --sizes 8 10 12 --out orb.jsonl
"""
import argparse, json, time
import numpy as np
from ced import Sector, measure, swap
import models


def build(a):
    """Return (d, hs, coeff, charges, gens, target_fn, label)."""
    if a.model == 'uls':
        m = models.uls(a.d)
        d = a.d
        tgt = lambda L: [L // d + (1 if i < L % d else 0) for i in range(d - 1)]
        return d, m['hs'], [1.0], m['charges'], None, tgt
    if a.model == 'sp4':
        m = models.sp4()
        coeff = a.coeff if a.coeff else [np.cos(a.theta), np.sin(a.theta)]
        return 4, m['hs'], coeff, m['charges'], m['gens'], lambda L: [L % 2, 0]
    if a.model == 'orbital':
        m = models.spin_orbital()
        return 4, m['hs'], a.coeff, m['charges'], m['gens'], lambda L: [L % 2, L % 2]
    if a.model.startswith('spin'):
        s = {'spin1': 1, 'spin32': 1.5, 'spin12': 0.5, 'spin2': 2}[a.model]
        m = models.spin_s(s)
        d = m['d']
        if a.basis == 'proj':
            hs = m['P']
        else:   # powers of X = S.S
            hs = [np.linalg.matrix_power(m['X'], p) for p in range(1, d)]
        extra = []
        sz = m['sz']
        if a.anis:
            extra = [np.kron(sz, sz), 0.5 * (np.kron(sz @ sz, np.eye(d)) + np.kron(np.eye(d), sz @ sz))]
        return d, hs + extra, a.coeff, m['charges'], m['gens'], lambda L: [a.sz2 if a.sz2 is not None else (L * (d - 1)) % 2]
    if a.model == 'son':
        m = models.son_vector(a.d)
        coeff = a.coeff if a.coeff else [np.cos(a.theta), np.sin(a.theta)]
        r = a.d // 2
        return a.d, m['hs'], coeff, m['charges'], m['gens'], lambda L: [0] * r if (a.d % 2 or L % 2 == 0) else [1] + [0] * (r - 1)
    if a.model == 'g2':
        m = models.g2()
        return 7, m['hs'], a.coeff, m['charges'], m['gens'], lambda L: [0, 0]
    if a.model == 'decoupled':   # spin-orbital with only A and B: two decoupled Heisenberg chains
        m = models.spin_orbital()
        return 4, m['hs'], [0.0, 1.0, a.ratio], m['charges'], m['gens'], lambda L: [L % 2, L % 2]
    raise ValueError(a.model)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--model', required=True)
    p.add_argument('--d', type=int, default=3)
    p.add_argument('--theta', type=float, default=None)
    p.add_argument('--coeff', type=float, nargs='+', default=None)
    p.add_argument('--ratio', type=float, default=3.0)
    p.add_argument('--basis', default='proj')
    p.add_argument('--anis', action='store_true')
    p.add_argument('--sz2', type=int, default=None)
    p.add_argument('--target', type=int, nargs='+', default=None)
    p.add_argument('--sizes', type=int, nargs='+', required=True)
    p.add_argument('--nlev', type=int, default=4)
    p.add_argument('--cas', type=int, default=0, help='Casimir labels for lowest n states per momentum')
    p.add_argument('--noent', action='store_true')
    p.add_argument('--nocache', action='store_true')
    p.add_argument('--label', default='')
    p.add_argument('--momenta', type=int, nargs='+', default=None)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    d, hs, coeff, charges, gens, tgt = build(a)
    for L in a.sizes:
        t = time.time()
        target = a.target if a.target is not None else tgt(L)
        sec = Sector(d, L, hs, charges, target, cache=not a.nocache)
        tb = time.time() - t
        r = measure(sec, coeff, nlev=a.nlev, momenta=a.momenta, entropy=not a.noent, gens=gens if a.cas else None, ncas=a.cas)
        r.update(model=a.model, d=d, coeff=[float(c) for c in coeff], theta=a.theta, target=[int(x) for x in target],
                 label=a.label, seconds=time.time() - t, build_seconds=tb)
        with open(a.out, 'a') as f:
            f.write(json.dumps(r) + '\n')
        print(json.dumps({k: r[k] for k in ('model', 'L', 'k0', 'E0', 'dim', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
