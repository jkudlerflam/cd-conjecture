#!/usr/bin/env python3
"""Finite-entanglement scaling of the spin-3/2 chart with a 4-site (or n-site) VUMPS cell.
Uses vumps.py from ../d2_followup_220_1282 (multi-site VUMPS, numpy/scipy only).

Per bond dimension it records the energy per site, the entanglement entropy on every bond of
the cell, every bond energy (so dimerization shows up as alternating bond energies), the
leading transfer-matrix eigenvalues and xi = -n / ln|l1/l0|.
c_fes is the slope of S against ln(xi) times 6, over successive chi.

usage: vumps_run.py --theta T --eta E --n 4 --chis 16 32 64 --out f.json   (T in units of pi)
       vumps_run.py --a a0 a1 a2 a3 ...
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('../d2_followup_220_1282', '../counterexample_search'):
    sys.path.insert(0, os.path.join(HERE, p))
import numpy as np
import vumps as V
import models


def bond_energies(psi, h, d):
    out = []
    n = psi.n
    for k in range(n):
        AC = np.tensordot(psi.AL[k], psi.C[(k + 1) % n], ([2], [0]))
        th = np.tensordot(AC, psi.AR[(k + 1) % n], ([2], [0]))
        t = th.transpose(1, 2, 0, 3).reshape(d * d, -1)
        out.append(float(np.real(np.vdot(t, h @ t)) / np.real(np.vdot(t, t))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--theta', type=float, default=None)
    ap.add_argument('--eta', type=float, default=0.0)
    ap.add_argument('--a', type=float, nargs=4, default=None)
    ap.add_argument('--n', type=int, default=4)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 32, 64])
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    if a.a is not None:
        coef = list(a.a)
    else:
        t = a.theta * np.pi
        coef = [0.0, np.sin(t) + a.eta, np.cos(t), np.sin(t) - a.eta]
    m = models.spin_s(1.5)
    h = sum(c * P for c, P in zip(coef, m['P']))
    h = (h + h.T) / 2
    W = V.mpo_from_term(h, 2, d=4)
    data = dict(theta=a.theta, eta=a.eta, a=coef, n=a.n, mode=a.mode, results=[])
    t0 = time.time()
    sz = np.real(m['sz'])

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['ebond'] = bond_energies(psi, h, 4)
        rec['sz'] = V.local_expect(psi, [sz])[0].tolist()
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        print('chi=%d e=%.10f err=%.1e it=%d S=%s xi=%.2f eb=%s t=%.0fs' % (
            rec['chi'], rec['e'], rec['err'], rec['iters'], np.round(rec['S'], 4).tolist(), rec['xi'],
            np.round(rec['ebond'], 5).tolist(), rec['time_s']), flush=True)

    V.run([W] * a.n, a.n, a.chis, d=4, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode,
          log=lambda rec: None, on_chi=save)


if __name__ == '__main__':
    main()
