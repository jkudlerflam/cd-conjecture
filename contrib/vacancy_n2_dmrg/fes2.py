#!/usr/bin/env python3
"""Finite-entanglement scaling (dense one-site VUMPS of ../d2_followup_220_1282/vumps.py) for LV_2,
its control, and the three-mode free reference, at large chi.

usage: fes2.py NAME --model lv  --taus 1 0.5 --eps E0 E1 [--vz 0 3] --J 0.5 --mu_v MU --chis ...
       fes2.py NAME --model kxx --vels 1 r1 r2 --fills f0 f1 f2 --chis ...
  lv  : logical vacancy model LV_n of ../vacancy_mobile_numerics/logical_vac.py (d = 1 + 2^n)
  kxx : k decoupled XX (free-fermion) chains with hoppings vels and fillings fills (chemical potential
        mu_i = -2 v_i cos(pi f_i)), c = k, one velocity per chain
Per chi it stores energy, gradient, S, the 24 leading transfer-matrix eigenvalues, and local densities.
Output data/fes/NAME.json (no MPS is written).
Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-28.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.join(HERE, '..', 'd2_followup_220_1282'), os.path.join(HERE, '..', 'vacancy_mobile_numerics'),
          os.path.join(HERE, 'lib')):
    sys.path.insert(0, p)
import numpy as np
import vumps as V
import logical_vac as LV


def model(a):
    sp = np.array([[0., 1.], [0., 0.]])
    if a.model == 'lv':
        nf = len(a.taus)
        h = LV.model_bond(nf, a.t, a.J, a.mu_v, a.taus, a.eps, a.vz)
        d = 1 + 2 ** nf
        n = np.diag(np.r_[0, np.ones(d - 1)])
        fl = [np.diag(np.r_[0, [(b >> (nf - 1 - f)) & 1 for b in range(d - 1)]]).astype(float) for f in range(nf)]
        return h, d, [n] + fl
    if a.model == 'kxx':
        k = len(a.vels)
        I = np.eye(2)

        def emb(o, i):
            m = np.ones((1, 1))
            for j in range(k):
                m = np.kron(m, o if j == i else I)
            return m
        h = 0
        nops = []
        for i, (v, f) in enumerate(zip(a.vels, a.fills)):
            A = np.kron(emb(sp, i), emb(sp.T, i))
            h = h - v * (A + A.T)
            ni = emb(np.diag([0., 1.]), i)        # state 1 = occupied (sp annihilates it)
            mu = -2 * v * np.cos(np.pi * f)
            h = h - 0.5 * mu * (np.kron(ni, np.eye(2 ** k)) + np.kron(np.eye(2 ** k), ni))
            nops.append(ni)
        return h, 2 ** k, nops
    raise ValueError(a.model)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='lv')
    for k, v in dict(t=1.0, J=0.5, mu_v=0.3).items():
        ap.add_argument('--' + k, type=float, default=v)
    ap.add_argument('--taus', type=float, nargs='+', default=[1.0, 0.5])
    ap.add_argument('--eps', type=float, nargs='+', default=[0.0, 0.0])
    ap.add_argument('--vz', type=float, nargs='+', default=None)
    ap.add_argument('--vels', type=float, nargs='+', default=[1.0, 0.45, 0.19])
    ap.add_argument('--fills', type=float, nargs='+', default=[0.5, 0.5, 0.5])
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32])
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--stol', type=float, default=1e-8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--nev', type=int, default=24)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    h, d, ops = model(a)
    W = V.mpo_from_term(h, 2, d=d)
    out = a.out or os.path.join(HERE, 'data', 'fes', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), results=[])
    t0 = time.time()

    def on_chi(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['local'] = V.local_expect(psi, ops).tolist()
        try:
            lam = V.transfer_spectrum(psi, nev=a.nev)
            rec['lamN'] = [[float(z.real), float(z.imag)] for z in lam]
        except Exception as e:  # noqa: BLE001
            rec['lamN_error'] = str(e)
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} local={np.round(rec['local'], 5).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)

    V.run([W], 1, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode='damped',
          log=lambda rec: None, on_chi=on_chi, stol=a.stol)


if __name__ == '__main__':
    main()
