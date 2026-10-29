#!/usr/bin/env python3
"""Open-chain U(1)^3 DMRG runs (qdmrg.py) for LV_2, its Neel control, and the free reference.

usage: run_obc.py NAME --model {lv,kxx} --L 120 --M 72 --N 24 18 [--vz 0 3] [--J 0.5 --taus 1 0.5]
                  --chis 64 128 256 512 [--sweeps 2 --maxsweeps 6]
  lv  : LV_2 open chain in the sector (M cells, N_f fermions of flavour f); mu_v and eps drop out
  kxx : three decoupled XX chains (d = 8) with particle numbers (M, N_0, N_1), the exact reference
Writes data/obc/NAME.json after every chi (S on every bond, densities, energy, truncation).
Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-28.
"""
import os, sys, json, argparse, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.join(HERE, '..', 'vacancy_mobile_numerics'), HERE):
    sys.path.insert(0, p)
import numpy as np
import qdmrg as Q
import models as MD


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='lv')
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--M', type=int, required=True)
    ap.add_argument('--N', type=int, nargs='+', required=True)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.5)
    ap.add_argument('--taus', type=float, nargs='+', default=[1.0, 0.5])
    ap.add_argument('--vz', type=float, nargs='+', default=None)
    ap.add_argument('--chis', type=int, nargs='+', default=[64, 128])
    ap.add_argument('--sweeps', type=int, default=2)
    ap.add_argument('--maxsweeps', type=int, default=6)
    ap.add_argument('--etol', type=float, default=1e-9)
    ap.add_argument('--noise', type=float, default=1e-3)
    ap.add_argument('--noise_sweeps', type=int, default=2)
    ap.add_argument('--out', default=None)
    ap.add_argument('--vels', type=float, nargs='+', default=None, help='kxx: chain hoppings (charge chain first)')
    ap.add_argument('--state', default=None, help='checkpoint file (scratch only; resumes above the saved chi)')
    a = ap.parse_args()
    if a.model == 'lv':
        import logical_vac as LV
        n = len(a.taus)
        h = LV.model_bond(n, a.t, a.J, 0.0, a.taus, [0.0] * n, a.vz)
        qloc = MD.lv_charges(n); cfg = MD.lv_config(a.L, a.M, a.N)
    else:
        h, qloc = MD.kxx_bond(1 + len(a.N), a.vels)
        cfg = MD.kxx_config(a.L, [a.M] + a.N)
    out = a.out or os.path.join(HERE, 'data', 'obc', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), results=[])
    state = None; chis = a.chis
    if a.state and os.path.exists(a.state):
        st = pickle.load(open(a.state, 'rb'))
        state = st['mps']; chis = [c for c in a.chis if c > st['chi']]
        if os.path.exists(out):
            data['results'] = json.load(open(out))['results']
        print(a.name, 'resuming above chi', st['chi'], flush=True)

    def on_chi(dm, rec):
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        if a.state:
            pickle.dump(dict(chi=rec['chi'], mps=dm.get_state()), open(a.state + '.tmp', 'wb'))
            os.replace(a.state + '.tmp', a.state)
    Q.run(h, qloc, a.L, cfg, chis, sweeps=a.sweeps, maxsweeps=a.maxsweeps, etol=a.etol,
          noise0=a.noise, noise_sweeps=a.noise_sweeps, state=state,
          log=lambda s: print(a.name, s, flush=True), on_chi=on_chi)


if __name__ == '__main__':
    main()
