#!/usr/bin/env python3
"""Open-chain U(1)^3 DMRG of the full microscopic d = 13 bond (codex/work_C88/bond13.py, read-only).

Conserved: Q (bond13 charge, twice the logical fermion number on the code), N (occupied sites) and
Mk (marker number).  The sector is fixed by a legal product state of Mk cells of length b = p + 1,
cell j = [R_0, R_{a_0}, ..., R_{a_{p-2}}, E] with Nlog bits set, spread evenly, and L - N vacancies
V_0 placed after markers.  mu and mu_sel drop out of a fixed (N, Mk) sector.

usage: run_mobile.py NAME --p 2 --L 48 --N 30 --Nlog 4 --chis 64 128 256 [couplings]
  --N L gives the full-occupancy (vacancy frozen out) control.
  --model ref : exact-reference stand-in, free fermions run through the same code (not used here).
Writes data/obc/NAME.json after every chi.
Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-29.
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'codex', 'work_C88'))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import numpy as np
import bond13 as B
import cqdmrg as Q


def charges():
    return np.array([[B.charge(a), int(a[0] != 'V'), int(B.marker(a))] for a in B.LET])


def config(p, L, N, Nlog):
    ix = {a: i for i, a in enumerate(B.LET)}
    b = p + 1
    assert N % b == 0
    Mk = N // b
    nb = Mk * (p - 1)
    bits = np.zeros(nb, int)
    if Nlog:
        bits[np.floor((np.arange(Nlog) + 0.5) * nb / Nlog).astype(int)] = 1
    nv = L - N
    vac_after = np.diff(np.floor(np.arange(Mk + 1) * nv / Mk).astype(int))
    cfg = []
    for j in range(Mk):
        cell = [ix[('R', 0)]] + [ix[('R', int(x))] for x in bits[j * (p - 1):(j + 1) * (p - 1)]] + [ix[('e', 0)]]
        cfg += cell + [ix[('V', 0)]] * int(vac_after[j])
    assert len(cfg) == L
    return cfg, Mk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--p', type=int, default=2)
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--N', type=int, required=True)
    ap.add_argument('--Nlog', type=int, required=True)
    ap.add_argument('--D', type=float, default=0.0)
    ap.add_argument('--Delta', type=float, default=0.5)
    ap.add_argument('--g', type=float, default=0.3)
    ap.add_argument('--lam', type=float, default=0.3)
    ap.add_argument('--td', type=float, default=0.2)
    ap.add_argument('--K', type=float, default=1e4)
    ap.add_argument('--J', type=float, default=0.5)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--chis', type=int, nargs='+', default=[64, 128])
    ap.add_argument('--sweeps', type=int, default=3)
    ap.add_argument('--maxsweeps', type=int, default=8)
    ap.add_argument('--etol', type=float, default=1e-11)
    ap.add_argument('--noise', type=float, default=1e-3)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    P = B.params(n=a.p, D=a.D, Delta=a.Delta, g=a.g, td=a.td, K=a.K, musel=0.0)
    h = B.bond(P, lam=a.lam, J=a.J, t=a.t, mu=0.0)
    assert np.abs(h - h.conj().T).max() == 0
    cfg, Mk = config(a.p, a.L, a.N, a.Nlog)
    out = a.out or os.path.join(HERE, 'data', 'obc', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), Mk=Mk, results=[])

    def on_chi(dm, rec):
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
    Q.run(h, charges(), a.L, cfg, a.chis, sweeps=a.sweeps, maxsweeps=a.maxsweeps, etol=a.etol,
          noise0=a.noise, log=lambda s: print(a.name, s, flush=True), on_chi=on_chi)


if __name__ == '__main__':
    main()
