#!/usr/bin/env python3
"""Finite-entanglement scaling off half filling with VUMPS (one-cell unit cell, infinite chain).

micro     : the period-3 type sector of the full chain as 32-dim cells P B P (exact, as in
            ../routeA_demonstrator/run_fes.py --model cell), plus a one-site field -hf Z_c on the B letter.
logical_x : logical XXZ chain with the exact two-cell couplings at this lam, plus -hf Z.
logical   : second-order logical chain, plus -hf Z.
The field hf is set so that the lam-dressed logical chain (logical_x) sits at down-spin density 1/3 by the
Bethe ansatz; the filling actually reached is measured (<Q_cell> = 1/2 - n).
Uses ../d2_followup_220_1282/vumps.py.
usage: vumps_fes.py NAME --model micro --lam 0.1 --rho 2 --chis 8 12 16 24 32 48 64
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
sys.path.insert(0, HERE)
import numpy as np
import vumps as V
import chain as CH
import bethe as BA
C = CH.C


def field_for(lam, rho):
    Jp, Jz = CH.couplings(rho)
    Jx, Jzx, _ = CH.exact_two_cell(lam, Jp, Jz)
    Jeff = 2 * abs(Jx)
    r = BA.at_density(2 * Jzx / abs(Jx), 1 / 3)
    return r['h'] * Jeff / 2          # -h Sz = -(h/2) Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='micro')
    ap.add_argument('--lam', type=float, default=0.1)
    ap.add_argument('--rho', type=float, default=2.0)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48])
    ap.add_argument('--maxit', type=int, default=1500)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--stol', type=float, default=1e-8)
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--tau_max', type=float, default=500.0)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    Jp, Jz = CH.couplings(a.rho)
    hf = field_for(a.lam, a.rho)
    if a.model == 'micro':
        H0, ports, code, qc = C.cell_ops(1.0, 1.0, 0.6)
        R, L = ports['R'], ports['L']
        terms = [(a.lam * Jp, R['sp'], L['sm']), (a.lam * Jp, R['sm'], L['sp']), (a.lam * Jz, R['z'], L['z'])]
        Zc = C.kron(C.I2, C.I2, C.Z, C.I2, C.I2)
        D = 2 + len(terms)
        W = np.zeros((D, D, 32, 32))
        W[0, 0] = np.eye(32); W[D - 1, D - 1] = np.eye(32)
        W[0, D - 1] = H0 - hf * Zc
        for i, (cf, A, B) in enumerate(terms):
            W[0, 1 + i] = cf * A
            W[1 + i, D - 1] = B
        d, ops = 32, [np.diag(qc)]
    else:
        if a.model == 'logical':
            Jxy, Jzz, _ = CH.second_order(Jp, Jz); Jxy, Jzz = Jxy * a.lam ** 2, Jzz * a.lam ** 2
        else:
            Jxy, Jzz, _ = CH.exact_two_cell(a.lam, Jp, Jz)
        h = (Jxy * (np.kron(C.SP, C.SM) + np.kron(C.SM, C.SP)) + Jzz * np.kron(C.Z, C.Z)
             - 0.5 * hf * (np.kron(C.Z, np.eye(2)) + np.kron(np.eye(2), C.Z)))
        W = V.mpo_from_term(h, 2, d=2)
        d, ops = 2, [0.5 * C.Z]
    out = a.out or os.path.join(HERE, 'data', 'fes', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), hf=hf, results=[])
    t0 = time.time()

    def on_chi(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['Qcell'] = float(V.local_expect(psi, ops)[0, 0])
        rec['n'] = 0.5 - rec['Qcell']
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.13f} err={rec['err']:.1e} it={rec['iters']} S={rec['S'][0]:.5f} "
              f"xi={rec['xi']:.3f} n={rec['n']:.5f} t={rec['time_s']:.0f}s", flush=True)
    V.run([W], 1, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode, log=lambda r: None,
          on_chi=on_chi, stol=a.stol, tau_max=a.tau_max)


if __name__ == '__main__':
    main()
