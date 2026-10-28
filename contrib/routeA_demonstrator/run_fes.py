#!/usr/bin/env python3
"""Finite-entanglement scaling with VUMPS for the demonstrator chain and its calibrations.

usage: run_fes.py NAME --model {micro,eff,xx} [--lam 0.1 --mu 0.55 --K 2 --kap 0.6 --D0 1 --D1 1 --Jp 2 --Jz -1]
                  [--Jxy J --Jzz J]  --n 3 --chis 16 24 32 ... [--seed 0] [--out data/fes/NAME.json]
Models
  micro : the full d = 10 bond of cell.py, n-site unit cell, random initial state (types not imposed)
  eff   : logical XXZ chain  H = Jxy (s+s- + s-s+) + Jzz ZZ  (d = 2), the effective chain
  xx    : XX chain  H = -(s+s- + s-s+)  (c = 1 calibration)
  cell  : the same microscopic Hamiltonian restricted to its period-3 type sector (exact, types are
          conserved), written as a chain of 32-dim cells P B P; energies are per cell (3 sites)
Records per chi: energy per site, gradient error, entropies on every bond of the cell, transfer-matrix
spectrum and xi (sites), per-site density of B letters, density of BB bonds, charge per site.
Uses ../d2_followup_220_1282/vumps.py.  Claude (Opus 5.5), 2026-10-27.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
sys.path.insert(0, HERE)
import numpy as np
import vumps as V
import cell as C


def model(a):
    if a.model == 'micro':
        h = C.bond(D0=a.D0, D1=a.D1, kap=a.kap, lam=a.lam, Jp=a.Jp, Jz=a.Jz, mu=a.mu, K=a.K)
        d = 10
        nB = np.diag([0.] * 2 + [1.] * 8)
        q = np.diag(C.charge_site())
        return h, d, [nB, q]
    if a.model == 'cell':
        # the period-3 type sector as a chain of 32-dim cells (P B P): H0 on site, lam V across cells
        H0, ports, code, qc = C.cell_ops(a.D0, a.D1, a.kap)
        R, L = ports['R'], ports['L']
        terms = [(a.lam * a.Jp, R['sp'], L['sm']), (a.lam * a.Jp, R['sm'], L['sp']), (a.lam * a.Jz, R['z'], L['z'])]
        D = 2 + len(terms)
        Wc = np.zeros((D, D, 32, 32))
        Wc[0, 0] = np.eye(32); Wc[D - 1, D - 1] = np.eye(32)
        Wc[0, D - 1] = H0 - a.mu * np.eye(32)
        for i, (cf, A, B) in enumerate(terms):
            Wc[0, 1 + i] = cf * A
            Wc[1 + i, D - 1] = B
        return Wc, 32, [np.diag(qc), code @ C.Z @ code.T]
    sp, sm, z = C.SP, C.SM, C.Z
    if a.model == 'eff':
        h = a.Jxy * (np.kron(sp, sm) + np.kron(sm, sp)) + a.Jzz * np.kron(z, z)
        return h, 2, [z]
    if a.model == 'xx':
        h = -(np.kron(sp, sm) + np.kron(sm, sp))
        return h, 2, [z]
    raise ValueError(a.model)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='micro')
    for k, v in dict(lam=0.1, mu=0.55, K=2.0, kap=0.6, D0=1.0, D1=1.0, Jp=2.0, Jz=-1.0, Jxy=-1.0, Jzz=0.0).items():
        ap.add_argument('--' + k, type=float, default=v)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32])
    ap.add_argument('--n', type=int, default=3)
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--stol', type=float, default=1e-8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tau_max', type=float, default=50.0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--switch_chi', type=int, default=None,
                    help='use --mode up to this chi, then damped imaginary-time steps from that state on')
    ap.add_argument('--state_dir', default=None, help='save the MPS after every chi (npz)')
    a = ap.parse_args()
    h, d, ops = model(a)
    W = h if a.model == 'cell' else V.mpo_from_term(h, 2, d=d)
    out = a.out or os.path.join(HERE, 'data', 'fes', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), results=[])
    t0 = time.time()
    hBB = None
    if a.model == 'micro':
        nB = np.diag([0.] * 2 + [1.] * 8)
        hBB = np.kron(nB, nB)

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['local'] = V.local_expect(psi, ops).tolist()
        if hBB is not None:
            bb = []
            for k in range(psi.n):
                AC = np.tensordot(psi.AL[k], psi.C[(k + 1) % psi.n], ([2], [0]))
                th = np.tensordot(AC, psi.AR[(k + 1) % psi.n], ([2], [0])).reshape(AC.shape[0], d * d, -1)
                nrm = np.vdot(th, th).real
                bb.append(float(np.einsum('asb,s,asb->', th.conj(), np.diag(hBB), th).real / nrm))
            rec['BB'] = bb
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} local={np.round(rec['local'], 4).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)

    def save_state(psi, rec):
        save(psi, rec)
        if a.state_dir:
            os.makedirs(a.state_dir, exist_ok=True)
            np.savez_compressed(os.path.join(a.state_dir, f"{a.name}_chi{rec['chi']}.npz"),
                                **{f'AL{k}': psi.AL[k] for k in range(psi.n)}, **{f'AR{k}': psi.AR[k] for k in range(psi.n)},
                                **{f'C{k}': psi.C[k] for k in range(psi.n)})

    if a.switch_chi is None:
        V.run([W] * a.n, a.n, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode,
              log=lambda rec: None, on_chi=save_state, stol=a.stol, tau_max=a.tau_max)
    else:
        lo = [c for c in a.chis if c <= a.switch_chi]
        hi = [c for c in a.chis if c > a.switch_chi]
        psi, _ = V.run([W] * a.n, a.n, lo, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode,
                       log=lambda rec: None, on_chi=save_state, stol=a.stol, tau_max=a.tau_max)
        # damped continuation: the first entry re-converges at the current chi, then chi grows
        V.run([W] * a.n, a.n, [lo[-1]] + hi, d=d, tol=a.tol, maxit=a.maxit, mode='damped', psi=psi,
              log=lambda rec: None, on_chi=lambda p, r: save_state(p, r) if r['chi'] in hi else None,
              stol=a.stol, tau_max=a.tau_max)


if __name__ == '__main__':
    main()
