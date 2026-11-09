#!/usr/bin/env python3
"""C177 (Claude Opus 5.5 agent, 2026-11-09): VUMPS (C166's vumps.py, unchanged module) for the full qubit chain at
the Fibonacci cusp, n-site cell, chi grown through a list, with checkpoint/restart as in C166's vumps_cells.py.

Changes relative to vumps_cells.py, all outside vumps.py:
- Bond from tune.json (or explicit generic terms), MPO fesx.mpo_cell with J = [[Jx, D, Jxz], [-D, Jy, 0], [Jxz, 0, Jz]],
  B = (hx, 0, hz), and a constant -shift per site (shift = Jz + hz) added to the on-site slot so that tolerances act on
  the low-energy scale.  Reported energies have the shift added back.
- mode 'vumps' (default here): the A_C and C eigenproblems are solved by Davidson with the diagonal of H_AC / H_C as
  preconditioner (dmrg2.davidson), because ARPACK stalls on the V/t2 scale separation.  mode 'damped' is C166's.
Per chi: e, S, xi, transfer spectrum, density <(1 - Z)/2> per cell site, bond energies, and the vector chirality.
usage: vumps_fc.py TAG --r 8 --n 2 --chis 16 32 ... --out f.json [--npz f.npz] [--jxpy x --jxz x --hx x]"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
import fesx as F
from dmrg2 import davidson

_orig_step = V.vumps_step


def _dav(apply, x0, diag, tol):
    shp = x0.shape
    x = x0.ravel()
    if not np.any(x):
        x = np.random.default_rng(0).normal(size=x.size) + 0j
    e, v, nmv, rn = davidson(lambda y: apply(y.reshape(shp)).ravel(), x.astype(complex), diag,
                             res_tol=max(tol, 1e-12), maxit=600)
    return e, v.reshape(shp)


def vumps_step_dav(psi, W, tol_eig, tol_env, cache, tau=None, krylov=16):
    if tau is not None:
        return _orig_step(psi, W, tol_eig, tol_env, cache, tau=tau, krylov=krylov)
    n = psi.n
    Ls, eL, _ = V.left_envs(psi.AL, psi.C, W, tol_env, cache.get('L'))
    Rs, eR, _ = V.right_envs(psi.AR, psi.C, W, tol_env, cache.get('R'))
    D = W[0].shape[0]
    cache['L'] = Ls[0][:, D - 1, :]
    cache['R'] = Rs[0][:, 0, :]
    AC, C = [], []
    for k in range(n):
        Rn = Rs[(k + 1) % n]
        Lk = Ls[k]
        la = np.einsum('awa->aw', Lk); rb = np.einsum('bvb->bv', Rn)
        dac = np.einsum('aw,wvs,bv->asb', la, np.einsum('wvss->wvs', W[k]), rb).ravel()
        _, ac = _dav(lambda x: V.hac(Lk, W[k], Rn, x), psi.AC[k], dac, tol_eig)
        rk = np.einsum('bvb->bv', Rs[k])
        dc = np.einsum('aw,bw->ab', la, rk).ravel()
        _, c = _dav(lambda x: V.hc(Lk, Rs[k], x), psi.C[k], dc, tol_eig)
        AC.append(ac / np.linalg.norm(ac)); C.append(c / np.linalg.norm(c))
    AL, AR = [], []
    err = 0.0
    for k in range(n):
        chil, d, chir = AC[k].shape
        Uac = V.polar_u(AC[k].reshape(chil * d, chir))
        Uc = V.polar_u(C[(k + 1) % n])
        AL.append((Uac @ Uc.conj().T).reshape(chil, d, chir))
        Uac2 = V.polar_u(AC[k].reshape(chil, d * chir))
        Uc2 = V.polar_u(C[k])
        AR.append((Uc2.conj().T @ Uac2).reshape(chil, d, chir))
        x1 = np.tensordot(AL[k], C[(k + 1) % n], ([2], [0]))
        x2 = np.tensordot(C[k], AR[k], ([1], [0]))
        for x in (x1, x2):
            ov = np.vdot(x, AC[k])
            ph = ov / abs(ov) if abs(ov) > 0 else 1.0
            err = max(err, float(np.linalg.norm(AC[k] - ph * x)))
    new = V.UMPS(AL, AR, C)
    new.AC = AC
    return new, 0.5 * (eL + eR) / n, err, (Ls, Rs)


def two_site(psi, k, O):
    n = psi.n
    k1 = (k + 1) % n
    A2 = np.tensordot(psi.AL[k], psi.AC[k1], ([2], [0]))
    t = np.einsum('uvst,astb->auvb', O.reshape(2, 2, 2, 2), A2)
    return complex(np.vdot(A2, t) / np.vdot(A2, A2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tag')
    ap.add_argument('--r', type=float, default=8.0)
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 16, 32])
    ap.add_argument('--maxit', type=int, default=300)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--mode', default='vumps')
    ap.add_argument('--tau0', type=float, default=0.5)
    ap.add_argument('--taumax', type=float, default=50.0, help='damped mode: largest imaginary-time step (in units of 1/energy)')
    ap.add_argument('--krylov', type=int, default=16)
    ap.add_argument('--hz', type=float, default=None)
    ap.add_argument('--mu', type=float, default=None, help='dressed one-particle mu_eff (tune.json linear map); overrides the tuned -2.2')
    ap.add_argument('--jxpy', type=float, default=0.0)
    ap.add_argument('--jxz', type=float, default=0.0)
    ap.add_argument('--hx', type=float, default=0.0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--npz', default=None)
    ap.add_argument('--init', default=None, help='npz of another run used as the starting state (its chi must equal chis[0])')
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    cp = [d for d in json.load(open(os.path.join(here, 'tune.json'))) if abs(d['r'] - a.r) < 1e-9][0]['couplings']
    Tn = [d for d in json.load(open(os.path.join(here, 'tune.json'))) if abs(d['r'] - a.r) < 1e-9][0]
    if a.mu is not None:
        cp['hz'] = cp['hz'] - (a.mu - Tn['target_mu']) * Tn['t2'] / Tn['dc0_dhz']
    if a.hz is not None:
        cp['hz'] = a.hz
    Jx = (1 + a.jxpy) / 2; Jy = (a.jxpy - 1) / 2
    J = np.array([[Jx, cp['D'], a.jxz], [-cp['D'], Jy, 0.0], [a.jxz, 0.0, cp['Jz']]])
    Bf = np.array([a.hx, 0.0, cp['hz']])
    shift = cp['Jz'] + cp['hz']
    W = F.mpo_cell([J] * a.n, [Bf] * a.n)
    D = W[0].shape[0]
    for k in range(a.n):
        W[k][0, D - 1] = W[k][0, D - 1] - shift * np.eye(2)
    V.vumps_step = vumps_step_dav
    P = V.PAULI
    hb = sum(J[p, q] * np.kron(P[p], P[q]) for p in range(3) for q in range(3)) + sum(Bf[p] * np.kron(P[p], np.eye(2)) for p in range(3))
    chir = np.kron(P[0], P[1]) - np.kron(P[1], P[0])
    xxyy = np.kron(P[0], P[0]) + np.kron(P[1], P[1])   # valley-polarization order parameter (odd under the J_x + J_y = 0 symmetry)
    data = dict(args=vars(a), couplings=dict(cp, Jx=Jx, Jy=Jy, Jxz=a.jxz, hx=a.hx), J=J.tolist(), B=Bf.tolist(),
                shift=shift, results=[])
    t0 = time.time()
    psi0 = None
    chis = list(a.chis)
    if a.npz and os.path.exists(a.npz):
        z = np.load(a.npz)
        nn = int(z['n'])
        psi0 = V.UMPS([z[f'AL{k}'] for k in range(nn)], [z[f'AR{k}'] for k in range(nn)], [z[f'C{k}'] for k in range(nn)])
        if os.path.exists(a.out):
            data['results'] = json.load(open(a.out))['results']
        done = {r['chi'] for r in data['results']}
        chis = [c for c in chis if c not in done and c > psi0.AL[0].shape[0]]
        if not chis:
            print('nothing to do', flush=True); return
        print('restart from', a.npz, 'chi', psi0.AL[0].shape[0], '->', chis, flush=True)
        _, _, _, envs = vumps_step_dav(psi0, W, 1e-6, 1e-8, {})
        psi0 = V.expand(psi0, W, envs, chis[0])

    if psi0 is None and a.init:
        z = np.load(a.init)
        nn = int(z['n'])
        assert nn == a.n
        psi0 = V.UMPS([z[f'AL{k}'] for k in range(nn)], [z[f'AR{k}'] for k in range(nn)], [z[f'C{k}'] for k in range(nn)])
        print('initial state from', a.init, 'chi', psi0.AL[0].shape[0], flush=True)

    def on_chi(psi, rec):
        rec['e'] = rec['e'] + shift
        rec['time_s'] = time.time() - t0
        lam = V.transfer_spectrum(psi, nev=12)
        rec['lam'] = [[float(z.real), float(z.imag)] for z in lam]
        rec['schmidt0'] = [float(s) for s in V.schmidt(psi, 0)[:40]]
        mag = V.local_expect(psi, [P[0], P[1], P[2]])
        rec['mag'] = mag.tolist()
        rec['n'] = [float((1 - mag[2, k]) / 2) for k in range(psi.n)]
        rec['bond_e'] = [two_site(psi, k, hb).real for k in range(psi.n)]
        rec['chir_z'] = [two_site(psi, k, chir).real for k in range(psi.n)]
        rec['xxyy'] = [two_site(psi, k, xxyy).real for k in range(psi.n)]
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'), default=float); os.replace(a.out + '.tmp', a.out)
        print(f"{a.tag} n={a.n} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 4).tolist()} xi={rec['xi']:.2f} n={np.round(rec['n'], 4).tolist()} "
              f"chir={np.round(rec['chir_z'], 4).tolist()} P={np.round(rec['xxyy'], 5).tolist()} arg1/pi={np.angle(complex(*rec['lam'][1])) / np.pi:.4f} "
              f"t={rec['time_s']:.0f}s", flush=True)
        if a.npz:
            dd = dict(n=psi.n)
            for k in range(psi.n):
                dd[f'AL{k}'] = psi.AL[k]; dd[f'AR{k}'] = psi.AR[k]; dd[f'C{k}'] = psi.C[k]
            np.savez(a.npz + '.tmp.npz', **dd); os.replace(a.npz + '.tmp.npz', a.npz)

    V.run(W, a.n, chis, d=2, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode, on_chi=on_chi, psi=psi0,
          log=lambda r: None, tau0=a.tau0, tau_max=a.taumax, krylov=a.krylov)
    print('done', a.out, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
