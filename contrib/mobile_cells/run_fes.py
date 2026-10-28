#!/usr/bin/env python3
"""Finite-entanglement scaling (VUMPS, one-site cell) for the spinful mobile selector and controls.

usage: run_fes.py NAME --model {cells,xx,xx2,txx} [--t 1 --u 1 --mu -1 --J 1 --Jz 0 --K 8] [--n 1]
                  --chis 16 24 32 48 ... [--mode damped] [--maxit 400] [--out data/fes/NAME.json]
Models
  cells : the d = 6 chain of ed.py (letters L, R, #up, #dn, Xup, Xdn), bond ed.dense_bond
  xx    : spin-1/2 XX chain (c = 1), calibration
  xx2   : two decoupled XX chains with velocities 1 and 0.4, as one d = 4 site (c = 2), calibration
  cells2: two-flavour payload of ed2.py (d = 10; --J is J_a, plus --Jb --V --hb), c = 3 expected
  txx   : d = 3 t-XX chain (hole, up, dn): hard-core hopping t, XX exchange J between adjacent
          particles, chemical potential mu per particle (c = 2 Luttinger liquid, Ogata-Shiba setting)
Writes after every chi: energy, gradient error, entropy, transfer-matrix eigenvalues (xi), local
densities (markers, heads, illegal pairs), <S^z>, <S^x>, run time.
Uses ../d2_followup_220_1282/vumps.py (numpy + scipy).  Claude (Opus 5.5), 2026-10-27.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
sys.path.insert(0, HERE)                 # a vumps.py copied next to this file (cluster) wins
import numpy as np
import vumps as V
from ed import dense_bond


def model(a):
    if a.model == 'cells':
        h = dense_bond(t=a.t, u=a.u, mu=a.mu, J=a.J, Jz=a.Jz, K=a.K, dJ=a.dJ)
        d = 6
        n = np.zeros((4, d, d))
        n[0][[2, 3], [2, 3]] = 1          # marker
        n[1][[4, 5], [4, 5]] = 1          # head
        n[2][[2, 4], [2, 4]] = 0.5; n[2][[3, 5], [3, 5]] = -0.5     # S^z
        sx = np.zeros((d, d)); sx[2, 3] = sx[3, 2] = sx[4, 5] = sx[5, 4] = 0.5
        n[3] = sx
        return h, d, n, 1.0 + 1.0
    if a.model == 'cells2':
        from ed2 import dense_bond2
        h = dense_bond2(t=a.t, u=a.u, mu=a.mu, Ja=a.J, Jb=a.Jb, V=a.V, ha=0.0, hb=a.hb, K=a.K)
        d = 10
        n = np.zeros((4, d, d))
        for x in range(2, 10):
            sa, sb = ((x - 2) % 4) >> 1, ((x - 2) % 4) & 1
            n[0][x, x] = 1.0 if x < 6 else 0.0            # marker
            n[1][x, x] = 1.0                              # wall
            n[2][x, x] = sa - 0.5                         # S^z_a
            n[3][x, x] = sb - 0.5                         # S^z_b
        return h, d, n, 3.0
    if a.model == 'xx':
        sp = np.array([[0, 1], [0, 0]]); sm = sp.T
        h = -0.5 * (np.kron(sp, sm) + np.kron(sm, sp))       # ferro sign: one-site cell suffices
        sz = np.diag([0.5, -0.5]); sx = 0.5 * (sp + sm)
        return h, 2, np.array([sz, sx]), 1.0
    if a.model == 'xx2':
        sp = np.array([[0, 1], [0, 0]]); sm = sp.T; I = np.eye(2)
        hx = -0.5 * (np.kron(sp, sm) + np.kron(sm, sp))
        # site = (a, b); bond = hx on the a chain + 0.4 hx on the b chain.  Build on the
        # ordering (a1 a2 b1 b2), then permute to (a1 b1 a2 b2).
        h4 = np.kron(hx, np.eye(4)) + 0.4 * np.kron(np.eye(4), hx)
        T = h4.reshape([2] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7)
        h = T.reshape(16, 16)
        sz = np.kron(np.diag([0.5, -0.5]), I)
        return h, 4, np.array([sz]), 2.0
    if a.model == 'txx':
        d = 3   # 0 hole, 1 up, 2 dn
        h = np.zeros((9, 9))
        ix = lambda x, y: x * 3 + y
        for s in (1, 2):
            h[ix(0, s), ix(s, 0)] += -a.t; h[ix(s, 0), ix(0, s)] += -a.t
        h[ix(1, 2), ix(2, 1)] += a.J / 2; h[ix(2, 1), ix(1, 2)] += a.J / 2
        for x in range(3):
            for y in range(3):
                h[ix(x, y), ix(x, y)] += 0.5 * a.mu * ((x > 0) + (y > 0))
        n = np.zeros((2, 3, 3)); n[0][[1, 2], [1, 2]] = 1; n[1][1, 1] = .5; n[1][2, 2] = -.5
        return h, d, n, 2.0
    raise ValueError(a.model)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='cells')
    for k, v in dict(t=1.0, u=1.0, mu=-1.0, J=1.0, Jz=0.0, K=8.0, dJ=0.0, Jb=0.35, V=0.3, hb=0.05).items():
        ap.add_argument('--' + k, type=float, default=v)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64])
    ap.add_argument('--n', type=int, default=1, help='unit cell (3 for the static period-3 crystal)')
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--stol', type=float, default=1e-7)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--state_dir', default=None, help='save the MPS after every chi here (npz, not committed)')
    ap.add_argument('--load', default=None, help='npz state to start from; its chi must equal the first --chis entry')
    a = ap.parse_args()
    h, d, ops, c_exp = model(a)
    W = V.mpo_from_term(h, 2, d=d)
    out = a.out or os.path.join(HERE, 'data', 'fes', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), c_expected=c_exp, results=[])
    t0 = time.time()

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['local'] = V.local_expect(psi, list(ops)).ravel().tolist()
        if a.model == 'cells':
            # density of illegal pairs (should be ~0): <Q> on one bond, from A_L C A_R
            n = psi.n
            AC = np.tensordot(psi.AL[0], psi.C[1 % n], ([2], [0]))
            th = np.tensordot(AC, psi.AR[1 % n], ([2], [0])).reshape(AC.shape[0], d * d, -1)
            nrm = np.vdot(th, th).real
            Qd = np.diag(dense_bond(t=0, u=0, mu=0, J=0, Jz=0, K=1.0))
            rec['illegal'] = float(np.einsum('asb,s,asb->', th.conj(), Qd, th).real / nrm)
        data['results'].append(rec)
        if a.state_dir:
            os.makedirs(a.state_dir, exist_ok=True)
            np.savez_compressed(os.path.join(a.state_dir, f"{a.name}_chi{rec['chi']}.npz"),
                                **{f'AL{k}': psi.AL[k] for k in range(psi.n)}, **{f'AR{k}': psi.AR[k] for k in range(psi.n)},
                                **{f'C{k}': psi.C[k] for k in range(psi.n)})
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} local={np.round(rec['local'], 5).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)

    psi0 = None
    if a.load:
        z = np.load(a.load)
        psi0 = V.UMPS([z[f'AL{k}'] for k in range(a.n)], [z[f'AR{k}'] for k in range(a.n)], [z[f'C{k}'] for k in range(a.n)])
        if os.path.exists(out):
            data['results'] = json.load(open(out))['results']
    V.run([W] * a.n, a.n, a.chis, d=d, psi=psi0, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode,
          log=lambda rec: None, on_chi=save, stol=a.stol)


if __name__ == '__main__':
    main()
