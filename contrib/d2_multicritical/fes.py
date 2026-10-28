#!/usr/bin/env python3
"""Finite-entanglement scaling with VUMPS (vumps.py, copied unchanged from
../d2_followup_220_1282/).  c_fit = 6 x slope of the cell-averaged S against
ln xi, xi = -n / ln|l1/l0| from the n-site cell transfer matrix.

usage: fes.py MODEL --n N --chis 8 12 16 24 32 48 64 [--mode damped] [--out ...]
MODEL is a control name (see controls()) or 'FAM:k=v,k=v,...' with FAM in
scan.params (e.g. 'kit:tau=1,delta=1,V=3,h=1.2').
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
import ring
from scan import bond_of


def embed_pair(hA, hB):
    """d = 4 site = (qubit A, qubit B); two-site term hA on (A1,A2) plus hB on (B1,B2).
    Ordering of the 16-dim space: (A1, B1, A2, B2), leftmost most significant."""
    I = np.eye(4)
    tA = np.kron(hA, I).reshape([2] * 8)   # (A1,A2,B1,B2 | A1',A2',B1',B2')
    tB = np.kron(I, hB).reshape([2] * 8)
    perm = [0, 2, 1, 3, 4, 6, 5, 7]       # -> (A1,B1,A2,B2 | ...)
    return (tA + tB).transpose(perm).reshape(16, 16)


def controls(name):
    X, Y, Z, I2 = ring.X, ring.Y, ring.Z, ring.I2
    ising = ring.bond(np.diag([-1., 0, 0]), [0, 0, -1.])           # v = 2, c = 1/2
    xx = ring.bond(np.diag([1., 1., 0]), [0, 0, 0])                  # v = 4, c = 1
    if name == 'ising':
        return ising, 2, 2, 0.5
    if name == 'xx':
        return xx, 2, 2, 1.0
    if name == 'cluster_ising':
        return -np.kron(np.kron(X, Z), X) + np.kron(np.kron(Y, Y), I2), 3, 2, 1.5
    if name == 'xx_pair':      # two decoupled XX chains as one d = 4 site
        return embed_pair(xx, xx), 2, 4, 2.0
    if name == 'ising_xxz':    # Ising (v = 2) x XXZ Delta = 0.5 (v = 2 pi sin(g)/g ~ 2.6): two velocities
        return embed_pair(ising, ring.bond(np.diag([0.5, 0.5, 0.25]), [0, 0, 0])), 2, 4, 1.5
    raise KeyError(name)


def model(name):
    if ':' in name:
        fam, rest = name.split(':', 1)
        kw = {s.split('=')[0]: float(s.split('=')[1]) for s in rest.split(',') if s}
        return bond_of(fam, **kw), 2, 2, None
    return controls(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64])
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    h, r, d, c = model(a.model)
    W = V.mpo_from_term(h, r, d=d)
    data = dict(model=a.model, c_exact=c, n=a.n, mode=a.mode, results=[])
    t0 = time.time()

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['S'] = [float(x) for x in rec['S']]
        if d == 2:
            rec['mag'] = V.local_expect(psi, list(V.PAULI)).tolist()
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        print(f"{a.model} n={a.n} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.mean(rec['S']):.5f} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s", flush=True)

    V.run([W] * a.n, a.n, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode,
          log=lambda rec: None, on_chi=save)
    S = np.array([np.mean(x['S']) for x in data['results']]); xi = np.array([x['xi'] for x in data['results']])
    if len(S) >= 4:
        print('c_fit (last 4):', 6 * np.polyfit(np.log(xi[-4:]), S[-4:], 1)[0])


if __name__ == '__main__':
    main()
