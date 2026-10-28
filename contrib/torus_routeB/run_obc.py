#!/usr/bin/env python3
"""Open-chain two-site DMRG (second estimator) for the task-1 chains and the task-2 effective
chains.  Reuses FDMRG of ../encoded_chain_numerics/fdmrg_md.py unchanged.

usage:
  run_obc.py model --case ii --s 1 --g 0.3 --L 32 --chis 32 64 128 192 --out f.json
  run_obc.py eff   --tag tuned --lam 0.2 --mu0 0.1 --L 32 --chis ... --out f.json
Saves the energy, every bond entropy S(l), the truncation error and <n_a(x)> per chi.

Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'encoded_chain_numerics'))
import numpy as np
import cells as C              # noqa: E402
import fdmrg_md                # noqa: E402
from fdmrg_md import FDMRG     # noqa: E402
from scipy import linalg as _la


def lanczos_lowest(op, k=1, which='SA', v0=None, tol=1e-7, ncv=24, maxiter=None, nrestart=6):
    """Restarted Lanczos (full reorthogonalization, warm start) replacing ARPACK in FDMRG.sweep,
    which needed thousands of matrix-vector products per update.  Same return signature as eigsh."""
    n = op.shape[0]
    v = np.asarray(v0, complex).copy() if v0 is not None else np.random.default_rng(0).normal(size=n) + 0j
    if np.linalg.norm(v) < 1e-14:
        v = np.random.default_rng(0).normal(size=n) + 0j
    v /= np.linalg.norm(v)
    m = 24
    for _ in range(nrestart):
        V = [v]
        al, be = [], []
        w = op.matvec(v)
        for j in range(m):
            a = np.vdot(V[j], w).real
            al.append(a)
            w = w - a * V[j] - (be[-1] * V[j - 1] if j > 0 else 0)
            for u in V:
                w = w - np.vdot(u, w) * u
            b = np.linalg.norm(w)
            if b < 1e-12 or j == m - 1:
                break
            be.append(b)
            V.append(w / b)
            w = op.matvec(V[-1])
        T = np.diag(al) + np.diag(be[:len(al) - 1], 1) + np.diag(be[:len(al) - 1], -1)
        e, s = _la.eigh(T)
        x = sum(s[i, 0] * V[i] for i in range(len(al)))
        x /= np.linalg.norm(x)
        res = np.linalg.norm(op.matvec(x) - e[0] * x)
        v = x
        if res < tol:
            break
    return np.array([e[0]]), x[:, None]


fdmrg_md.eigsh = lanczos_lowest
import fermions as F           # noqa: E402
import tcases as M             # noqa: E402


def build(a):
    if a.mode == 'model':
        k, h, inter = M.two_site_term(a.case, a.s, a.g, a.seed, no=bool(a.no))
        W = C.mpo_from_terms([16], [(0, 2, h)])
        return [16] * a.L, W * a.L
    from run_vumps import load_eff
    T1, T2, T3, Q2, al2 = load_eff(a.tag, HERE)
    T2, T3 = a.kappa * T2, a.kappa * T3                  # kappa: amplification of the lam^2 terms
    mu0 = a.mu0 * 4 * a.lam * float(np.mean(al2))
    sc = 1.0 / (a.lam * float(np.mean(al2)))             # energy unit: the mean first-order hopping
    if a.order >= 2:
        terms = [(0, 2, sc * (a.lam * T1 + a.lam ** 2 * T2 + mu0 * Q2)), (0, 3, sc * a.lam ** 2 * T3)]
    else:
        terms = [(0, 2, sc * (a.lam * T1 + mu0 * Q2))]
    W = C.mpo_from_terms([16], terms)
    return [16] * a.L, W * a.L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['model', 'eff'])
    ap.add_argument('--case', default='ii')
    ap.add_argument('--s', type=float, default=1.0)
    ap.add_argument('--g', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--no', type=int, default=0, help='model: 1 = Hartree-subtracted densities')
    ap.add_argument('--tag', default='tuned')
    ap.add_argument('--lam', type=float, default=0.2)
    ap.add_argument('--mu0', type=float, default=0.0)
    ap.add_argument('--order', type=int, default=2)
    ap.add_argument('--kappa', type=float, default=1.0)
    ap.add_argument('--L', type=int, default=32)
    ap.add_argument('--chis', type=int, nargs='+', default=[32, 64, 128])
    ap.add_argument('--sweeps', type=int, default=3)
    ap.add_argument('--maxsweeps', type=int, default=8)
    ap.add_argument('--etol', type=float, default=1e-7)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    dims, Ws = build(a)
    c, P = F.site_ops(4)
    nops = [x.conj().T @ x for x in c]
    data = dict(args=vars(a), recs=[])
    t0 = time.time()
    dm = FDMRG(Ws, dims, seed=a.seed + 1, chi0=a.chis[0])     # random start at the first chi (no empty-flavour trap)
    for chi in a.chis:
        Eprev = None
        for sw in range(a.maxsweeps):
            E, tr = dm.sweep(chi)
            if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < a.etol * a.L:
                break
            Eprev = E
        dens = []
        for x in range(a.L - 1):
            dens.append([dm.two_site(x, np.kron(n, np.eye(16))) for n in nops])
        rec = dict(chi=chi, E=dm.energy(), sweeps=sw + 1, trunc=tr, S=dm.S.tolist(), n=dens, time_s=time.time() - t0)
        data['recs'].append(rec)
        print(f"L={a.L} chi={chi} E={rec['E']:.10f} sweeps={sw+1} trunc={tr:.1e} Smid={dm.S[a.L//2-1]:.5f} "
              f"n={np.round(np.mean(dens, 0), 4).tolist()} t={rec['time_s']:.0f}s", flush=True)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)


if __name__ == '__main__':
    main()
