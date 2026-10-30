"""Exact ring ground energies (Lanczos) and the feasibility check of the bootstrap:
the ED ground-space state, with all its marginals, must satisfy every SDP constraint."""
import os
for _n in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_n, '1')
import argparse, json, time
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

from models import bonds


def ring_op(h, d, L):
    h4 = h.reshape(d, d, d, d)

    def mv(v):
        psi = v.reshape((d,) * L)
        out = np.zeros_like(psi)
        for i in range(L):
            j = (i + 1) % L
            t = np.moveaxis(psi, (i, j), (0, 1))
            t = np.tensordot(h4, t, axes=([2, 3], [0, 1]))
            out += np.moveaxis(t, (0, 1), (i, j))
        return out.ravel()
    D = d ** L
    return LinearOperator((D, D), matvec=mv, dtype=float)


def ground(h, d, L, nev=4):
    H = ring_op(h, d, L)
    nev = min(nev, d ** L - 2)
    w, V = eigsh(H, k=nev, which='SA', tol=1e-12)
    o = np.argsort(w); w = w[o]; V = V[:, o]
    deg = np.sum(w < w[0] + 1e-8)
    return w, V[:, :deg]


def marginal(V, d, L, sites):
    m = len(sites)
    rho = 0
    for v in V.T:
        psi = v.reshape((d,) * L)
        rest = [i for i in range(L) if i not in sites]
        t = np.transpose(psi, list(sites) + rest).reshape(d ** m, -1)
        rho = rho + t @ t.T
    return rho / V.shape[1]


def feasibility(model, L, k, s, mloc, momenta, pair_psd=True, opset='full', pair_s=None):
    """Evaluate every SDP constraint on the exact ring ground-space state."""
    from ring_boot import Ring, momenta_list, site_basis
    import cvxpy as cp
    B = {nm: (dd, hh) for nm, dd, hh in bonds()}
    d, h = B[model]
    w, V = ground(h, d, L)
    R = Ring(h, d, L, k, s, mloc, momenta_list(L, momenta), use_pair_psd=pair_psd, opset=opset, pair_s=pair_s)
    R.build()
    mats, _ = site_basis(d)
    cache = {}
    x = np.zeros(R.nvar)
    for c, i in R.reg.items():
        sites = [p for p, _ in c]
        key = tuple(sites)
        if key not in cache:
            cache[key] = marginal(V, d, L, sites)
        op = np.array([[1.0]])
        for _, m in c:
            op = np.kron(op, mats[m])
        x[i] = np.trace(cache[key] @ op)
    R.x.value = x
    res = dict(model=model, L=L, k=k, s=s, E0_per_site=float(w[0] / L),
               energy_from_x=float(np.ravel(R.energy.value)[0]), ground_degeneracy=int(V.shape[1]))
    worst_eq = 0.0; min_eig = np.inf
    for c in R.cons:
        if isinstance(c, cp.constraints.PSD):
            M = c.args[0].value
            min_eig = min(min_eig, float(np.linalg.eigvalsh((M + M.T) / 2)[0]))
        else:
            worst_eq = max(worst_eq, float(np.max(np.abs(c.args[0].value - c.args[1].value))) if hasattr(c, 'args') else 0.0)
    res['worst_equality'] = worst_eq
    res['min_psd_eigenvalue'] = min_eig
    res['nvar'] = R.nvar; res['n_eq'] = getattr(R, 'n_eq', 0)
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['energies', 'feasible'], required=True)
    p.add_argument('--model'); p.add_argument('--L', type=int)
    p.add_argument('--Ls', default='')
    p.add_argument('--k', type=int); p.add_argument('--s', type=int, default=1)
    p.add_argument('--mloc', type=int, default=None)
    p.add_argument('--momenta', default='all')
    p.add_argument('--opset', default='full')
    p.add_argument('--pair-s', type=int, default=None)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if 'SLURM_JOB_ID' not in os.environ:
        p.error('run inside a SLURM allocation')
    if a.mode == 'energies':
        B = {nm: (dd, hh) for nm, dd, hh in bonds()}
        d, h = B[a.model]
        rows = []
        for L in [int(x) for x in a.Ls.split(',')]:
            t = time.time(); w, V = ground(h, d, L, nev=4)
            rows.append(dict(model=a.model, L=L, E0=float(w[0]), gaps=[float(x - w[0]) for x in w],
                             deg=int(V.shape[1]), seconds=time.time() - t))
            print(json.dumps(rows[-1]), flush=True)
        json.dump(rows, open(a.out, 'w'))
    else:
        mloc = a.mloc if a.mloc is not None else a.k - 2
        r = feasibility(a.model, a.L, a.k, a.s, mloc, a.momenta, opset=a.opset, pair_s=a.pair_s)
        print(json.dumps(r), flush=True)
        json.dump(r, open(a.out, 'w'))


if __name__ == '__main__':
    main()
