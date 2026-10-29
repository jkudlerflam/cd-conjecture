"""Small, local verification only. Large finite chains require --cluster.

Models: XX, Heisenberg, XXZ(1/2), SU(3) swap. Bond energies are computed
in balanced charge sectors. These are numerical checks, not SDP certificates.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def basis(d, length):
    counts = [length // d + (i < length % d) for i in range(d)]
    states = []
    def rec(remaining, color, row):
        if color == d - 1:
            q = row.copy()
            q[list(remaining)] = color
            states.append(q)
            return
        for sites in itertools.combinations(remaining, counts[color]):
            q = row.copy()
            q[list(sites)] = color
            rec(tuple(i for i in remaining if i not in sites), color+1, q)
    rec(tuple(range(length)), 0, np.zeros(length, dtype=np.int8))
    states = np.array(states)
    powers = d ** np.arange(length, dtype=np.int64)
    codes = states @ powers
    ix = np.argsort(codes)
    return states[ix], powers, codes[ix]


def bonds(d, length, anisotropy):
    states, powers, codes = basis(d, length)
    dim = len(codes)
    ix = np.arange(dim)
    ops = []
    for j in range(length-1):
        a, b = states[:, j], states[:, j+1]
        unequal = a != b
        targets = codes + (b-a)*powers[j] + (a-b)*powers[j+1]
        rows = np.searchsorted(codes, targets[unequal])
        if d == 2:
            diag = anisotropy*np.where(unequal, -.25, .25)
            off = .5
        else:
            diag = (~unequal).astype(float)
            off = 1.0
        op = coo_matrix((np.r_[diag, np.full(unequal.sum(), off)],
                         (np.r_[ix, rows], np.r_[ix, ix[unequal]])),
                        shape=(dim, dim)).tocsr()
        ops.append(op)
    return ops


def ground(ops, w):
    mat = sum(a*b for a, b in zip(w, ops))
    dim = mat.shape[0]
    if dim <= 8:
        ev, vec = np.linalg.eigh(mat.toarray())
        e, state = ev[0], vec[:, 0]
    else:
        v0 = np.cos(np.arange(dim)*np.sqrt(2.0))
        ev, vec = eigsh(mat, k=1, which='SA', tol=1e-12, v0=v0)
        e, state = ev[0], vec[:, 0]
    energies = np.array([state @ (op @ state) for op in ops])
    return float(e), energies, float(np.linalg.norm(mat@state-e*state))


def positive_interval_lp(w, g):
    n = len(w)
    columns, costs = [], []
    for left in range(n):
        for right in range(left+1, n+1):
            col = np.zeros(n)
            col[left:right] = 1
            columns.append(col)
            costs.append(g[right-left])
    res = linprog(costs, A_eq=np.array(columns).T, b_eq=w,
                  bounds=(0, None), method='highs')
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.fun), float(np.sum(res.x))


def small_checks():
    models = [('XX', 2, 0., -1/np.pi, 1.),
              ('Heisenberg', 2, 1., .25-np.log(2), np.pi/2),
              ('XXZ_half', 2, .5, -.375, 3*np.sqrt(3)/4),
              ('SU3', 3, None, 1-np.log(3)-np.pi/(3*np.sqrt(3)), 4*np.pi/3)]
    out = []
    for name, d, anisotropy, e0, A in models:
        max_n = 12 if d == 2 else 9
        open_defects = {}
        cached = {}
        for n in range(2, max_n+1):
            ops = bonds(d, n, anisotropy)
            en, _, residual = ground(ops, np.ones(n-1))
            open_defects[n-1] = (n-1)*e0-en
            cached[n] = ops
        for n in range(4 if d == 2 else 6, max_n+1, d):
            x = np.arange(1, n, dtype=float)
            for shape in ['quadratic', 'ramp']:
                w = x*(n-x)/n if shape == 'quadratic' else x
                en, local, residual = ground(cached[n], w)
                deficit = e0*w.sum()-en
                vals = np.unique(w)
                prev, lc = 0., 0.
                for t in vals:
                    lc += (t-prev)*open_defects[int(np.sum(w>=t))]
                    prev = t
                lp, mass = positive_interval_lp(w, open_defects)
                factor = (x*(n-x)/n)**2 if shape == 'quadratic' else x*x
                raw = factor*(e0-local)
                kernel = np.array([1, 2, 1])/4 if d==2 else np.array([1, 2, 3, 2, 1])/9
                r = len(kernel)//2
                smooth = np.convolve(e0-local, kernel, mode='valid')*factor[r:-r]
                M = .5 if d == 2 else 1.
                out.append(dict(model=name, n=n, dimension=cached[n][0].shape[0],
                                shape=shape, e0=e0, A=A, M=M, energy=en,
                                defect=deficit, layer_cake_allowance=lc,
                                optimized_interval_allowance=lp,
                                interval_coefficient_mass=mass, max_weight=float(w.max()),
                                raw_scaled=raw.tolist(), smoothed_scaled=smooth.tolist(),
                                predicted_smooth=A/(24*np.pi),
                                sharp_conjectured_smooth=M*np.log2(d)/24,
                                energy_sum_residual=float(abs(w@local-en)),
                                eigenvector_residual=residual,
                                min_uniform_block_defect=float(min(open_defects.values()))))
                assert deficit <= lp + 1e-8
                assert lp <= lc + 1e-8
                assert mass >= w.max()-1e-8
    return out


def large_gaussian():
    if 'SLURM_JOB_ID' not in os.environ:
        raise RuntimeError('Large finite-chain computations require a SLURM job.')
    out = []
    for n in [512, 1024, 2048, 4096, 8192]:
        x = np.arange(1, n, dtype=float)
        for shape in ['quadratic', 'ramp']:
            w = x*(n-x)/n if shape=='quadratic' else x
            vals, vecs = eigh_tridiagonal(np.zeros(n), w/2,
                                         select='i', select_range=(0, n//2-1))
            local = np.sum(vecs[:-1]*vecs[1:], axis=1)
            factor = (x*(n-x)/n)**2 if shape=='quadratic' else x*x
            smooth = np.convolve(-1/np.pi-local, [1/4, 1/2, 1/4], mode='valid')
            out.append(dict(n=n, shape=shape, defect=float(-w.sum()/np.pi-vals.sum()),
                            local_energy_sum_residual=float(abs(w@local-vals.sum())),
                            sample=[dict(x=int(j), raw_scaled=float(factor[j-1]*(-1/np.pi-local[j-1])),
                                         smoothed_scaled=float(factor[j-1]*smooth[j-2]))
                                    for j in [17, 65, n//8, n//4, n//2] if 2<=j<n-1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cluster', action='store_true')
    ap.add_argument('--output', default='small_checks_results.json')
    args = ap.parse_args()
    start = time.perf_counter()
    rows = large_gaussian() if args.cluster else small_checks()
    data = dict(mode='cluster' if args.cluster else 'small',
                elapsed_seconds=time.perf_counter()-start, rows=rows)
    Path(args.output).write_text(json.dumps(data, indent=2)+'\n')
    for row in rows:
        print({k: v for k,v in row.items() if k not in ['raw_scaled', 'smoothed_scaled']})
    print('elapsed_seconds', data['elapsed_seconds'])


if __name__ == '__main__':
    main()
