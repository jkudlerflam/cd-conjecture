"""Small numerical checks. Run in this directory with OPENBLAS_NUM_THREADS=1.

All numerical results are diagnostics, not certified spectral enclosures.
"""
import argparse
import itertools
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


@lru_cache(None)
def structure(L, n):
    states = [sum(1 << x for x in c) for c in itertools.combinations(range(L), n)]
    index = {s: i for i, s in enumerate(states)}
    rows, cols, dirs, diag = [], [], [], []
    for i, s in enumerate(states):
        zz = 0.0
        for x in range(L):
            y = (x + 1) % L
            bx, by = (s >> x) & 1, (s >> y) & 1
            zz += 0.25 if bx == by else -0.25
            if bx != by:
                rows.append(index[s ^ (1 << x) ^ (1 << y)])
                cols.append(i)
                dirs.append(1 if bx else -1)
        diag.append(zz)
    return np.array(rows), np.array(cols), np.array(dirs), np.array(diag)


def matrix(L, n, delta, a=np.pi, J=1.0):
    rows, cols, dirs, diag = structure(L, n)
    size = len(diag)
    d = np.arange(size)
    data = np.r_[J * 0.5 * np.exp(1j * a * dirs), J * delta * diag]
    H = coo_matrix((data, (np.r_[rows, d], np.r_[cols, d])), shape=(size, size)).tocsr()
    if abs(np.sin(a)) < 1e-14:
        H = H.real
    return H


@lru_cache(10000)
def energy(L, n, delta, a=np.pi, J=1.0):
    H = matrix(L, n, delta, a, J)
    if H.shape[0] <= 100:
        return float(np.linalg.eigvalsh(H.toarray())[0])
    w, v = eigsh(H, k=1, which='SA', tol=2e-12,
                 v0=np.linspace(1.0, 1.1, H.shape[0]), maxiter=10000)
    resid = float(np.linalg.norm(H @ v[:, 0] - w[0] * v[:, 0]))
    if resid > 2e-9:
        raise RuntimeError(f'ED residual {resid} for {L,n,delta,a}')
    return float(w[0])


def bond_M(J, jz, B=0.0):
    ev = np.array([jz / 4 - B / 2, jz / 4 + B / 2,
                   -jz / 4 - J / 2, -jz / 4 + J / 2])
    return float(np.ptp(ev) / 2)


def twist_check(maxL=10, points=25):
    worst = (-1.0, None)
    minimum_error = 0.0
    count = 0
    for L in range(2, maxL + 1):
        for n in range(1, L // 2 + 1):
            for delta in [-1.001, -1, -0.999999, -0.999, -0.99, -0.9, 0, .5, 1, 2, 3]:
                base = energy(L, n, delta)
                for fluxdiff in np.linspace(0, np.pi, points):
                    val = energy(L, n, delta, np.pi + float(fluxdiff) / L)
                    ratio = (val - base) / (np.pi**2 / (4 * L))
                    minimum_error = min(minimum_error, val - base)
                    exact_bound = L * (1 - np.cos(fluxdiff / L)) / 2
                    assert val - base <= exact_bound + 2e-9, (L, n, delta, fluxdiff)
                    if ratio > worst[0]:
                        worst = (ratio, [L, n, delta, float(fluxdiff)])
                    count += 1
        print(f'twist L={L} done', flush=True)
    return dict(samples=count, worst_ratio=worst[0], at=worst[1],
                minimum_energy_difference=minimum_error, maxL=maxL, points=points)


def sanity_check(maxL=8):
    rows = []
    for jz, dm in itertools.product([-.95, -.5, 0, .5, 1, 1.5], [0, .3]):
        J, theta = float(np.hypot(1, dm)), float(np.arctan2(dm, 1))
        delta = jz / J
        for L in range(2, maxL + 1, 2):
            small = np.array([energy(L, n, delta, theta, J) for n in range(L // 2 + 1)])
            large = np.array([energy(2 * L, n, delta, theta, J) for n in range(L + 1)])
            for B in [0, .4, .8]:
                n = int(np.argmin(small - B * (L / 2 - np.arange(len(small)))))
                EL = small[n] - B * (L / 2 - n)
                E2 = np.min(large - B * (L - np.arange(len(large))))
                defect = E2 - 2 * EL
                raw = large[2 * n] - 2 * small[n]
                opt = energy(2 * L, 2 * n, delta, np.pi, J) - 2 * energy(L, n, delta, np.pi, J)
                M = bond_M(J, jz, B)
                assert defect <= raw + 2e-9
                assert defect <= opt + np.pi**2 * M / (4 * L) + 2e-9
                rows.append(dict(Jperp=1, DM=dm, Jz=jz, B=B, L=L, n=n, M=M,
                                 defect=float(defect), raw_sector=float(raw), optimized_sector=float(opt),
                                 field_slack=float(raw-defect),
                                 combined_slack=float(opt+np.pi**2*M/(4*L)-defect),
                                 Q=float(L*defect/M)))
        print(f'sanity Jz={jz}, DM={dm} done', flush=True)
    return rows


def free_check():
    err = 0.0
    for L in range(2, 11):
        for n in range(L // 2 + 1):
            expected = -np.sin(np.pi*n/L)/np.sin(np.pi/L)
            err = max(err, abs(energy(L, n, 0.0)-expected))
    assert err < 2e-10
    return dict(max_energy_error=err, maxL=10)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['free', 'twist', 'sanity'], required=True)
    p.add_argument('--maxL', type=int)
    p.add_argument('--points', type=int, default=25)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    if args.mode == 'free':
        result = free_check()
    elif args.mode == 'twist':
        result = twist_check(args.maxL or 10, args.points)
    else:
        result = sanity_check(args.maxL or 8)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result if isinstance(result, dict) else dict(
        samples=len(result), max_Q=max(x['Q'] for x in result),
        min_field_slack=min(x['field_slack'] for x in result),
        min_combined_slack=min(x['combined_slack'] for x in result)), indent=2))
