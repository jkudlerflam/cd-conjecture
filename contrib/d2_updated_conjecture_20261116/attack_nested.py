"""Finite-ring attack on a proposed c <= 2 bound; no central-charge inference.

Requires numpy only. Occupation n=(1-Z)/2, periodic spin boundary conditions.
H(eps)=sum nn + eps sum (n_left-n_right) X + eps**2(-T1+2N).
The last term cancels the second-order hopping and self-energy, leaving
2 eps**2 sum n_j n_(j+2). We inspect the next low-energy manifold, P2,
whose particles have separation at least three. All number sectors are kept.
The optional --pair adds adjacent-pair creation of size pair*eps**2 and
retains all of its effects during the same full-space extraction.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np


def operators(L):
    dim = 1 << L
    E = np.zeros(dim)
    D2 = np.zeros(dim)
    F = np.zeros((dim, dim))
    W = np.zeros_like(F)
    for s in range(dim):
        n = [(s >> j) & 1 for j in range(L)]
        E[s] = sum(n[j] * n[(j+1) % L] for j in range(L))
        D2[s] = sum(n[j] * n[(j+2) % L] for j in range(L))
        W[s, s] = 2 * sum(n)
        for j in range(L):
            F[s ^ (1 << j), s] += n[(j-1) % L] - n[(j+1) % L]
            k = (j+1) % L
            if n[j] != n[k]:
                W[s ^ (1 << j) ^ (1 << k), s] -= 1
    assert np.max(abs(F-F.T)) == 0
    assert np.max(abs(W-W.T)) == 0
    return E, D2, F, W


def run(L, eps_values, pair=0.0):
    E, D2, F, W = operators(L)
    # A second invisible direction of size eps**2: adjacent-pair creation.
    # It leaves the order-eps**2 blockade unchanged, but can affect order eps**4.
    if pair != 0:
        for s in range(1 << L):
            for j in range(L):
                k = (j+1)%L
                if ((s >> j)&1) == ((s >> k)&1):
                    W[s^(1 << j)^(1 << k),s] += pair
    P = np.flatnonzero(E == 0)
    Q = np.flatnonzero(E > 0)
    P2 = np.flatnonzero((E == 0) & (D2 == 0))
    A = F[np.ix_(Q, P)]
    K2 = W[np.ix_(P, P)] - A.T @ (A / E[Q, None])
    err2 = float(np.max(abs(K2 - np.diag(2*D2[P]))))
    assert err2 < 1e-12
    A2 = F[np.ix_(Q, P2)] / E[Q, None]
    B2 = W[np.ix_(Q, P2)]
    K3 = A2.T @ F[np.ix_(Q, Q)] @ A2 - A2.T @ B2 - B2.T @ A2
    err3 = float(np.max(abs(K3)))
    pos = {s: i for i, s in enumerate(P2)}
    N = np.array([int(s).bit_count() for s in P2])
    rows = []
    for eps in eps_values:
        H = np.diag(E) + eps*F + eps**2*W
        evals, vecs = np.linalg.eigh(H)
        m = len(P2)
        A = vecs[P2, :m]
        u, sv, vh = np.linalg.svd(A)
        polar = u @ vh
        He = (polar * evals[:m]) @ polar.T
        norm = eps**4
        vals = {'eps': eps, 'low_band_gap_over_eps2': float((evals[m]-evals[m-1])/eps**2),
                'min_overlap': float(min(sv)), 'norm_over_eps4': float(np.max(abs(He))/norm),
                'number_mixing_over_eps4': float(np.max(abs(He * (N[:,None] != N[None,:])))/norm),
                'onsite_over_eps4': float((He[pos[1],pos[1]]-He[pos[0],pos[0]])/norm)}
        vals['hops_over_eps4'] = {str(r): float(He[pos[1 << r],pos[1]]/norm) for r in range(1,L//2+1)}
        vals['pair_interactions_over_eps4'] = {
            str(r): float((He[pos[1+(1 << r)],pos[1+(1 << r)]]-He[pos[1],pos[1]]
                          -He[pos[1 << r],pos[1 << r]]+He[pos[0],pos[0]])/norm)
            for r in range(3,L//2+1)}
        largest = np.unravel_index(np.argmax(abs(He * (N[:,None] != N[None,:]))),He.shape)
        vals['largest_mixing_states'] = [int(P2[i]) for i in largest]
        rows.append(vals)
        print(json.dumps({'L':L, 'pair':pair, **vals}), flush=True)
    return {'L':L, 'pair':pair, 'dim_P':len(P), 'dim_P2':len(P2), 'second_order_error':err2,
            'third_order_P2_norm':err3, 'rows':rows}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--sizes', type=int, nargs='+', default=[8,9,10])
    p.add_argument('--eps', type=float, nargs='+', default=[0.10,0.075,0.05,0.035])
    p.add_argument('--pair', type=float, nargs='+', default=[0.0])
    p.add_argument('--output', default=str(Path(__file__).with_name('nested_results.json')))
    args = p.parse_args()
    start = time.perf_counter()
    results = [run(L,args.eps,pair) for pair in args.pair for L in args.sizes]
    out = {'purpose':'Nested-blockade diagnostic, not a CFT measurement',
           'numpy_version':np.__version__, 'results':results, 'runtime_seconds':time.perf_counter()-start}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
