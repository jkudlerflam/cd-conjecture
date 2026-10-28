#!/usr/bin/env python3
"""Rank-k symmetric SU(N) chains (h = sum_j H_{k-j} P_j, ||h|| = H_k/2, the C14 bond divided by 2
and shifted): ED data against the Bethe/WZW prediction c = k(N^2-1)/(N+k), v = pi/N,
e_0 = H_k - (1/N) sum_{j=1}^k [psi(1+(j-1)/N) - psi(j/N)] (C14 eq. 5 in this normalization).

Per ring length: A_e0(L) = -6 L (E_0 - e_0 L)/pi (exact e_0), v(L) = L gap(k0 +- 2pi/L)/(2pi),
c_half between consecutive vacuum lengths. Fits: A_cas over all vacuum lengths (e_0 free) and
the linear-in-1/L^2 extrapolation of A_e0(L) and v(L) over the two largest lengths.
  python3 analyze_symrep.py data/symrep_local.jsonl data/symrep_big.jsonl
"""
import sys
import numpy as np
from scipy.special import digamma
from analyze import load, vel
import core


def theory(N, k):
    Hk = sum(1.0 / j for j in range(1, k + 1))
    c = k * (N * N - 1) / (N + k)
    v = np.pi / N
    e0 = Hk - sum(digamma(1 + (j - 1) / N) - digamma(j / N) for j in range(1, k + 1)) / N
    return dict(c=c, v=v, A=c * v, hn=Hk / 2, e0=e0, R=c * v / (2 * np.pi * Hk / 2))


def extrap(Ls, ys):
    x = 1 / np.asarray(Ls[-2:], float) ** 2
    y = np.asarray(ys[-2:], float)
    return y[1] - (y[1] - y[0]) / (x[1] - x[0]) * x[1]


def main(paths):
    recs = load(paths)
    rows = []
    for (model, par), rs in sorted(recs.items()):
        if model != 'symrep':
            continue
        N, k = int(par[0]), int(par[1])
        th = theory(N, k)
        Ls = sorted(L for L in rs if L % N == 0 and L >= 2 * (N == 2) + 2)
        if not Ls:
            continue
        hn = rs[Ls[0]]['hnorm']
        A_e0 = [-6 * L * (rs[L]['E0'] - th['e0'] * L) / np.pi for L in Ls]
        v = [vel(rs[L]) for L in Ls]
        ch = [core.c_half(rs[a]['S_half'], a, rs[b]['S_half'], b) for a, b in zip(Ls[:-1], Ls[1:])]
        print('N=%d k=%d d=%d  ||h||=%.4f (H_k/2=%.4f)  theory: c=%.4f v=%.4f A=%.4f e0=%.6f R=%.4f'
              % (N, k, round(np.prod([N + i for i in range(k)]) / np.prod(range(1, k + 1))), hn, th['hn'],
                 th['c'], th['v'], th['A'], th['e0'], th['R']))
        for i, L in enumerate(Ls):
            print('   L=%2d E0/L=%.6f k0=%d  A_e0=%.4f (%.3f of theory)  v=%.4f (%.3f)  c_half=%s'
                  % (L, rs[L]['E0'] / L, rs[L]['k0'], A_e0[i], A_e0[i] / th['A'], v[i], v[i] / th['v'],
                     ('%.3f' % ch[i - 1]) if i else '-'))
        out = dict(N=N, k=k, Ls=Ls, A_e0=A_e0[-1], v=v[-1])
        if len(Ls) >= 2:
            out.update(A_ex=extrap(Ls, A_e0), v_ex=extrap(Ls, v), c_half=ch[-1])
            print('   extrapolated (two largest L, linear in 1/L^2): A=%.4f (%.3f)  v=%.4f (%.3f)  c_half*v_ex=%.4f'
                  % (out['A_ex'], out['A_ex'] / th['A'], out['v_ex'], out['v_ex'] / th['v'], ch[-1] * out['v_ex']))
        if len(Ls) >= 3:
            A, e0, b = core.fit_cas(Ls, [rs[L]['E0'] for L in Ls])
            out['A_cas'] = A
            print('   free-e0 Casimir fit over %s: A=%.4f (%.3f) e0=%.6f (exact %.6f)' % (Ls, A, A / th['A'], e0, th['e0']))
        rows.append(out)
    return rows


if __name__ == '__main__':
    main(sys.argv[1:])
