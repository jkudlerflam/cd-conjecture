"""Exact Gaussian half-line XX profile by scalar spectral quadrature.

No finite-chain diagonalization. The Jacobi matrix has J[n,n+1]=(n+1)/2.
P_n are orthonormal for sech(pi*l) dl. Subtract the uniform vacuum using
C_uniform=(1+exp(2*pi*J))**(-1), before integrating, to avoid cancellation.
This is a lightweight recurrence, including at large bond coordinates.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre


def profile(max_x, order=400, cutoff=12.0):
    z, wt = roots_legendre(order)
    lam = (z + 1) * cutoff / 2
    wt = wt * cutoff / 2
    kernel = 2 * wt / (np.cosh(np.pi * lam) * (1 + np.exp(2 * np.pi * lam)))
    pm = np.zeros_like(lam)
    p = np.ones_like(lam)
    deficit = np.empty(max_x)
    for n in range(max_x):
        pp = (2 * lam * p - n * pm) / (n + 1)
        deficit[n] = np.dot(kernel, p * pp)
        pm, p = p, pp
    return deficit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-x', type=int, default=8193)
    ap.add_argument('--order', type=int, default=400)
    ap.add_argument('--output', default='xx_halfline_results.json')
    args = ap.parse_args()
    start = time.perf_counter()
    a = profile(args.max_x, args.order)
    b = profile(args.max_x, args.order * 2)
    rows = []
    for x in [3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, 2049, 4097, 8191]:
        if x + 1 > args.max_x:
            continue
        raw = a[x - 1]
        smooth = (a[x - 2] + 2 * a[x - 1] + a[x]) / 4
        rows.append(dict(x=x, deficit=raw, raw_scaled=x*x*raw,
                         binomial_scaled=x*x*smooth,
                         alternating_normalized=2*np.pi*x*np.log(x)*raw,
                         quadrature_difference=abs(a[x-1]-b[x-1])))
    result = dict(method='half-line spectral quadrature', max_x=args.max_x,
                  orders=[args.order, 2*args.order], cutoff=12.0,
                  absolute_tail_bound=float(np.exp(-24*np.pi)),
                  max_quadrature_difference=float(np.max(np.abs(a-b))),
                  smooth_prediction=1/(24*np.pi), rows=rows,
                  elapsed_seconds=time.perf_counter()-start)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
