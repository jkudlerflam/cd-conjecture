#!/usr/bin/env python3
"""A/(2 pi ||h||) against d for the rank-k symmetric SU(N) family (Bethe/WZW formula),
SU(d) swaps, and the envelope R_max(d) = max over (N, k) with binom(N+k-1, k) <= d.
Compares with R = (1/2) log2 d, i.e. A = pi ||h|| log2 d."""
import numpy as np
from math import comb, log


def R(N, k):
    Hk = sum(1.0 / j for j in range(1, k + 1))
    return k * (N * N - 1) / (N * (N + k) * Hk)


pts = []
for N in range(2, 400):
    for k in range(1, 400):
        d = comb(N + k - 1, k)
        if d > 10 ** 40:
            break
        pts.append((d, R(N, k), N, k))
pts.sort()
print('first (N,k) with R > 1 by d:', [(d, N, k, round(r, 4)) for d, r, N, k in pts if r > 1][:6])
best = -1
print('%12s %8s %6s %6s %9s %9s %12s' % ('d', 'R_max', 'N', 'k', 'R/ln d', '0.5log2d', 'R*lnln d/ln d'))
env = []
for d, r, N, k in pts:
    if r > best:
        best = r
        env.append((d, r, N, k))
for d, r, N, k in env:
    if d in (2, 3, 6, 10, 15, 21, 28, 36) or (d > 36 and log(d, 10) % 1 < 0.05) or d == env[-1][0]:
        print('%12.3g %8.4f %6d %6d %9.4f %9.3f %12.4f' % (d, r, N, k, r / log(d), 0.5 * log(d, 2), r * log(log(d)) / log(d) if d > 3 else np.nan))
viol = [(d, r) for d, r, N, k in pts if r > 0.5 * log(d, 2) + 1e-12]
print('violations of R <= (1/2)log2 d in this family:', viol[:5])
print('max R/ln d over family: %.4f; SU(2) Heisenberg gives %.4f' % (max(r / log(d) for d, r, N, k in pts), 0.5 / log(2)))
