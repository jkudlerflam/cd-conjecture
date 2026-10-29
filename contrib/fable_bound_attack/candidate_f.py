#!/usr/bin/env python3
"""Envelope of the known chains against candidate functions f(d) in A <= pi M f(d).
Known analytic values of A/(pi M): Heisenberg 1 (d=2), SU(N) fundamental 2(1-1/N) (d=N), the
integrable SU(N) rank-k symmetric chains 2 k (N^2-1)/(N (N+k) H_k) (d = binom(N+k-1, k)), the
free-fermion value 2/pi per species, and the surveyed non-integrable points (all below the
SU(d) value at their d).  Candidates: log2 d, the harmonic number H_{d-1}, 2(1-1/d)."""
import json
from math import comb, log, log2
import numpy as np

def harmonic(n):
    return sum(1.0 / j for j in range(1, n + 1))

def R_sym(N, k):
    return 2 * k * (N * N - 1) / (N * (N + k) * harmonic(k))

rows = []
dmax = 200
env = {}
for N in range(2, dmax + 1):
    for k in range(1, 60):
        d = comb(N + k - 1, k)
        if d > dmax:
            break
        env[d] = max(env.get(d, 0), R_sym(N, k))
# running envelope over all d' <= d (a bound must be monotone in d by stacking with a trivial factor: see README Lemma 1)
best = 0; out = []
for d in range(2, dmax + 1):
    best = max(best, env.get(d, 0))
    out.append({'d': d, 'envelope': best, 'log2': log2(d), 'H_{d-1}': harmonic(d - 1), '2(1-1/d)': 2 * (1 - 1 / d),
                'env/log2': best / log2(d), 'env/H': best / harmonic(d - 1), 'env/2(1-1/d)': best / (2 * (1 - 1 / d))})
# large d asymptotics of the symmetric-family envelope
big = []
for D in [1e3, 1e6, 1e9, 1e12, 1e20, 1e40]:
    b = 0
    for N in range(2, 20000):
        for k in range(1, 400):
            d = comb(N + k - 1, k)
            if d > D:
                break
            b = max(b, R_sym(N, k))
        if comb(N, 1) > D:
            break
    big.append({'d': D, 'envelope': b, 'log2': log2(D), 'H': log(D) + 0.5772, 'ratio_log2': b / log2(D), 'ratio_H': b / (log(D) + 0.5772)})
json.dump({'table': out, 'large_d': big}, open('data/candidate_f.json', 'w'), indent=1)
for r in out[:12] + [out[19], out[26], out[35], out[99], out[-1]]:
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
print('first d where 2(1-1/d) fails:', next(r['d'] for r in out if r['env/2(1-1/d)'] > 1 + 1e-12))
for b in big:
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in b.items()})
print('max env/H over d<=200:', max(r['env/H'] for r in out), 'max env/log2:', max(r['env/log2'] for r in out))
