#!/usr/bin/env python3
"""My refit of the exact finite-size Bethe energies in astra_bethe/bethe_exact_energies.json:
A(L) = -6L[E_0(L) - e_0 L]/pi with the Bethe e_0, fitted to A_inf + corrections over the upper
two thirds or upper half of the sizes. The 1/ln^3 L + 1/L^2 form reproduces the exact A at
N = 2..5 to 3e-4 in R, and brackets come from the other forms."""
import json
import numpy as np
from scipy.special import digamma

r = json.load(open('astra_bethe/bethe_exact_energies.json'))
forms = {'1/L^2': lambda L: [1 / L ** 2], '1/L^2+1/L^4': lambda L: [1 / L ** 2, 1 / L ** 4],
         '1/ln^3': lambda L: [1 / np.log(L) ** 3], '1/ln^3+1/L^2': lambda L: [1 / np.log(L) ** 3, 1 / L ** 2]}
for f in r['families']:
    N, k = f['N'], f['k']
    Hk = sum(1 / j for j in range(1, k + 1))
    A0 = k * (N * N - 1) / (N + k) * np.pi / N
    e0 = Hk - sum(digamma(1 + (j - 1) / N) - digamma(j / N) for j in range(1, k + 1)) / N
    d = [x for x in f['data'] if x['L'] >= 2 * N]
    L = np.array([x['L'] for x in d], float); E = np.array([x['E'] for x in d])
    AL = -6 * L * (E - e0 * L) / np.pi
    res = {}
    for name, cols in forms.items():
        for lo in (len(L) // 3, len(L) // 2):
            X = np.vstack([np.ones_like(L[lo:])] + cols(L[lo:])).T
            res[(name, int(L[lo]))] = np.linalg.lstsq(X, AL[lo:], rcond=None)[0][0] / (2 * np.pi * Hk / 2)
    best = [v for (n, l), v in res.items() if n == '1/ln^3+1/L^2']
    print('N=%d k=%d L=%d..%d  R_pred=%.4f  R(L_max)=%.4f  R_fit(1/ln^3+1/L^2)=%.4f  range over forms %.4f..%.4f'
          % (N, k, L[0], L[-1], A0 / (2 * np.pi * Hk / 2), AL[-1] / (2 * np.pi * Hk / 2), np.mean(best), min(res.values()), max(res.values())))
