#!/usr/bin/env python3
"""Additivity under stacking: Heisenberg (x) Ising as one d = 4 site.
h = h_H (x) 1 + 1 (x) h_I with the two-site factors reordered to (s1 t1)(s2 t2).
Checks ||h|| = ||h_H|| + ||h_I|| (also for random pairs of bonds) and A = pi + 1."""
import numpy as np, json
import core, chains as ch


def stack(h1, d1, h2, d2):
    A = np.kron(h1, np.eye(d2 * d2)).reshape([d1, d1, d2, d2] * 2)       # (s1 s2 t1 t2) x (...)
    B = np.kron(np.eye(d1 * d1), h2).reshape([d1, d1, d2, d2] * 2)
    perm = [0, 2, 1, 3, 4, 6, 5, 7]                                      # -> (s1 t1 s2 t2)
    D = d1 * d2
    return (A + B).transpose(perm).reshape(D * D, D * D)


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    for d1, d2 in ((2, 2), (2, 3), (3, 2)):
        g = lambda d: (lambda m: m + m.conj().T)(rng.normal(size=(d * d, d * d)) + 1j * rng.normal(size=(d * d, d * d)))
        h1, h2 = g(d1), g(d2)
        print('random d1=%d d2=%d: ||stack||=%.10f  ||h1||+||h2||=%.10f' % (d1, d2, core.hnorm(stack(h1, d1, h2, d2), d1 * d2), core.hnorm(h1, d1) + core.hnorm(h2, d2)))
    hH, hI = ch.models()['Heisenberg'][0], ch.models()['Ising'][0]
    h = stack(hH, 2, hI, 2)
    print('Heisenberg+Ising: ||h||=%.6f (1 + sqrt2 = %.6f)' % (core.hnorm(h, 4), 1 + np.sqrt(2)))
    q = np.array([[1], [1], [-1], [-1]])            # 2 Sz of the Heisenberg factor (state = 2*s + t)
    Ls, E = [], []
    for L in (6, 8, 10, 12):
        r = core.ring(h, 4, L, charges=q, targets=[[0]], ks=range(L // 2 + 1), nlev=1)
        Ls.append(L); E.append(r['E0'])
        print(L, r['E0'], r['k0'], flush=True)
    A, e0, b = core.fit_cas(Ls, E)
    print('A_cas=%.4f (pi+1 = %.4f), ratio=%.4f, weighted average of 0.5 and 0.1125 = %.4f'
          % (A, np.pi + 1, A / (2 * np.pi * core.hnorm(h, 4)), (np.pi + 1) / (2 * np.pi * (1 + np.sqrt(2)))))
