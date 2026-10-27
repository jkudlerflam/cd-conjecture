#!/usr/bin/env python3
"""Independent checks of the Astra A claims used in README.md (numpy, scipy; under a minute).

(a) delta_2 for the C8 n = 40 example (theta_a = sqrt(p_a)/100, k_a = pi/2 + sqrt(q_a)/100).
(b) delta_2 and delta_4 for Astra's Sidon scaffold at n = 40 (N = 80 Fermi points, q = 83).
    delta_4 is found from near-coincident signed pair sums with disjoint index pairs.
(c) One-loop flow of the so(4) current-current couplings g_AB J_R^A J_L^B,
    dg_EF/dl = (1/2) f^ACE f^BDF g_AB g_CD: isotropic g < 0 is irrelevant, and a small
    anisotropy on one current coupling runs away for one sign and flows to the Gaussian line
    for the other (the isotropic point is on the separatrix).
"""
import itertools
import numpy as np
from scipy.integrate import solve_ivp


def primes(k):
    ps, n = [], 2
    while len(ps) < k:
        if all(n % q for q in ps if q * q <= n):
            ps.append(n)
        n += 1
    return ps


def delta2(P):
    X = np.sort(np.concatenate([P, -P]) % (2 * np.pi))
    return float(np.min(np.diff(np.concatenate([X, [X[0] + 2 * np.pi]]))))


def delta4(P, window=400):
    N = len(P)
    idx, vals = [], []
    for i in range(N):
        for j in range(i + 1, N):
            for s, t in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                idx.append((i, j)); vals.append(s * P[i] + t * P[j])
    vals = np.array(vals) % (2 * np.pi)
    order = np.argsort(vals); vs = vals[order]
    vs_ext = np.concatenate([vs, vs + 2 * np.pi]); order_ext = np.concatenate([order, order])
    best = np.inf
    for a in range(len(vs)):
        ia = idx[order_ext[a]]
        for b in range(a + 1, a + window):
            ib = idx[order_ext[b]]
            if len(set(ia) | set(ib)) == 4:
                best = min(best, vs_ext[b] - vs_ext[a]); break
    return float(best)


def part_a():
    ps = primes(80)
    th = np.sqrt(ps[:40]) / 100; k = np.pi / 2 + np.sqrt(ps[40:]) / 100
    P = np.concatenate([th + k, th - k])
    print('(a) C8 n=40 example: delta_2 = %.6e' % delta2(P))


def part_b():
    N, q = 80, 83
    A = [2 * q * i + (i * i % q) for i in range(N)]
    dA = min(np.diff(A)); B = -(-dA // 2); T = 2 * A[-1] + 4 * B + 2
    P = np.array([2 * np.pi / T * (a + B + 1 / 6) for a in A])
    print('(b) Sidon scaffold: T = %d, delta_2 = %.10f (198 pi/T = %.10f), delta_4 = %.6e (2 pi/(3T) = %.6e)'
          % (T, delta2(P), 198 * np.pi / T, delta4(P), 2 * np.pi / (3 * T)))
    # velocities with the pairing p_{a,-} = P[a-1], p_{a,+} = P[a+39]
    kk = (P[40:] - P[:40]) / 2
    print('    min velocity 2 sin k_a = %.4f' % np.min(2 * np.sin(kk)))


def part_c(m=4):
    pairs = list(itertools.combinations(range(m), 2)); nA = len(pairs)
    L = []
    for (a, b) in pairs:
        M = np.zeros((m, m), complex); M[a, b] = -1j; M[b, a] = 1j; L.append(M)
    f = np.zeros((nA, nA, nA))
    for A_ in range(nA):
        for B_ in range(nA):
            C = L[A_] @ L[B_] - L[B_] @ L[A_]
            for E in range(nA):
                f[A_, B_, E] = np.real(np.trace(C @ L[E]) / 2j)
    print('(c) so(%d): sum_AC f^ACE f^ACE = %.3f (2(m-2) = %d)' % (m, np.einsum('ACE,ACE->E', f, f)[0], 2 * (m - 2)))

    def rhs(t, y):
        g = y.reshape(nA, nA)
        return (0.5 * np.einsum('ACE,BDF,AB,CD->EF', f, f, g, g)).reshape(-1)

    def blow(t, y):
        return np.max(np.abs(y)) - 1.0
    blow.terminal = True
    for pert in (0.0, +0.01, -0.01):
        g = -0.05 * np.eye(nA); g[0, 0] += pert
        sol = solve_ivp(rhs, [0, 3000], g.reshape(-1), rtol=1e-8, atol=1e-12, events=blow)
        gf = sol.y[:, -1].reshape(nA, nA)
        print('    g0 = -0.05, anisotropy %+.2f on J^{01}: stops at l = %.1f, max|g| = %.4g'
              % (pert, sol.t[-1], np.max(np.abs(gf))))


if __name__ == '__main__':
    part_a(); part_b(); part_c()
