"""Clock-source Gram sums for product cell ground states (C89/C94 vacuum-clock sources).

h_i is the coefficient of J A_i: mu/2 at the two marker bonds (m-1, m) of each cell, and the clock
hop -(|q><q+1| + h.c.) at bond start+q of a non-pinned cell. phi is the product of isolated cell
ground states w_q = sin(pi (q+1)/(b+1)) (normalized). Returns S(r) = sum_i <h_i h_{i+r}>, r mod N,
so that B = sum_r chi_N(r) S(r).
"""
import math
import numpy as np


def selector(b):
    def a(l): return 2 * (1 - math.cos(math.pi / (l + 1)))
    def d(l): return (l + 1) * a(l) - l * a(l + 1)
    return 2 - (d(b) + d(b - 1)) / 2


def cell_energy(l, mu):
    return mu - 2 * math.cos(math.pi / (l + 1))


def frozen(lengths, pins, mu):
    return sum(mu if p else cell_energy(b, mu) for b, p in zip(lengths, pins))


def gram(lengths, pins, mu):
    N = sum(lengths)
    scal = np.zeros(N); mean = np.zeros(N)
    blocks = []            # (bond indices, connected covariance of hop operators)
    start = 0
    for b, p in zip(lengths, pins):
        m = start + b - 1
        scal[(m - 1) % N] += mu / 2; scal[m % N] += mu / 2
        if not p and b >= 2:
            w = np.sin(np.pi * np.arange(1, b + 1) / (b + 1)); w /= np.linalg.norm(w)
            K = np.zeros((b - 1, b))          # K_q w as vectors
            for q in range(b - 1):
                K[q, q] = -w[q + 1]; K[q, q + 1] = -w[q]
            kmean = K @ w
            cov = K @ K.T - np.outer(kmean, kmean)
            idx = np.arange(start, start + b - 1)
            mean[idx] += kmean
            blocks.append((idx, cov))
        start += b
    hbar = scal + mean
    return hbar, blocks


def S_of_r(lengths, pins, mu):
    N = sum(lengths)
    hbar, blocks = gram(lengths, pins, mu)
    f = np.fft.fft(hbar)
    S = np.fft.ifft(np.conj(f) * f).real       # S(r) = sum_i hbar_i hbar_{i+r}
    for idx, cov in blocks:
        for a, i in enumerate(idx):
            r = (idx - i) % N
            np.add.at(S, r, cov[a])
    return S


def B_coeff(chi, lengths, pins, mu):
    return float(np.dot(chi, S_of_r(lengths, pins, mu)))
