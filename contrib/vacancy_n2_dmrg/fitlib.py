"""Fits of open-chain entanglement profiles.  Claude (Opus 5.5), 2026-10-28."""
import numpy as np


def chord(L, l):
    return np.log(2 * (L + 1) / np.pi * np.sin(np.pi * l / (L + 1)))


def profile_fit(S, L, qs, lmin=4, p=1.0, mid=None):
    """Least-squares fit of S(l), lmin <= l <= L - lmin, to
        (c/6) X + s0 + sum_q [a_q cos(q l) + b_q sin(q l)] e^{-p X},  X = ln[(2(L+1)/pi) sin(pi l/(L+1))].
    qs: oscillation wavevectors.  Returns c and the rms residual."""
    l = np.arange(lmin, L - lmin + 1)
    if mid is not None:     # restrict to the central window |l - L/2| <= mid * L / 2
        l = l[np.abs(l - L / 2) <= mid * L / 2]
    X = chord(L, l)
    cols = [X / 6, np.ones_like(X)]
    for q in qs:
        cols += [np.cos(q * l) * np.exp(-p * X), np.sin(q * l) * np.exp(-p * X)]
    A = np.column_stack(cols)
    y = np.asarray(S)[l]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    return float(coef[0]), float(np.sqrt(np.mean(r ** 2)))


def wavevectors(L, Ns, combos=False):
    q = [2 * np.pi * N / (L + 1) for N in Ns if 0 < N < L]
    if combos:
        q2 = []
        for i in range(len(q)):
            for j in range(i + 1, len(q)):
                q2 += [q[i] + q[j], abs(q[i] - q[j])]
        q = q + q2
    return q
