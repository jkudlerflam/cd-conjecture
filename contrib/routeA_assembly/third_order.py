#!/usr/bin/env python3
"""Third-order two-cell coefficient of C19's final reflection cell, by operator class.

Builds C19's cell with codex/work_C19/uniform_cell.py (imported read-only; run
with `python3 -B` so no bytecode is written there), computes

    K3 = P V Q H0^-1 Q V Q H0^-1 Q V P        (PVP = 0, C19 section 10)

for two cells by full local diagonalization (same method as C19's
completion_checks.third_order, rewritten to keep memory near 2 GB at n = 3),
and splits its traceless part into Majorana-bilinear classes:

    same-flavour number conserving   (renormalizes t_a, theta_a, mu_a)
    same-flavour pairing             (U(1)_a breaking, momentum 2 theta_a)
    cross-flavour conversion         (momentum p_b - p_a)
    cross-flavour pairing            (momentum p_a + p_b)
    non-Gaussian residual            (operator norm after removing all bilinears)

Everything is reported relative to rho^2, the coefficient scale of the
second-order target, so a value c means a correction of relative size
c * lambda / Delta.  Usage:  python3 -B third_order.py [n ...]   (default 1 2 3)
"""
import os
import sys
import json
import time
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parents[1] / 'codex' / 'work_C19'))
import numpy as np                      # noqa: E402
import scipy.linalg as la               # noqa: E402
import uniform_cell as UC               # noqa: E402  (C19, read-only)


def K3_two_cell(c):
    H, G = c['H'], c['G']
    D, N = c['D'], c['N']
    NN = N * N
    ev, U = la.eigh(H)
    L = [U.conj().T @ f @ U for f in c['FL']]
    R = [U.conj().T @ f @ U for f in c['FR']]
    gg = U.conj().T @ G
    den = (ev[:, None] + ev[None, :]).ravel()
    inv = np.zeros_like(den)
    inv[den > .5] = 1 / den[den > .5]
    w = np.zeros((D * D, NN), complex)
    for r, l in zip(R, L):
        w += np.kron(r @ gg, l @ gg)
    w *= inv[:, None]
    w3 = w.reshape(D, D, NN)
    out = np.zeros((NN, NN), complex)
    for r, l in zip(R, L):
        t1 = (r @ w3.reshape(D, D * NN)).reshape(D, D, NN)          # sum_b r[a,b] w[b,d,k]
        t2 = (l @ t1.transpose(1, 0, 2).reshape(D, D * NN))           # sum_d l[c,d] t1[a,d,k]
        t2 = t2.reshape(D, D, NN).transpose(1, 0, 2)                  # [a,c,k]
        out += w.conj().T @ t2.reshape(D * D, NN)
        del t1, t2
    return out


def classify(h, n):
    """Bilinear coefficients of h in the two-cell Majorana basis, by class."""
    gs, _, P = UC.modes(n)
    N = 2 ** n
    I = np.eye(N, dtype=complex)
    Gam = [np.kron(g, I) for g in gs] + [np.kron(P, g) for g in gs]   # cell 0, cell 1
    m = len(Gam)
    d = N * N
    h = h - np.trace(h) / d * np.eye(d)
    A = np.zeros((m, m))
    fit = np.zeros_like(h)
    for a in range(m):
        for b in range(a):
            op = 1j * Gam[a] @ Gam[b]
            coef = np.real(np.trace(op @ h)) / d
            A[a, b], A[b, a] = coef, -coef
            fit += coef * op
    resid = float(la.norm(h - fit, 2))
    flav = [(a % (2 * n)) // 2 for a in range(m)]
    eps = np.array([[0., 1.], [-1., 0.]])
    cls = dict(same_conserving=0., same_pairing=0., cross_conversion=0., cross_pairing=0.)
    for a in range(0, m, 2):
        for b in range(0, m, 2):
            B = A[a:a + 2, b:b + 2]
            comm = (B + eps @ B @ eps.T) / 2
            anti = (B - eps @ B @ eps.T) / 2
            same = flav[a] == flav[b]
            cls['same_conserving' if same else 'cross_conversion'] = max(
                cls['same_conserving' if same else 'cross_conversion'], float(np.abs(comm).max()))
            cls['same_pairing' if same else 'cross_pairing'] = max(
                cls['same_pairing' if same else 'cross_pairing'], float(np.abs(anti).max()))
    return cls, resid, float(la.norm(h, 2))


def main(ns):
    out = []
    for n in ns:
        t0 = time.time()
        c = UC.reflection_cell(n, True)
        rho2 = float(c['rho']) ** 2
        tgt, _, tgt_norm = classify(c['H2'], n)        # second-order bond, for normalization
        K3 = K3_two_cell(c)
        herm = float(la.norm(K3 - K3.conj().T))
        K3 = (K3 + K3.conj().T) / 2
        cls, resid, tl = classify(K3, n)
        rec = dict(n=n, D=c['D'], rho2=rho2, seconds=None, hermiticity_error=herm,
                   K3_norm=float(la.norm(K3, 2)), K3_traceless_norm=tl,
                   c3_traceless=tl / rho2, nonGaussian_over_rho2=resid / rho2,
                   nonGaussian_abs=resid,
                   target_max_coefficient=tgt['same_conserving'],
                   K3_classes_over_rho2={k: v / rho2 for k, v in cls.items()},
                   K3_classes_over_target={k: v / tgt['same_conserving'] for k, v in cls.items()},
                   target_classes_over_rho2={k: v / rho2 for k, v in tgt.items()})
        rec['seconds'] = time.time() - t0
        out.append(rec)
        print(json.dumps(rec), flush=True)
        (HERE / 'third_order_results.json').write_text(json.dumps(out, indent=1))
    return out


if __name__ == '__main__':
    main([int(x) for x in sys.argv[1:]] or [1, 2, 3])
