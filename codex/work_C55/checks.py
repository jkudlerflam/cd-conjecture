"""Small controls for C55, not an interacting RG calculation.

Run python3 -B checks.py in this folder. Only NumPy and the standard library
are used. Results are printed and saved locally.
"""

import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


def shell(p, energy, z, nr=72, na=256, linear=False):
    """One smooth annular shell, integrated in elliptic polar coordinates."""
    nodes, weights = leggauss(nr)
    radius = 0.75 + 0.25 * nodes
    weights = 0.25 * weights
    angle = 2 * np.pi * (np.arange(na) + 0.5) / na
    rr = radius[:, None]
    v = 2 * np.sin(p)
    q0 = energy * rr * np.sin(angle)[None, :]
    q = energy * rr * np.cos(angle)[None, :] / v
    # Stable difference of cosines, without small-energy subtraction loss.
    e = 4 * np.sin(p + q / 2) * np.sin(q / 2)
    if linear:
        e = v * q
    u = 2 * rr - 1
    cutoff = np.exp(4 - 1 / (u * (1 - u)))
    multiplier = cutoff / (-1j * q0 + e)
    phase = np.exp(1j * (q0 * z[0] + q * z[1]))
    measure = energy**2 * rr / v * weights[:, None]
    return np.sum(measure * multiplier * phase) / (2 * np.pi * na)


def shell_checks():
    rows = []
    max_linear_zero = 0.0
    max_scaled = 0.0
    for p in (0.9, 1.9):
        for z in ((0.0, 0.0), (0.3, 1.0), (1.2, -2.0)):
            vals = []
            for h in range(10):
                en = 0.125 * 2.0 ** (-h)
                value = shell(p, en, z)
                scaled = abs(value) / en**2
                vals.append(scaled)
                max_scaled = max(max_scaled, scaled / (1 + abs(z[0]) + abs(z[1])))
                linzero = abs(shell(p, en, (0.0, 0.0), linear=True))
                max_linear_zero = max(max_linear_zero, linzero)
            # The E^2 coefficient is finite and approaches a nonzero limit.
            assert max(vals) < 0.2
            assert abs(vals[-1] - vals[-2]) < 2e-5 * max(1.0, vals[-1])
            rows.append({"p": p, "z": z, "scaled_first": vals[0],
                         "scaled_last": vals[-1]})
    assert max_linear_zero < 1e-14
    # Reflection-symmetric half filling removes the curvature tadpole at z=0.
    half_zero = abs(shell(np.pi / 2, 0.125, (0.0, 0.0)))
    assert half_zero < 1e-14
    return {"linear_origin_residual": max_linear_zero,
            "half_filled_origin_residual": half_zero,
            "max_E2_weighted_ratio": max_scaled, "cases": rows}


def cofactor_checks():
    rng = np.random.default_rng(55)
    worst = 0.0
    for n in range(1, 8):
        a = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / 5
        b = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / 7
        nodes, weights = leggauss(n + 1)
        rhs = 0j
        for s, w in zip((nodes + 1) / 2, weights / 2):
            m = a + s * b
            derivative = 0j
            for i in range(n):
                for j in range(n):
                    minor = np.delete(np.delete(m, i, axis=0), j, axis=1)
                    derivative += b[i, j] * (-1) ** (i + j) * np.linalg.det(minor)
            rhs += w * derivative
        lhs = np.linalg.det(a + b) - np.linalg.det(a)
        err = abs(lhs - rhs)
        worst = max(worst, err)
        assert err < 1e-12
    return {"determinant_interpolation_residual": worst}


def sequence_checks():
    gamma, alpha, kappa = 2.0, 0.3, 0.9
    analytic = 1 / ((1 - gamma ** (-kappa)) * (gamma ** (kappa - alpha) - 1))
    direct = sum(gamma ** (-kappa * n) *
                 sum(gamma ** (alpha * r) for r in range(1, n + 1))
                 for n in range(1, 250))
    assert abs(analytic - direct) < 1e-12
    tail = sum(gamma ** (-j - 1) * gamma ** (-alpha * j) for j in range(250))
    d_inv = 1 / (gamma - gamma ** (-alpha))
    assert abs(tail - d_inv) < 1e-14

    # Independent finite relevant recursion and terminal-tail formula.
    lower = -70
    beta = {h: 0.01 * gamma ** (alpha * h) * (1 + 0.2 * math.sin(h))
            for h in range(lower + 1, 1)}
    nu = {lower: 0.0}
    for h in range(lower + 1, 1):
        nu[h] = (nu[h - 1] - beta[h]) / gamma
    worst = max(abs(nu[h] + sum(gamma ** (k - h - 1) * beta[k]
                               for k in range(lower + 1, h + 1)))
                for h in beta)
    assert worst < 1e-14
    return {"history_sum_residual": abs(analytic - direct),
            "relevant_tail_norm_residual": abs(tail - d_inv),
            "finite_recurrence_residual": worst}


def root_and_source_checks():
    p0, eps = 1.9, 0.004
    left, right = 1.8, 2.0
    vmin = min(2 * math.sin(left), 2 * math.sin(right))
    c = 0.8
    radius = 2 * c * eps / vmin
    def f(p):
        return 2 * (math.cos(p0) - math.cos(p)) + eps * (0.5 + 0.2 * math.cos(p))
    lo, hi = p0 - radius, p0 + radius
    assert f(lo) < 0 < f(hi)
    for _ in range(80):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    root = (lo + hi) / 2
    assert abs(root - p0) <= radius
    assert abs(f(root)) < 1e-14
    pair = np.exp(-1j * root) - np.exp(1j * root)
    pair_error = abs(pair + 2j * math.sin(root))
    assert pair_error < 1e-14
    return {"synthetic_root": root, "root_residual": abs(f(root)),
            "displacement": abs(root - p0), "proved_radius": radius,
            "pair_overlap_residual": pair_error}


if __name__ == "__main__":
    result = {"scope": "Numerical controls of C55 lemmas, not interacting RG",
              "shells": shell_checks(), "cofactors": cofactor_checks(),
              "sequences": sequence_checks(), "root_sources": root_and_source_checks()}
    Path("checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
