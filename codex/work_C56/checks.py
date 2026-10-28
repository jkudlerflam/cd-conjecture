"""Small checks of C56's algebra and conditional counterterm construction.

These do not test the missing renormalized-tree bounds or the physical chain.
"""
import json
import math
import random
from pathlib import Path

import numpy as np


def mst_weight(points, metric):
    reached = {0}
    cost = 0.0
    distances = [metric(points[0], p) for p in points]
    while len(reached) < len(points):
        j = min((j for j in range(len(points)) if j not in reached),
                key=lambda j: distances[j])
        cost += distances[j]
        reached.add(j)
        for k in range(len(points)):
            distances[k] = min(distances[k], metric(points[j], points[k]))
    return cost


def network_check():
    rng = random.Random(56028)
    worst_ratio = 0.0
    for _ in range(120):
        nodes = {(0, 0)}
        target = rng.randrange(2, 22)
        while len(nodes) < target:
            x, t = rng.choice(sorted(nodes))
            dx, dt = rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
            nodes.add((x + dx, t + dt))
        tau, a, kappa = 1.7, 0.3, 0.04
        points = []
        for x, j in sorted(nodes):
            # Many external times in the same node, including duplicated legs.
            for _ in range(rng.randrange(1, 7)):
                p = (x, 2 * tau * (j + rng.random()))
                points.extend([p, p])
        metric = lambda p, q: a * abs(p[0] - q[0]) + kappa * abs(p[1] - q[1])
        measured = mst_weight(points, metric)
        bound = 2 * a * len(nodes) + 4 * kappa * tau * len(nodes)
        assert measured <= bound + 1e-12
        worst_ratio = max(worst_ratio, measured / bound)
    return worst_ratio


def diameter_counterexample():
    # Verify the terminal-tree formula by an independent MST for small grids.
    for s in range(2, 9):
        points = [(s * i, s * j) for i in range(s) for j in range(s)]
        q = 2.0
        metric = lambda p, r: q * math.log1p(abs(p[0] - r[0]) + abs(p[1] - r[1]))
        measured = mst_weight(points, metric)
        exact = q * (s * s - 1) * math.log1p(s)
        assert abs(measured - exact) < 1e-10
    return {
        str(s): (-s * math.log(2) - 2 * s * s * math.log(2)
                 - 2 * s * (s - 1) + 2 * (s * s - 1) * math.log1p(s))
        for s in (8, 16, 32, 64)
    }


def relevant_recursion():
    theta, eps, ell = 0.23, 0.007, 0.4
    d = 2 + ell * eps
    q = 2 ** (-theta) / d
    denominator = d - 2 ** (-theta)
    worst = 0.0
    errors = {}
    for N in (0, 2, 8, 20):
        max_error = 0.0
        for p in np.linspace(-0.1, 0.1, 101):
            f = 1 + 0.3 * math.sin(p)
            a = 0.0  # a_(-N-1)
            for h in range(-N, 1):
                a = (a - eps * f * 2 ** (theta * h)) / d
            exact_infinite = -eps * f / denominator
            exact_finite = exact_infinite * (1 - q ** (N + 1))
            worst = max(worst, abs(a - exact_finite))
            max_error = max(max_error, abs(a - exact_infinite))
        assert max_error <= eps * 1.03 / denominator * q ** (N + 1) + 1e-16
        errors[str(N)] = max_error
    assert worst < 1e-15
    return worst, errors


def two_channel_root_equation():
    # Four chiral branches, widely unequal signed speeds. The RG forcing is
    # a solvable test recursion, not a beta function derived for C39.
    v = np.array([1.0, -0.001, 3.0, -0.02])
    eps, theta, ell = 1e-7, 0.23, 0.4
    den = 2 + ell * eps - 2 ** (-theta)
    matrix = np.array([[0.2, 0.1, 0.0, -0.1],
                       [0.1, -0.1, 0.2, 0.0],
                       [0.0, 0.1, 0.1, 0.1],
                       [-0.1, 0.0, 0.1, 0.2]])
    source = np.array([1.0, 0.4, -0.3, 0.7])
    y = np.linalg.solve(np.diag(v) + eps / den * matrix, -eps / den * source)
    exact = np.arcsin(y)
    p = np.zeros(4)
    for _ in range(15):
        nu = eps / den * (source + matrix @ np.sin(p))
        p = np.arcsin(-nu / v)
    residual = v * np.sin(p) + eps / den * (source + matrix @ np.sin(p))
    assert np.max(np.abs(p - exact)) < 1e-15
    assert np.max(np.abs(residual)) < 1e-20
    Cnu = 1.1 / den
    bound = 2 * Cnu * eps / np.min(np.abs(v))
    assert np.max(np.abs(p)) <= bound
    return {"max_root_error": float(np.max(np.abs(p - exact))),
            "max_equation_residual": float(np.max(np.abs(residual))),
            "max_root_shift": float(np.max(np.abs(p))),
            "proved_general_shift_bound": float(bound)}


def main():
    residual, errors = relevant_recursion()
    results = {
        "status": "Small algebra and geometry controls only; physical RG remains unproved.",
        "connected_network_max_ratio": network_check(),
        "diameter_example_log_terms": diameter_counterexample(),
        "relevant_recursion_max_formula_error": residual,
        "counterterm_uniform_grid_errors": errors,
        "four_branch_toy_root_equation": two_channel_root_equation(),
    }
    Path("checks.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
