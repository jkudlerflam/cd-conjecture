"""[Proved] Finite algebra checks supporting RESULT.md, not a dynamics simulation."""

import math
import numpy as np


def symmetric_function(matrix, function):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * function(values)) @ vectors.T


def scaling_metric(stiffness, exchange):
    root = symmetric_function(stiffness, np.sqrt)
    inverse_root = symmetric_function(stiffness, lambda values: 1 / np.sqrt(values))
    absolute = symmetric_function(root @ exchange @ root, np.abs)
    return inverse_root @ absolute @ inverse_root


def check_counts():
    pairs = [
        (a, b)
        for a in range(-2, 3)
        for b in range(-2, 3)
        if (a - b) % 2 == 0 and a * a + b * b <= 8
    ]
    for n in range(1, 9):
        counts = np.zeros((9, 2), dtype=np.int64)
        counts[0, 0] = 1
        for _ in range(n):
            updated = np.zeros_like(counts)
            for a, b in pairs:
                norm = a * a + b * b
                for weight in range(9 - norm):
                    for parity in range(2):
                        updated[weight + norm, (parity + b) % 2] += counts[weight, parity]
            counts = updated
        shell4 = int(counts[4, 0] + 2 * counts[0, 0])
        shell8 = int(counts[8, 0] + 2 * counts[4, 0])
        choose = lambda j: math.comb(n, j) if n >= j else 0
        formula = 1 + 8 * n + 32 * choose(2) + 96 * choose(3) + 128 * choose(4)
        assert (shell4 + shell8) // 2 == formula
        print(f"[Proved by finite enumeration] n={n}, shell4={shell4}, shell8={shell8}, M={formula}: match")


def check_metric():
    size = 4
    exchange = np.block([[np.zeros((size, size)), np.eye(size)], [np.eye(size), np.zeros((size, size))]])
    raw = np.arange(1, (2 * size) ** 2 + 1, dtype=float).reshape(2 * size, 2 * size)
    raw = np.sin(raw) / 5
    stiffness = np.eye(2 * size) + raw.T @ raw
    metric = scaling_metric(stiffness, exchange)
    error = np.max(np.abs(metric @ exchange @ metric - exchange))
    assert error < 1e-12
    print(f"[Derived with standard bosonization input, numerical check] General metric identity residual={error:.3g}")

    velocities = np.array([0.8, 1.0, 1.3, 1.7])
    parameters = np.array([0.7, 0.9, 1.1, 1.4])
    stiffness = np.diag(np.concatenate([velocities / parameters, velocities * parameters]))
    expected = np.diag(np.concatenate([parameters, 1 / parameters]))
    error = np.max(np.abs(scaling_metric(stiffness, exchange) - expected))
    assert error < 1e-12
    print(f"[Derived with standard bosonization input, numerical check] Decoupled metric residual={error:.3g}")


if __name__ == "__main__":
    check_counts()
    check_metric()
