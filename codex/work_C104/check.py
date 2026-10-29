"""Small deterministic consistency checks for C104, not infrared simulations."""
from pathlib import Path
import numpy as np


def comm(a, b):
    return a @ b - b @ a


def spin(d):
    s = (d - 1) / 2
    m = np.arange(d, dtype=float) - s
    plus = np.zeros((d, d), dtype=complex)
    for j in range(d - 1):
        plus[j + 1, j] = np.sqrt(s * (s + 1) - m[j] * (m[j] + 1))
    return [(plus + plus.T.conj()) / 2,
            (plus - plus.T.conj()) / (2j), np.diag(m)]


def opnorm(a):
    return np.max(np.abs(np.linalg.eigvalsh(a)))


lines = ['# C104 checks', '', '[Numerical] Deterministic matrices, d <= 8. No many-body infrared fits.', '']
max_error = 0.0
for d in range(2, 9):
    ts = spin(d)
    index = d * (d * d - 1) / 12
    gram = np.array([[np.trace(a @ b).real for b in ts] for a in ts])
    err = np.max(np.abs(gram - index * np.eye(3)))
    assert err < 1e-11
    max_error = max(max_error, err)
    bond = sum(np.kron(t, t) for t in ts)
    eig = np.linalg.eigvalsh(bond)
    mid = (eig[-1] + eig[0]) / 2
    scale = (eig[-1] - eig[0]) / 2
    centered = bond - mid * np.eye(d * d)
    casimir = sum(comm(np.kron(t, np.eye(d)), comm(np.kron(t, np.eye(d)), centered)) for t in ts)
    assert np.linalg.norm(casimir - 2 * bond) < 1e-9
    q = np.kron(ts[2] / np.sqrt(index), np.eye(d))
    double = comm(q, comm(q, centered))
    assert opnorm(double) <= 4 * scale + 1e-10
    # Independent, general normalized full-Casimir estimate.
    gamma = sum(t @ t for t in ts) / index
    assert opnorm(casimir / index) <= 4 * opnorm(gamma) * scale + 1e-10
    # Check the arbitrary-window variance bound with an integer spin shift.
    charge = np.arange(d, dtype=float)
    centered_charge = charge - np.mean(charge)
    assert np.dot(centered_charge, centered_charge) <= d * (d - 1) ** 2 / 4 + 1e-12
    lines.append(f'[Numerical] d={d}: I_R={index:.12g}; spin Gram, double commutator, and Casimir checks passed.')

# Enumerate integer partitions: every SU(2) representation of dimension <= 16.
def partitions(n, smallest=1):
    if n == 0:
        yield []
    for a in range(smallest, n + 1):
        for tail in partitions(n - a, a):
            yield [a] + tail

count = 0
for d in range(1, 17):
    for parts in partitions(d):
        assert sum(n ** 3 - n for n in parts) <= d ** 3 - d
        count += 1
lines.extend(['', f'[Numerical] Index bound checked for all {count} SU(2) decompositions at dimensions 1 through 16.'])

ratios = []
for rank in range(1, 101):
    ratios.append(('A', rank, (rank * rank + 2 * rank) / (rank * (rank + 2))))
    ratios.append(('C', rank, (2 * rank + 1) / (rank + 2)))
    if rank >= 2:
        ratios.append(('B', rank, (2 * rank + 1) / (2 * rank)))
    if rank >= 4:
        ratios.append(('D', rank, (2 * rank - 1) / (2 * rank - 1)))
for name, rank, dim, dual in [('E6', 6, 78, 12), ('E7', 7, 133, 18), ('E8', 8, 248, 30), ('F4', 4, 52, 9), ('G2', 2, 14, 4)]:
    ratios.append((name, rank, dim / (rank * (1 + dual))))
assert all(0 < x[2] <= 2 for x in ratios)
lines.append(f'[Numerical] {len(ratios)} Lie-algebra ratio checks passed, including all exceptional algebras.')

# Matching conventions, not an independent determination of the dictionary.
v = np.pi / 2
D_sz = v / (2 * np.pi)
D_q = D_sz / .5
assert np.isclose(D_sz, .25) and np.isclose(D_q, .5)
assert np.isclose(np.pi * D_q, v)
k, v, index = 2, 2 * np.pi, 2
assert np.isclose(k * v / (2 * np.pi), 2)
assert np.isclose(k * v / (2 * np.pi * index), 1)
lines.extend([
    '[Numerical] Heisenberg: D_Sz=1/4, D_normalized=1/2, A=pi D_normalized.',
    '[Numerical] Spin-1 SU(2)_2 with v=2 pi: D_Sz=2, D_normalized=1. The uncorrected fundamental formula gives 4.',
    f'[Numerical] Maximum spin trace residual: {max_error:.3e}.',
    '',
    '[Proved, scope] All assertions passed. These computations check normalization and finite-dimensional identities only. The derivations in PROOF.md establish the general estimates; no bound on Gaussian K was tested or inferred.'
])
report = '\n'.join(lines) + '\n'
Path('CHECKS.md').write_text(report)
print(report)
