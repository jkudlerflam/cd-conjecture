"""Finite checks for PROOF.md; numpy only."""
import json
from pathlib import Path
import numpy as np

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def tensor(factors):
    out = np.ones((1, 1), dtype=complex)
    for factor in factors:
        out = np.kron(out, factor)
    return out


results = []
for n in (1, 2, 3):
    gammas = []
    for k in range(n):
        for pauli in (X, Y):
            gammas.append(tensor([Z] * k + [pauli] + [I] * (n-k-1)))
    parity = tensor([Z] * n)
    rights = [1j * gamma @ parity for gamma in gammas]
    d = 2**n
    unit = np.eye(d)
    car_error = max(np.max(np.abs(a @ b + b @ a - (2*unit if i == j else 0)))
                    for i, a in enumerate(gammas) for j, b in enumerate(gammas))
    mixed_error = max(np.max(np.abs(a @ b - b @ a - (2j*parity if i == j else 0)))
                      for i, a in enumerate(gammas) for j, b in enumerate(rights))
    bond = sum(np.kron(r, g) for r, g in zip(rights, gammas))
    k12, k23 = np.kron(bond, unit), np.kron(unit, bond)
    comm = k12 @ k23 - k23 @ k12
    predicted = 2j * sum(tensor([r, parity, g]) for r, g in zip(rights, gammas))
    comm_error = float(np.max(np.abs(comm-predicted)))
    bond_norm = float(np.max(np.abs(np.linalg.eigvalsh(bond))))
    comm_norm = float(np.max(np.abs(np.linalg.eigvalsh(1j*comm))))
    assert car_error < 1e-12 and mixed_error < 1e-12 and comm_error < 1e-12
    assert abs(bond_norm-2*n) < 1e-10 and abs(comm_norm-4*n) < 1e-10
    results.append(dict(n=n, car_error=float(car_error), mixed_error=float(mixed_error),
                        commutator_identity_error=comm_error, bond_norm=bond_norm,
                        adjacent_commutator_norm=comm_norm,
                        minimum_bond_error_at_J1=float(np.sqrt((2*n)**2+2*n)-2*n)))

collar_checks = 0
for depth in range(8):
    for gate_range in range(1, 6):
        spread = depth*(gate_range-1)
        width = 2*spread+1
        crossing = list(range(-spread, spread+1))
        union_low = min(i-spread for i in crossing)
        union_high = max(i+1+spread for i in crossing)
        assert 1-union_low == width and union_high == width
        block_length = 2*width
        assert set(range(width)).isdisjoint(range(block_length-width, block_length))
        collar_checks += 1

report = dict(clifford_checks=results, collar_checks=collar_checks, passed=True)
path = Path(__file__).with_name('results.json')
path.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
