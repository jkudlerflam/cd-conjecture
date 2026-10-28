"""Finite physical-chain diagnostics for the density-window theorem.

Run from any directory with numpy/scipy. No campaign model is imported.
The internal bond can violate every internal symmetry; occupied number
remains exact. Output is saved next to this file.
"""
from itertools import product
from fractions import Fraction
from pathlib import Path
import json
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def physical_sector(length, number, bond):
    q = int(round(np.sqrt(len(bond))))
    words = [w for w in product(range(q + 1), repeat=length)
             if sum(x != 0 for x in w) == number]
    index = {w: i for i, w in enumerate(words)}
    tr, tc, tv, vr, vc, vv = [], [], [], [], [], []
    for col, word in enumerate(words):
        for x in range(length):
            y = (x + 1) % length
            a, b = word[x], word[y]
            if (a == 0) != (b == 0):
                other = list(word)
                other[x], other[y] = b, a
                tr.append(index[tuple(other)])
                tc.append(col)
                tv.append(-1.0)
            elif a and b:
                inp = q * (a - 1) + b - 1
                for out in range(q * q):
                    coefficient = bond[out, inp]
                    if coefficient:
                        c, d = divmod(out, q)
                        other = list(word)
                        other[x], other[y] = c + 1, d + 1
                        vr.append(index[tuple(other)])
                        vc.append(col)
                        vv.append(coefficient)
    shape = (len(words), len(words))
    kinetic = coo_matrix((tv, (tr, tc)), shape=shape).tocsr()
    interaction = coo_matrix((vv, (vr, vc)), shape=shape,
                             dtype=complex).tocsr()
    return kinetic, interaction


def ground(matrix):
    if matrix.nnz == 0 or not np.any(matrix.data):
        return 0.0, 0.0
    if matrix.shape[0] <= 4:
        values, vectors = np.linalg.eigh(matrix.toarray())
        energy, state = values[0], vectors[:, 0]
    else:
        # A deterministic nonuniform start avoids symmetry-orthogonality.
        rng = np.random.default_rng(20261028 + len(matrix.data))
        start = rng.normal(size=matrix.shape[0])
        if np.iscomplexobj(matrix.data):
            start = start + 1j * rng.normal(size=matrix.shape[0])
        values, vectors = eigsh(matrix, k=1, which='SA', tol=2e-12,
                                maxiter=20000, v0=start)
        energy, state = values[0], vectors[:, 0]
    residual = np.linalg.norm(matrix @ state - energy * state)
    assert residual < 2e-8, residual
    return float(energy), float(residual)


def main():
    started = time.monotonic()
    rng = np.random.default_rng(241023)
    z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    random_bond = (z + z.conj().T) / 2
    random_bond /= np.max(np.abs(np.linalg.eigvalsh(random_bond)))
    cases = []
    worst_free_error = worst_residual = 0.0
    for name, bond, lengths in [
            ('generic_complex', random_bond, range(3, 8)),
            ('attractive_scalar', -np.eye(4), (5, 6)),
            ('repulsive_scalar', np.eye(4), (5, 6))]:
        norm = float(np.max(np.abs(np.linalg.eigvalsh(bond))))
        for length in lengths:
            rows = []
            for number in range(length + 1):
                kinetic, interaction = physical_sector(length, number, bond)
                diff = interaction - interaction.conj().T
                assert not diff.nnz or np.max(np.abs(diff.data)) < 1e-14
                free, residual = ground(kinetic)
                reference = -2 * np.sin(np.pi * number / length) / np.sin(np.pi / length)
                free_error = abs(free - reference)
                assert free_error < 2e-9, (name, length, number, free, reference)
                worst_free_error = max(worst_free_error, free_error)
                worst_residual = max(worst_residual, residual)
                energies = {}
                for coupling in (0.001, 0.04, 0.6, 4.0):
                    energy, residual = ground(kinetic + coupling * interaction)
                    worst_residual = max(worst_residual, residual)
                    lower = reference - coupling * norm * number
                    upper = reference + coupling * norm * number
                    assert lower - 2e-8 <= energy <= upper + 2e-8
                    energies[str(coupling)] = energy
                rows.append(dict(number=number, dimension=kinetic.shape[0],
                                 free_energy=free, energies=energies))
            minima = []
            for coupling in (0.001, 0.04, 0.6, 4.0):
                minimum = min(row['energies'][str(coupling)] for row in rows)
                minimizing = [row['number'] for row in rows
                              if row['energies'][str(coupling)] < minimum + 2e-8]
                bound = 3 * np.pi * coupling * norm / 16 + np.pi**2 / (32 * length**2)
                for number in minimizing:
                    assert (number / length - 0.5)**2 <= bound + 1e-12
                minima.append(dict(coupling=coupling, numbers=minimizing,
                                   squared_density_bound=bound))
            cases.append(dict(bond=name, length=length, norm=norm,
                              sectors=rows, minima=minima))

    # Exact rational checks of the thermodynamic window and C33 transfer.
    assert Fraction(3 * 22, 7 * 16000) < Fraction(1, 1600)
    assert 3 * Fraction(21, 40)**2 < 1
    assert Fraction(1600, 441) > 3
    margin = Fraction(19, 40 * 24) * Fraction(31558, 10**11)
    assert margin > Fraction(624585, 10**14)
    reciprocal_margin = 2 * np.pi * (1 - 21 * np.sqrt(3) / 40)
    assert reciprocal_margin > 0.5
    gaussian_lower_dimension = 3 * np.sqrt(3 / 5)
    assert gaussian_lower_dimension > 2
    output = dict(
        status='passed', scope='finite density/energy bounds, not CFT data',
        runtime_seconds=time.monotonic() - started,
        kinetic_max_error=worst_free_error, max_eigenvector_residual=worst_residual,
        safe_density_window=['19/40', '21/40'],
        reference_physical_momentum_lower_bound=float(margin),
        nonzero_reciprocal_margin=reciprocal_margin,
        gaussian_lower_dimension=gaussian_lower_dimension,
        cases=cases)
    Path(__file__).with_name('density_checks.json').write_text(
        json.dumps(output, indent=2) + '\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'cases'}, indent=2))


if __name__ == '__main__':
    main()
