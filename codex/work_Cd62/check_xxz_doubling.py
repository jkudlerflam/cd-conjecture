"""Check the matched-sector proof against spin matrices and free references.

Run with python3 -B check_xxz_doubling.py. Only the standard library is used.
All energies use J = 1 and the negative real spin-exchange matrix element.
"""

from itertools import combinations
from math import atan, atanh, cos, cosh, pi, sin, sinh, sqrt, tan, tanh


def momentum(x, eta):
    return 2 * atan(tan(eta / 2) * tanh(x))


def phase(x, eta):
    return 2 * atan(cos(eta) / sin(eta) * tanh(x))


def bare_density(x, eta):
    return sin(eta) / (pi * (cosh(2 * x) + cos(eta)))


def scattering(x, eta):
    return sin(2 * eta) / (2 * pi * (sinh(x) ** 2 + sin(eta) ** 2))


def levels(size, count):
    return [(j - (count - 1) / 2) / size for j in range(count)]


def split_levels(old_levels, length, s):
    return [t + sign * s / (4 * length) for t in old_levels for sign in (-1, 1)]


def counting_residual(roots, size, targets, eta):
    return [
        (momentum(x, eta) + sum(phase(x - y, eta) for y in roots) / size)
        / (2 * pi) - t
        for x, t in zip(roots, targets)
    ]


def gaudin(roots, size, eta):
    matrix = []
    for i, x in enumerate(roots):
        row = [-scattering(x - y, eta) / size for y in roots]
        row[i] += bare_density(x, eta) - sum(row)
        matrix.append(row)
    return matrix


def positive_solve(matrix, rhs):
    """Cholesky solve for the positive real Gaudin matrix."""
    n = len(rhs)
    lower = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            entry = matrix[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            lower[i][j] = sqrt(entry) if i == j else entry / lower[j][j]
    solution = list(rhs)
    for i in range(n):
        solution[i] = (solution[i] - sum(lower[i][j] * solution[j] for j in range(i))) / lower[i][i]
    for i in reversed(range(n)):
        solution[i] = (solution[i] - sum(lower[j][i] * solution[j] for j in range(i + 1, n))) / lower[i][i]
    return solution


def solve_roots(size, targets, eta, start):
    roots = list(start)
    if not roots:
        return roots
    for _ in range(80):
        residual = counting_residual(roots, size, targets, eta)
        error = max(map(abs, residual))
        if error < 3e-15:
            return roots
        correction = positive_solve(gaudin(roots, size, eta), residual)
        step = 1.0
        for _ in range(45):
            trial = [x - step * dx for x, dx in zip(roots, correction)]
            if max(map(abs, trial)) < 100:
                new_error = max(map(abs, counting_residual(trial, size, targets, eta)))
                if new_error < error:
                    roots = trial
                    break
            step *= 0.5
        else:
            raise ArithmeticError("Newton line search stalled")
    raise ArithmeticError("Bethe residual did not converge")


def ground_roots(size, count, eta):
    targets = levels(size, count)
    roots = [atanh(tan(pi * t)) for t in targets]
    angle = pi / 2
    while angle > eta:
        angle = max(0.65 * angle, eta)
        roots = solve_roots(size, targets, angle, roots)
    return roots


def energy(size, roots, eta):
    return -size * cos(eta) / 4 - pi * sin(eta) * sum(bare_density(x, eta) for x in roots)


def dressed_velocity(roots, size, eta):
    bare = [2 * pi * sin(momentum(x, eta)) for x in roots]
    rhs = [bare_density(x, eta) * w for x, w in zip(roots, bare)]
    return positive_solve(gaudin(roots, size, eta), rhs)


def spin_energy(size, count, eta):
    """Perron iteration in the spin basis, stopped by two-sided row ratios."""
    masks = [sum(1 << site for site in occupied) for occupied in combinations(range(size), count)]
    indices = {mask: i for i, mask in enumerate(masks)}
    shifted_diagonal = []
    neighbours = []
    shift = size
    for mask in masks:
        adjacent = 0
        hops = []
        for site in range(size):
            following = (site + 1) % size
            a, b = (mask >> site) & 1, (mask >> following) & 1
            adjacent += a * b
            if a != b:
                hops.append(indices[mask ^ (1 << site) ^ (1 << following)])
        diagonal = cos(eta) * (-size / 4 + count - adjacent)
        shifted_diagonal.append(shift - diagonal)
        neighbours.append(hops)
    state = [1.0] * len(masks)
    for _ in range(20000):
        new = [d * v + 0.5 * sum(state[j] for j in hops)
               for d, v, hops in zip(shifted_diagonal, state, neighbours)]
        ratios = [b / a for a, b in zip(state, new)]
        low, high = min(ratios), max(ratios)
        if high - low < 2e-12:
            return shift - (high + low) / 2
        scale = max(new)
        state = [v / scale for v in new]
    raise ArithmeticError("Perron bounds did not converge")


def check_references():
    errors = []
    for size in range(2, 11):
        for count in range(size // 2 + 1):
            free_energy = -sin(pi * count / size) / sin(pi / size)
            errors.append(abs(spin_energy(size, count, pi / 2) - free_energy))
            errors.append(abs(spin_energy(size, count, 0) + size / 4))
    assert max(errors) < 3e-12
    print(f"XX and isotropic ferromagnet references: {len(errors)}, max error {max(errors):.3g}")


def check_spin_endpoints():
    errors, ratios = [], []
    for eta in (0.02, 0.2, pi / 4):
        for length in (2, 3, 4, 5):
            for n in range(length // 2 + 1):
                energies = []
                for size, count in ((length, n), (2 * length, 2 * n)):
                    roots = ground_roots(size, count, eta)
                    exact_spin = spin_energy(size, count, eta)
                    errors.append(abs(energy(size, roots, eta) - exact_spin))
                    energies.append(exact_spin)
                ratios.append(length * abs(energies[1] - 2 * energies[0]) / (pi * sin(eta)))
    assert max(errors) < 3e-12 and max(ratios) <= 1 + 3e-10
    print(f"Bethe vs spin sector energies: {len(errors)}, max error {max(errors):.3g}")
    print(f"Spin doubling comparisons: {len(ratios)}, largest fraction of bound {max(ratios):.6f}")


def check_interpolation():
    ratios, differences = [], []
    for eta in (0.01, 0.15, pi / 4):
        for length, n in ((5, 1), (6, 3), (9, 3), (12, 6), (18, 7)):
            old = ground_roots(length, n, eta)
            roots = [x for x in old for _ in range(2)]
            old_levels = levels(length, n)
            for s in (0.0, 0.38, 0.8, 1.0):
                roots = solve_roots(2 * length, split_levels(old_levels, length, s), eta, roots)
                v = dressed_velocity(roots, 2 * length, eta)
                variation = sum(abs(b - a) for a, b in zip(v, v[1:]))
                ratios.append(variation / (4 * pi * sin(eta)))
                if s == 0.38:
                    step = 2e-4
                    two_energies = []
                    for offset in (-step, step):
                        targets = split_levels(old_levels, length, s + offset)
                        shifted = solve_roots(2 * length, targets, eta, roots)
                        two_energies.append(energy(2 * length, shifted, eta))
                    numerical = (two_energies[1] - two_energies[0]) / (2 * step)
                    analytic = sum(b - a for a, b in zip(v[::2], v[1::2])) / (4 * length)
                    differences.append(abs(numerical - analytic))
    assert max(ratios) <= 1 + 1e-12 and max(differences) < 2e-8
    print(f"Sampled velocity variations: {len(ratios)}, largest fraction of bound {max(ratios):.6f}")
    print(f"Energy derivative comparisons: {len(differences)}, max error {max(differences):.3g}")


if __name__ == "__main__":
    check_references()
    check_spin_endpoints()
    check_interpolation()
