#!/usr/bin/env python3
"""Closed-book F1 checks. Run with python3 f1.py.

Uses only the standard library. No input files are read. These checks verify
finite-dimensional constraints and central-charge arithmetic, not an IR limit.
See f1_results.md for proofs and the assumptions needed for continuum claims.
"""

from fractions import Fraction as F
from itertools import product
from math import comb


def rank(matrix):
    """Exact rational row rank, including empty matrices."""
    a = [[F(x) for x in row] for row in matrix]
    if not a:
        return 0
    pivot = 0
    for col in range(len(a[0])):
        hit = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if hit is None:
            continue
        a[pivot], a[hit] = a[hit], a[pivot]
        scale = a[pivot][col]
        a[pivot] = [x / scale for x in a[pivot]]
        for i in range(len(a)):
            if i != pivot and a[i][col]:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def su_c(n, k):
    return F(k * (n * n - 1), k + n) if n > 1 else F(0)


def rectangular_dimension(n, m, q):
    """SU(n) rectangle with q rows and m columns, from Weyl's formula."""
    ans = F(1)
    for i in range(1, q + 1):
        for j in range(q + 1, n + 1):
            ans *= F(m + j - i, j - i)
    assert ans.denominator == 1
    return ans.numerator


def occupation_sector(charges, target):
    return [s for s in product((0, 1), repeat=len(charges))
            if dot(charges, s) == target]


def projected_hop_count(states):
    """Count nonzero allowed single-flavor hops on two code sites.

    Fermion signs are irrelevant to whether a matrix element vanishes.
    """
    allowed = set(states)
    count = 0
    for left, right in product(states, repeat=2):
        for a in range(len(left)):
            if left[a] == right[a]:
                continue
            new_left, new_right = list(left), list(right)
            new_left[a], new_right[a] = right[a], left[a]
            count += (tuple(new_left) in allowed and tuple(new_right) in allowed)
    return count


def bilinear_action(bits, creation, annihilation):
    """Return (target bits, fermion sign) for f^dagger_i f_j, or None."""
    if not bits & (1 << annihilation):
        return None
    sign = (-1) ** (bits & ((1 << annihilation) - 1)).bit_count()
    after = bits ^ (1 << annihilation)
    if after & (1 << creation):
        return None
    sign *= (-1) ** (after & ((1 << creation) - 1)).bit_count()
    return after | (1 << creation), sign


def common_kernel_dimension(basis, generators):
    """Exact common kernel of bilinears, including images outside the basis."""
    rows = []
    for generator in generators:
        image_rows = {}
        for col, bits in enumerate(basis):
            for coefficient, creation, annihilation in generator:
                result = bilinear_action(bits, creation, annihilation)
                if result is not None:
                    target, sign = result
                    row = image_rows.setdefault(target, [0] * len(basis))
                    row[col] += coefficient * sign
        rows.extend(row for row in image_rows.values() if any(row))
    return len(basis) - rank(rows)


def check_explicit_singlet_spaces():
    # These checks construct exterior-Fock actions rather than use dimensions.
    for m in range(1, 6):
        # Spin index 0,1; color index a; orbital index spin*m+a.
        # Vanishing color Cartans at total occupation m require one per color.
        color_basis = [sum(1 << (spin[a] * m + a) for a in range(m))
                       for spin in product((0, 1), repeat=m)]
        color_generators = [
            [(1, a, b), (1, m + a, m + b)]
            for a in range(m) for b in range(m) if a != b
        ]
        assert common_kernel_dimension(color_basis, color_generators) == m + 1

        # Symplectic Cartans force each (a_i,b_i) pair empty or full.
        sp_basis = [sum((3 << (2 * a)) * occupied[a] for a in range(m))
                    for occupied in product((0, 1), repeat=m)]
        sp_generators = []
        # Matrices [[A,B],[C,-A^t]], with B and C symmetric.
        for a in range(m):
            for b in range(m):
                sp_generators.append([(1, 2*a, 2*b), (-1, 2*b+1, 2*a+1)])
            sp_generators.append([(1, 2*a, 2*a+1)])
            sp_generators.append([(1, 2*a+1, 2*a)])
            for b in range(a + 1, m):
                sp_generators.append([(1, 2*a, 2*b+1), (1, 2*b, 2*a+1)])
                sp_generators.append([(1, 2*a+1, 2*b), (1, 2*b+1, 2*a)])
        assert common_kernel_dimension(sp_basis, sp_generators) == m + 1

    for n in (4, 6):
        so_generators = [[(1, a, b), (-1, b, a)]
                         for a in range(n) for b in range(a + 1, n)]
        assert common_kernel_dimension(list(range(1 << n)), so_generators) == 2
    print("Explicit Fock-space singlet kernels passed: U(m) color and "
          "Sp(2m) for 1<=m<=5, SO(n) for n=4,6.")


def check_weighted_family():
    print("Weighted U(1) codes. Formal cosets are NOT established projected IRs.")
    for n in (4, 6, 8, 10, 12):
        charges = [2**j for j in range(n - 1)] + [2**(n - 1) - 1]
        target = charges[-1]
        states = occupation_sector(charges, target)
        expected = [(0,) * (n - 1) + (1,), (1,) * (n - 1) + (0,)]
        assert states == expected
        assert {sum(s) % 2 for s in states} == {1}
        assert sum(charges) == 2 * target  # Half-filled mean field is compatible.
        assert projected_hop_count(states) == 0

        # All additional linear occupation constraints, including backgrounds.
        constraints = []
        for j in range(1, n - 1):
            row = [0] * n
            row[0], row[j] = -1, 1
            constraints.append(row)
        constraints.append([1] + [0] * (n - 2) + [1])
        backgrounds = [0] * (n - 2) + [1]
        assert all([dot(row, s) for row in constraints] == backgrounds
                   for s in states)
        assert rank(constraints) == n - 1
        gram = [[dot(a, b) for b in constraints] for a in constraints]
        assert rank(gram) == n - 1
        assert rank(constraints + [charges]) == n - 1

        # An explicit current killed by the code but retained by the one-U(1) coset.
        killed_current = [2, -3, 1] + [0] * (n - 3)
        assert dot(killed_current, charges) == 0
        assert all(dot(killed_current, s) == 0 for s in states)
        assert dot(killed_current, killed_current) == 14
        kappa = dot(charges, charges)
        predicted_kappa = F(4**(n - 1) - 1, 3) + (2**(n - 1) - 1)**2
        assert kappa == predicted_kappa
        print(f"  n_c={n:2d}, d=2, level={kappa:8d}, "
              f"one-U(1) c={n-1:2d}, completed-torus c=1, "
              "killed surviving current level=14")


def check_rectangular_family():
    cases = 0
    for n in range(2, 11):
        for m in range(1, 11):
            for q in range(1, n):
                d = rectangular_dimension(n, m, q)
                c = su_c(n, m)
                gauge_c = su_c(m, n) + 1
                assert c + gauge_c == n * m
                assert d >= comb(n + m - 1, m)
                assert d >= 1 + m * (n - 1)
                assert c <= m * (n - 1) <= d - 1
                if q < n - 1:
                    ratio = F(rectangular_dimension(n, m, q + 1), d)
                    expected = F(comb(m + n - q - 1, n - q - 1),
                                 comb(m + q, q))
                    assert ratio == expected
                cases += 1
    print(f"Rectangular SU(N) representations: {cases} exact checks passed "
          "(2<=N<=10, 1<=m<=10, 1<=q<N).")


def check_other_families():
    for m in range(1, 101):
        spin_c = F(3 * m, m + 2)
        sp_c = F(m * (2 * m + 1), m + 2)
        color_c = su_c(m, 2) + 1
        assert spin_c + sp_c == 2 * m
        assert spin_c + color_c == 2 * m
        assert spin_c <= m  # d=m+1.
    for n in range(4, 102, 2):
        so_dim = n * (n - 1) // 2
        so_c = F(2 * so_dim, 2 + n - 2)
        assert n - so_c == 1  # Two local SO(n)-invariant states.
    for n in range(2, 31):
        for q in range(1, n):
            assert n - 1 <= comb(n, q) - 1
    for n in range(2, 101):
        para_c = su_c(2, n) - 1
        assert para_c == F(2 * (n - 1), n + 2)
        assert para_c <= n - 1
    print("Sp(2m), SO(n), fixed-charge exterior powers, and parafermion "
          "central-charge checks passed.")


def check_local_projection_and_examples():
    for d in range(2, 9):
        states = occupation_sector([1] * d, 1)
        assert len(states) == d
        assert projected_hop_count(states) == 0
        assert rank([[s[j] - states[0][j] for j in range(d)]
                     for s in states[1:]]) == d - 1
    even_two = [s for s in product((0, 1), repeat=2) if sum(s) % 2 == 0]
    assert even_two == [(0, 0), (1, 1)]
    assert len(even_two) == 2
    assert projected_hop_count(even_two) == 0
    examples = [
        ("Ising", 2, F(1, 2)),
        ("three-state Potts", 3, F(4, 5)),
        ("four-state Potts", 4, F(1)),
        ("spin-1 TB", 3, su_c(2, 2)),
        ("spin-3/2 TB", 4, su_c(2, 3)),
        ("gapless t-J", 3, F(2)),
        ("two-mode Hubbard regime", 4, F(2)),
        ("SU(4)_1 on spinful site", 4, su_c(4, 1)),
    ]
    print("Example arithmetic (IR identifications are inputs, not numeric results):")
    for name, d, c in examples:
        assert c <= d - 1
        print(f"  {name:28s} d={d}, c={str(c):>3s}, d-1={d-1}")
    print("All tested hard-projected single-parton hopping operators vanish.")


def main():
    check_weighted_family()
    check_rectangular_family()
    check_other_families()
    check_explicit_singlet_spaces()
    check_local_projection_and_examples()
    print("PASS. No continuum limit, parent Hamiltonian, or universal "
          "interacting-chain bound was numerically established.")


if __name__ == "__main__":
    main()
