"""Finite exact algebra checks for C72. This is not a GN tree construction."""

from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
import json


def clean(p):
    return {e: Q(c) for e, c in p.items() if c}


def add(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def scale(p, a):
    return clean({e: a * c for e, c in p.items()})


def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            ef = tuple(x + y for x, y in zip(e, f))
            out[ef] = out.get(ef, Q(0)) + c * d
    return clean(out)


def derivative(p, beta):
    out = {}
    for e, c in p.items():
        if all(x >= y for x, y in zip(e, beta)):
            f = tuple(x - y for x, y in zip(e, beta))
            a = c
            for x, y in zip(e, beta):
                a *= factorial(x) // factorial(x - y)
            out[f] = a
    return clean(out)


def at_zero(p, indices):
    return {e: c for e, c in p.items() if all(e[i] == 0 for i in indices)}


def remainder(p, r, indices):
    return {e: c for e, c in p.items() if sum(e[i] for i in indices) >= r}


def multiindices(n, order):
    return [e for e in product(range(order + 1), repeat=n) if sum(e) == order]


def multi_factorial(e):
    ans = 1
    for x in e:
        ans *= factorial(x)
    return ans


def monomial(e):
    return {tuple(e): Q(1)}


def flag(p, r, s_indices, n):
    if r == 0:
        return p
    c_indices = tuple(i for i in range(n) if i not in s_indices)
    ans = remainder(p, r, s_indices)
    for order in range(r):
        for alpha in multiindices(n, order):
            if any(alpha[i] for i in c_indices):
                continue
            coefficient = at_zero(derivative(p, alpha), s_indices)
            coefficient = remainder(coefficient, r - order, c_indices)
            ans = add(ans, scale(mul(monomial(alpha), coefficient),
                                 Q(1, multi_factorial(alpha))))
    return ans


def integral_remainder(p, r, indices, n):
    """Evaluate (2) monomial by monomial, including the beta integral."""
    if r == 0:
        return p
    ans = {}
    other = tuple(i for i in range(n) if i not in indices)
    for gamma in multiindices(n, r):
        if any(gamma[i] for i in other):
            continue
        dp = derivative(p, gamma)
        integrated = {}
        for e, c in dp.items():
            k = sum(e[i] for i in indices)
            # Integral t^k (1-t)^(r-1) dt.
            beta_integral = Q(factorial(k) * factorial(r - 1),
                              factorial(k + r))
            integrated[e] = c * beta_integral * Q(r, multi_factorial(gamma))
        ans = add(ans, mul(monomial(gamma), integrated))
    return ans


counts = dict(flag_identities=0, integral_remainders=0,
              moment_first_identities=0, factor_normal_forms=0,
              partial_commutations=0, local_counter_cases=0,
              cross_component_cancellations=0, obstruction_sequence=0,
              nested_flag_reductions=0, nested_flag_moment_commutations=0,
              retained_pair_interpolations=0)

for n in (1, 2, 3):
    p = {e: Q((-1) ** sum(e) * (sum(e) + 1),
              1 + sum((i + 1) * x for i, x in enumerate(e)))
         for e in product(range(4), repeat=n)}
    all_indices = tuple(range(n))
    betas = [e for b in range(3) for e in multiindices(n, b)]
    powers = [e for e in product(range(3), repeat=n) if sum(e) <= 4]
    for mask in product((False, True), repeat=n):
        s_indices = tuple(i for i, yes in enumerate(mask) if yes)
        for r in range(3):
            assert flag(p, r, s_indices, n) == remainder(p, r, all_indices)
            counts['flag_identities'] += 1
            assert integral_remainder(p, r, s_indices, n) == remainder(p, r, s_indices)
            counts['integral_remainders'] += 1
        for beta in betas:
            dp = derivative(p, beta)
            lhs = derivative(add(p, scale(at_zero(p, s_indices), -1)), beta)
            rhs = dp if any(beta[i] for i in s_indices) else add(
                dp, scale(at_zero(dp, s_indices), -1))
            assert lhs == rhs
            counts['partial_commutations'] += 1
            b = sum(beta)
            for m in range(3):
                lhs = derivative(remainder(p, m, all_indices), beta)
                rhs = flag(dp, max(0, m - b), s_indices, n)
                assert lhs == rhs
                counts['moment_first_identities'] += 1
                for power in powers:
                    pa = monomial(power)
                    a = sum(power)
                    lhs = derivative(remainder(mul(pa, p), m, all_indices), beta)
                    rhs = {}
                    for delta in product(*(range(x + 1) for x in beta)):
                        pd = derivative(pa, delta)
                        if not pd:
                            continue
                        d = sum(delta)
                        remaining_beta = tuple(x - y for x, y in zip(beta, delta))
                        r = max(0, m - b - a + d)
                        weight = 1
                        for x, y in zip(beta, delta):
                            weight *= comb(x, y)
                        h = derivative(p, remaining_beta)
                        rhs = add(rhs, scale(mul(pd, flag(h, r, s_indices, n)), weight))
                    assert lhs == rhs, (n, mask, beta, m, power)
                    counts['factor_normal_forms'] += 1

# Field-load check with four scalar labeled ports. This tests the local
# counter only, without treating one-step closure as closure under contraction.
for powers in product(range(3), repeat=4):
    a = sum(powers)
    for b in range(3):
        for m in range(3):
            for delta in product(range(3), repeat=4):
                d = sum(delta)
                if d > b or any(x > y for x, y in zip(delta, powers)):
                    continue
                r = max(0, m - b - a + d)
                B, A, J = b - d + r, a - d + r, a + b - d + r
                assert B <= 2 and A <= max(a, 2) and J >= m
                for new_degree in multiindices(4, r):
                    final = tuple(x - y + z for x, y, z in zip(powers, delta, new_degree))
                    assert max(final) <= 2
                counts['local_counter_cases'] += 1

for degree in range(2, 11):
    p = {(j, k): Q((-1) ** j * (j + k + 1), j + 2 * k + 1)
         for j in range(degree + 1) for k in range(degree + 1)}
    for b in (1, 2):
        pc = derivative(p, (b, 0))
        first = remainder(pc, 2, (1,))
        second = remainder(at_zero(pc, (1,)), max(0, 2 - b), (0,))
        third = mul(monomial((0, 1)), remainder(
            at_zero(derivative(p, (b, 1)), (1,)), max(0, 1 - b), (0,)))
        full = add(first, second, third)
        expected = remainder(pc, max(0, 2 - b), (0, 1))
        assert full == expected
        if b == 1:
            assert add(first, third) == remainder(pc, 1, (1,))
        else:
            assert full == pc
        counts['cross_component_cancellations'] += 1

# For K_N=N^-2 exp(iN(c+s)), each derivative has absolute magnitude N^(k-2).
# The coefficient in (22) has magnitude N. No floating-point exponential needed.
for N in range(1, 101):
    for order in range(3):
        assert Q(N ** order, N ** 2) <= 1
    assert Q(N ** 3, N ** 2) == N
    counts['obstruction_sequence'] += 1


def reduce_flags(flags):
    kept = {}
    for subset, order in flags:
        if order and (order not in kept or len(subset) < len(kept[order])):
            kept[order] = subset
    if 1 in kept and 2 in kept and set(kept[2]).issubset(kept[1]):
        del kept[1]
    return [(subset, order) for order, subset in kept.items()]


def apply_flags(p, flags):
    for subset, order in flags:
        p = remainder(p, order, subset)
    return p


for n in (1, 2, 3):
    p = {e: Q(sum(e) + 1, multi_factorial(e))
         for e in product(range(4), repeat=n)}
    betas = [e for b in range(3) for e in multiindices(n, b)]
    for length in range(1, 41):
        for seed in range(4):
            # Every subset is an initial segment, hence all are nested.
            flags = [(tuple(range(1 + (j + seed) % n)),
                      1 + ((j // n + seed) % 2)) for j in range(length)]
            full = apply_flags(p, flags)
            reduced = reduce_flags(flags)
            assert len(reduced) <= 2
            assert full == apply_flags(p, reduced)
            counts['nested_flag_reductions'] += 1
            for beta in betas:
                lhs = derivative(full, beta)
                changed = [(subset, max(0, m - sum(beta[i] for i in subset)))
                           for subset, m in flags]
                reduced_changed = reduce_flags(changed)
                rhs = apply_flags(derivative(p, beta), reduced_changed)
                assert lhs == rhs
                assert len(reduced_changed) <= 2
                counts['nested_flag_moment_commutations'] += 1
    for size_b in range(1, n):
        for size_a in range(size_b + 1, n + 1):
            B = tuple(range(size_b))
            C = tuple(range(size_b, size_a))
            A = B + C
            lhs = apply_flags(p, [(A, 2), (B, 1)])
            rhs = remainder(p, 2, B)
            for j in B:
                ej = tuple(int(i == j) for i in range(n))
                piece = remainder(at_zero(derivative(p, ej), B), 1, C)
                rhs = add(rhs, mul(monomial(ej), piece))
            assert lhs == rhs
            counts['retained_pair_interpolations'] += 1

result = dict(status='passed', arithmetic='exact rational and integer', **counts,
              scope='Local projector and interpolation identities; no global GN routing test.')
Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
