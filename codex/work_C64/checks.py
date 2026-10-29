"""Finite algebra controls for C64. No full GN routing or RG calculation."""

from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import json


def trim(p):
    p = list(map(Q, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    out = [Q(0)] * max(len(p), len(q))
    for i, x in enumerate(p):
        out[i] += x
    for i, x in enumerate(q):
        out[i] += x
    return trim(out)


def scale(p, a):
    return trim([a * x for x in p])


def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)


def deriv(p, order=1):
    for _ in range(order):
        p = [i * x for i, x in enumerate(p)][1:] or [Q(0)]
    return trim(p)


def evaluate(p, x):
    result = Q(0)
    for a in reversed(p):
        result = result * x + a
    return result


def affine(p, offset, slope):
    result = [Q(0)]
    power = [Q(1)]
    for a in p:
        result = add(result, scale(power, a))
        power = mul(power, [offset, slope])
    return result


def integrate_unit(p):
    return sum((a / (i + 1) for i, a in enumerate(p)), Q(0))


leibniz_cases = 0
taylor_cases = 0
error_cases = 0
anchors = [(Q(0), Q(1)), (Q(1), Q(1)),
           (Q(-2, 3), Q(4, 5)), (Q(3, 2), Q(-7, 4)),
           (Q(1, 2), Q(0))]
for b in range(3):
    monomial = [Q(0)] * b + [Q(1)]
    for degree in range(13):
        g = [Q((-1) ** j * (j + 2), j + 1)
             for j in range(degree + 1)]
        weighted = mul(monomial, g)
        for c in range(3):
            lhs = deriv(weighted, c)
            rhs = [Q(0)]
            for d in range(min(b, c) + 1):
                coefficient = comb(c, d) * factorial(b) // factorial(b - d)
                residual_weight = [Q(0)] * (b - d) + [Q(1)]
                rhs = add(rhs, scale(mul(residual_weight, deriv(g, c - d)),
                                     coefficient))
            assert lhs == rhs, (b, c, degree)
            leibniz_cases += 1
        for z0, r in anchors:
            for c in (1, 2):
                lhs = evaluate(weighted, z0 + r)
                lhs -= sum((r ** j * evaluate(deriv(weighted, j), z0)
                            / factorial(j) for j in range(c)), Q(0))
                integrand = affine(deriv(weighted, c), z0, r)
                if c == 2:
                    integrand = mul(integrand, [Q(1), Q(-1)])
                rhs = r ** c * integrate_unit(integrand) / factorial(c - 1)
                assert lhs == rhs, (b, c, degree, z0, r)
                taylor_cases += 1

            z = z0 + r
            covariance_only = z ** b * (evaluate(g, z) - evaluate(g, z0))
            whole_weight = evaluate(weighted, z) - evaluate(weighted, z0)
            expected_error = (z0 ** b - z ** b) * evaluate(g, z0)
            assert covariance_only - whole_weight == expected_error
            error_cases += 1


cut_cases = 0
maximum_coefficient_mass = 0
for b in range(3):
    for c in (1, 2):
        coefficient_mass = sum(comb(c, d) * factorial(b) // factorial(b - d)
                               for d in range(min(b, c) + 1))
        maximum_coefficient_mass = max(maximum_coefficient_mass, coefficient_mass)
        assert coefficient_mass <= 7
        for child_gap in range(6):
            for parent_gap in range(6):
                h = -20
                s = h + parent_gap
                t = s + child_gap
                for d in range(min(b, c) + 1):
                    paths = [('t', 's')] * (c - d)
                    paths += [('t', 'h')] * d
                    paths += [('s', 'h')] * (b - d)
                    scales = dict(t=t, s=s, h=h)
                    direct = sum(scales[sink] - scales[source]
                                 for source, sink in paths)
                    child_count = sum(source == 't' for source, sink in paths)
                    parent_count = sum(sink == 'h' for source, sink in paths)
                    assert child_count == c and parent_count == b
                    assert direct == c * (s - t) + b * (h - s)
                    assert direct == -child_count * child_gap - parent_count * parent_gap
                    assert (c - d) - (b - d) == c - b
                    cut_cases += 1


# The full three-contraction matrix has determinant
# lambda*u*a*d - lambda^3*u*b*c. Differentiate by each cross edge.
bbf_cases = 0
for j in range(1, 25):
    u, a, b, c, d = (Q(j, j + 1), Q(2, j + 2), Q(-3, j + 3),
                     Q(4, j + 4), Q(-5, j + 5))
    full = [Q(0), u * a * d, Q(0), -u * b * c]
    r_tree = [u * a * d, Q(0), -u * b * c]
    left_tree_one = [Q(0), Q(0), -u * b * c]
    left_tree_two = list(left_tree_one)
    assert deriv(full) == add(add(r_tree, left_tree_one), left_tree_two)
    assert integrate_unit(r_tree) == u * a * d - u * b * c / 3
    assert sum(map(integrate_unit, [r_tree, left_tree_one, left_tree_two])) == u * (a * d - b * c)
    bbf_cases += 1


g = [Q(1), Q(1)]
weighted = mul([Q(0), Q(0), Q(1)], g)
bad_value = 4 * (evaluate(g, Q(2)) - evaluate(g, Q(1)))
correct_value = evaluate(weighted, Q(2)) - evaluate(weighted, Q(1))
assert (bad_value, correct_value, bad_value - correct_value) == (4, 10, -6)

result = {
    'status': 'passed',
    'arithmetic': 'exact rational and integer',
    'weighted_leibniz_identities': leibniz_cases,
    'anchored_taylor_identities': taylor_cases,
    'covariance_only_error_identities': error_cases,
    'two_cut_incidence_cases': cut_cases,
    'bbf_three_contraction_cases': bbf_cases,
    'maximum_local_leibniz_coefficient_mass': maximum_coefficient_mass,
    'nonzero_operator_error': {
        'covariance_only': int(bad_value),
        'whole_weight': int(correct_value),
        'difference': int(bad_value - correct_value)
    },
    'scope': 'Local identities and two-cut records, not a global GN routing theorem.'
}
Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
