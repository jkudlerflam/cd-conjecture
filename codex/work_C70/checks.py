"""Exact finite algebra controls for C70, not a GN routing search."""
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import json


def clean(p):
    return {k: Q(v) for k, v in p.items() if v}


def add(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def scale(p, c):
    return clean({e: c * v for e, v in p.items()})


def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            key = tuple(a + b for a, b in zip(e, f))
            out[key] = out.get(key, Q(0)) + c * d
    return clean(out)


def prod(ps, n):
    out = {(0,) * n: Q(1)}
    for p in ps:
        out = mul(out, p)
    return out


def power(p, k, n):
    return prod([p] * k, n)


def deriv(p, i, order=1):
    for _ in range(order):
        out = {}
        for e, c in p.items():
            if e[i]:
                f = list(e)
                f[i] -= 1
                out[tuple(f)] = c * e[i]
        p = out
    return clean(p)


def affine(offset, slopes):
    n = len(slopes)
    out = {(0,) * n: Q(offset)}
    for i, c in enumerate(slopes):
        e = [0] * n
        e[i] = 1
        out[tuple(e)] = Q(c)
    return clean(out)


def evaluate(p, x):
    return sum((c * scalar_product(v ** k for v, k in zip(x, e))
                for e, c in p.items()), Q(0))


def scalar_product(xs):
    ans = Q(1)
    for x in xs:
        ans *= x
    return ans


def radial_integral(p, r):
    # x_i=u*r_i. Integrate u from zero to one.
    return sum((c * scalar_product(v ** k for v, k in zip(r, e))
                / (sum(e) + 1) for e, c in p.items()), Q(0))


def restrict_zero(p, indices):
    return {e: c for e, c in p.items() if all(e[i] == 0 for i in indices)}


def remainder(p, m):
    return {e: c for e, c in p.items() if sum(e) >= m}


product_cases = integral_cases = exterior_nonzero = 0
for n in (1, 2, 3):
    for b in range(3):
        for q in range(1, 5):
            zslopes = [Q((-1) ** i, i + 1) for i in range(n)]
            z = affine(Q(2, 3), zslopes)
            g = add(affine(1, [0] * n), z, scale(power(z, 2, n), Q(1, 3)))
            gp = add(affine(1, [0] * n), scale(z, Q(2, 3)))
            weight = mul(power(z, b, n), g)
            ys = [affine(Q(j + 1, j + 2),
                         [Q((i + j) % 2 + 1, 3 * (j + 1)) for i in range(n)])
                  for j in range(q)]
            # These can stand for already differentiated exterior fields.
            fields = [add(affine(j + 1, [0] * n), power(y, 2, n))
                      for j, y in enumerate(ys)]
            aa = affine(1, [Q(1, i + 2) for i in range(n)])
            ab = affine(2, [Q(-1, i + 3) for i in range(n)])
            ba = affine(-1, [Q(2, i + 4) for i in range(n)])
            bb = affine(3, [Q(-2, i + 5) for i in range(n)])
            determinant = add(mul(aa, bb), scale(mul(ab, ba), -1))
            whole = prod([weight, determinant] + fields, n)
            r = [Q((-1) ** i * (i + 2), i + 3) for i in range(n)]
            integral = Q(0)
            for i in range(n):
                covterm = prod([scale(mul(power(z, b, n), gp), zslopes[i]),
                                determinant] + fields, n)
                consumed = {} if b == 0 else prod([
                    scale(mul(power(z, b - 1, n), g), b * zslopes[i]),
                    determinant] + fields, n)
                # Determinant derivative kept as four product terms.
                detprime = add(mul(deriv(aa, i), bb), mul(aa, deriv(bb, i)),
                               scale(mul(deriv(ab, i), ba), -1),
                               scale(mul(ab, deriv(ba, i)), -1))
                detterm = prod([weight, detprime] + fields, n)
                extterms = []
                for j in range(q):
                    changed = list(fields)
                    changed[j] = deriv(changed[j], i)
                    extterms.append(prod([weight, determinant] + changed, n))
                rhs = add(covterm, consumed, detterm, *extterms)
                assert rhs == deriv(whole, i)
                if add(*extterms):
                    exterior_nonzero += 1
                integral += r[i] * radial_integral(rhs, r)
                product_cases += 1
            assert integral == evaluate(whole, r) - evaluate(whole, [Q(0)] * n)
            integral_cases += 1


# A transparent missing-exterior-term test: x^2(1+x)*(1+v*x)^2.
n = 1
x = affine(0, [1])
v = Q(1, 2)
w = mul(power(x, 2, n), add(affine(1, [0]), x))
phi = power(affine(1, [v]), 2, n)
full = mul(w, phi)
omitted = radial_integral(mul(deriv(w, 0), phi), [Q(1)])
correct = evaluate(full, [Q(1)]) - evaluate(full, [Q(0)])
missing = radial_integral(mul(w, deriv(phi, 0)), [Q(1)])
assert correct - omitted == missing and missing != 0


survival_cases = homogeneous_cases = 0
for n in (1, 2, 3):
    p = {e: Q((-1) ** sum(e) * (sum(e) + 1), 1 + sum((i + 1) * k for i, k in enumerate(e)))
         for e in product(range(4), repeat=n)}
    for mask in product((False, True), repeat=n):
        surv = [i for i, yes in enumerate(mask) if yes]
        contracted = [i for i, yes in enumerate(mask) if not yes]
        k0s = restrict_zero(p, surv)
        survpart = add(p, scale(k0s, -1))
        const = restrict_zero(p, range(n))
        contractpart = add(k0s, scale(const, -1))
        assert add(survpart, contractpart) == remainder(p, 1)
        # Integrating q_j partial_j K(q_C,u*q_S) divides by S-degree.
        factored = {}
        for e, c in p.items():
            degree = sum(e[j] for j in surv)
            if degree:
                factored[e] = sum((c * Q(e[j], degree) for j in surv), Q(0))
        assert clean(factored) == survpart
        assert restrict_zero(survpart, surv) == {}
        survival_cases += 1
    for a in range(5):
        pa = {(a,) + (0,) * (n - 1): Q(1)}
        for m in (1, 2):
            assert remainder(mul(pa, p), m) == mul(pa, remainder(p, max(0, m - a)))
            homogeneous_cases += 1


# Levels increase toward the root; source=0, first exterior level=1.
def flow(paths, edges):
    return [sum(source <= e < sink for source, sink in paths) for e in range(edges)]


cut_cases = 0
for span in product(range(1, 4), repeat=3):
    levels = [0, -span[0], -span[0] - span[1], -sum(span)]
    for b in range(3):
        for sinks in product((2, 3), repeat=b):
            old = [(1, h) for h in sinks]
            oldflow = flow(old, 3)
            # Covariance and exterior-field derivatives, including later sinks.
            for k in (1, 2, 3):
                paths = old + [(0, k)]
                f = flow(paths, 3)
                assert f[0] == 1
                assert all(f[e] >= oldflow[e] for e in range(1, 3))
                exponent = sum(levels[h] - levels[t] for t, h in paths)
                assert exponent == -sum(f[e] * span[e] for e in range(3))
                cut_cases += 1
            for consumed_index in range(b):
                h = sinks[consumed_index]
                paths = [p for j, p in enumerate(old) if j != consumed_index] + [(0, h)]
                f = flow(paths, 3)
                assert f == flow(old + [(0, 1)], 3)
                exponent = sum(levels[k] - levels[t] for t, k in paths)
                assert exponent == -sum(f[e] * span[e] for e in range(3))
                cut_cases += 1


cancel_cases = 0
for degree in range(2, 41):
    term1 = Q(2 * degree * (degree - 1), degree)
    term2 = Q(degree * (degree - 1) * (degree - 2), degree)
    assert term1 + term2 == degree * (degree - 1)
    cancel_cases += 1


out = dict(
    status='passed', arithmetic='exact rational and integer',
    affine_product_rule_cases=product_cases,
    radial_transposition_cases=integral_cases,
    nonzero_exterior_derivative_cases=exterior_nonzero,
    survival_split_cases=survival_cases,
    homogeneous_remainder_cases=homogeneous_cases,
    cut_record_cases=cut_cases,
    taylor_cancellation_cases=cancel_cases,
    missing_exterior_term_example=dict(correct=str(correct),
                                      omitted_exterior=str(omitted),
                                      difference=str(missing)),
    scope='Finite local algebra controls. No global routing construction or projected-tree counterexample.'
)
Path(__file__).with_name('checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
