"""Exact finite controls for C74. No general GN routing is tested."""
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import json


def clean(p):
    return {e: Q(a) for e, a in p.items() if a}


def add(*ps):
    ans = {}
    for p in ps:
        for e, a in p.items():
            ans[e] = ans.get(e, Q(0)) + a
    return clean(ans)


def scale(p, a):
    return clean({e: a * b for e, b in p.items()})


def mul(p, q):
    ans = {}
    for e, a in p.items():
        for f, b in q.items():
            ef = tuple(x + y for x, y in zip(e, f))
            ans[ef] = ans.get(ef, Q(0)) + a * b
    return clean(ans)


def deriv(p, i, order=1):
    for _ in range(order):
        ans = {}
        for e, a in p.items():
            if e[i]:
                f = list(e)
                f[i] -= 1
                ans[tuple(f)] = a * e[i]
        p = ans
    return clean(p)


def jet(p, m, indices):
    return {e: a for e, a in p.items() if sum(e[i] for i in indices) < m}


def rem(p, m, indices):
    return add(p, scale(jet(p, m, indices), -1))


def restrict(p, indices):
    return {e: a for e, a in p.items() if all(e[i] == 0 for i in indices)}


def contract(w, f):
    # Integrate the first coordinate over [-1,1] with measure dc/2.
    # The second coordinate is the surviving external momentum.
    ans = {}
    for (c, q), a in mul(w, f).items():
        if c % 2 == 0:
            e = (0, q)
            ans[e] = ans.get(e, Q(0)) + a / (c + 1)
    return clean(ans)


def polynomial(seed, degree):
    return {(c, q): Q((-1) ** (seed + c) * (1 + c + q),
                       1 + seed + 2 * c + 3 * q)
            for c in range(degree + 1) for q in range(degree + 1)}


counts = dict(contraction_jet_cases=0, child_parent_cases=0,
              moment_differentiated_cases=0, spatial_split_cases=0,
              field_counter_cases=0, ladder_multiplier_cases=0,
              cut_cases=0, bbf_cases=0)

for seed in range(1, 9):
    for degree in range(1, 5):
        f = polynomial(seed, degree)
        w = polynomial(seed + 2, degree + 1)
        for m in (1, 2):
            lhs = rem(contract(w, f), m, (1,))
            rhs = contract(w, rem(f, m, (1,)))
            for alpha in range(m):
                fj = restrict(deriv(f, 1, alpha), (1,))
                term = rem(contract(w, fj), m - alpha, (1,))
                rhs = add(rhs, scale(mul({(0, alpha): Q(1)}, term),
                                     Q(1, factorial(alpha))))
            assert lhs == rhs
            counts['contraction_jet_cases'] += 1
            for b in range(3):
                normalized = rem(deriv(contract(w, f), 1, b),
                                 max(0, m - b), (1,))
                assert deriv(rhs, 1, b) == normalized
                counts['moment_differentiated_cases'] += 1

        k0q = restrict(f, (1,))
        k00 = restrict(f, (0, 1))
        lhs = rem(contract(w, add(f, scale(k00, -1))), 1, (1,))
        surv = contract(w, add(f, scale(k0q, -1)))
        reset = contract(add(w, scale(restrict(w, (1,)), -1)),
                         add(k0q, scale(k00, -1)))
        assert lhs == add(surv, reset)
        counts['child_parent_cases'] += 1

# Explicit false commute control, with a nonzero integrated child jet.
f = {(2, 0): Q(1)}
w = {(0, 0): Q(1), (0, 1): Q(1)}
bad_commute = contract(w, rem(f, 1, (1,)))
correct = rem(contract(w, rem(f, 1, (0, 1))), 1, (1,))
assert bad_commute == {} and correct == {(0, 1): Q(1, 3)}

# Position-space identity (16)-(17), independently evaluated as scalar tests.
def aval(x, degree):
    return sum(Q((-1) ** j * (j + 1), j + 2) * x ** j
               for j in range(degree + 1))


for b, z, r, u in product((Q(-2, 3), Q(0), Q(4, 5)), repeat=4):
    for degree in range(1, 6):
        gr = 1 + z + z ** 2
        gl = lambda x: 2 - x + x ** 3
        a = lambda x: aval(x, degree)
        raw = gr * (gl(z + r) * a(b + z + u) - gl(z) * a(b + z))
        contact = gr * (gl(z + r) - gl(z)) * a(b)
        surviving = gr * gl(z + r) * (a(b + z + u) - a(b + z))
        reset = gr * (gl(z + r) - gl(z)) * (a(b + z) - a(b))
        assert raw - contact == surviving + reset
        counts['spatial_split_cases'] += 1

# Enumerate all retained field degree lists through six fields.
# Contracted labels can be discarded from any longer incoming list, so these
# include the possible outputs of a contraction with the induction hypothesis.
for p in (2, 4, 6):
    m = {2: 2, 4: 1, 6: 0}[p]
    for degrees in product(range(3), repeat=p):
        a = sum(degrees)
        r = max(0, m - a)
        for new in product(range(r + 1), repeat=p):
            if sum(new) != r:
                continue
            final = tuple(x + y for x, y in zip(degrees, new))
            assert max(final) <= 2
            # Independent polynomial check of homogeneity at this vertex.
            mon = {degrees: Q(1)}
            kernel = {tuple(0 for _ in range(p)): Q(2),
                      tuple(1 for _ in range(p)): Q(-3)}
            assert rem(mul(mon, kernel), m, tuple(range(p))) == mul(
                mon, rem(kernel, r, tuple(range(p))))
            counts['field_counter_cases'] += 1

# Repeated actual scalar multiplication contractions. Projection each time
# agrees with projecting once, at every depth, without evaluating derivatives.
for seed in range(1, 7):
    initial = {(0, k): Q(seed + k, k + 1) for k in range(5)}
    original = rem(initial, 1, (1,))
    iterative = original
    accumulated = {(0, 0): Q(1)}
    for depth in range(1, 41):
        bubble = {(0, 0): Q(1, 2),
                  (0, 1): Q((-1) ** depth, seed + depth + 1)}
        iterative = rem(mul(bubble, iterative), 1, (1,))
        accumulated = mul(bubble, accumulated)
        assert iterative == mul(original, accumulated)
        assert restrict(iterative, (1,)) == {}
        counts['ladder_multiplier_cases'] += 1

for length in range(1, 65):
    for seed in range(1, 9):
        gaps = [1 + (seed * (j + 1)) % 17 for j in range(length + 1)]
        scales = [0]
        for gap in gaps:
            scales.append(scales[-1] - gap)
        # Node 0 is the initial endpoint; node 1 is the first contraction.
        # The last node is the external sink.
        for paths in ([(0, length + 1)], [(0, 1), (1, length + 1)]):
            flow = [sum(start <= edge < end for start, end in paths)
                    for edge in range(length + 1)]
            assert flow == [1] * (length + 1)
            exponent = sum(scales[end] - scales[start] for start, end in paths)
            assert exponent == scales[-1] == -sum(gaps)
            assert exponent == -sum(f * gap for f, gap in zip(flow, gaps))
            counts['cut_cases'] += 1

# For two cross contractions, det(lambda G)=lambda^2 u v.
# Each of the two choices of tree line contributes integral_0^1 lambda u v.
for j in range(1, 41):
    u, v = Q(j, j + 1), Q(2 - j, j + 2)
    one_tree = u * v / 2
    assert one_tree + one_tree == u * v
    counts['bbf_cases'] += 1

out = dict(status='passed', arithmetic='exact rational and integer', **counts,
           false_commutation_example=dict(correct='q/3', omitted_jet='0'),
           maximum_checked_ladder_depth=64,
           scope='Finite algebra and ladder controls, not a general GN routing theorem.')
Path(__file__).with_name('checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
