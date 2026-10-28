from fractions import Fraction as Q
from math import comb, factorial
import json

# Ascending coefficients, with exact arithmetic independent of routing labels.
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

def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)

def deriv(p, b):
    for _ in range(b):
        p = [i * x for i, x in enumerate(p)][1:] or [Q(0)]
    return trim(p)

def rem(p, m):
    return trim([Q(0) if i < m else x for i, x in enumerate(p)])

normal_forms = 0
for m in range(3):
    for b in range(5):
        for a in range(13):
            pa = [0] * a + [1]
            for degree in range(9):
                f = [Q((-1) ** i * (i + 1), i + 2)
                     for i in range(degree + 1)]
                lhs = deriv(rem(mul(pa, f), m), b)
                rhs = [Q(0)]
                for d in range(min(a, b) + 1):
                    r = max(0, m - b - a + d)
                    term = mul(deriv(pa, d), rem(deriv(f, b - d), r))
                    rhs = add(rhs, [comb(b, d) * x for x in term])
                assert lhs == rhs, (m, b, a, degree, lhs, rhs)
                normal_forms += 1

local_counts = 0
for m in range(3):
    for b in range(3):
        for a in range(17):
            for d in range(min(a, b) + 1):
                r = max(0, m - b - a + d)
                assert a + b - d + r == max(m, a + b - d)
                assert a + b - d + r >= m
                assert 0 <= b - d + r <= 2
                local_counts += 1

# Non-chain tree. A source may end at any ancestor; endpoint sources have scale 0.
parent = {'u': 'root', 'v': 'root', 'w': 'u',
          'e1': 'w', 'e2': 'w', 'e3': 'u', 'e4': 'v'}
base = {'root': -12, 'u': -8, 'v': -6, 'w': -3,
        'e1': 0, 'e2': 0, 'e3': 0, 'e4': 0}
paths = [('e1', 'u'), ('e2', 'root'), ('e3', 'root'),
         ('e4', 'root'), ('w', 'u'), ('u', 'root')]
flow = {v: 0 for v in parent}
for source, sink in paths:
    node = source
    while node != sink:
        flow[node] += 1
        node = parent[node]
cut_cases = 0
for stretch in range(1, 101):
    scale = {k: v * stretch for k, v in base.items()}
    direct = sum(scale[sink] - scale[source] for source, sink in paths)
    cuts = -sum(flow[v] * (scale[v] - scale[parent[v]]) for v in parent)
    assert direct == cuts
    cut_cases += 1

# Verify the cancellation between the two integrals in PROOF (8).
# For a monomial q^k, coefficients of q^(k-2) are integrated exactly.
integral_cases = 0
for k in range(2, 33):
    first = Q(2 * k * (k - 1), k)
    second = Q(k * (k - 1) * (k - 2), k)
    assert first + second == k * (k - 1)
    integral_cases += 1

out = dict(normal_form_polynomial_cases=normal_forms,
           local_count_cases=local_counts,
           branched_cut_cases=cut_cases,
           taylor_integral_cases=integral_cases,
           arithmetic='exact rational and integer',
           status='passed',
           scope='finite algebra checks, not a routing construction or RG theorem')
Path = __import__('pathlib').Path
Path('checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
