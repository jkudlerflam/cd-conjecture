"""Independent determinant/root/pair check of C39 exact certificates.

The original exact Gram evaluator supplies A and B; everything after that
uses an independent algorithm. In particular: permutation determinants,
endpoint signs (no Sturm implementation), and rational Cayley bounds
(no trigonometric series). Python standard library only.
"""
from fractions import Fraction as F
from itertools import permutations, combinations_with_replacement
from pathlib import Path
import importlib.util
import json
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "codex/work_C39/certify.py"
spec = importlib.util.spec_from_file_location("c39_exact_input", SOURCE)
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def neg(a):
    return (-a[1], -a[0])


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    vals = [x*y for x in a for y in b]
    return min(vals), max(vals)


def abs_lower(a):
    return F(0) if a[0] <= 0 <= a[1] else min(abs(x) for x in a)


def poly_mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return out


def eval_poly(p, x):
    total = F(0)
    for a in p[::-1]:
        total = total*x+a
    return total


def direct_polynomial(A, B, E, eta, mu):
    n = len(A)
    out = [F(0)]*(2*n+1)
    for perm in permutations(range(n)):
        inversions = sum(perm[i] > perm[j]
                         for i in range(n) for j in range(i+1, n))
        p = [F((-1)**inversions)]
        for i, j in enumerate(perm):
            p = poly_mul(p, [-2*A[i][j]+mu*E[i][j],
                             4*eta*B[i][j],
                             2*A[i][j]+mu*E[i][j]])
        out = [x+y for x, y in zip(out, p)]
    return out


def pair_numerator(x, y, z, w):
    one = (F(1), F(1))
    return sub(mul(add(x, y), sub(one, mul(z, w))),
               mul(add(z, w), sub(one, mul(x, y))))


def verify(n):
    start = time.monotonic()
    data = json.loads((HERE/f"certificate_n{n}.json").read_text())
    w, v, u, eta, mu = (F(data[s]) for s in ("w", "v", "u", "eta", "mu"))
    GA, GB = source.gram(n, w), source.gram(n, w*v)
    A = [[GA[i][j]/w**(2*i+2*j) for j in range(n)] for i in range(n)]
    B = [[GB[i][j]/w**(2*i+2*j) for j in range(n)] for i in range(n)]
    E = [[F(i == j)*u**(4*i) for j in range(n)] for i in range(n)]
    p = direct_polynomial(A, B, E, eta, mu)
    assert p == list(map(F, data["polynomial"]))
    assert p[-1] > 0  # k=pi, the Cayley point at infinity, is not a root.
    intervals = [tuple(map(F, pair)) for pair in data["root_intervals"]]
    assert len(intervals) == 2*n
    assert all(intervals[i][1] < intervals[i+1][0] for i in range(2*n-1))
    assert all(eval_poly(p, lo)*eval_poly(p, hi) < 0 for lo, hi in intervals)
    # 2n disjoint sign-changing intervals for degree 2n: one simple root
    # in each, and no other real or complex roots, counting multiplicity.
    labels = list(combinations_with_replacement(range(2*n), 2))
    lower = None
    count = 0
    for i, pair in enumerate(labels):
        for other in labels[:i]:
            x, y, z, w_ = (intervals[j] for j in (*pair, *other))
            numerator = abs_lower(pair_numerator(x, y, z, w_))
            assert numerator > 0
            # |sin((p_x+p_y-p_z-p_w)/2)| = |numerator| / sqrt(product(1+x^2)).
            # Bound sqrt(product) above by product, since every factor >=1.
            denominator_upper = F(1)
            for interval in (x, y, z, w_):
                denominator_upper *= 1+max(abs(t) for t in interval)**2
            bound = 2*numerator/denominator_upper
            lower = bound if lower is None else min(lower, bound)
            count += 1
    return {
        "n": n, "simple_fermi_points": 2*n,
        "unordered_pairs_with_replacement": len(labels),
        "independently_checked_pair_comparisons": count,
        "delta4_rational_lower": str(lower),
        "delta4_display": float(lower),
        "seconds": time.monotonic()-start,
        "checks": "exact determinant, endpoint signs, rational Cayley pair separation",
        "input_scope": "Gram matrices inherited from C39 exact evaluator",
    }


if __name__ == "__main__":
    results = []
    for n in range(1, 5):
        result = verify(n)
        results.append(result)
        print({k: v for k, v in result.items() if k != "delta4_rational_lower"}, flush=True)
    (HERE/"verification.json").write_text(json.dumps(results, indent=2)+"\n")
