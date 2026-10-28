"""Exact checks of the Laurent-polynomial source construction, standard library only."""
from fractions import Fraction as Q
from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class G:
    re: Q = Q(0)
    im: Q = Q(0)

    def __add__(self, other):
        return G(self.re + other.re, self.im + other.im)

    def __mul__(self, other):
        return G(self.re * other.re - self.im * other.im,
                 self.re * other.im + self.im * other.re)

    def scale(self, scalar):
        return G(self.re * scalar, self.im * scalar)

    def conj(self):
        return G(self.re, -self.im)


ZERO, ONE = G(), G(Q(1))


def multiply(p, q):
    out = {}
    for r, a in p.items():
        for s, b in q.items():
            out[r + s] = out.get(r + s, ZERO) + a * b
    return out


def evaluate(p, z):
    radius = max(abs(r) for r in p)
    positive = [ONE]
    negative = [ONE]
    for _ in range(radius):
        positive.append(positive[-1] * z)
        negative.append(negative[-1] * z.conj())
    out = ZERO
    for r, c in p.items():
        out = out + c * (positive[r] if r >= 0 else negative[-r])
    return out


def point(t):
    return G(Q(1 - t * t, 1 + t * t), Q(2 * t, 1 + t * t))


def check(n):
    # Distinct rational points on the unit circle, with increasingly close
    # spacings. All operations and cardinality assertions are exact.
    points = [point(t) for t in range(-(n // 2), n - n // 2)]
    for z in points:
        assert z * z.conj() == ONE
    for a, za in enumerate(points):
        polynomial = {0: ONE}
        denominator = Q(1)
        for b, zb in enumerate(points):
            if a == b:
                continue
            # sin^2((k-k_b)/2) = 1/2 - (z/z_b + z_b/z)/4.
            factor = {-1: zb.scale(Q(-1, 4)),
                      0: G(Q(1, 2)),
                      1: zb.conj().scale(Q(-1, 4))}
            polynomial = multiply(polynomial, factor)
            value = evaluate(factor, za)
            assert value.im == 0 and value.re > 0
            denominator *= value.re
        polynomial = {r: c.scale(1 / denominator)
                      for r, c in polynomial.items()}
        assert max(abs(r) for r in polynomial) == n - 1
        assert all(polynomial[-r] == c.conj()
                   for r, c in polynomial.items())
        for b, zb in enumerate(points):
            assert evaluate(polynomial, zb) == (ONE if a == b else ZERO)
    return {"branches": n, "range": n - 1,
            "cardinal_entries_checked": n * n,
            "arithmetic": "exact Gaussian rational", "passed": True}


if __name__ == "__main__":
    results = [check(n) for n in (2, 4, 8, 12, 24, 46)]
    text = json.dumps(results, indent=2) + "\n"
    Path(__file__).with_name("interpolation_checks.json").write_text(text)
    print(text, end="")
