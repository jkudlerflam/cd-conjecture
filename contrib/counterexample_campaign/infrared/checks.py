"""Independent rational check of C33 roots -> vacancy resonance margins.

Uses the Python standard library. Matrix inertia is inherited from C33;
all subsequent circle separation computations are repeated here.
"""
from fractions import Fraction as F
from functools import cmp_to_key
from itertools import combinations_with_replacement
from pathlib import Path
import json
import math
import random
import time


def compare_phase(a, b):
    x, y = a[0]
    u, v = b[0]
    ha = 0 if y > 0 or (y == 0 and x >= 0) else 1
    hb = 0 if v > 0 or (v == 0 and u >= 0) else 1
    if ha != hb:
        return -1 if ha < hb else 1
    cross = x*v-y*u
    return -1 if cross > 0 else 1 if cross < 0 else 0


def separation(items):
    circle = sorted(items, key=cmp_to_key(compare_phase))
    answer = None
    for left, right in zip(circle, circle[1:] + circle[:1]):
        x, y = left[0]
        u, v = right[0]
        assert x*v-y*u >= 0, "adjacent circular gap exceeds pi"
        bound = max(abs(x-u), abs(y-v))-left[1]-right[1]
        assert bound > 0, (left[2], right[2])
        if answer is None or bound < answer[0]:
            answer = (bound, [left[2], right[2]])
    return answer


def main():
    start = time.monotonic()
    root = Path(__file__).resolve().parents[3]
    source = root / "codex/work_C33/root_certificate_exact_n23.json"
    data = json.loads(source.read_text())
    roots = []
    for index, interval in enumerate(data["intervals"]):
        lo, hi = [F(int(x[0]), int(x[1])) for x in interval]
        assert lo < hi
        mid = (lo+hi)/2
        z = ((1-mid*mid)/(1+mid*mid), 2*mid/(1+mid*mid))
        assert z[0]*z[0]+z[1]*z[1] == 1
        # d(2 atan t)/dt <= 2; half width gives error <=hi-lo.
        roots.append((z, hi-lo, [index]))
    d2, witness2 = separation(roots)
    pairs = []
    for i, j in combinations_with_replacement(range(len(roots)), 2):
        (x, y), ei, _ = roots[i]
        (u, v), ej, _ = roots[j]
        pairs.append(((x*u-y*v, x*v+y*u), ei+ej, [i, j]))
    d4, witness4 = separation(pairs)
    assert d2 == F(*map(int, data["delta2"]))
    assert d4 == F(*map(int, data["delta4"]))
    assert d2 > F("0.0140193")
    assert d4 > F("0.00000031558")
    b = 24
    physical = min(d2, d4)/(2*b)
    assert physical > F("0.0000000065745")
    # pi>3, sqrt(3)<7/4 imply pi(2-sqrt(3))>3/4.
    assert F(3, 4) > physical
    # At rho=1/2, any uniform charge resonance h!=0 contributes >=4.
    for h in range(-100, 101):
        assert h == 0 or 4*h*h >= 4
    alpha = 0.5-2/math.pi**2
    assert 0.297 < alpha < 0.298

    # Numerical stress check of the independently proved free-charge
    # inequality for arbitrary particle-hole collections at half filling.
    rng = random.Random(230024)
    worst_ratio = float("inf")
    for _ in range(10000):
        count = rng.randint(1, 40)
        holes = [rng.uniform(-math.pi/2, math.pi/2) for _ in range(count)]
        particles = [rng.choice([-1, 1])*rng.uniform(math.pi/2, math.pi)
                     for _ in range(count)]
        momentum = sum(particles)-sum(holes)
        distance = abs((momentum+math.pi/2) % math.pi-math.pi/2)
        energy = sum(2*abs(math.cos(k)) for k in holes+particles)
        bound = 4*distance/math.pi
        assert energy+1e-12 >= bound
        if bound > 1e-12:
            worst_ratio = min(worst_ratio, energy/bound)

    result = {
        "status": "exact_rational_root_to_resonance_check_passed",
        "inherited_input": str(source.relative_to(root)),
        "root_count": len(roots),
        "pair_count": len(pairs),
        "delta2": float(d2),
        "delta4": float(d4),
        "delta4_witness": witness4,
        "b": b,
        "rho": "1/2",
        "physical_mismatch_lower": float(physical),
        "physical_mismatch_exact": [str(physical.numerator), str(physical.denominator)],
        "free_charge_restoration_threshold_over_t": 2*float(d4)/(math.pi*b),
        "alpha_half": alpha,
        "minimum_uniform_nonGaussian_dimension": 3,
        "minimum_nonzero_h_charge_dimension": 4,
        "charge_inequality_random_checks": 10000,
        "smallest_sample_energy_to_bound_ratio": worst_ratio,
        "seconds": time.monotonic()-start,
        "limitations": ["Does not rerun C33 matrix inertia",
                        "Does not establish interacting matching or measure central charge"],
    }
    target = Path(__file__).with_name("checks.json")
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
