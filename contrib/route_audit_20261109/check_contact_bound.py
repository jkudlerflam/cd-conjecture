"""Exact arithmetic checks for CONTACT_BOUND.md; standard library only.

Run from the repository root with:
    python3 -B contrib/route_audit_20261109/check_contact_bound.py

This checks the proof's elementary inequalities, not a many-body model.
The saved decimal values are diagnostic; all assertions use exact rationals
or integers except the bisection of the improved implicit scalar bound.
"""

from fractions import Fraction as F
import json
import math
from pathlib import Path


def binary_entropy(a):
    if a <= 0.0 or a >= 1.0:
        return 0.0
    return -a * math.log(a) - (1.0 - a) * math.log1p(-a)


def entropy_residual(a, c=0.25, t=1.0, temperature=None):
    if temperature is None:
        temperature = t / (16.0 * math.log(2.0))
    return (
        -c * (1.0 - a)
        + 4.0 * t * math.sqrt(a * (1.0 - a))
        + temperature * binary_entropy(a)
    )


def main():
    checks = {}

    # log(264)<6 follows from a positive Taylor lower bound on e^6.
    exp6_lower = sum((F(6**r, math.factorial(r)) for r in range(8)), F(0))
    assert exp6_lower > 264
    checks["exp6_lower"] = str(exp6_lower)
    checks["exp6_lower_minus_264"] = str(exp6_lower - 264)

    # log(2)=2*sum_{r>=0} (1/3)^(2r+1)/(2r+1)>2/3.
    log2_lower = 2 * (F(1, 3) + F(1, 81))
    assert log2_lower > F(2, 3)
    checks["log2_lower"] = str(log2_lower)

    # Upper numerator in (7) is 128 sqrt(263)-2083, strictly negative.
    square_margin = 2083**2 - 128**2 * 263
    assert square_margin == 29897
    assert square_margin > 0
    checks["contact_squared_integer_margin"] = square_margin

    # Zero-temperature inequality at c=t/4 gives exactly 1/257.
    c = F(1, 4)
    zero_temperature_bound = c * c / (c * c + 16)
    assert zero_temperature_bound == F(1, 257)
    checks["zero_temperature_contact_bound"] = str(zero_temperature_bound)

    # cos(pi/16)>=1-pi^2/512>1-(22/7)^2/512>49/50.
    cosine_lower = 1 - F(22, 7) ** 2 / 512
    assert cosine_lower > F(49, 50)
    checks["cosine_lower"] = str(cosine_lower)
    pin_margin = F(49, 6600) - F(1, 135)
    assert pin_margin == F(1, 59400)
    assert pin_margin > F(6, 2**608)
    checks["pin_margin_before_port"] = str(pin_margin)
    checks["port_upper_bound"] = "6 * 2^-608"
    checks["pin_activation_improvement"] = str(F(4096, 135))
    checks["contact_bound_improvement"] = str(F(4096, 264))

    # Useful decimal diagnostics, not used as certificates in the proof.
    lo, hi = 0.0, 0.5
    for _ in range(80):
        mid = (lo + hi) / 2
        if entropy_residual(mid) < 0:
            lo = mid
        else:
            hi = mid
    checks["diagnostic"] = {
        "residual_at_1_over_264": entropy_residual(1 / 264),
        "coarse_exact_upper_bound_at_1_over_264": (
            128 * math.sqrt(263) - 2083
        ) / (32 * 264),
        "sharp_implicit_contact_bound_in_window": (lo + hi) / 2,
        "reciprocal_sharp_implicit_bound": 2 / (lo + hi),
        "certified_contact_lower_bound": 1 / 264,
        "old_contact_lower_bound": 1 / 4096,
        "new_pin_activation_lower_bound_over_J": 1 / 135,
        "old_pin_activation_lower_bound_over_J": 1 / 4096,
    }

    result = {
        "status": "all exact arithmetic assertions passed",
        "scope": "contact and pin inequalities; no CFT or finite-K claim",
        "checks": checks,
    }
    text = json.dumps(result, indent=2) + "\n"
    (Path(__file__).parent / "contact_checks.json").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
