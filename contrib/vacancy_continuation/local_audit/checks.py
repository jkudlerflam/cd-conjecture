"""Independent normalization diagnostics, not a lattice CFT test."""
import json
import math
from pathlib import Path


def main():
    t = 1.0
    rho = 0.5
    gap_density = 1 / rho - 1
    v = 2 * t * rho * math.sin(math.pi * rho)
    K = 1 / rho**2
    epp = 2 * math.pi * t * rho**3 * math.sin(math.pi * rho)
    alpha = lambda r: (1 / (1 + r)
                       - math.sin(math.pi / (1 + r))**2 * (1 + r) / math.pi**2)
    a = -(0.25 + 1 / math.pi**2)
    step = 1e-5
    numerical_a = (alpha(gap_density + step) - alpha(gap_density - step)) / (2 * step)
    assert abs(numerical_a - a) < 1e-9
    assert abs(K / (math.pi * v) - 1 / epp) < 1e-14
    maximum = 0.0
    for k in [10.0**i for i in range(-10, 3)]:
        for omega in [0.0] + [10.0**i for i in range(-10, 3)]:
            normalized = v**2 * k**2 / (omega**2 + v**2 * k**2)
            assert 0 <= normalized <= 1
            maximum = max(maximum, normalized)
    C = a*a/(2*epp)
    canonical_errors = []
    for n in [2, 3, 7, 24, 101]:
        source = [((j*7) % 11 - 5) / 5 for j in range(n)]
        lhs = -sum(source[i]*C*((1 if i == j else 0)-1/n)*source[j]
                   for i in range(n) for j in range(n))
        rhs = -C*sum(s*s for s in source) + C*sum(source)**2/n
        canonical_errors.append(abs(lhs-rhs))
    assert max(canonical_errors) < 1e-12
    out = {"status": "passed diagnostics; continuum assumptions not tested",
           "K": K, "v": v, "e_gap_second_derivative": epp,
           "alpha": alpha(gap_density), "alpha_derivative_exact": a,
           "alpha_derivative_finite_difference": numerical_a,
           "C_times_t": C*t, "max_normalized_dynamic_response": maximum,
           "max_canonical_quadratic_form_error": max(canonical_errors)}
    Path(__file__).with_name("checks.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
