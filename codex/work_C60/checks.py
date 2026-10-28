"""Small controls for C60's exact free-pair and degree-count formulas."""
import itertools
import json
import math
from pathlib import Path


def fermi_kernel(x, p):
    return p / math.pi if x == 0 else math.sin(p * x) / (math.pi * x)


def wick_pair(x, p):
    return fermi_kernel(x, p) ** 2 - fermi_kernel(x - 1, p) * fermi_kernel(x + 1, p)


def closed_pair(x, p):
    return (math.sin(p) ** 2 / (x * x - 1)
            - math.sin(p * x) ** 2 / (x * x * (x * x - 1))) / math.pi ** 2


fillings = [math.pi / 3, math.pi / 5, 0.83]
errors = [abs(wick_pair(x, p) - closed_pair(x, p))
          for p in fillings for x in range(2, 201)]
assert max(errors) < 1e-14

p = math.pi / 3
asymptotic = 3 / (4 * math.pi ** 2)
tail = []
for x in [3, 30, 300, 3000]:
    # Evaluate the exact remainder without subtracting two leading terms.
    remainder = (math.sin(p) ** 2 - math.sin(p * x) ** 2) / (
        math.pi ** 2 * x * x * (x * x - 1))
    scaled = x ** 4 * remainder
    predicted = asymptotic * x * x / (x * x - 1)
    assert abs(scaled - predicted) < 1e-13
    tail.append({"x": x, "x4_remainder": scaled, "nonzero_limit": asymptotic})

degree_cases = 0
for n in range(1, 5):
    for half_degrees in itertools.product(range(1, 8), repeat=n):
        total = sum(half_degrees)
        exponent = sum(max(1, m - 1) for m in half_degrees)
        for contractions in range(n - 1, total):
            output_half_degree = total - contractions
            assert exponent >= max(1, output_half_degree - 1)
            degree_cases += 1

# Direct rational identity for one fixed-cutoff absorption step.
shell_errors = []
for outer, inner in [(1.0, 0.2), (1.0, 1.0), (0.7, 0.0), (0.0, 0.0)]:
    for z in [0.0, 0.02, -0.03, 0.01 + 0.02j]:
        denominator = 0.31 - 0.29j
        difference = ((1 + z) * outer / (denominator * (1 + z * outer))
                      - inner / denominator)
        shell = (1 + z) * (outer - inner) / (denominator * (1 + z * outer))
        shell_errors.append(abs(difference - shell))
assert max(shell_errors) < 1e-14

result = {
    "status": "passed",
    "wick_closed_max_absolute_error": max(errors),
    "free_pair_tail": tail,
    "degree_inequality_cases": degree_cases,
    "fixed_cutoff_shell_identity_max_absolute_error": max(shell_errors),
    "scope": "Exact-formula controls only. No interacting RG or routing test."
}
Path(__file__).with_name("checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
