"""Independent free kinetic/twist and harmonic checks of DENSITY_WINDOW."""
import itertools
import json
import math
from pathlib import Path


def circle(x):
    return abs((x + math.pi) % (2 * math.pi) - math.pi)


def main():
    kinetic_worst_error = 0.0
    kinetic_cases = 0
    for L in range(2, 65):
        for N in range(L + 1):
            exact = -2 * math.sin(math.pi * N / L) / math.sin(math.pi / L)
            parity = math.pi if N % 2 == 0 else 0.0
            for internal_twist in [0.0, .13, math.pi / 3, math.pi, 5.9]:
                levels = sorted(-2 * math.cos((2 * math.pi * m + parity + internal_twist) / L)
                                for m in range(L))
                E = sum(levels[:N])
                assert E >= exact - 2e-12
                if internal_twist == 0:
                    kinetic_worst_error = max(kinetic_worst_error, abs(E-exact))
                kinetic_cases += 1
    assert kinetic_worst_error < 2e-12
    # Formula checks over a deliberately nonuniform, nonresonant branch set.
    points = [math.sqrt(p) for p in [2, 3, 5, 7, 11, 13]]
    pairs = [points[i] + points[j]
             for i in range(len(points)) for j in range(i, len(points))]
    delta = min(circle(a-b) for a, b in itertools.combinations(points, 2))
    delta = min(delta, min(circle(a-b) for a,b in itertools.combinations(pairs, 2)))
    assert delta > 1e-6
    payloads = [(a-b, 1) for a,b in itertools.permutations(points, 2)]
    payloads += [(a-b, 2) for a,b in itertools.permutations(pairs, 2)]
    payloads += [(0, d) for d in [0, 1, 2]]
    margins = []
    harmonic_cases = 0
    for b in [1, 2, 5, 24]:
        for rho in [19/40, .5, 21/40]:
            bound = min(19*delta/(40*b), 2*math.pi*(1-21*math.sqrt(3)/40))
            for q, dim in payloads:
                for ell in range(-2*b-5, 2*b+6):
                    Q = q + 2*math.pi*ell
                    dimension = dim + (Q/(2*math.pi*b))**2
                    if dimension >= 3 or abs(Q) < 1e-14:
                        continue
                    mismatch = circle(rho*Q/b)
                    assert mismatch >= bound - 1e-13
                    margins.append(mismatch-bound)
                    harmonic_cases += 1
    out = {"status": "passed independent numerical diagnostics",
           "kinetic_twist_cases": kinetic_cases,
           "kinetic_equality_max_error": kinetic_worst_error,
           "harmonic_cases": harmonic_cases,
           "harmonic_bound_min_slack": min(margins),
           "thermodynamic_max_density_deviation_at_budget": math.sqrt(3*math.pi/16000),
           "window_half_width": 1/40}
    Path(__file__).with_name('density_review_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
