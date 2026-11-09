"""Exact rational checks for the half-line binding theorem; stdlib only."""
from fractions import Fraction as F
import json


WITNESS = [10000, -6789, -5763, 8811, 975, -8294, 3047, 5872,
           -5424, -2571, 5885, -535, -4742, 2634, 2711, -3374,
           -631, 2896, -809, -1716, 1276, 507, -860, 157]


def forms(x, u, v, g):
    x = [F(a) for a in x]
    norm = sum(a*a for a in x)
    eth = -4-u*u/2
    energy = (4*u*sum(a*b for a, b in zip(x, x[1:]))
              + 4*sum(a*b for a, b in zip(x, x[2:]))
              + (-2+v)*x[0]**2-g*x[1]**2)
    # Include the first zero site after the support: A f is nonzero there.
    ext = [F(0)] + x + [F(0), F(0)]
    residual = [ext[n-1]+u/2*ext[n]+ext[n+1]
                for n in range(1, len(x)+2)]
    factored = 2*sum(a*a for a in residual)+v*x[0]**2-g*x[1]**2
    vc = g*u*u/(2*(2-g))
    completed = (2*sum(a*a for a in residual[1:])
                 + (2-g)*(x[1]+u/(2-g)*x[0])**2
                 + (v-vc)*x[0]**2)
    assert energy-eth*norm == factored == completed
    return energy/norm, eth, vc


def main():
    cases = 0
    vectors = [WITNESS, [1, 0], [0, 1], [3, -1, 4, 2, -7],
               [1, 2, 3, 4, 5, 6, 7]]
    for u in [F(1, 3), F(1), F(2), F(7, 2)]:
        for g in [F(1, 5), F(2, 3), F(3, 2)]:
            vc = g*u*u/(2*(2-g))
            for v in [F(0), vc, vc+1]:
                for x in vectors:
                    e, eth, _ = forms(x, u, v, g)
                    if v >= vc:
                        assert e >= eth
                    cases += 1
    e, eth, vc = forms(WITNESS, F(1), F(0), F(2, 3))
    gap = eth-e
    assert gap > F(1, 100)
    print(json.dumps({
        "status": "all exact checks passed",
        "rational_identity_cases": cases,
        "support_relative_separations": [2, len(WITNESS)+1],
        "norm_squared": sum(x*x for x in WITNESS),
        "rayleigh_quotient_exact": str(e),
        "rayleigh_quotient_decimal": float(e),
        "threshold_exact": str(eth),
        "binding_lower_bound_exact": str(gap),
        "binding_lower_bound_decimal": float(gap),
        "certified_binding_greater_than": "1/100",
        "critical_v_exact": str(vc),
        "proof_scope": "two-body second-order effective model at fixed total momentum"
    }, indent=2))


if __name__ == "__main__":
    main()
