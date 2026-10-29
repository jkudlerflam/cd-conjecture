"""Exact finite-g,K transport and all-length midpoint selector at n=16."""
import json
import math
from fractions import Fraction as F
import numpy as np
from rational_bands import integer_ldl


def atan_bounds(x, terms=60):
    s = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    t = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(s, s+t), max(s, s+t)


def cos_bounds(x, terms=30):
    s = t = F(1)
    for j in range(1, terms):
        t = -t*x*x/F((2*j-1)*(2*j))
        s += t
    t = -t*x*x/F((2*terms-1)*(2*terms))
    return min(s, s+t), max(s, s+t)


class Interval:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)

    def __add__(self, x):
        if not isinstance(x, Interval):
            x = Interval(x)
        return Interval(self.lo+x.lo, self.hi+x.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, x):
        return self+-x

    def __rsub__(self, x):
        return -self+x

    def __mul__(self, x):
        assert x >= 0
        return Interval(self.lo*x, self.hi*x)
    __rmul__ = __mul__

    def __truediv__(self, x):
        return self*F(1, x)


def main():
    n, gp, kp, lp = 16, 288, 1280, 704
    a = np.load('derivative_n16.npz', allow_pickle=True)
    R, ds = integer_ldl(a['F_numerator']+a['F_numerator'].T)
    bits = int(a['bits'])
    powers = [(int(d).bit_length()-1)//2 for d in ds]
    trace = sum((F((1 << (2*bits))*sum(int(x)**2 for x in R[i]),
                   1 << (2*powers[i])) for i in range(n)), F(0))
    det = int(a['detC'])
    assert F(64, 18**6) > F(1, 1 << 20)
    unit = F(1, (1 << 20)*det*det*(1 << (2*gp)))
    conversion = 4*n*det*det*trace/F(1 << (2*gp))
    penalty = 16*n*trace*F(1 << 64, 1 << kp)/unit
    budget = F(1, 1 << 35)
    assert conversion+penalty < budget/8
    # Only a finite two-cell diagnostic, not a uniform matching theorem.
    fourth = trace*F(1 << 50, 1 << (2*lp))/unit
    ab, cd = atan_bounds(F(1, 5)), atan_bounds(F(1, 239))
    pi = Interval(16*ab[0]-4*cd[1], 16*ab[1]-4*cd[0])
    def cospi(d, p=1):
        return Interval(cos_bounds(pi.hi*p/d)[0], cos_bounds(pi.lo*p/d)[1])
    aa = lambda b: 2*(1-cospi(b+1))
    dd = lambda b: (b+1)*aa(b)-b*aa(b+1)
    mu = 2-(dd(17)+dd(16))/2
    ee = lambda b: (mu-2*cospi(b+1))/b
    delta = min((ee(16)-ee(17)).lo, (ee(18)-ee(17)).lo)
    gap = 2*(cospi(18)-cospi(18, 2))
    assert delta > F(1, 1 << 17)
    assert ee(17).hi < -F(1, 1 << 9)
    assert gap.lo > F(1, 16)
    def lg(x):
        return math.log10(x.numerator)-math.log10(x.denominator)
    out = dict(n=n, g_power=-gp, K_min_power=kp, lambda_power=-lp,
               field_normalized_power=-149, beta_denominator=str(det),
               trace_log10=lg(trace), Eunit_lower_log10=lg(unit),
               conversion_error_log10=lg(conversion), penalty_error_log10=lg(penalty),
               budget_power=-35, combined_below_one_eighth_budget=True,
               finite_two_cell_fourth_order_ratio_log10=lg(fourth),
               selector_density_gap_lower=float(delta),
               selector_density_bound_power=-17, selected_density_upper=float(ee(17).hi),
               selected_density_bound_power=-9, mu_sel=[float(mu.lo), float(mu.hi)],
               vacuum_gap_lower=float(gap.lo))
    exact = {key: [str(val.numerator), str(val.denominator)]
             for key, val in [('trace', trace), ('conversion', conversion),
                              ('penalty', penalty), ('budget', budget),
                              ('selector_gap_certified_lower', F(1, 1 << 17))]}
    with open('budget_n16.json', 'w') as f:
        json.dump(out, f, indent=2)
    with open('budget_exact_n16.json', 'w') as f:
        json.dump(exact, f, indent=2)
    print(json.dumps(out))


if __name__ == '__main__':
    main()
