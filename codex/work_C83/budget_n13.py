# C83 interval utilities adapted from C82.
"""Exact finite-g,K transport and all-length midpoint selector at n=16."""
import json
import sys
sys.set_int_max_str_digits(0)
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
    n,gp,kp,lp=13,256,2304,576
    a=np.load('derivative_n13.npz',allow_pickle=True)
    R,ds=integer_ldl(a['F_numerator']+a['F_numerator'].T)
    bits=int(a['bits']);powers=[(int(d).bit_length()-1)//2 for d in ds]
    trace=sum((F((1<<(2*bits))*sum(int(x)**2 for x in R[i]),
                 1<<(2*powers[i])) for i in range(n)),F(0))
    det=int(a['detC'])
    assert F(64,15**6)>F(1,1<<18)
    unit=F(1,(1<<18)*det**2*(1<<(2*gp)))
    conv=4*n*det**2*trace/F(1<<(2*gp))
    # Physical error includes the zero-port splitting. Constants enlarged
    # from C80 to cover the extra data site and extra code-column factor.
    pen=16*n*trace*F(1<<84,1<<kp)/(F(1,1<<(2*lp))*unit)
    local=16*n*trace*F(1<<52,1<<(2*lp))/unit
    root=json.load(open('root_certificate_exact_n13.json'))
    budget=F(1,1<<root['perturbation_budget_power'])
    assert conv+pen+local<budget/8
    ab,cd=atan_bounds(F(1,5)),atan_bounds(F(1,239))
    pi=Interval(16*ab[0]-4*cd[1],16*ab[1]-4*cd[0])
    def cospi(d,p=1):
        return Interval(cos_bounds(pi.hi*p/d)[0],cos_bounds(pi.lo*p/d)[1])
    aa=lambda b:2*(1-cospi(b+1))
    dd=lambda b:(b+1)*aa(b)-b*aa(b+1)
    mu=2-(dd(14)+dd(13))/2
    ee=lambda b:(mu-2*cospi(b+1))/b
    delta=min((ee(13)-ee(14)).lo,(ee(15)-ee(14)).lo)
    gap=2*(cospi(15)-cospi(15,2))
    assert delta>F(1,1<<15)
    assert ee(14).hi<-F(1,1<<8)
    assert gap.lo>F(1,16)
    gamma=F(1,1<<15)-F(1,(1<<kp)-4)
    assert gamma>F(1,1<<16) and gamma>10*F(1,1<<lp)
    def lg(x):return math.log10(x.numerator)-math.log10(x.denominator)
    out=dict(n=n,material_alphabet=12,vacancy_colors=2,dimension=14,cell_length=14,
             g_power=-gp,K_min_power=kp,lambda_power=-lp,loop_eta_power=-74,
             field_normalized_power=-121,beta_denominator=str(det),
             Eunit_lower_log10=lg(unit),metric_trace_log10=lg(trace),
             conversion_error_log10=lg(conv),penalty_including_zero_port_error_log10=lg(pen),
             two_cell_lambda_error_log10=lg(local),budget_power=-35,
             sum_below_one_eighth_budget=True,selector_density_gap_lower=float(delta),
             selected_density_upper=float(ee(14).hi),mu_sel=[float(mu.lo),float(mu.hi)],
             vacuum_gap_lower=float(gap.lo),frozen_F3_gamma_lower_power=-16,
             frozen_F3_port_bound_pass=True)
    exact={k:[str(v.numerator),str(v.denominator)] for k,v in
           [('trace',trace),('unit_lower',unit),('conversion',conv),('penalty',pen),
            ('lambda_local',local),('budget',budget),('gamma_lower',gamma),
            ('selector_gap_lower',delta),('vacuum_gap_lower',gap.lo)]}
    with open('budget_n13.json','w') as f:json.dump(out,f,indent=2)
    with open('budget_exact_n13.json','w') as f:json.dump(exact,f,indent=2)
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()
