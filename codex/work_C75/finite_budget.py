"""Exact rational transport of the n=15 band certificate to finite K,g."""
import json,math
import numpy as np
from fractions import Fraction as Q
from rational_bands import integer_ldl

def main():
    n=15;gpow=256;Kpow=1152;lampow=640
    a=np.load('derivative_n15.npz',allow_pickle=True);R,ds=integer_ldl(a['F_numerator']+a['F_numerator'].T)
    bits=int(a['bits']);powers=[(int(d).bit_length()-1)//2 for d in ds]
    trace=sum((Q((1<<(2*bits))*sum(int(x)**2 for x in R[i]),1<<(2*powers[i])) for i in range(n)),Q(0))
    det=int(a['detC']);unit_lower=Q(1,(1<<19)*det*det*(1<<(2*gpow)))
    # alpha^4 >=64/17^6 >2^-19, using sin(pi/17)>=2/17.
    assert Q(64,17**6)>Q(1,1<<19)
    conversion=4*n*det*det*trace/Q(1<<(2*gpow))
    # PROOF: ||H2(K,g)-H2(infinity,g)|| <= 2^64/K.
    penalty=16*n*trace*Q(1<<64,1<<Kpow)/unit_lower
    budget=Q(1,1<<32)
    assert conversion+penalty<budget/8
    # A separate two-cell fourth-order comparison. This is NOT a
    # volume-uniform infrared estimate.
    lambda_local=trace*Q(1<<50,1<<(2*lampow))/unit_lower
    def lg(x):return math.log10(x.numerator)-math.log10(x.denominator)
    rec=dict(n=n,alphabet=15,cell_length=16,g_power=-gpow,K_power=Kpow,lambda_power=-lampow,
             Eunit_lower_log10=lg(unit_lower),trace_metric_log10=lg(trace),
             conversion_whitened_error_log10=lg(conversion),penalty_whitened_error_log10=lg(penalty),
             combined_below_one_eighth_budget=True,budget_power=-32,local_lambda_ratio_log10=lg(lambda_local),
             delta2_lower=0.014554025253406,delta4_lower=1.1444044416118114e-7)
    print(json.dumps(rec),flush=True);json.dump(rec,open('finite_budget.json','w'),indent=2)
    exact={k:[str(v.numerator),str(v.denominator)] for k,v in [('trace',trace),('conversion',conversion),('penalty',penalty),('budget',budget)]}
    json.dump(exact,open('finite_budget_exact.json','w'),indent=2)
if __name__=='__main__':main()
