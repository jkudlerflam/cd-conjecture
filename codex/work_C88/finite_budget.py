# C84 local source copied from C80. See PROOF.md for changes.
"""Rational transport to finite conversion, penalty and two-cell port strength."""
import json,math
from fractions import Fraction as Q
import numpy as np
from rational_bands import integer_ldl

def main():
 n=13;gp=272;kp=2560;lp=608
 a=np.load(f'derivative_n{n}.npz',allow_pickle=True)
 R,ds=integer_ldl(a['F_numerator']+a['F_numerator'].T)
 bits=int(a['bits']);powers=[(int(d).bit_length()-1)//2 for d in ds]
 trace=sum((Q((1<<(2*bits))*sum(int(x)**2 for x in R[i]),1<<(2*powers[i])) for i in range(n)),Q(0))
 det=int(a['detC'])
 assert Q(64,16**6)==Q(1,1<<18)
 unit=Q(1,(1<<18)*det**2*(1<<(2*gp)))
 conv=4*n*det**2*trace/Q(1<<(2*gp))
 # Includes the zero-port splitting in PHYSICAL units, not just H2.
 pen=16*n*trace*Q(1<<84,1<<kp)/(Q(1,1<<(2*lp))*unit)
 lam=16*n*trace*Q(1<<52,1<<(2*lp))/unit
 root=json.load(open(f'root_certificate_exact_n{n}.json'))
 budget=Q(1,1<<root['perturbation_budget_power'])
 assert conv+pen+lam<budget/8
 def lg(x):return math.log10(x.numerator)-math.log10(x.denominator)
 out=dict(n=n,alphabet=11,cell_length=15,g_power=-gp,K_power=kp,lambda_power=-lp,
  Eunit_lower_log10=lg(unit),metric_trace_log10=lg(trace),conversion_error_log10=lg(conv),
  penalty_including_zero_port_error_log10=lg(pen),two_cell_lambda_error_log10=lg(lam),
  budget_power=-root['perturbation_budget_power'],sum_below_one_eighth_budget=True)
 exact={k:[str(v.numerator),str(v.denominator)] for k,v in [('trace',trace),('unit_lower',unit),('conversion',conv),('penalty',pen),('lambda_local',lam),('budget',budget)]}
 json.dump(out,open('finite_budget.json','w'),indent=2)
 json.dump(exact,open('finite_budget_exact.json','w'),indent=2)
 print(json.dumps(out),flush=True)
if __name__=='__main__':main()
