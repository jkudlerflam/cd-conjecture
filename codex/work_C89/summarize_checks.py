"""Read only this folder's checks and write compact schedule arithmetic."""
from pathlib import Path
import json
import math
from check_defects import selector,cell_energy

HERE=Path(__file__).resolve().parent
mu=selector(15)
e=cell_energy(15,mu)/15
a=0.5-2/math.pi**2
pin=2*math.cos(math.pi/16)
wrong={str(l):dict(excess=cell_energy(l,mu)-l*e,
                  leading_cost_over_J=a*(cell_energy(l,mu)-l*e),
                  density_margin=cell_energy(l,mu)/l-e)
       for l in [1,2,14,16,30,60,1000]}
d0,d1=1,1004
for j in range(1,15): d0,d1=d1,1004*d1-d0
alpha=math.sqrt(2/16)*math.sin(math.pi/16)
log2unit=4*math.log2(alpha)-2*math.log2(d1)-544
ed=json.loads((HERE/'ed_results.json').read_text())
residuals=[c[k]['eigenpair_residual'] for c in ed['clock'] for k in ['reference','defect']]
residuals += [s[k] for s in ed['scalar'] for k in ['eigenpair_residual','reference_residual']]
out=dict(tag='numerical_arithmetic_and_finite_ED_summary',mu_sel=mu,e15=e,
         alpha_half=a,pin_bare=pin,pin_leading_over_J=a*pin,wrong=wrong,
         weighted_margin_conservative=2**-17,
         alpha_weighted_margin=a*2**-17,
         log2_Eunit=log2unit,log2_lambda2=-1216,
         log2_lambda2_Eunit=log2unit-1216,
         log2_field=log2unit-1216-120,
         one_marker_asymptotic_length_threshold=4/(a*abs(e)),
         clock_comparisons=len(ed['clock']),scalar_comparisons=len(ed['scalar']),
         largest_eigenpair_residual=max(residuals),
         all_clock_costs_positive=all(c['cost']>0 for c in ed['clock']),
         all_scalar_costs_positive=all(c['cost']>0 for c in ed['scalar']),
         max_clock_dimension=max(c[k]['dimension'] for c in ed['clock'] for k in ['reference','defect']),
         max_scalar_dimension=max(c['dimension'] for c in ed['scalar']))
(HERE/'checks_summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
