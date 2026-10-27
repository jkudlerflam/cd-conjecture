"""Finite-size diagnostics. Velocity-based c estimates are explicitly conditional."""
import json, math
from pathlib import Path
import numpy as np
r=json.loads(Path('d1_momenta.json').read_text())
results={}
for model in ['ik','tl']:
 rows=[]
 for a in r:
  if a['model']!=model or a['momentum_integer']!=0:continue
  L=a['L'];b=next(t for t in r if t['model']==model and t['L']==L and t['momentum_integer']==1)
  ev=np.array(a['eigenvalues'])[:,0]; e0=ev[0]
  gaps=L*(ev-e0)
  nonzero=gaps[gaps>1e-7]
  vp=L*(b['eigenvalues'][0][0]-e0)/(2*math.pi)
  rows.append({'L':L,'E0':e0,'L_neutral_gap':float(nonzero[0]),
    'neutral_gap_ratios':list(nonzero[:4]/nonzero[0]),'velocity_descendant_assumption':vp})
 rows.sort(key=lambda a:a['L'])
 v_recalled=math.pi/(math.pi-6*.19) if model=='ik' else 2*math.pi*math.sin(.8)/(math.pi-.8)
 for i,a in enumerate(rows):
  if i==0:continue
  b=rows[i-1]
  A=-(a['E0']/a['L']-b['E0']/b['L'])/(1/a['L']**2-1/b['L']**2)
  a['pair_cv_eff']=6*A/math.pi
  a['pair_c_eff_using_descendant_velocity']=6*A/(math.pi*a['velocity_descendant_assumption'])
  a['pair_c_eff_using_unverified_recalled_velocity']=6*A/(math.pi*v_recalled)
 # Test linearity of 1/sqrt(Lgap) against log L; it cannot establish asymptotics.
 selected=rows[1:]
 x=np.log([a['L'] for a in selected]);y=1/np.sqrt([a['L_neutral_gap'] for a in selected])
 coef=np.polyfit(x,y,1);pred=np.polyval(coef,x)
 logfit={'slope':coef[0],'intercept':coef[1],'max_y_residual':max(abs(y-pred)),
   'warning':'Short-size consistency test only; power-law alternatives are not excluded.'}
 results[model]={'velocity_formula_status':'recalled/guessed, not derived or externally verified',
                  'velocity_candidate':v_recalled,'rows':rows,'log_gap_test':logfit}
Path('d1_diagnostics.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
