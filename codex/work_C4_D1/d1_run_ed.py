"""Small-sector ED, largest sector dimension 184756."""
import json
from pathlib import Path
from d1 import spectrum
out=[]
for model,param,lengths in [('ik',.19,[4,6,8,10,12]),('tl',.4,[4,6,8,10,12,14,16,18,20]),('tl_af',0.,[4,8,12,16])]:
 for L in lengths:
  for charge in ([0,2] if model in ['ik','tl'] else [0]):
   result=spectrum(model,L,param,charge,count=18)
   out.append(result)
   Path('d1_spectra.json').write_text(json.dumps(out,indent=2)+'\n')
   print(model,L,charge,'dim',result['dimension'],'E0',result['eigenvalues'][0],
         'gap',result['eigenvalues'][1][0]-result['eigenvalues'][0][0],
         'residual',result['eigen_residual_max'],'seconds',round(result['seconds'],2),flush=True)
