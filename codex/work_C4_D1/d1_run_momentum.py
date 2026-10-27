import json
from pathlib import Path
from d1 import momentum_spectrum
out=[]
for model,param,lengths in [('ik',.19,[4,6,8,10,12]),('tl',.4,[4,8,12,16,20])]:
 for L in lengths:
  result=momentum_spectrum(model,L,param,count=10)
  out.extend(result)
  Path('d1_momenta.json').write_text(json.dumps(out,indent=2)+'\n')
  print(model,L,[(r['momentum_integer'],r['eigenvalues'][0],r['full_sector_residual_max']) for r in result],flush=True)
