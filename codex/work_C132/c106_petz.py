import json
import fr_lib as F
h=F.c106_bond(1,.03,2);L=6;_,v,e=F.ground_states(h,L,k=8);M=F.M_of_bond(h)
vals=[]
for l in [2,3,4]:
    x=L*L*(F.petz_defect(h,v,L,l)-e/L)/M
    vals.append(dict(L=L,l=l,petz=x));print(vals[-1],flush=True)
from pathlib import Path
Path('c106_petz_results.json').write_text(json.dumps(vals))
