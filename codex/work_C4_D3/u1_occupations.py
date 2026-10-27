import json,numpy as np
from d3 import Sector,measure
from analyze import cf
rows=[json.loads(x) for x in open('refine_u1.jsonl')]
for r in rows:
 if r['L']!=14 or r['coeff'][1]!=1:continue
 sec=Sector('spin1',14,r['charge'],cache_blocks=False)
 e,v,_=sec.solve(r['coeff'],r['ground_k'],3);psi=sec.expand(r['ground_k'],v[:,0]);sec._blocks.clear()
 code=sec.basis.copy();K=np.zeros(len(code),dtype=int)
 for i in range(14):K+=(code%3)!=1;code//=3
 prob=abs(psi)**2;mean=float(prob@K);var=float(prob@(K-mean)**2)
 out={'label':r['label'],'L':14,'M':r['charge'][0],'K_nonzero_mean':mean,'K_variance':var,'E0':float(e[0])}
 with open('u1_occupations.jsonl','a') as f:f.write(json.dumps(out)+'\n')
 print(out,flush=True)
