"""Check the remaining conserved-charge sectors at the largest ED sizes."""
import json,time
import numpy as np
from scipy.sparse.linalg import eigs
from pathlib import Path
from d1 import sector_hamiltonian,ik_density,tl_three_density
out=[]
for model,d,L,g,local in [('ik',3,12,.19,ik_density(.19)),('tl',2,20,.4,tl_three_density(.4))]:
 for charge in range(0,(d-1)*L+1,2):
  start=time.time();h,_,_=sector_hamiltonian(local,d,L,charge)
  if h.shape[0]<=20: vals,vecs=np.linalg.eig(h.toarray())
  else:
   rng=np.random.default_rng(1907)
   v0=rng.normal(size=h.shape[0])+1j*rng.normal(size=h.shape[0])
   vals,vecs=eigs(h,k=2,which='SR',tol=1e-10,maxiter=10000,v0=v0)
  a=int(np.argmin(vals.real));v=vecs[:,a];z=vals[a]
  r={'model':model,'L':L,'charge_2Sz':charge,'dimension':h.shape[0],
     'lowest_real_part_eigenvalue':[float(z.real),float(z.imag)],
     'residual':float(np.linalg.norm(h@v-z*v)/np.linalg.norm(v)), 'seconds':time.time()-start}
  out.append(r);Path('d1_ground_audit.json').write_text(json.dumps(out,indent=2)+'\n')
  print(r,flush=True)
