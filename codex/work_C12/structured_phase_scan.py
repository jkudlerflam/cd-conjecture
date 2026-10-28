import json,time,numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from momentum_ed import OrbitHamiltonian
from qutrit_observables import entropy_profile,cfit
results=[]
for L,N in [(12,4),(15,5),(18,6)]:
 ob=OrbitHamiltonian(L,N,U=-14/9+.9)
 doublons=sum((((ob.reps>>(2*j))&3)==2).astype(float) for j in range(L))
 sectors={k:ob.sector(k) for k in range(L//2+1)}
 warm={}
 for delta in [.3,.45,.6,.75,.9,1.05,1.2]:
  best=None;levels=[];start=time.time()
  for k,(H,valid) in sectors.items():
   hh=H+diags((delta-.9)*doublons[valid])
   rng=np.random.default_rng(100+k)
   v0=rng.normal(size=hh.shape[0])+1j*rng.normal(size=hh.shape[0])
   if k in warm:v0=warm[k]+1e-3*v0/np.linalg.norm(v0)
   es,vs=eigsh(hh,k=1,which='SA',tol=2e-10,v0=v0,ncv=28)
   warm[k]=vs[:,0]
   levels.append(dict(k=k,E=float(es[0])))
   if best is None or es[0]<best[0]:best=(float(es[0]),k,vs[:,0])
  states,v=ob.reconstruct(best[1],best[2]);prof=entropy_profile(states,v,L)
  row=dict(L=L,N=N,delta=delta,E=best[0],momentum=best[1],levels=levels,entropy=prof,fit=cfit(prof,L),half=prof[-1],seconds=time.time()-start)
  print(json.dumps({k:row[k] for k in ['L','N','delta','E','momentum','fit','half','seconds']}),flush=True)
  results.append(row)
  with open('structured_phase_scan_results.json','w') as f:json.dump(results,f,indent=2)
