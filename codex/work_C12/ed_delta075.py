import json,time
from momentum_ed import OrbitHamiltonian
from qutrit_observables import entropy_profile,cfit
out=[]
for L,N in [(21,7),(24,8)]:
 ob=OrbitHamiltonian(L,N,U=-14/9+.75)
 print('Built',L,len(ob.reps),flush=True)
 ks=range(L//2+1) if L==21 else [0,1,4,8,12]
 levels=[];best=None
 for k in ks:
  es,vs,r=ob.solve(k,nev=2,tol=3e-10);levels.append(r)
  print(L,k,es.tolist(),flush=True)
  if best is None or es[0]<best[0]:best=(float(es[0]),k,vs[:,0])
 states,v=ob.reconstruct(best[1],best[2]);prof=entropy_profile(states,v,L)
 row=dict(L=L,N=N,delta=.75,levels=levels,ground_k=best[1],entropy=prof,fit=cfit(prof,L),half=prof[-1])
 out.append(row);print('FIT',L,row['fit'],row['half'],flush=True)
 with open('ed_delta075_results.json','w') as f:json.dump(out,f,indent=2)
