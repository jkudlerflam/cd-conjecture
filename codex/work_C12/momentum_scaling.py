import json,time
from momentum_ed import OrbitHamiltonian
results=[]
for L,N in [(9,3),(12,4),(15,5),(18,6),(21,7)]:
 ob=OrbitHamiltonian(L,N,U=-14/9+.9)
 print('Built',L,N,ob.full_dim,len(ob.reps),ob.build_seconds,flush=True)
 levels=[]
 for k in range(L//2+1):
  es,vs,r=ob.solve(k,nev=3,tol=1e-10)
  levels.append(r)
  print(L,k,es.tolist(),r['seconds'],flush=True)
 rr=dict(L=L,N=N,levels=levels)
 results.append(rr)
 with open('momentum_scaling_results.json','w') as f:json.dump(results,f,indent=2)
