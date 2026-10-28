import json,time
from momentum_ed import OrbitHamiltonian
results=[]
for delta,L,ns in [(.9,12,range(2,9)),(.9,18,range(4,9)),(.75,18,range(4,9))]:
 row=dict(delta=delta,L=L,sectors=[])
 for N in ns:
  ob=OrbitHamiltonian(L,N,U=-14/9+delta)
  best=None
  for k in range(L//2+1):
   es,vs,r=ob.solve(k,nev=1,tol=5e-10)
   if best is None or es[0]<best['E']:best=dict(N=N,E=float(es[0]),k=k)
  row['sectors'].append(best)
  print(delta,L,best,flush=True)
  with open('compressibility_check_results.json','w') as f:json.dump(results+[row],f,indent=2)
 results.append(row)
