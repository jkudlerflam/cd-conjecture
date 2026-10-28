import itertools,numpy as np,json,time,signal
from momentum_ed import OrbitHamiltonian
signal.alarm(1200)
results=[]
for L in [36,48,60]:
 start=time.time();codes=[]
 for sites in itertools.combinations_with_replacement(range(L),4):
  if sites[0]==sites[2] or sites[1]==sites[3]:continue
  codes.append(sum(1<<(2*j) for j in sites))
 states=np.array(sorted(codes),dtype=object)
 ob=OrbitHamiltonian(L,4,custom_basis=states)
 row=dict(L=L,build_seconds=time.time()-start,full_dim=len(states),orbit_dim=len(ob.reps))
 for k,key in [(0,'opposite'),(L//3,'same')]:
  es,vs,r=ob.solve(k,nev=1,tol=1e-11);row[key]=float(es[0])
 row['difference']=row['opposite']-row['same'];row['scaled_difference']=L**3*row['difference']
 print(json.dumps(row),flush=True);results.append(row)
 with open('long_scattering_results.json','w') as f:json.dump(results,f,indent=2)
