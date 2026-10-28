import json,time
from momentum_ed import OrbitHamiltonian
from qutrit_observables import entropy_profile,cfit
L,N=24,8
ob=OrbitHamiltonian(L,N,U=-14/9+.9)
print('Built L24',ob.full_dim,len(ob.reps),ob.build_seconds,flush=True)
rows=[];best=None
for k in [0,1,4,8,12]:
 es,vs,row=ob.solve(k,nev=2,tol=3e-10)
 print(json.dumps(row),flush=True);rows.append(row)
 if k!=1 and (best is None or es[0]<best[0]):best=(es[0],k,vs[:,0].copy())
 with open('ed_L24_results.json','w') as f:json.dump(dict(L=L,N=N,levels=rows),f,indent=2)
e,k,v=best
states,full=ob.reconstruct(k,v)
profile=entropy_profile(states,full,L)
out=dict(L=L,N=N,levels=rows,selected_ground_k=k,entropy=profile,fit=cfit(profile,L),norm=float(sum(abs(full)**2)))
print('Entropy fit',out['fit'],flush=True)
with open('ed_L24_results.json','w') as f:json.dump(out,f,indent=2)
