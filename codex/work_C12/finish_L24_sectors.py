import json,signal
from momentum_ed import OrbitHamiltonian
from qutrit_observables import entropy_profile,cfit
signal.alarm(900)
ob=OrbitHamiltonian(24,8,U=-14/9+.9)
base=json.load(open('ed_L24_results.json'))
try:levels=json.load(open('ed_L24_remaining_sectors.json'))
except FileNotFoundError:levels=[]
best=None
for k in [2,3,5,6,7,9,10,11]:
 if k in {r['k'] for r in levels}:continue
 es,vs,r=ob.solve(k,nev=1,tol=5e-10);levels.append(r)
 print(json.dumps(r),flush=True)
 if es[0]<min(x['energies'][0] for x in base['levels']):
  if best is None or es[0]<best[0]:best=(es[0],k,vs[:,0])
 with open('ed_L24_remaining_sectors.json','w') as f:json.dump(levels,f,indent=2)
if best:
 states,v=ob.reconstruct(best[1],best[2]);prof=entropy_profile(states,v,24)
 with open('ed_L24_new_ground.json','w') as f:json.dump(dict(k=best[1],E=float(best[0]),entropy=prof,fit=cfit(prof,24)),f,indent=2)
 print('New lower sector',best[1],best[0],flush=True)
