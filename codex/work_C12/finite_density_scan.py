import json,time
from qutrit_ed import solve
from qutrit_observables import entropy_profile,cfit,describe
results=[]
for delta in [0,0.15,0.3,0.45,0.6,0.9,1.2,1.8]:
    es,vs,states,r=solve(16,6,k=6,U=-14/9+delta,tol=3e-9)
    prof=entropy_profile(states,vs[:,0],16)
    r.update(delta=delta,levels=describe(states,es,vs,16),entropy=prof,fit=cfit(prof,16))
    results.append(r)
    print(json.dumps({key:r[key] for key in ['L','N','delta','energies','seconds','fit','levels']}),flush=True)
    with open('finite_density_scan_results.json','w') as f:json.dump(results,f,indent=2)
