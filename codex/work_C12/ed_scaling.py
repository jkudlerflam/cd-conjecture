import json,sys
from qutrit_ed import solve
from qutrit_observables import entropy_profile,cfit,momentum_eigenstate,occupations
import numpy as np
results=[]
delta=float(sys.argv[1]) if len(sys.argv)>1 else .9
for L,N in [(9,3),(12,4),(15,5),(18,6),(21,7)]:
    es,vs,states,r=solve(L,N,k=2,U=-14/9+delta,tol=1e-9)
    v,mom=momentum_eigenstate(states,vs[:,0],L)
    prof=entropy_profile(states,v,L)
    r.update(delta=delta,momentum=mom,entropy=prof,fit=cfit(prof,L),doublons=float(np.sum(abs(v)**2*np.sum(occupations(states,L)==2,axis=1))))
    print(json.dumps({k:r[k] for k in ['L','N','energies','seconds','fit','momentum','doublons']}),flush=True)
    results.append(r)
    with open('ed_scaling_delta_'+str(delta)+'.json','w') as f:json.dump(results,f,indent=2)
