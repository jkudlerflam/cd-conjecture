import json
import fr_lib as F
from optimal_recovery import solve,gauge_scale
from scan_sdp import models
with open('targeted_results.jsonl','a') as f:
    for name,h in models():
        if name=='C106':continue
        d=F.bond_dim(h)
        cases=[(8,4)] if d==2 else [(6,3)]
        for L,l in cases:
            _,v,e=F.ground_states(h,L,k=8);M=gauge_scale(h)
            r=F.ti_marginal(v,d,L,l)
            z=solve(h,r,l,d,eps=2e-6,max_iters=2500)
            z.update(name=name,L=L,l=l,M=M,E=e,petz=L*L*(z['petz_energy']-e/L)/M,opt=L*L*(z['energy']-e/L)/M,dual_scaled=L*L*(z['corrected_dual']-e/L)/M)
            f.write(json.dumps(z)+'\n');f.flush();print(z,flush=True)
