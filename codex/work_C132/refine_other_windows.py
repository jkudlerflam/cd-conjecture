import json
import fr_lib as F
from optimal_recovery import solve,gauge_scale
from scan_sdp import models
with open('refined_other_windows.jsonl','a') as f:
    for name,h in models():
        if name not in ['XX','Heis','XXZ0.5','XXZ2','TFIM1.0','TFIM1.6']:continue
        L=6;d=2;_,v,e=F.ground_states(h,L);M=gauge_scale(h)
        for l in [2,4,5]:
            z=solve(h,F.ti_marginal(v,d,L,l),l,d,eps=1e-9,solver='CLARABEL')
            z.update(name=name,L=L,l=l,M=M,E=e,petz=L*L*(z['petz_energy']-e/L)/M,opt=L*L*(z['energy']-e/L)/M,dual_scaled=L*L*(z['corrected_dual']-e/L)/M)
            f.write(json.dumps(z)+'\n');f.flush();print(z,flush=True)
