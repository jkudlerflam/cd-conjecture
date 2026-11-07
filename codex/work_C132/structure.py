import json,numpy as np
from scan_sdp import models
from optimal_recovery import solve
import fr_lib as F
with open('structure_results.jsonl','w') as f:
    for name,h in models():
        if name not in ['XX','Heis','XXZ0.5','TFIM1.0','TFIM1.6']:continue
        _,v,e=F.ground_states(h,6);r=F.ti_marginal(v,2,6,3)
        z=solve(h,r,3,2,eps=1e-9,solver='CLARABEL',return_matrix=True)
        j=z['J'];eig=np.linalg.eigvalsh(j)
        a=dict(name=name,kraus_weights=eig[eig>1e-6].tolist(),distance_petz=z['distance_petz'],input_discarded_mass=z['input_discarded_mass'],output_discarded_mass=z['output_discarded_mass'])
        f.write(json.dumps(a)+'\n');f.flush();print(a,flush=True)
