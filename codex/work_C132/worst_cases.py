import json,numpy as np
import fr_lib as F
from optimal_recovery import solve,gauge_scale
bonds=np.load('random_bonds.npy')
with open('worst_results.jsonl','a') as f:
    for n in [5,13]:
        h=bonds[n];M=gauge_scale(h)
        for L in [4,6,8,10,12]:
            w,v,e=F.ground_states(h,L,k=4)
            ls=range(2,L) if L<=8 else [L//2,L//2+1]
            vals=[]
            for l in ls:
                val=F.petz_defect(h,v,L,l);z=dict(name='random'+str(n),L=L,l=l,petz=L*L*(val-e/L)/M)
                if L<=6:
                    r=F.ti_marginal(v,2,L,l);a=solve(h,r,l,2,eps=1e-9,solver='CLARABEL')
                    z.update(opt=L*L*(a['energy']-e/L)/M,dual_scaled=L*L*(a['corrected_dual']-e/L)/M,tp_res=a['tp_res'],rec_res=a['rec_res'],min_eig=a['min_eig'],gap=a['gap'],status=a['status'])
                f.write(json.dumps(z)+'\n');f.flush();print(z,flush=True)
