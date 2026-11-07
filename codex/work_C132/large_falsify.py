"""Cluster continuation. Run only on Typhon; ED and dense marginals become large."""
import json,numpy as np
import fr_lib as F
from falsify import gs
from optimal_recovery import solve,gauge_scale,support
from pathlib import Path
P=[2*F.sx,2*F.sy,2*F.sz];I=np.eye(2)
bonds=np.load('random_bonds.npy')
with open('large_falsify_results.jsonl','a') as out:
    for L in [14,16,18,20]:
        cases=[('XXZ'+str(delta),F.xxz_bond(delta)) for delta in [-.9,-.99,-.9999,-1.]]
        for c in [-2.,-.5,.5,2.,5.]:
            g=1+c/L;cases.append(('TFIMc'+str(c),(-np.kron(P[2],P[2])-g*(np.kron(P[0],I)+np.kron(I,P[0]))/2).real))
        cases.extend(('random'+str(i),bonds[i]) for i in [5,13])
        for name,h in cases:
            if name=='XXZ-1.0':
                z=dict(name=name,L=L,petz_opt= L/(L-1),method='proved formula, full ground-space mixture')
                out.write(json.dumps(z)+'\n');out.flush();continue
            _,v,e=F.ground_states(h,L,k=8);M=gauge_scale(h)
            vals=[L*L*(F.petz_defect(h,v,L,l)-e/L)/M for l in [L//2,L//2+1]]
            z=dict(name=name,L=L,petz=min(vals))
            out.write(json.dumps(z)+'\n');out.flush();print(z,flush=True)
