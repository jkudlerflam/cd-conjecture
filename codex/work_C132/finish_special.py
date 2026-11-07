import json
import fr_lib as F
from optimal_recovery import solve,gauge_scale,support
with open('special_results.jsonl','a') as f:
    for name,h,L,l in [('SU3',F.swap_bond(3),6,4),('C106',F.c106_bond(1,.03,2),6,3),('C106',F.c106_bond(1,.03,2),6,4)]:
        d=F.bond_dim(h);_,v,e=F.ground_states(h,L,k=8);M=gauge_scale(h);r=F.ti_marginal(v,d,L,l)
        b=len(support(F.ptrace_first(r,d,l))[0]);q=len(support(r)[0]);print(name,L,l,b,q,flush=True)
        if b*q>256:continue
        z=solve(h,r,l,d,eps=1e-6,max_iters=3000)
        z.update(name=name,L=L,l=l,M=M,E=e,petz=L*L*(z['petz_energy']-e/L)/M,opt=L*L*(z['energy']-e/L)/M,dual_scaled=L*L*(z['corrected_dual']-e/L)/M)
        f.write(json.dumps(z)+'\n');f.flush();print(z,flush=True)
