import json,sys,time
import numpy as np
import fr_lib as F
from optimal_recovery import solve,gauge_scale
from pathlib import Path

def models():
    out=[('XX',F.xxz_bond(0)),('Heis',F.heis_bond()),('XXZ0.5',F.xxz_bond(.5)),('XXZ2',F.xxz_bond(2))]
    X=2*F.sx; Z=2*F.sz; I=np.eye(2)
    for g in [1.,1.6]:out.append(('TFIM'+str(g),(-np.kron(Z,Z)-g*(np.kron(X,I)+np.kron(I,X))/2).real))
    out.append(('SU3',F.swap_bond(3)))
    out.append(('C106',F.c106_bond(1,.03,2)))
    return out

if __name__=='__main__':
    cap=int(sys.argv[1]) if len(sys.argv)>1 else 128
    done={(z['name'],z['L'],z['l']) for z in map(json.loads,Path('sdp_results.jsonl').read_text().splitlines()) if z['status']=='optimal' and z['min_eig']>-1e-7} if Path('sdp_results.jsonl').exists() else set()
    with open('sdp_results.jsonl','a') as f:
        for name,h in models():
            d=F.bond_dim(h); M=gauge_scale(h)
            Ls=([4,6,8,10,12] if '--long' in sys.argv else [4,6,8]) if d==2 else ([3,6,9] if '--long' in sys.argv and d==3 else ([3,6] if d==3 else [4,6]))
            for L in Ls:
                _,v,e=F.ground_states(h,L,k=8)
                for l in range(2,L):
                    if (name,L,l) in done:continue
                    if d**l>2048:continue
                    r=F.ti_marginal(v,d,L,l)
                    from optimal_recovery import support
                    b=len(support(F.ptrace_first(r,d,l))[0]); q=len(support(r)[0])
                    if b*q>cap:continue
                    z=solve(h,r,l,d)
                    z.update(name=name,L=L,l=l,M=M,E=e,deg=v.shape[1],petz=L*L*(z['petz_energy']-e/L)/M,opt=L*L*(z['energy']-e/L)/M)
                    f.write(json.dumps(z)+'\n');f.flush()
                    print(name,L,l,'Petz %.6f opt %.6f'%(z['petz'],z['opt']),z['status'],flush=True)
