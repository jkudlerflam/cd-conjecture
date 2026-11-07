import json,sys
from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse.linalg as spla
import fr_lib as F
from optimal_recovery import gauge_scale,solve

def gs(h,L):
    if np.max(abs(h-F.xxz_bond(-1)))<1e-12 and L%2==0:
        inds=np.arange(2**L); nums=np.array([int(i).bit_count() for i in inds])
        phase=np.ones(2**L)
        for site in range(0,L,2):phase*=1-2*((inds>>(L-1-site))&1)
        v=np.zeros((2**L,L+1))
        for n in range(L+1):
            ids=nums==n;v[ids,n]=phase[ids]/np.sqrt(np.sum(ids))
        return v,-L/4
    H=F.ring_hamiltonian(h,L); N=H.shape[0]
    if N<=256:w,v=la.eigh(H.toarray())
    else:
        k=min(24,N-2); w,v=spla.eigsh(H,k=k,which='SA',tol=1e-11)
        z=np.argsort(w);w=w[z];v=v[:,z]
    ids=w-w[0]<1e-9
    if ids.all():raise RuntimeError('ground space not exhausted')
    # ARPACK degeneracies are orthonormalized before forming the equal mixture.
    vv=la.qr(v[:,ids],mode='economic')[0]
    return vv,float(w[0])

def run(name,h,L,f,do_sdp=False):
    if (name,L) in DONE:return DONE[name,L]
    d=2; M=gauge_scale(h); v,e=gs(h,L)
    ls=sorted(set([max(2,L//2),min(L-1,L//2+1)]))
    vals=[]
    for l in ls:
        pet=F.petz_defect(h,v,L,l)
        vals.append((L*L*(pet-e/L)/M,l))
    best,l=min(vals)
    z=dict(name=name,L=L,M=M,E=e,deg=v.shape[1],petz=best,l=l)
    if do_sdp:
        r=F.ti_marginal(v,d,L,l); a=solve(h,r,l,d,eps=1e-8,solver="CLARABEL")
        z.update(opt=L*L*(a['energy']-e/L)/M,status=a['status'],tp_res=a['tp_res'],rec_res=a['rec_res'],min_eig=a['min_eig'])
    f.write(json.dumps(z)+'\n');f.flush(); print(z,flush=True)
    return z

DONE={}
if __name__=='__main__':
    if Path('falsify_results.jsonl').exists():
        DONE={(z['name'],z['L']):z for z in map(json.loads,Path('falsify_results.jsonl').read_text().splitlines())}
    rng=np.random.default_rng(132)
    paulis=[2*F.sx,2*F.sy,2*F.sz]; I=np.eye(2)
    with open('falsify_results.jsonl','a') as f:
        for delta in [-.9,-.99,-.9999,-1.]:
            for L in [4,6,8,10,12]:run('XXZ'+str(delta),F.xxz_bond(delta),L,f,L<=6)
        for c in [-2.,-.5,.5,2.,5.]:
            for L in [4,6,8,10,12]:
                g=1+c/L; h=(-np.kron(paulis[2],paulis[2])-g*(np.kron(paulis[0],I)+np.kron(I,paulis[0]))/2).real
                run('TFIMc'+str(c),h,L,f,L<=6)
        bonds=[]
        for k in range(16):
            coeff=rng.normal(size=12)
            h=sum(coeff[3*a+b]*np.kron(paulis[a],paulis[b]) for a in range(3) for b in range(3))
            h+=sum(coeff[9+a]*(np.kron(paulis[a],I)+np.kron(I,paulis[a]))/2 for a in range(3))
            h/=gauge_scale(h);bonds.append(h)
            for L in [4,6,8,10]:run('random'+str(k),h,L,f,L<=6)
        np.save('random_bonds.npy',np.array(bonds))
