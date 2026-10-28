"""Sparse fixed-charge ED. Explicit local bond conventions are in local_h()."""
import itertools, math, json, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

def local_h(U=-14/9,g=2,tp=10/9,V11=0,V12=20,V22=20,t21=0,mu=0):
    h=np.zeros((9,9))
    for a in range(3):
        for b in range(3):
            h[3*a+b,3*a+b]=U*((a==2)+(b==2))/2-mu*(a+b)/2
    for a,b,v in [(1,1,V11),(1,2,V12),(2,1,V12),(2,2,V22)]:h[3*a+b,3*a+b]+=v
    for a,b,c,d,v in [(1,0,0,1,-1),(2,0,0,2,-tp),(1,1,2,0,-g),(1,1,0,2,g),(2,1,1,2,-t21)]:
        h[3*a+b,3*c+d]=h[3*c+d,3*a+b]=v
    return h

def basis(L,N):
    if L>31:raise ValueError('uint64 basis requires L<=31')
    codes=[]
    for sites in itertools.combinations_with_replacement(range(L),N):
        if any(sites[j]==sites[j+2] for j in range(N-2)):continue
        codes.append(sum(1<<(2*j) for j in sites))
    return np.sort(np.asarray(codes,dtype=np.uint64))

def matrix(L,N,pbc=True,**pars):
    states=basis(L,N);dim=len(states)
    rows=[];cols=[];vals=[]
    h=local_h(**pars)
    diag=np.zeros(dim)
    for i in range(L if pbc else L-1):
        j=(i+1)%L
        a=((states>>np.uint64(2*i))&np.uint64(3)).astype(np.int8)
        b=((states>>np.uint64(2*j))&np.uint64(3)).astype(np.int8)
        local=3*a+b
        diag+=np.diag(h)[local]
        for old in range(9):
            for new in range(9):
                v=h[new,old]
                if old==new or v==0:continue
                ix=np.flatnonzero(local==old)
                if not len(ix):continue
                aa,bb=divmod(new,3);ao,bo=divmod(old,3)
                delta=(aa-ao)*(1<<(2*i))+(bb-bo)*(1<<(2*j))
                targets=states[ix]+np.uint64(delta) if delta>=0 else states[ix]-np.uint64(-delta)
                dest=np.searchsorted(states,targets)
                if not np.all(states[dest]==targets):raise ValueError('basis mismatch')
                rows.append(dest);cols.append(ix);vals.append(np.full(len(ix),v))
    rows.append(np.arange(dim));cols.append(np.arange(dim));vals.append(diag)
    H=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(dim,dim)).tocsr()
    return H,states

def solve(L,N,k=2,pbc=True,tol=1e-9,**pars):
    start=time.time();H,states=matrix(L,N,pbc,**pars)
    rng=np.random.default_rng(17+L+N)
    if H.shape[0]<=k+1:
        es,vs=np.linalg.eigh(H.toarray());es,vs=es[:k],vs[:,:k]
    else:
        es,vs=eigsh(H,k=k,which='SA',tol=tol,v0=rng.normal(size=H.shape[0]),maxiter=20000)
        ix=np.argsort(es);es,vs=es[ix],vs[:,ix]
    return es,vs,states,dict(L=L,N=N,dim=len(states),energies=es.tolist(),binding=float(es[0]+2*N),seconds=time.time()-start,residual=float(np.linalg.norm(H@vs[:,0]-es[0]*vs[:,0])))

if __name__=='__main__':
    results=[]
    for L,N in [(12,2),(18,2),(24,2),(12,3),(18,3),(24,3),(12,4),(18,4),(24,4),(12,5),(18,5),(24,5),(12,6),(18,6)]:
        e,v,s,r=solve(L,N,k=2)
        print(json.dumps(r),flush=True);results.append(r)
        with open('fewbody_results.json','w') as f:json.dump(results,f,indent=2)
