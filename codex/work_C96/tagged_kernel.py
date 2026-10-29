"""Finite tagged-contact checks. Run only from work_C96."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from itertools import combinations
import math
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, cg
HERE=Path(__file__).resolve().parent
assert Path.cwd().resolve()==HERE

def compositions(v,n):
    for sep in combinations(range(v+n-1),n-1):
        a=(-1,)+sep+(v+n-1,)
        yield tuple(a[i+1]-a[i]-1 for i in range(n))

def bath(n,v):
    gs=np.array(list(compositions(v,n)),dtype=int)
    idx={tuple(g):i for i,g in enumerate(gs)}
    rows=[]; cols=[]
    for i,g in enumerate(gs):
        for j in range(n):
            if g[j]>0:
                for dest in [(j+1)%n,(j-1)%n]:
                    out=g.copy();out[j]-=1;out[dest]+=1
                    rows.append(idx[tuple(out)]);cols.append(i)
    h=coo_matrix((-np.ones(len(rows)),(rows,cols)),shape=(len(gs),len(gs))).tocsr()
    xs=np.zeros_like(gs); xs[:,1:]=np.cumsum(gs[:,:-1]+1,axis=1)
    l=n+v
    vals=np.ones(len(gs))
    for i in range(n):
        for j in range(i+1,n):
            vals*=2*np.sin(math.pi*(xs[:,j]-xs[:,i])/l)
    vals/=np.linalg.norm(vals)
    en=-2*math.sin(math.pi*n/l)/math.sin(math.pi/l)
    aa=(gs==0).astype(float)
    a=float(vals@(aa[:,0]*vals))
    return h,xs,vals,en,aa,a

def momenta(n,v):
    l=n+v
    parity=0 if n%2 else .5
    for m in range(n):
        for inds in combinations(range(l),n):
            if (sum(inds)+round(n*parity)+m)%l:
                continue
            p=2*math.pi*(np.array(inds)+parity+m/n)/l
            yield m,inds,p

def contact_ff(p,k,n,l):
    x=np.arange(l)
    up=np.exp(1j*np.outer(x,p));uk=np.exp(1j*np.outer(x,k))
    rest=up[2:].conj().T@uk[2:]
    z=np.outer(up[0].conj(),uk[0]);w=np.outer(up[1].conj(),uk[1])
    coeff=(np.linalg.det(rest+z+w)-np.linalg.det(rest+z)
           -np.linalg.det(rest+w)+np.linalg.det(rest))
    return coeff/(n*l**(n-1))

def spectral_row(n,v,verify_vectors=False):
    l=n+v
    k=2*math.pi*(np.arange(n)-(n-1)/2)/l
    e0=-2*np.cos(k).sum()
    raw=[]; sums=np.zeros(n); weights=np.zeros(n)
    xi=0 if n%2 else .5
    grid=2*math.pi*(np.arange(l)+xi)/l
    filled=set(np.argsort(-2*np.cos(grid))[:n])
    uniform_error=0.
    if verify_vectors:
        h,xs,psi,_,aa,a=bath(n,v)
        vg=np.linalg.det(np.exp(1j*xs[:,:,None]*k[None,None,:]))/math.sqrt(n*l**(n-1))
        vmax=0.;fmax=0.;vs=[]
    for m,inds,p in momenta(n,v):
        e=-2*np.cos(p).sum()
        f=contact_ff(p,k,n,l)
        if m==0 and e-e0>1e-10:
            particles=sorted(set(inds)-filled)
            holes=sorted(filled-set(inds))
            expected=0.
            if len(particles)==2:
                p1,p2=grid[particles];h1=grid[holes[0]]
                expected=abs(2/(n*l)*(math.cos(p1-h1)-math.cos(p2-h1)))
            uniform_error=max(uniform_error,abs(abs(f)-expected))
        if verify_vectors:
            vp=np.linalg.det(np.exp(1j*xs[:,:,None]*p[None,None,:]))/math.sqrt(n*l**(n-1))
            vmax=max(vmax,np.linalg.norm(h@vp-e*vp))
            fmax=max(fmax,abs(np.vdot(vp,aa[:,0]*vg)-f))
            vs.append(vp)
        weights[m]+=abs(f)**2
        if e-e0>1e-10:
            sums[m]+=abs(f)**2/(e-e0)
        raw.append((m,e,f))
    row=np.fft.fft(sums).real
    mean=n/l-math.sin(math.pi*n/l)**2/(n*l*math.sin(math.pi/l)**2)
    assert abs(weights.sum()-mean)<2e-11
    assert uniform_error<2e-11
    ans=dict(N=n,V=v,dimension=len(raw),row=row.tolist(),eigenvalues=(n*sums).tolist(),spectral_weight=float(weights.sum()),projector_error=float(abs(weights.sum()-mean)),uniform_form_factor_error=float(uniform_error))
    if verify_vectors:
        mat=np.array(vs).T
        ans.update(eigen_residual=float(vmax),ff_residual=float(fmax),orthogonality_error=float(np.max(abs(mat.conj().T@mat-np.eye(len(vs))))))
    return ans

def resolvent_row(n,v):
    h,xs,psi,en,aa,a=bath(n,v)
    def mv(z):return h@z-en*z+psi*np.dot(psi,z)
    op=LinearOperator(h.shape,matvec=mv,dtype=float)
    src=(aa[:,0]-a)*psi
    z,info=cg(op,src,rtol=1e-12,atol=1e-14)
    assert info==0
    row=(aa-a).T@(psi*z)
    return row,float(np.linalg.norm(mv(z)-src))

if __name__=='__main__':
    import json
    results=[]
    for n,v in [(3,2),(3,3),(4,4),(5,5),(6,6),(7,7),(8,8)]:
        ans=spectral_row(n,v,verify_vectors=n<=4)
        row,err=resolvent_row(n,v)
        ans['resolvent_error']=err
        ans['row_error']=float(np.max(abs(row-ans['row'])))
        rho=n/(n+v)
        aq=-rho*rho+rho/math.pi*math.sin(2*math.pi*rho)-math.sin(math.pi*rho)**2/math.pi**2
        epp=2*math.pi*rho**3*math.sin(math.pi*rho)
        c=aq*aq/(2*epp)
        reg=row+c/n
        dist=np.minimum(np.arange(n),n-np.arange(n))
        ans.update(c_density=c,K=float((1+dist)@abs(reg)),regular_row=reg.tolist())
        assert ans['row_error']<2e-10,ans
        for field in ['eigen_residual','ff_residual','orthogonality_error']:
            if field in ans:assert ans[field]<2e-10,(field,ans)
        results.append(ans)
        print(n,v,ans,flush=True)
    (HERE/'tagged_checks.json').write_text(json.dumps(results,indent=2)+'\n')
