"""ED of the uniform reflection completion in an EXACT endpoint sector.

Each unused endpoint Bell pair stays in Omega. This sector is invariant,
contains the full ground code, and omits no low state for small lambda.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh
from uniform_cell import *
from checks import logical_chain
from completion_checks import gaussian_residual


def endpoint(c,used):
    D,N,k,K=c['D'],c['N'],c['k'],c['K']
    cols=[]
    for flag in range(k+1):
        ge=np.zeros((D,N),complex)
        for l in range(2):
            for r in range(2):
                for z in range(N): ge[(2*l+r)*K+flag*4*N+(2*l+r)*N+z,z]=.5
        cols.append(ge)
        cols.extend([f@ge for f in c['FR' if used=='R' else 'FL']])
    T=np.concatenate(cols,axis=1)
    ops=[c['H'],c['field_phys']]+c['FR' if used=='R' else 'FL']
    leakage=max(la.norm(op@T-T@(T.conj().T@op@T)) for op in ops)
    assert leakage<1e-11
    return dict(D=T.shape[1],H=T.conj().T@c['H']@T,G=T.conj().T@c['G'],
                field_phys=T.conj().T@c['field_phys']@T,
                F=[T.conj().T@f@T for f in c['FR' if used=='R' else 'FL']],leakage=float(leakage))

def run(c,M,lambdas):
    left=endpoint(c,'R');right=endpoint(c,'L')
    cells=[left]+[c]*(M-2)+[right]
    dims=[x['D'] for x in cells];dim=math.prod(dims);low=c['N']**M
    def axis(a,x,j):
        y=np.tensordot(a,x.reshape(dims),axes=(1,j))
        return np.moveaxis(y,0,j).reshape(-1)
    GG=kron_all([x['G'] for x in cells]);target=logical_chain(c,M)
    rng=np.random.default_rng(919+c['n']*10+M);v0=rng.normal(size=dim)+1j*rng.normal(size=dim)
    ans=[]
    for lam in lambdas:
        start=time.time()
        local=[x['H']+lam**2*c['rho']**2*x['field_phys'] for x in cells]
        def mv(v):
            w=sum((axis(a,v,j) for j,a in enumerate(local)),np.zeros(dim,complex))
            for j in range(M-1):
                Rs=left['F'] if j==0 else c['FR']
                Ls=right['F'] if j+1==M-1 else c['FL']
                for r,l in zip(Rs,Ls):w+=lam*axis(r,axis(l,v,j+1),j)
            return w
        H=LinearOperator((dim,dim),matvec=mv,dtype=complex)
        ee,U=eigsh(H,k=low+1,which='SA',tol=1e-12,ncv=max(2*(low+1)+1,32),v0=v0)
        Q,_=la.qr(U,mode='economic');HQ=np.column_stack([H@Q[:,a] for a in range(Q.shape[1])])
        ev,V=la.eigh(Q.conj().T@HQ);U=Q@V
        er=max(la.norm(H@U[:,a]-ev[a]*U[:,a]) for a in range(low+1))
        a,s,b=la.svd(GG.conj().T@U[:,:low]);polar=a@b
        eff=(polar*ev[:low])@polar.conj().T;rem=eff-lam**2*target
        trrem=rem-np.trace(rem)/low*np.eye(low)
        rec=dict(n=c['n'],cells=M,dimension=dim,full_dimension=c['D']**M,local_dimensions=dims,
                 lambda_=lam,seconds=time.time()-start,endpoint_leakage=max(left['leakage'],right['leakage']),
                 eigenpair_residual=float(er),gap=float(ev[low]-ev[low-1]),
                 error=float(la.norm(rem,2)),traceless_error=float(la.norm(trrem,2)),
                 gaussian_residual=gaussian_residual(eff,c['n'],M),
                 low_eigenvalues=ev[:low].tolist(),heff_real=eff.real.tolist(),heff_imag=eff.imag.tolist())
        print(json.dumps({k:v for k,v in rec.items() if not k.startswith('heff')}),flush=True)
        ans.append(rec)
    return ans

if __name__=='__main__':
    out=[]
    for n,M in [(1,2),(1,3),(2,2)]:
        out+=run(reflection_cell(n),M,[.04,.02])
        with open('reflection_chain_results.json','w') as f:json.dump(out,f,indent=2)
