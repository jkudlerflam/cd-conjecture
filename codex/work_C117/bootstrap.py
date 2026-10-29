"""Small real-bond local SDPs, scipy dual search and exact integer certificates.
All optimization output is a lower bound, not a claimed SDP optimum.
Run --small locally; larger experiments belong in SLURM.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize


def sbasis(n):
    out=[]
    for i in range(n):
        for j in range(i,n):
            a=np.zeros((n,n)); a[i,j]=a[j,i]=1
            out.append(a)
    return np.array(out)


def bonds():
    x=np.array([[0,1],[1,0.]])
    z=np.diag([1.,-1.]); y=np.array([[0,-1j],[1j,0]])
    xy=(np.kron(x,x)+np.kron(y,y)).real/4
    out=[('Heisenberg',2,xy+np.kron(z,z)/4),('XXZ_half',2,xy+np.kron(z,z)/8)]
    p=np.zeros((9,9))
    for i in range(3):
        for j in range(3): p[3*i+j,3*j+i]=1
    out.append(('SU3_swap',3,p))
    rng=np.random.default_rng(117)
    for d in (2,3):
        a=rng.integers(-8,9,size=(d*d,d*d))
        out.append(('random_real_d'+str(d),d,(a+a.T)/16))
    return out


def setup(h,d,k,ground):
    n=d**(k-1); D=d**k
    C=np.kron(h,np.eye(d**(k-2)))
    bas=sbasis(n)
    DX=np.array([np.kron(b,np.eye(d))-np.kron(np.eye(d),b) for b in bas])
    # Delete one diagonal coefficient to remove the scalar gauge.
    DX=DX[:-1]; bas=bas[:-1]
    S=np.zeros((0,D,D)); G=np.zeros((0,0,D,D))
    if ground:
        assert k>=3
        H3=np.kron(h,np.eye(d))+np.kron(np.eye(d),h)
        ops=[]
        for a in range(d):
            for b in range(d):
                o=np.zeros((d,d)); o[a,b]=1
                ops.append(np.kron(np.eye(d),np.kron(o,np.eye(d))))
        ss=[]
        for a in range(d):
            for b in range(a+1,d):
                o=ops[a*d+b]-ops[b*d+a]
                ss.append(np.kron(H3@o-o@H3,np.eye(d**(k-3))))
        S=np.array(ss)
        G=np.array([[np.kron((oa.T@(H3@ob-ob@H3)+(oa.T@H3-H3@oa.T)@ob)/2,
                                 np.eye(d**(k-3))) for ob in ops] for oa in ops])
    return C,DX,S,G


def unpack(z,nx,ns,m):
    return z[:nx],z[nx:nx+ns],z[nx+ns:].reshape(m,m)


def matrix(C,DX,S,G,x,y,B):
    K=C+np.einsum('a,aij->ij',x,DX)
    if len(y): K=K+np.einsum('a,aij->ij',y,S)
    if B.size: K=K-np.einsum('ab,abij->ij',B@B.T,G)
    return (K+K.T)/2


def certify(C,DX,S,G,x,y,B):
    # All problem coefficients are dyadic. These scale checks are exact here.
    def integers(a,p):
        v=np.rint(a*(2**p))
        assert np.array_equal(v/(2**p),a)
        return np.array([int(a) for a in v.ravel()],dtype=object).reshape(v.shape)
    xi=np.rint(x*2**16).astype(object); yi=np.rint(y*2**16).astype(object)
    bi=np.rint(B*2**12).astype(object)
    xi=np.array([int(a) for a in xi],dtype=object)
    yi=np.array([int(a) for a in yi],dtype=object)
    bi=np.array([[int(a) for a in row] for row in bi],dtype=object).reshape(B.shape)
    ci=integers(C,32)
    # Convert integral floats to Python integers before exact operations.
    ci=np.array([[int(a) for a in row] for row in ci],dtype=object)
    for a,t in zip(xi,DX): ci+=int(a)*integers(t,16)
    for a,t in zip(yi,S): ci+=int(a)*integers(t,16)
    if B.size:
        wi=bi@bi.T
        for a in range(len(B)):
            for b in range(len(B)): ci-=wi[a,b]*integers(G[a,b],8)
    ci=np.array([[int(a) for a in row] for row in ci],dtype=object)
    assert np.array_equal(ci,ci.T)
    K=np.asarray(ci,dtype=float)/2**32
    lam=eigh(K,subset_by_index=[0,0],eigvals_only=True)[0]
    ei=int(np.floor((lam-1e-7)*2**32))
    Q=np.linalg.cholesky(K-ei/2**32*np.eye(len(K)))
    qi=np.array([[int(a) for a in row] for row in np.rint(Q*2**24)],dtype=object)
    R=ci*2**16-ei*2**16*np.eye(len(K),dtype=object)-qi@qi.T
    delta=max(sum(abs(int(a)) for a in row) for row in R)
    numerator=ei*2**16-delta
    cert={'K_integer':ci.tolist(),'K_denominator':2**32,
          'X_coeff_integer':xi.tolist(),'X_denominator':2**16,
          'stationarity_integer':yi.tolist(),'stationarity_denominator':2**16,
          'B_integer':bi.tolist(),'B_denominator':2**12,
          'Q_integer':qi.tolist(),'Q_denominator':2**24,
          'e_trial_integer':ei,'residual_row_bound_integer':int(delta),
          'bound_numerator':int(numerator),'bound_denominator':2**48}
    return numerator/2**48,cert


def solve(h,d,k,ground,maxiter=500):
    C,DX,S,G=setup(h,d,k,ground); nx=len(DX); ns=len(S); m=len(G)
    z=np.zeros(nx+ns+m*m)
    if m: z[nx+ns:]=(np.eye(m)*.03).ravel()
    records=[]
    for tau in (.03,.006,.0012,.00024):
        def fun(z):
            x,y,B=unpack(z,nx,ns,m)
            K=matrix(C,DX,S,G,x,y,B)
            e,U=eigh(K); weights=np.exp(-(e-e[0])/tau); weights/=weights.sum()
            rho=(U*weights)@U.T
            f=-e[0]+tau*np.log(np.exp(-(e-e[0])/tau).sum())
            gx=-np.einsum('ij,aji->a',rho,DX)
            gy=-np.einsum('ij,aji->a',rho,S)
            gb=np.zeros((m,m))
            if m:
                GR=np.einsum('ij,abji->ab',rho,G)
                gb=(GR+GR.T)@B
            return f,np.r_[gx,gy,gb.ravel()]
        opt=minimize(fun,z,jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-13,'gtol':1e-8,'maxls':30})
        z=opt.x; records.append({'tau':tau,'iterations':int(opt.nit),'success':bool(opt.success),'message':str(opt.message)})
    x,y,B=unpack(z,nx,ns,m)
    lb,cert=certify(C,DX,S,G,x,y,B)
    return lb,cert,records


def ring(h,d,L):
    D=d**L; H=np.zeros((D,D))
    digits=np.array(np.unravel_index(np.arange(D),(d,)*L)).T
    powers=d**np.arange(L-1,-1,-1)
    for col,word in enumerate(digits):
        for i in range(L):
            j=(i+1)%L; ab=d*word[i]+word[j]
            for cd in range(d*d):
                c,b=divmod(cd,d)
                row=col+(c-word[i])*powers[i]+(b-word[j])*powers[j]
                H[row,col]+=h[cd,ab]
    return H


def main():
    p=argparse.ArgumentParser(); p.add_argument('--small',action='store_true'); p.add_argument('--extended',action='store_true')
    args=p.parse_args()
    if not args.small and not args.extended: p.error('select --small or --extended')
    if args.extended and 'SLURM_JOB_ID' not in os.environ: p.error('extended mode requires a SLURM allocation')
    out=Path('certificates'); out.mkdir(exist_ok=True)
    rows=[]
    for name,d,h in bonds():
        ks=[2,3] if not args.extended else ([2,3,4,5] if d==2 else [2,3,4])
        for k in ks:
            for gs in ([False,True] if k>=3 else [False]):
                t=time.time(); lb,cert,rec=solve(h,d,k,gs)
                cert.update({'model':name,'d':d,'k':k,'ground_constraints':gs,'bond':h.tolist()})
                path=out/(name+'_k'+str(k)+('_ground' if gs else '_marginal')+'.json')
                path.write_text(json.dumps(cert,indent=2)+'\n')
                row={'model':name,'d':d,'k':k,'ground':gs,'lower_per_site':lb,'seconds':time.time()-t,'optimizer':rec,'certificate':str(path)}
                rows.append(row); print(json.dumps(row),flush=True)
    # Tiny rings only. These floating eigenvalues are diagnostics, not certified endpoints.
    energies=[]
    for name,d,h in bonds():
        for L in ([4,6] if d==2 else [3,4]):
            E=float(eigh(ring(h,d,L),subset_by_index=[0,0],eigvals_only=True)[0])
            energies.append({'model':name,'d':d,'L':L,'energy':E,'per_site':E/L})
    Path('numerics.json').write_text(json.dumps({'bounds':rows,'ring_diagnostics':energies},indent=2)+'\n')

if __name__=='__main__': main()
