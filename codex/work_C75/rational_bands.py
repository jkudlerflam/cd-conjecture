# C75 reuse of C33 scalar-certificate implementation. See PROOF.md.
"""Fraction-free whitening and a certified-entry loop derivative.
No floating point is used until the final small Bloch matrices are formed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,time,math,json
from scipy.linalg import eig
from loop_scan import margins

def integer_ldl(N):
    A=N.copy();n=len(A);R=np.eye(n,dtype=object);prev=1;ds=[]
    for k in range(n):
        p=int(A[k,k]);assert p>0
        ds.append(prev*p)
        if k==n-1:break
        col=A[:,k].copy()
        for i in range(k+1,n):
            for j in range(n):
                v=p*R[i,j]-col[i]*R[k,j]
                assert v%prev==0;R[i,j]=v//prev
        for i in range(k+1,n):
            for j in range(k+1,n):
                v=p*A[i,j]-col[i]*col[j]
                assert v%prev==0;A[i,j]=v//prev
        prev=p
    # All congruence operations are exact integer identities.
    assert np.array_equal(R@N@R.T,np.diag(np.array(ds,dtype=object)))
    return R,ds

def congruence(B,R,ds,extra=1):
    C=R@B@R.T;n=len(ds)
    powers=[(int(d).bit_length()-1)//2 for d in ds]
    diag=np.array([int(d)/(1<<(2*b)) for d,b in zip(ds,powers)])
    out=np.array([[int(C[i,j])/(int(extra)*(1<<(powers[i]+powers[j]))) for j in range(n)] for i in range(n)])
    return out/np.sqrt(diag[:,None]*diag[None,:])

def derivative(n=23,D=1000,bits=256,max_seconds=175):
    start=time.monotonic();m=n+1;S=2*n+1;scale=1<<bits
    d0,d1=1,D+4
    for k in range(2,m+1):d0,d1=d1,(D+4)*d1-d0
    detC=d1;steps=detC.bit_length()+bits+10
    den=np.array([D+4]*m+[4]*n,dtype=object);den=den[:,None]+den[None,:]
    def solve(B):
        X=np.zeros(B.shape,dtype=object);div=den if B.ndim==2 else den[:,:,None]
        for k in range(steps):
            Y=B.copy();Y[1:,...]+=X[:-1,...];Y[:-1,...]+=X[1:,...]
            Y[:,1:,...]+=X[:,:-1,...];Y[:,:-1,...]+=X[:,1:,...]
            X=Y//div
            if k%40==0 and time.monotonic()-start>max_seconds:raise TimeoutError((n,k,B.shape))
        return X
    B=np.zeros((S,S),dtype=object);B[0,m]=detC*scale;X=solve(B)
    U=np.zeros((S,S,2*n),dtype=object)
    for j in range(n):U[m+j,:,j]=X[m+j,:];U[:,m+j,n+j]=X[:,m+j]
    Y=solve(U);RU=Y[:,:,:n];RV=Y[:,:,n:]
    def K(X):
        # K[M,end]=-1, K[end,M]=+1, H_loop'=i K.
        Z=np.zeros_like(X)
        Z[0,...]-=X[-1,...];Z[-1,...]+=X[0,...]
        Z[:,0,...]-=X[:,-1,...];Z[:,-1,...]+=X[:,0,...]
        return Z
    Z=solve(-K(X))
    A=np.zeros((n,n),dtype=object);F=np.zeros((n,n),dtype=object)
    for i in range(n):
        for j in range(n):
            F[i,j]=2*sum(X[m+i,s]*RV[m+i,s,j] for s in range(S))
            A[i,j]=sum(RU[m+j,s,i]*Z[m+j,s]+RV[s,m+j,i]*Z[s,m+j] for s in range(S))
    A=A-A.T
    for Y in [RU,RV]:
        # Only boundary rows/columns of K Y contribute.
        for i in range(n):
            for j in range(n):
                A[i,j]-=sum(-Y[0,s,i]*Y[-1,s,j]+Y[-1,s,i]*Y[0,s,j]
                             -Y[s,0,i]*Y[s,-1,j]+Y[s,-1,i]*Y[s,0,j] for s in range(S))
    assert np.array_equal(A,-A.T)
    np.savez(f'derivative_n{n}.npz',F_numerator=F,derivative_numerator=A,detC=np.array(detC,dtype=object),bits=bits,D=D)
    print(json.dumps(dict(stage='integer_derivative',n=n,seconds=time.monotonic()-start)),flush=True)
    return F,A,detC,bits

def bands_from(F,A,detC,bits,n):
    start=time.monotonic();N=F+F.T
    # S=(F+F.T)/2/scale^2. Work with N, giving symmetric hopping 1/2.
    R,ds=integer_ldl(N)
    Ft=congruence(F,R,ds)
    At=congruence(A,R,ds,extra=detC)
    norm=np.linalg.norm(2*At,2)
    eta_power=max(0,math.ceil(math.log2(norm)+10))
    eta=2.**(-eta_power)
    # Use the certified tree inverse-trace lower eigenvalue as field scale.
    cert=json.load(open(f'certificate_n{n}.json'))
    kappa=10**cert['log10_lambda_min_lower']
    scale=1<<bits
    field=.3*kappa*congruence(np.eye(n,dtype=object)*(scale*scale),R,ds)
    # In units 2 alpha^4 beta^2: T=-F, M=-2 i eta A/beta-denominator.
    Mt=-2j*eta*At+field
    Tt=-Ft
    I=np.eye(n);Z=np.zeros((n,n))
    zs=eig(np.block([[-Mt,-Tt.T],[I,Z]]),np.block([[Tt,Z],[Z,I]]),right=False)
    roots=np.sort(np.angle(zs));d2,d4,wi=margins(roots)
    # Symmetric tree roots identify partner branches under this small perturbation.
    centres=(roots[:n]+roots[n:][::-1])/2
    rec=dict(n=n,status='numerical_linear_loop_response',eta_power_of_two=-eta_power,
        whitened_loop_norm=float(eta*norm),unit_circle_error=float(max(abs(abs(zs)-1))),
        roots=roots.tolist(),centres=centres.tolist(),delta2=d2,delta4=d4,witness=wi,
        whitening_hermitian_residual=float(np.max(abs((Ft+Ft.T)/2-.5*I))),seconds=time.monotonic()-start)
    print(json.dumps(rec),flush=True)
    json.dump(rec,open(f'rational_bands_n{n}.json','w'),indent=2)
    np.savez(f'whitened_n{n}.npz',T=Tt,M=Mt,loop_derivative=At,R=R,ds=np.array(ds,dtype=object),field=field)
    return rec
if __name__=='__main__':
    import sys
    n=int(sys.argv[1]) if len(sys.argv)>1 else 23
    F,A,detC,bits=derivative(n)
    bands_from(F,A,detC,bits,n)
