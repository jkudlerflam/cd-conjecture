# C80 local copy of C75 source, used as a reference or scalar certificate.
# C75 reuse of C33 scalar-certificate implementation. See PROOF.md.
"""Integer interval certificate for the normalized finite-D g^2 hopping.
All positivity decisions use Python integers, not floating-point eigenvalues.
The scalar product-path inverse is evaluated by floor-rounded Jacobi steps.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import time,json,math,sys

def bareiss(A,positive=False):
    A=[list(map(int,row)) for row in A];n=len(A)
    if n==0:return 1
    prev=1;sgn=1
    for k in range(n-1):
        if A[k][k]==0:
            if positive:return None
            row=next((i for i in range(k+1,n) if A[i][k]),None)
            if row is None:return 0
            A[k],A[row]=A[row],A[k];sgn=-sgn
        pivot=A[k][k]
        if positive and pivot<=0:return None
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=A[i][j]*pivot-A[i][k]*A[k][j]
                assert val%prev==0
                A[i][j]=val//prev
            A[i][k]=0
        prev=pivot
    if positive and A[-1][-1]<=0:return None
    return sgn*A[-1][-1]

def log10int(i):
    return math.log10(i)

def certificate(n,D=1000,bits=256,max_seconds=175):
    start=time.monotonic();m=n+1;S=2*n+1
    # det(C), C has diagonal D+4 and off-diagonal -1.
    p0,p1=1,D+4
    for k in range(2,m+1):p0,p1=p1,(D+4)*p1-p0
    detC=p1;scale=1<<bits;N=detC.bit_length()+bits+10
    diag=np.array([D+4]*m+[4]*n,dtype=object)
    den=diag[:,None]+diag[None,:]
    def solve(B):
        X=np.zeros(B.shape,dtype=object)
        div=den if B.ndim==2 else den[:,:,None]
        for k in range(N):
            Y=B.copy()
            Y[1:,...]+=X[:-1,...];Y[:-1,...]+=X[1:,...]
            Y[:,1:,...]+=X[:,:-1,...];Y[:,:-1,...]+=X[:,1:,...]
            X=Y//div
            if k%40==0 and time.monotonic()-start>max_seconds:raise TimeoutError((n,k,B.shape))
        return X
    B=np.zeros((S,S),dtype=object);B[0,m]=detC*scale
    X=solve(B)
    V=np.zeros((S,S,n),dtype=object)
    for j in range(n):V[:,m+j,j]=X[:,m+j]
    Y=solve(V)
    # x error <=3/scale, y error <=4/scale, including truncation and rounding.
    Q=np.empty((n,n),dtype=object);U=np.empty((n,n),dtype=object)
    for i in range(n):
        for j in range(n):
            Q[i,j]=sum(X[m+i,s]*Y[m+i,s,j] for s in range(S))
            U[i,j]=sum((X[m+i,s]+3)*(Y[m+i,s,j]+4) for s in range(S))
    L=Q+Q.T;E=U+U.T-L
    err=max(sum(row) for row in E)
    det=bareiss(L,positive=True)
    if det is None:
        result=dict(n=n,D=D,bits=bits,status='lower_matrix_not_positive',seconds=time.monotonic()-start)
    else:
        minors=0
        for k in range(n):
            inds=[i for i in range(n) if i!=k]
            minors+=bareiss(L[np.ix_(inds,inds)])
        # lambda_min(L/scale^2)>=det/minors/scale^2.
        good=det>err*minors
        lower_num=det-err*minors
        result=dict(n=n,D=D,bits=bits,iterations=N,product_dimension=S*S,
            status='proved_positive' if good else 'interval_too_wide',
            log10_lambda_min_lower=(log10int(lower_num)-log10int(minors)-2*bits*math.log10(2)) if good else None,
            log10_lambda_min_upper=log10int(n*det+err*minors)-log10int(minors)-2*bits*math.log10(2),
            log10_lambda_max_lower=log10int(max(L[i,i] for i in range(n)))-2*bits*math.log10(2),
            log10_lambda_max_upper=log10int(sum(U[i,i]*2 for i in range(n)))-2*bits*math.log10(2),
            log10_entry_radius=log10int(err)-2*bits*math.log10(2),
            log10_beta=-log10int(detC),seconds=time.monotonic()-start)
        np.savez(f'certificate_n{n}_D{D}.npz',lower_numerators=L,radius_numerators=E,F_half_numerators=Q,
                 detC=np.array(detC,dtype=object),bits=bits)
    return result

if __name__=='__main__':
    ns=list(map(int,sys.argv[1:])) or [2,4,8,12,16,23]
    for n in ns:
        r=certificate(n);print(json.dumps(r),flush=True)
        json.dump(r,open(f'certificate_n{n}.json','w'),indent=2)
