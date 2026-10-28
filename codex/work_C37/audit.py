"""Independent finite checks for the modular proof components."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import scipy.sparse as sp
from modular import SparseLU,inv,dot,P,matrices
from certify_pairs import pair_polynomial,gcd,poly,det
import json,time,math

def test_solver():
    records=[]
    for n in range(1,6):
        s=2*n+1
        H=np.diag([4]*n+[1004]*(n+1)).astype(np.int64)
        for j in range(s-1):H[j,j+1]=H[j+1,j]=-1
        H[0,-1]=-24297;H[-1,0]=24297
        A=(np.kron(H,np.eye(s,dtype=np.int64))+np.kron(np.eye(s,dtype=np.int64),H))%P
        B=np.arange(s*s*7,dtype=np.int64).reshape(s*s,7)%P
        X=SparseLU(sp.csc_matrix(A)).solve(B)
        Y=dot(inv(A),B)
        assert np.array_equal(X,Y)
        assert np.array_equal(dot(A,X),B)
        records.append(dict(n=n,dimension=s*s,rhs=7,match=True))
    return records

def test_symmetric_square():
    records=[]
    for m in [2,4,6,10,20,46]:
        roots=list(range(2,m+2));f=np.array([1],dtype=np.int64)
        for z in roots:f=np.convolve(f,[-z,1])%P
        q=pair_polynomial(f,P);ref=np.array([1],dtype=np.int64)
        for i,x in enumerate(roots):
            for y in roots[i:]:ref=np.convolve(ref,[-x*y,1])%P
        assert np.array_equal(q,ref)
        records.append(dict(roots=m,pairs=m*(m+1)//2,match=True))
    return records

def direct_leading(n,D=2,a=3,L=1,E=2):
    p=P;ii=24297;td=pow(pow(10,L,p),-1,p);s=2*n+1
    def h(typ,sg):
        H=np.diag([a]*n+[a+D]*(n+1)).astype(np.int64)
        for j,k in [(j,j+1) for j in range(n-1)]+[(n+j,n+j+1) for j in range(n)]+[(0,n)]:H[j,k]=H[k,j]=-1
        phase=sg*ii if typ==0 else 1
        con=-sg*ii if typ==0 else 1
        H[-1,n-1]-=td*phase;H[n-1,-1]-=td*con
        return H%p
    rs={};lus={};source=np.zeros(s*s,dtype=np.int64);source[(s-1)*s]=1
    for sign in [1,-1]:
        hs=[h(t,sign) for t in range(3)]
        for x in range(3):
            for y in range(3):
                A=(np.kron(hs[x],np.eye(s,dtype=np.int64))+np.kron(np.eye(s,dtype=np.int64),hs[y]))%p
                lus[x,y,sign]=inv(A)
                rs[x,y,sign]=dot(lus[x,y,sign],source).reshape(s,s)
    def Y(bits,sign):
        out={}
        for x,y,c in [(0,0,-sign*ii*(-1)**((bits&((1<<n)-1)).bit_count())),(1,2,sign*ii),(2,1,-sign*ii)]:
            R=rs[x,y,sign]*c%p
            for cell,typ in enumerate([x,y]):
                local=(bits>>(cell*n))&((1<<n)-1)
                for j in range(n):
                    bit=(local>>j)&1
                    if typ==0:ns=1 if bit else 2
                    elif typ==1:
                        if bit:continue
                        ns=0
                    else:
                        if not bit:continue
                        ns=0
                    sg=(-1)**((local&((1<<j)-1)).bit_count())
                    key=(ns,y,bits^(1<<(cell*n+j))) if cell==0 else (x,ns,bits^(1<<(cell*n+j)))
                    if key not in out:out[key]=np.zeros((s,s),dtype=np.int64)
                    if cell==0:out[key][j,:]=(out[key][j,:]+sg*R[j,:])%p
                    else:out[key][:,j]=(out[key][:,j]+sg*R[:,j])%p
        return out
    bits=[0]+[1<<j for j in range(2*n)]
    plus=[Y(b,1) for b in bits];minus=[Y(b,-1) for b in bits]
    out=np.zeros((2*n+1,2*n+1),dtype=np.int64)
    for j,Yj in enumerate(plus):
        for key,v in Yj.items():
            w=dot(lus[key[0],key[1],1],v.ravel())
            for i,Yi in enumerate(minus):
                if key in Yi:out[i,j]=(out[i,j]-int(dot(Yi[key].ravel(),w)))%p
    H=(out[1:,1:]-out[0,0]*np.eye(2*n,dtype=np.int64))%p
    T,Tc,M,eps,im=matrices(n,D,a,L,E)
    assert np.array_equal(H[:n,n:],T)
    assert np.array_equal(H[n:,:n],Tc)
    assert np.array_equal((H[:n,:n]+H[n:,n:])%p,M)
    return dict(n=n,pair_dimension=s*s,logical_inputs=2*n+1,match=True)

if __name__=='__main__':
    start=time.monotonic()
    assert all(P%d for d in range(2,math.isqrt(P)+1))
    assert 24297**2%P==P-1
    data=dict(prime=P,prime_trial_division=True,solver=test_solver(),
      symmetric_square=test_symmetric_square(),leading=[direct_leading(n) for n in [1,2,3]],seconds=time.monotonic()-start)
    json.dump(data,open('audit_results.json','w'),indent=2)
    print(json.dumps(data),flush=True)
