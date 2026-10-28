"""Finite-field specialization of the polynomial-size leading g^2 bond.
Integer operations modulo p. NumPy int64 products are used only when
all intermediate integer sums are below 2**63.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
import time

P=65521

def dot(a,b,p=P):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    assert a.shape[-1]*(p-1)**2 < 2**63
    return (a@b)%p

def inv(a,p=P):
    a=np.array(a,dtype=np.int64,copy=True); n=len(a); b=np.eye(n,dtype=np.int64)
    for k in range(n):
        row=next((j for j in range(k,n) if a[j,k]),None)
        if row is None:raise ArithmeticError('singular')
        if row!=k:a[[row,k]]=a[[k,row]];b[[row,k]]=b[[k,row]]
        r=pow(int(a[k,k]),-1,p);a[k]=a[k]*r%p;b[k]=b[k]*r%p
        c=a[:,k].copy();c[k]=0
        a=(a-c[:,None]*a[k])%p;b=(b-c[:,None]*b[k])%p
    return b

class SparseLU:
    def __init__(self,A,p=P):
        self.p=p;self.N=A.shape[0]
        graph=A.copy().astype(float);graph.data[:]=1
        graph=graph+sp.eye(self.N,format='csc')*10
        lu=splu(graph,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0)
        self.order=np.argsort(lu.perm_c)
        B=A[self.order,:][:,self.order].tocsr()
        rows=[{int(j):int(v)%p for j,v in zip(B.indices[B.indptr[i]:B.indptr[i+1]],B.data[B.indptr[i]:B.indptr[i+1]]) if int(v)%p} for i in range(self.N)]
        cols=[set() for _ in rows]
        for i,row in enumerate(rows):
            for j in row:
                if i>j:cols[j].add(i)
        self.L=[[] for _ in rows];self.U=[];self.di=[]
        for k in range(self.N):
            pivot=rows[k].get(k,0)
            if not pivot:raise ArithmeticError(('zero pivot',k))
            di=pow(pivot,-1,p);self.di.append(di)
            rr=[(j,v) for j,v in rows[k].items() if j>k and v]
            self.U.append((np.array([j for j,v in rr],dtype=int),np.array([v for j,v in rr],dtype=np.int64)))
            for i in list(cols[k]):
                val=rows[i].pop(k,0)
                if not val:continue
                fac=val*di%p;self.L[i].append((k,fac))
                for j,v in rr:
                    old=rows[i].get(j,0);new=(old-fac*v)%p
                    if new:
                        rows[i][j]=new
                        if i>j:cols[j].add(i)
                    elif old:
                        del rows[i][j]
                        if i>j:cols[j].discard(i)
            rows[k].clear();cols[k].clear()
        self.L=[(np.array([j for j,v in r],dtype=int),np.array([v for j,v in r],dtype=np.int64)) for r in self.L]
    def solve(self,B):
        vec=B.ndim==1
        if vec:B=B[:,None]
        p=self.p;X=np.array(B[self.order],dtype=np.int64,copy=True)
        for i,(js,vs) in enumerate(self.L):
            if len(js):X[i]=(X[i]-dot(vs,X[js],p))%p
        for i in range(self.N-1,-1,-1):
            js,vs=self.U[i]
            if len(js):X[i]=(X[i]-dot(vs,X[js],p))%p
            X[i]=X[i]*self.di[i]%p
        out=np.empty_like(X);out[self.order]=X
        return out[:,0] if vec else out

def matrices(n,D=1000,a=4,td_exp=2000,eps_exp=1000,p=P,td_zero=False,eps_zero=False):
    ii=next(j for j in range(2,p) if j*j%p==p-1)
    td=pow(pow(10,td_exp,p),-1,p);eps=pow(pow(10,eps_exp,p),-1,p)
    if td_zero:td=0
    if eps_zero:eps=0
    s=2*n+1
    def carrier(typ,sgn):
        H=np.diag([a]*n+[a+D]*(n+1)).astype(np.int64)
        for j,k in [(j,j+1) for j in range(n-1)]+[(n+j,n+j+1) for j in range(n)]+[(0,n)]:H[j,k]=H[k,j]=-1
        c=td*(sgn*ii if typ==0 else 1)%p
        cc=td*(-sgn*ii if typ==0 else 1)%p
        H[s-1,n-1]-=c;H[n-1,s-1]-=cc
        return sp.csc_matrix(H%p)
    Hp=[carrier(j,1) for j in [0,1]];Hm=[carrier(j,-1) for j in [0,1]]
    I=sp.eye(s,dtype=np.int64,format='csc')
    def pair(x,y,minus=False):
        H=Hm if minus else Hp
        return (sp.kron(H[x!=0],I)+sp.kron(I,H[y!=0])).tocsc()
    b=np.zeros(s*s,dtype=np.int64);b[(s-1)*s]=1
    sources={}
    for key in [(0,0),(1,2),(2,1)]:
        x,y=key
        for sign in [1,-1]:
            lu=SparseLU(pair(x,y,sign==-1),p)
            sol=lu.solve(b)
            assert np.array_equal((pair(x,y,sign==-1)@sol)%p,b)
            sources[key,sign]=sol.reshape(s,s)
    # Conversion source vectors, before their bit-dependent signs and port phases.
    basis={};groups={}
    for x,y in [(0,0),(1,2),(2,1)]:
        for cell,typ in enumerate([x,y]):
            for ns in ([1,2] if typ==0 else [0]):
                target=(ns,y) if cell==0 else (x,ns)
                for j in range(n):
                    key=(x,y,cell,ns,j);idx=len(groups.setdefault(target,[]));groups[target].append(key);basis[key]=(target,idx)
    gram={}
    for target,keys in groups.items():
        vp=np.zeros((s*s,len(keys)),dtype=np.int64);vm=np.zeros_like(vp)
        for col,(x,y,cell,ns,j) in enumerate(keys):
            for sign,V in [(1,vp),(-1,vm)]:
                R=sources[(x,y),sign];z=V[:,col].reshape(s,s)
                if cell==0:z[j,:]=R[j,:]
                else:z[:,j]=R[:,j]
        lu=SparseLU(pair(*target),p);w=lu.solve(vp)
        gram[target]=-dot(vm.T,w,p)%p
        assert np.array_equal((pair(*target)@w)%p,vp%p)
    def terms(bits):
        out=[]
        for x,y in [(0,0),(1,2),(2,1)]:
            coeff=(-ii*((-1)**((bits&((1<<n)-1)).bit_count()))) if x==0 else (ii if x==1 else -ii)
            # All source coefficients are purely imaginary, so conjugation is minus.
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
                    target,idx=basis[(x,y,cell,ns,j)]
                    out.append(((target,bits^(1<<(cell*n+j))),idx,coeff*sg%p))
        return out
    ts=[terms(b) for b in [0]+[1<<j for j in range(2*n)]]
    out=np.zeros((2*n+1,2*n+1),dtype=np.int64)
    for i,xx in enumerate(ts):
        for j,yy in enumerate(ts):
            total=0
            for (key,k,c) in xx:
                for (key2,l,d) in yy:
                    if key==key2:total-=c*d*int(gram[key[0]][k,l])
            out[i,j]=total%p
    h=(out[1:,1:]-out[0,0]*np.eye(2*n,dtype=np.int64))%p
    T=h[:n,n:];Tc=h[n:,:n];M=(h[:n,:n]+h[n:,n:])%p
    return T,Tc,M,eps,ii
