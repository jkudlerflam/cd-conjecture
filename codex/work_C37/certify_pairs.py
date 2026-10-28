"""Exact finite-field pair-product certificate, with replayable polynomial data."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import json,time,sys
from modular import matrices,P

def det(A,p):
    A=np.array(A,dtype=np.int64,copy=True);n=len(A);v=1
    for k in range(n):
        r=next((j for j in range(k,n) if A[j,k]),None)
        if r is None:return 0
        if r!=k:A[[r,k]]=A[[k,r]];v=-v
        d=int(A[k,k]);v=v*d%p
        if k+1<n:
            f=A[k+1:,k]*pow(d,-1,p)%p
            A[k+1:,k+1:]=(A[k+1:,k+1:]-f[:,None]*A[k,k+1:])%p
    return v%p

def poly(T,Tc,M,eps,p):
    n=len(T);B=(M+eps*np.eye(n,dtype=np.int64))%p
    vals=np.array([det((T*z*z+B*z+Tc)%p,p) for z in range(2*n+1)],dtype=np.int64)
    out=np.zeros(2*n+1,dtype=np.int64);fall=np.array([1],dtype=np.int64);fact=1
    for k in range(2*n+1):
        out[:k+1]=(out[:k+1]+fall*int(vals[0])*pow(fact,-1,p))%p
        vals=np.diff(vals)%p
        if k<2*n:
            fall=np.convolve(fall,np.array([-k,1],dtype=np.int64))%p;fact=fact*(k+1)%p
    return out

def trim(a):
    a=list(map(int,a))
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def gcd(a,b,p):
    a=trim(a);b=trim(b)
    while any(b):
        inv=pow(b[-1],-1,p)
        while len(a)>=len(b) and any(a):
            c=a[-1]*inv%p;shift=len(a)-len(b)
            for j,x in enumerate(b):a[j+shift]=(a[j+shift]-c*x)%p
            a=trim(a)
        a,b=b,a
    return [(x*pow(a[-1],-1,p))%p for x in a]

def pair_polynomial(f,p):
    m=len(f)-1;N=m*(m+1)//2
    if not f[-1] or not f[0]:raise ArithmeticError('degree drop or zero root')
    c=(np.array(f,dtype=np.int64)*pow(int(f[-1]),-1,p))%p
    s=[m]
    for k in range(1,2*N+1):
        val=sum(int(c[m-j])*s[k-j] for j in range(1,min(k-1,m)+1))
        if k<=m:val+=k*int(c[m-k])
        s.append(-val%p)
    t=[0]+[((s[k]*s[k]+s[2*k])*pow(2,-1,p))%p for k in range(1,N+1)]
    q=[1]
    for k in range(1,N+1):q.append(-sum(q[k-j]*t[j] for j in range(1,k+1))*pow(k,-1,p)%p)
    return list(reversed(q))

def run(n,D=1000,a=4,td_exp=2000,eps_exp=1000,p=P):
    start=time.monotonic();T,Tc,M,eps,ii=matrices(n,D,a,td_exp,eps_exp,p)
    f=poly(T,Tc,M,eps,p)
    checks=[]
    for z in [2*n+1,2*n+3,2*n+7]:
        value=sum(int(c)*pow(z,j,p) for j,c in enumerate(f))%p
        direct=det((T*z*z+(M+eps*np.eye(n,dtype=np.int64))*z+Tc)%p,p)
        assert value==direct
        checks.append([z,value])
    q=pair_polynomial(f,p)
    g=gcd(q,[(j*q[j])%p for j in range(1,len(q))],p)
    r=dict(n=n,D=D,a=a,td_power=-td_exp,field_power=-eps_exp,p=p,image_of_i=ii,
      carrier_dimension=2*n+1,pair_resolvent_dimension=(2*n+1)**2,
      f=f.tolist(),pair_product_polynomial=q,gcd=g,extra_determinant_checks=checks,
      status='proved_no_pair_product_collision' if len(g)==1 else 'inconclusive_modular_collision',
      seconds=time.monotonic()-start)
    json.dump(r,open(f'pairs_n{n}.json','w'),indent=2)
    print(json.dumps({k:v for k,v in r.items() if k not in ['f','pair_product_polynomial','gcd']}),flush=True)
    return r

if __name__=='__main__':
    for n in (list(map(int,sys.argv[1:])) or [1,2,3,4]):run(n)
