"""C38 exact arithmetic checks and small numerical diagnostics.
Run with python3 -B checks.py. No imports from other work folders.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
from fractions import Fraction as F
from math import factorial
import json
import numpy as np
from scipy.linalg import solve_continuous_lyapunov


def circ_gap(points):
    points = sorted(p % 1 for p in points)
    gaps = [b-a for a,b in zip(points,points[1:])]
    gaps.append(1+points[0]-points[-1])
    assert min(gaps) > 0
    return min(gaps)


def target(n):
    r,b=F(1,16),F(1,64)
    ps=[(-F(1,4) if i%2==0 else F(1,4))+b*r**i for i in range(2*n)]
    pair=[(ps[i]+ps[j])%1 for i in range(2*n) for j in range(i,2*n)]
    d2,d4=circ_gap(ps),circ_gap(pair)
    lower=(1-2*r)*b*r**(2*n-1)
    assert d2 >= lower and d4 >= lower
    centers=[(ps[2*a]+ps[2*a+1])/2 for a in range(n)]
    assert all(centers[a+1]==r*r*centers[a] for a in range(n-1))
    return dict(n=n,points=2*n,pair_sums=len(pair),
                lower_normalized=str(lower),delta2_over_lower=str(d2/lower),
                delta4_over_lower=str(d4/lower),tag="proved by rational enumeration")


def C(i,j):
    return F(factorial(2*i+2*j+2),factorial(2*i+1)*factorial(2*j+1)*2**(2*i+2*j+3))


def gram(a,mobile=False):
    n=2
    A=a*np.eye(n)+np.diag([3.,1.]) if mobile else a*np.eye(n)
    A[0,1]=A[1,0]=-1.
    E=np.diag([1.,0.])
    B=np.kron(np.eye(n),A)
    Z=np.zeros_like(B)
    L=np.block([[B,-np.kron(np.eye(n),E)],[Z,B]])
    v=np.r_[np.zeros(n*n),np.eye(n).reshape(-1,order="F")]
    P=solve_continuous_lyapunov(L,np.outer(v,v))
    G=P[np.ix_([0,3],[0,3])]
    S=np.diag([1.,a*a])
    normalized=a**3*S@G@S
    ref=np.array([[float(C(i,j)) for j in range(n)] for i in range(n)])
    return dict(a=a,mobile=mobile,lyapunov_dimension=8,
                squared_correlation=float(G[0,1]**2/(G[0,0]*G[1,1])),
                normalized_relative_error=float(np.linalg.norm(normalized-ref)/np.linalg.norm(ref)),
                tag="numerical, double precision, n=2")


def mixing():
    ans=[]
    for q in [10,100,1000,10000]:
        t2=F(1,q**4); eps=F(1,10)
        # h11=1, h22=t2*f2, offdiag=eps/q^2.
        off=eps/F(q*q)
        correction=off*off
        assert correction/t2==eps*eps
        ans.append(dict(q=q,relative_schur_correction=str(correction/t2)))
    return ans

if __name__=="__main__":
    assert C(0,1)**2/(C(0,0)*C(1,1))==F(2,5)
    correlations=[]
    for j in range(64):
        ratio=C(j,j+1)**2/(C(j,j)*C(j+1,j+1))
        assert ratio==F((4*j+4)*(4*j+3),(4*j+6)*(4*j+5))
        if j in [0,1,2,7,31,63]: correlations.append(dict(j=j,ratio=str(ratio)))
    out=dict(targets=[target(n) for n in [1,2,3,4,8,16,23,32,64,128]],
             adjacent_correlations=correlations,
             first_schur_pivot=str(C(1,1)-C(0,1)**2/C(0,0)),
             mixing=mixing(),
             gram_diagnostics=[gram(a,m) for m in [False,True] for a in [8.,16.,32.,64.,128.]])
    with open("checks_results.json","w") as f: json.dump(out,f,indent=2)
    for x in out["targets"]: print("exact target",x["n"],x["pair_sums"],"pair sums; margin ratio",x["delta4_over_lower"])
    print("exact first squared correlation",F(2,5),"Schur pivot",out["first_schur_pivot"])
    for x in out["gram_diagnostics"]: print("Gram",x)
