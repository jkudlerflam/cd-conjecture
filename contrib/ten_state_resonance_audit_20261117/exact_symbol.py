"""Exact rational scalar inverse and the two-flavour four-root identity."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np


def inverse(A):
    n=len(A);B=[[F(int(A[i,j])) for j in range(n)]+[F(int(i==j)) for j in range(n)] for i in range(n)]
    for i in range(n):
        pivot=next(j for j in range(i,n) if B[j][i])
        B[i],B[pivot]=B[pivot],B[i];v=B[i][i];B[i]=[x/v for x in B[i]]
        for j in range(n):
            if j!=i:
                v=B[j][i]
                if v:B[j]=[x-v*y for x,y in zip(B[j],B[i])]
    return np.array([row[n:] for row in B],dtype=object)


def main():
    n=2;s=n+2;H=2*np.eye(s,dtype=int)-np.eye(s,k=1,dtype=int)-np.eye(s,k=-1,dtype=int)
    D=np.kron(H,np.eye(s,dtype=int))+np.kron(np.eye(s,dtype=int),H)
    R=inverse(D);Q=np.array([[F(int(x)) for x in row] for row in D],dtype=object)
    assert np.array_equal(Q@R,np.eye(s*s,dtype=int))
    source=np.array([F(int(i==(s-1)*s)) for i in range(s*s)],dtype=object)
    x=(R@source).reshape(s,s);U=np.full((s*s,n),F(0),dtype=object);V=U.copy()
    for i in range(n):U[i*s:(i+1)*s,i]=x[i,:];V[i::s,i]=x[:,i]
    T=-2*(U.T@R@V);S=(T+T.T)/2;A=(T-T.T)/2
    determinant=lambda M:M[0,0]*M[1,1]-M[0,1]*M[1,0]
    detS=determinant(S);detT=determinant(T);rho2=A[0,1]**2/(-detS)
    assert detS<0 and detT<0 and 0<rho2<1
    center=F(1,5);radius=F(2,3);x0=(1-rho2)*center
    eps2=4*(rho2+(1-rho2)*(radius**2-rho2*center**2))
    trS=S[0,0]+S[1,1]
    assert trS<0 and x0>0 and eps2*(trS*trS-4*detS)>4*trS*trS
    for root in (center-radius,center+radius):
        assert eps2-4*((root-x0)**2+rho2*(1-root**2))==0
    out=dict(T=[[str(x) for x in row] for row in T],detS=str(detS),detT=str(detT),rho_squared=str(rho2),x0=str(x0),epsilon_squared=str(eps2),cos_roots=[str(center-radius),str(center+radius)],
             identity='det h(k) = 4*(-detS)*(1-rho_squared)*(radius²-(cos(k)-center)²)',exact_inverse_and_root_checks=True)
    (Path(__file__).parent/'exact_symbol.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
