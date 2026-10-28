"""Exact g^2 hopping and onsite matrix using only 4n source solves."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from structured import carrier,visibility,Leading

def fast(n,D=100.,a=2.2,td=0.,loopflux=np.pi):
    S=2*n+1
    H0,Hc,_=carrier(n,D,a,td,0.,loopflux-np.pi/2)
    def factor(A,B):return splu((sp.kron(A,sp.eye(S))+sp.kron(sp.eye(S),B)).tocsc())
    lu00=factor(H0,H0);lucc=factor(Hc,Hc);lu0c=factor(H0,Hc);luc0=factor(Hc,H0)
    b=np.zeros(S*S,complex);b[(S-1)*S]=1
    x=lu00.solve(b).reshape(S,S);y=lucc.solve(b).reshape(S,S)
    def restrict(v):
        U=np.zeros((S*S,n),complex);V=U.copy()
        for j in range(n):U[j*S:(j+1)*S,j]=v[j,:];V[j::S,j]=v[:,j]
        return U,V
    U,V=restrict(x);W,Z=restrict(y)
    RU=luc0.solve(U);RV=lu0c.solve(V);RW=lu0c.solve(W);RZ=luc0.solve(Z)
    # First cross matrix above has (j,i), second has (i,j).
    T=-2*visibility(n)*np.real((Z.conj().T@RU).T+W.conj().T@RV)
    M=-2j*visibility(n)*np.imag(U.conj().T@RU+W.conj().T@RW+V.conj().T@RV+Z.conj().T@RZ)
    return T,M
if __name__=='__main__':
 for n in [2,3,4]:
    P=dict(n=n,D=10.,a=3.,td=1e-5)
    T,M=fast(**P);R,N,_=Leading(**P,loopflux=np.pi).calculate()
    print(n,'T relative',np.max(abs(T-R))/np.max(abs(T)),'M imaginary relative',np.max(abs(M-1j*N.imag))/np.max(abs(M)),flush=True)
