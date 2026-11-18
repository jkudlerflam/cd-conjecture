"""Actual C84 loop derivative, with n active sites and n+1 data vertices.

All matrices omit the common physical factor lambda^2 g^2 |v0 vp|^2.
The imaginary loop is H[M,D_n]=-i t_d, on both conversion species.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import splu


def carrier(n,D,w,a=4.,td=0.):
    p=n+1;s=2*p+1
    H=np.diag([a]*p+[a+D]*(p+1)).astype(complex if td else float)
    for j in range(p-1):H[j,j+1]=H[j+1,j]=-w
    for j in range(p,2*p):H[j,j+1]=H[j+1,j]=-1
    H[0,p]=H[p,0]=-1
    J=np.zeros((s,s));J[-1,p-1]=-1;J[p-1,-1]=1
    if td:H+=1j*td*J
    return H,J


def restrictions(x,n,s):
    U=np.zeros((s*s,n),dtype=x.dtype);V=U.copy()
    x=x.reshape(s,s)
    for j in range(n):U[j*s:(j+1)*s,j]=x[j,:];V[j::s,j]=x[:,j]
    return U,V


def coefficients(n,D,w,a=4.,td=0.):
    H,J=carrier(n,D,w,a,td);s=len(H)
    L=sp.kronsum(sp.csc_matrix(H),sp.csc_matrix(H),format='csc')
    lu=splu(L);b=np.zeros(s*s);b[(s-1)*s]=1
    x=lu.solve(b);U,V=restrictions(x,n,s)
    RU,RV=lu.solve(U),lu.solve(V)
    C=U.conj().T@RV;A=U.conj().T@RU+V.conj().T@RV
    T=-2*C.real;M=-2j*A.imag
    if td:return T,M
    J2=sp.kronsum(sp.csc_matrix(J),sp.csc_matrix(J),format='csc')
    y=-lu.solve(J2@x)  # R'(0)b=i y
    U1,V1=restrictions(y,n,s)
    B=U.T@lu.solve(U1)-U1.T@RU+V.T@lu.solve(V1)-V1.T@RV
    B-=RU.T@(J2@RU)+RV.T@(J2@RV)
    assert la.norm(B+B.T)<1e-9*max(la.norm(B),1e-100)
    return T,-2j*B


def low_limit(n,w,a=4.):
    """Leading large-D tensors for eta=t_d D^p and the delay skew."""
    p=n+1;A=a*np.eye(p)-w*(np.eye(p,k=1)+np.eye(p,k=-1))
    R=la.inv(np.kron(A,np.eye(p))+np.kron(np.eye(p),A))
    x0=R[:,0];xN=R[:,n*p]
    U,V=restrictions(x0,n,p);UN,VN=restrictions(xN,n,p)
    G=2*U.T@R@V
    B=U.T@R@UN-UN.T@R@U+V.T@R@VN-VN.T@R@V
    TN=-2*UN.T@R@VN
    # Differentiate one return filter in its two time orderings.
    y=R@np.kron(A,np.eye(p))[:,0]
    YU,YV=restrictions(y,n,p)
    Q=-YU.T@R@V-U.T@R@YV
    Q[0,:]+=x0@V
    assert la.norm(Q+Q.T)<1e-10*la.norm(G)
    return G,B,TN,Q


def roots_and_response(T,M1,F):
    n=len(T);S=-(T+T.T)/2;L=la.cholesky(S,lower=True);C=la.inv(L)
    T=C@T@C.T;M1=C@M1@C.T;F=C@F@C.T
    I=np.eye(n);Z=np.zeros((n,n))
    zs=la.eigvals(np.block([[-F,-T.T],[I,Z]]),np.block([[T,Z],[Z,I]]))
    assert np.max(abs(abs(zs)-1))<1e-7
    positive=sorted(float(np.angle(z)) for z in zs if np.angle(z)>0)
    result=[]
    for k in positive:
        H=T*np.exp(1j*k)+T.T*np.exp(-1j*k)+F
        e,U=la.eigh(H);j=np.argmin(abs(e));u=U[:,j]
        speed=float(np.vdot(u,(1j*T*np.exp(1j*k)-1j*T.T*np.exp(-1j*k))@u).real)
        response=-2*float(np.vdot(u,M1@u).real)/speed
        result.append((k,response))
    return result
