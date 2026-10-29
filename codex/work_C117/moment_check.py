"""Small finite-ring check of a massive contribution to the q^4 moment."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from scipy.linalg import eigh
from bootstrap import ring

def comm(a,b):return a@b-b@a

def embed_bond(h,L,i):
    D=2**L; out=np.zeros((D,D)); j=(i+1)%L
    for col in range(D):
        word=[(col//2**(L-1-k))%2 for k in range(L)]
        for a in range(2):
            for b in range(2):
                row=col+(a-word[i])*2**(L-1-i)+(b-word[j])*2**(L-1-j)
                out[row,col]=h[2*a+b,2*word[i]+word[j]]
    return out

rows=[]
X=np.array([[0.,1],[1,0]]); Z=np.diag([1.,-1]); I=np.eye(2)
for J in (.025,.05,.1):
    g=1.; L=6; h=-J*np.kron(X,X)-g/2*(np.kron(Z,I)+np.kron(I,Z))
    H=ring(h,2,L); en,U=eigh(H); psi=U[:,0]; hs=[embed_bond(h,L,i) for i in range(L)]
    b=[float(psi@comm(hs[0],comm(H,hs[r]))@psi) for r in range(3)]
    a2=-b[1]/2-2*b[2]; a4=b[1]/24+2*b[2]/3
    q=2*np.pi/L; hq=sum(np.exp(-1j*q*i)*hs[i] for i in range(L))
    moment=float((psi@comm(hq.conj().T,comm(H,hq))@psi).real/(2*L))
    poly=b[1]*(np.cos(q)-1)+b[2]*(np.cos(2*q)-1)
    assert abs(moment-poly)<1e-12
    assert abs(b[0]+2*b[1]+2*b[2])<1e-12
    rows.append({'J':J,'g':g,'L':L,'a2':a2,'a4':a4,'a4_over_gJ2':a4/(g*J*J),
                 'moment_polynomial_residual':abs(moment-poly),'M':float(np.sqrt(J*J+g*g))})
print(json.dumps(rows,indent=2))
