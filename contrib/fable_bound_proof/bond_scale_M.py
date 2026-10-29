"""Gauge-optimized bond scale M = min_{o,b} ||h + o x 1 - 1 x o - b|| for XXZ (d=2) and the SU(3) swap (d=3).
The spectral norm is convex in (o,b), so a local minimum is global; we run several random starts of Nelder-Mead/Powell."""
import numpy as np
from scipy.optimize import minimize
rng=np.random.default_rng(0)
def herm_basis(d):
    B=[]
    for i in range(d):
        for j in range(i,d):
            if i==j:
                if i<d-1:
                    m=np.zeros((d,d),complex); m[i,i]=1; m[i+1,i+1]=-1; B.append(m)
            else:
                m=np.zeros((d,d),complex); m[i,j]=m[j,i]=1; B.append(m)
                m=np.zeros((d,d),complex); m[i,j]=-1j; m[j,i]=1j; B.append(m)
    return B
def M_opt(h,d):
    B=herm_basis(d); I=np.eye(d)
    def spread(x):
        o=sum(xi*b for xi,b in zip(x,B))
        g=h+np.kron(o,I)-np.kron(I,o)
        w=np.linalg.eigvalsh(g); return (w[-1]-w[0])/2
    best=spread(np.zeros(len(B)))
    for s in range(6):
        x0=rng.normal(size=len(B))*0.3
        r=minimize(spread,x0,method='Powell',options={'xtol':1e-10,'ftol':1e-12,'maxiter':20000})
        best=min(best,r.fun)
    return spread(np.zeros(len(B))),best
sx=np.array([[0,1],[1,0]])/2; sy=np.array([[0,-1j],[1j,0]])/2; sz=np.diag([.5,-.5])
for D in [0.0,0.5,1.0,-0.5]:
    h=np.kron(sx,sx)+np.kron(sy,sy)+D*np.kron(sz,sz)
    print(f"XXZ Delta={D}: half-spread without gauge {M_opt(h,2)[0]:.6f}, optimized M {M_opt(h,2)[1]:.6f}")
P=np.zeros((9,9))
for a in range(3):
    for b in range(3): P[3*a+b,3*b+a]=1
print("SU(3) swap: ", M_opt(P,3))
P2=np.zeros((4,4))
for a in range(2):
    for b in range(2): P2[2*a+b,2*b+a]=1
print("SU(2) swap: ", M_opt(P2,2))
