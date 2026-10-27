import numpy as np, scipy.linalg as la
from scipy.optimize import least_squares
exec(open('lm_fp.py').read().split('out=[]')[0])
# Hermitian basis (16 real params)
B=[]
for i in range(4):
    E=np.zeros((4,4),complex); E[i,i]=1; B.append(E)
for i in range(4):
    for j in range(i+1,4):
        E=np.zeros((4,4),complex); E[i,j]=E[j,i]=1; B.append(E)
        E=np.zeros((4,4),complex); E[i,j]=-1j; E[j,i]=1j; B.append(E)
def sig(p, real=False):
    G=sum(pi*Bi for pi,Bi in zip(p,B))
    if real: G=G.real
    w,v=la.eigh(G); E=(v*np.exp(w-w.max()))@v.conj().T; return E/np.trace(E).real
def resid(p,real=False):
    s=sig(p,real); r3,a,b,lam=step(s)
    d=np.concatenate([(a-s).ravel(),(b-s).ravel()])
    return np.concatenate([d.real,d.imag])
rng=np.random.default_rng(7); sols=[]
for t in range(300):
    scale=[0.3,1,3][t%3]
    p0=rng.normal(size=16)*scale
    r=least_squares(resid,p0,bounds=(-12,12),xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    s=sig(r.x); ev=np.sort(la.eigvalsh(s))
    if np.abs(r.fun).max()<1e-10:
        lam=step(s)[3]; sols.append((3*lam/h14,ev,np.abs(r.fun).max()))
print("roots found:",len(sols),"of 300")
if sols:
    cs=np.array([x[0] for x in sols]); print("c_Delta values:",np.unique(np.round(cs,6)))
    print("min eig over roots:",min(x[1][0] for x in sols))
    for x in sols[:6]: print(round(x[0],8),np.round(x[1],8),x[2])
# Jacobian at tracial point: rank of linearized residual (isolated up to symmetries?)
p0=np.zeros(16)
J=np.array([(resid(p0+1e-6*e)-resid(p0-1e-6*e))/2e-6 for e in np.eye(16)]).T
sv=la.svdvals(J); print("Jacobian singular values at tracial:",np.round(sv,6))
