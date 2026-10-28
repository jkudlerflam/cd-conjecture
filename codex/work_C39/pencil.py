"""Gram pencil and its Fermi points, independently of finite-bus ED."""
import numpy as np
import scipy.linalg as la
from scipy.optimize import brentq
from itertools import combinations_with_replacement

def gram(n,w,a=4.):
 A=a*np.eye(n)-w*(np.eye(n,k=1)+np.eye(n,k=-1))
 E=np.zeros((n,n));E[0,0]=1
 AA=np.kron(np.eye(n),A)
 BB=np.block([[AA,-np.kron(np.eye(n),E)],[np.zeros_like(AA),AA]])
 v=np.r_[np.zeros(n*n),np.eye(n).reshape(-1)]
 P=la.solve_continuous_lyapunov(BB,np.outer(v,v))
 ix=np.arange(n)*(n+1)
 return P[np.ix_(ix,ix)]

def normalized(n,w):
 d=w**(2*np.arange(n));return gram(n,w)/np.outer(d,d)

def target(n,w=.25,v=.1,u=.1,eta=.05,mu=1e-10):
 A=normalized(n,w); dv=v**(2*np.arange(n))
 B=normalized(n,w*v)*np.outer(dv,dv)
 E=np.diag(u**(4*np.arange(n)))
 r,X=la.eigh(B,A);r=r[::-1];X=X[:,::-1]
 theta=-np.arctan(eta*r)
 slopes=np.einsum('ij,ik,kj->j',X,E,X)/(2*np.sqrt(1+(eta*r)**2))
 # Congruence by an A-whitener keeps scalar eigenvalues well resolved.
 L=la.cholesky(A,lower=True);Li=la.solve_triangular(L,np.eye(n),lower=True)
 BB=Li@B@Li.T;EE=Li@E@Li.T
 def h(k):return -2*np.eye(n)*np.cos(k)+2*eta*BB*np.sin(k)+mu*EE
 roots=[]
 for j,t in enumerate(theta):
  for sign in [-1,1]:
   p=t+sign*np.pi/2; ev=lambda k:la.eigvalsh(h(k))[j if sign==1 else n-1-j]
   # Sorting index follows descending pencil eigenvalues, positive sine.
   width=min(.1, max(mu*max(slopes)*4,1e-9))
   # use determinant sign near a predicted isolated root
   eig=la.eigvalsh(h(p));idx=int(np.argmin(abs(eig)))
   ev=lambda k:la.eigvalsh(h(k))[idx]
   roots.append(brentq(ev,p-width,p+width,xtol=1e-15))
 pairs=list(combinations_with_replacement(range(2*n),2))
 circ=lambda x:abs((x+np.pi)%(2*np.pi)-np.pi)
 d2=min((circ(roots[i]-roots[j]) for i in range(2*n) for j in range(i)),default=np.inf)
 d4=min((circ(roots[i]+roots[j]-roots[k]-roots[l]) for z,(i,j) in enumerate(pairs) for k,l in pairs[:z]),default=np.inf)
 return dict(r=r.tolist(),theta=theta.tolist(),slopes=slopes.tolist(),roots=roots,delta2=d2,delta4=d4)

if __name__=='__main__':
 import json
 for n in range(1,4): print(n,json.dumps(target(n,v=.3,u=.3,mu=1e-7),indent=2))
