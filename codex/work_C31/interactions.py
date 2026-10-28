"""Leading non-Gaussian coefficient of the second-order junction, tau^4."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.linalg as la
from itertools import combinations
from reduced import Reduced
from model import Chain

def inner(A,B):return sum(np.vdot(X,B[k]) for k,X in A.items() if k in B)
def coeff(n,xL=1.,phi=0.):
 r=Reduced(n,xL,phi);Ys=[]
 for b in range(1<<(2*n)):
  X=r.resolve(r.source(b))
  Y=r.convert(X)
  Y2=r.convert(r.resolve(Y));Ys.append(Y2)
 RYs=[r.resolve(Y) for Y in Ys]
 K=-np.array([[inner(A,B) for B in RYs] for A in Ys])
 return K

def normal_order(K,m):
 vac=K[0,0];h=K[np.ix_([1<<i for i in range(m)],[1<<i for i in range(m)])]-vac*np.eye(m)
 G=np.eye(1<<m,dtype=complex)*vac
 for st in range(1<<m):
  for i in range(m):
   for j in range(m):
    z=Chain.apply([(i,1),(j,0)],[],st)
    if z is not None:sg,new=z;G[new,st]+=h[i,j]*sg
 V=K-G;pairs=list(combinations(range(m),2));ids=[(1<<i)+(1<<j) for i,j in pairs]
 return V[np.ix_(ids,ids)],h,V,pairs

def vertex(V,pairs,ua,ub,uc,ud):
 left=np.array([ua[i]*ub[j]-ub[i]*ua[j] for i,j in pairs])
 right=np.array([uc[i]*ud[j]-ud[i]*uc[j] for i,j in pairs])
 return np.vdot(left,V@right)

if __name__=='__main__':
 import json,time
 n=2;t=__import__('time').time();T,M=Reduced(n).leading();s,U=la.eigh(-T)
 K=coeff(n);V,h,NG,pairs=normal_order(K,2*n)
 out=[]
 for mf in [.1,.3,1.]:
  mu=mf*s[0];p=np.arccos(mu/(2*s));v=2*s*np.sin(p)
  states={}
  for a in range(n):
   for chir in [-1,1]:states[a,chir]=np.concatenate([U[:,a],np.exp(1j*chir*p[a])*U[:,a]])
  g=np.array([[vertex(V,pairs,states[a,1],states[b,-1],states[a,1],states[b,-1]).real for b in range(n)] for a in range(n)])
  J=vertex(V,pairs,states[0,1],states[0,-1],states[1,1],states[1,-1])
  dKi=-g/(np.pi*(v[:,None]+v[None,:]));l=np.array([1.,-1.]);x=-l@dKi@l
  row=dict(mu_factor=mf,p=p.tolist(),v=v.tolist(),g2=g.tolist(),J_real=J.real,J_imag=J.imag,x_per_tau2=x,J_over_pi_sqrtvv_per_tau2=abs(J)/(np.pi*np.sqrt(v.prod())))
  out.append(row);print(json.dumps(row))
 print('time',time.time()-t,'nonGaussian norm',la.norm(NG,2),'quartic norm',la.norm(V,2))
 np.savez('centered_quartic_n2.npz',K4tau=K,V=V,h=h,T=T)
 json.dump(dict(rows=out,nonGaussian_norm=la.norm(NG,2),seconds=time.time()-t),open('centered_results.json','w'),indent=2)
