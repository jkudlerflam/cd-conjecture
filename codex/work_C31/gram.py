import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.integrate import quad_vec
from reduced import Reduced

def gram(n):
 r=Reduced(n);et,U=r.eigs['T'];eh,V=r.eigs['H'];p=r.ports[1]
 coeff=np.array([V[p,:,None]*V[a,:,None].conj()*U[a,None,:]*U[p,None,:].conj() for a in r.ds])
 de=et[None,:]-eh[:,None]
 def f(t):
  mn=np.minimum(eh[:,None],et[None,:]);ad=abs(de)
  kernel=np.exp(-mn*t)*t*np.where(ad*t<1e-7,1-ad*t/2+(ad*t)**2/6,-np.expm1(-ad*t)/np.where(ad==0,1,ad*t))
  F=np.sum(coeff*kernel,axis=(1,2)).real
  return np.outer(F,F)
 G,err=quad_vec(f,0,np.inf,epsabs=1e-35,epsrel=1e-10)
 T,M=r.leading();ratio=-T/G
 print(n,ratio, 'error',err,flush=True)
 return G
if __name__=='__main__':
 for n in [1,2,3]:gram(n)
