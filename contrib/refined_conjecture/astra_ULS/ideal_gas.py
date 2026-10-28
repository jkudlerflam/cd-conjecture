import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from scipy.optimize import brentq

def gas(t):
 def density(mu): return quad(lambda q:expit((mu-q*q)/t),0,max(20,np.sqrt(t)*12),epsabs=1e-10)[0]/np.pi
 mu=brentq(lambda mu:density(mu)-1,-max(100,t*20),20+t)
 def entropy(q):
  x=(mu-q*q)/t;p=expit(x)
  return np.logaddexp(0,x)-p*x
 s=quad(entropy,0,max(20,np.sqrt(t)*12),epsabs=1e-10)[0]/np.pi
 return mu,s,6*s/t
if __name__=='__main__':
 for t in [.5,1,2,3,4,5,6,8,10,20,100]:print(t,gas(t))
 root=brentq(lambda t:gas(t)[2]-.9,6,8)
 print('first_0.9',root,gas(root),'ell/N',2*np.pi/root)
