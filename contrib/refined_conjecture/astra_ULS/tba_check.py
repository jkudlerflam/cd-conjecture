import numpy as np
from scipy.special import digamma, expit
from scipy.optimize import anderson
from scipy.fft import rfft, irfft, rfftfreq
import time, argparse

def e0(N): return 1+2/N*(digamma(1/N)+np.euler_gamma)
class TBA:
 def __init__(self,N,M=32,K=1024,R=None):
  self.N,self.M,self.K=N,M,K
  self.R=R or max(16,4*N)
  self.u=(np.arange(K)-K//2)*(2*self.R/K)
  self.dx=2*self.R/K
  self.sh=1/(2*np.cosh(np.pi*rfftfreq(K,self.dx)))
  r=np.arange(1,N)[:,None,None]; n=np.arange(1,M+2)[None,:,None]
  self.xinf=np.log(n*(n+N)/(r*(N-r)))
  self.xbase=np.broadcast_to(self.xinf[:,:M],(N-1,M,K)).copy()
  self.Linf=np.logaddexp(0,self.xinf)
  self.drive=np.zeros_like(self.xbase)
  self.drive[0,0]=-np.pi/np.cosh(np.pi*self.u)
  rr=np.arange(1,N)[:,None]
  self.rho=np.sin(np.pi*rr/N)/(N*(np.cosh(2*np.pi*self.u/N)-np.cos(np.pi*rr/N)))
 def conv(self,z): return irfft(rfft(z,axis=-1)*self.sh,n=self.K,axis=-1)
 def mapping(self,x,T):
  L=np.logaddexp(0,x); l=L-x
  b=np.zeros_like(x)
  b[:,1:]+=L[:,:-1]; b[:,:-1]+=L[:,1:]; b[:,-1:]+=self.Linf[:,-1:]
  b[1:]-=l[:-1]; b[:-1]-=l[1:]
  return self.drive/T+self.conv(b)
 def solve(self,T,x=None,tol=1e-9):
  if x is None: x=self.xbase+self.drive/T
  shape=x.shape
  def res(z):
   xx=z.reshape(shape); return (self.mapping(xx,T)-xx).ravel()
  t0=time.time()
  try: sol=anderson(res,x.ravel(),M=8,alpha=.7,w0=.01,f_tol=tol,maxiter=1200)
  except Exception as ex:
   if isinstance(ex.args[0],np.ndarray): sol=ex.args[0]
   else: raise
  x=sol.reshape(shape)
  err=np.max(abs(res(sol)))
  # T derivative z = d x / d log(T)
  p=expit(x)
  def tangent(z):
   zz=z.reshape(shape); b=np.zeros_like(zz)
   b[:,1:]+=p[:,:-1]*zz[:,:-1]; b[:,:-1]+=p[:,1:]*zz[:,1:]
   b[1:]+=(1-p[:-1])*zz[:-1]; b[:-1]+=(1-p[1:])*zz[1:]
   return (-self.drive/T+self.conv(b)-zz).ravel()
  try: z=anderson(tangent,(-self.drive/T).ravel(),M=8,alpha=.7,w0=.01,f_tol=tol,maxiter=1200).reshape(shape)
  except Exception as ex: z=ex.args[0].reshape(shape)
  a=np.sum(self.rho*np.logaddexp(0,x[:,0]))*self.dx
  s=a+np.sum(self.rho*p[:,0]*z[:,0])*self.dx
  return dict(N=self.N,T=T,f=e0(self.N)-T*a,s=s,ratio=s/(self.N*(self.N-1)*T/6),residual=err,tangent_residual=np.max(abs(tangent(z.ravel()))),seconds=time.time()-t0),x
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=2);ap.add_argument('--M',type=int,default=24);ap.add_argument('--K',type=int,default=512);ap.add_argument('--R',type=float,default=None);ap.add_argument('--T',type=float,nargs='+',default=[.2,.4,.6,1.,3.]);a=ap.parse_args()
 obj=TBA(a.N,a.M,a.K,a.R); x=None
 for T in a.T:
  out,x=obj.solve(T,x);print(out,flush=True)
