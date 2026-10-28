"""SU(N) symmetric-rank-k rational chain TBA.
Adapted, inside C14 only, from astra_ULS/tba_check.py. The source and
free-energy probe move from string 1 to string k. h has ferromagnetic
reference energy zero and two-site eigenvalues 2*(H_{k-j}-H_k).
python3 fusion_tba.py --N 8 --k 2 --T 0.04 0.02 0.01
"""
import numpy as np
from scipy.special import digamma, expit
from scipy.optimize import anderson
from scipy.fft import rfft, irfft, rfftfreq
import time, argparse

def ground_energy(N,k):
 j=np.arange(1,k+1,dtype=float)
 return -2/N*np.sum(digamma(1+(j-1)/N)-digamma(j/N))
class TBA:
 def __init__(self,N,k=1,M=32,K=1024,R=None):
  self.k=k
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
  self.drive[0,k-1]=-np.pi/np.cosh(np.pi*self.u)
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
  a=np.sum(self.rho*np.logaddexp(0,x[:,self.k-1]))*self.dx
  s=a+np.sum(self.rho*p[:,self.k-1]*z[:,self.k-1])*self.dx
  return dict(N=self.N,T=T,f=ground_energy(self.N,self.k)-T*a,s=s,c_entropy=6*s/(self.N*T), c_free=12*a/(self.N*T), k=self.k, M=self.M, K=self.K, R=self.R,residual=err,tangent_residual=np.max(abs(tangent(z.ravel()))),seconds=time.time()-t0),x

if __name__ == '__main__':
 import json
 from pathlib import Path
 ap=argparse.ArgumentParser()
 for name,default in [('N',2),('k',1),('M',24),('K',512)]: ap.add_argument('--'+name,type=int,default=default)
 ap.add_argument('--R',type=float,default=None)
 ap.add_argument('--T',type=float,nargs='+',default=[.1,.05,.02])
 ap.add_argument('--output',default='tba_results.json')
 a=ap.parse_args(); obj=TBA(a.N,a.k,a.M,a.K,a.R); x=None
 path=Path(a.output); results=json.loads(path.read_text()) if path.exists() else []
 for T in a.T:
  row,x=obj.solve(T,x)
  row['c_exact']=a.k*(a.N*a.N-1)/(a.N+a.k)
  row['norm']=float(sum(1/j for j in range(1,a.k+1)))
  row['velocity']=2*np.pi/a.N
  row['ratio_entropy']=row['c_entropy']/(a.N*row['norm'])
  row['ratio_free']=row['c_free']/(a.N*row['norm'])
  results.append(row); path.write_text(json.dumps(results,indent=2)); print(json.dumps(row),flush=True)
