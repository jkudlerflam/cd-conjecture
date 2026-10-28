"""Exact coefficient of lambda^2 tau^2 using carrier-pair Sylvester solves.
Physical fermion signs are evaluated in a literal occupation basis.
No many-body basis or small-tau subtraction is used.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy.linalg import eigh
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from functools import lru_cache
from model import Chain

class Reduced:
 def __init__(self,n,xL=1.,phi=.0,w=.3,method="positive"):
  self.method=method
  self.n=n;self.S=3*n+3;self.xL=xL;self.phi=phi
  # Physical order differs from upstream only by putting data after dimers.
  self.data=list(range(4*self.S,4*self.S+2*n))
  self.vac=sum(1<<(2*i) for i in range(2*self.S))
  self.ports=[2*n+2,0]
  self.ds=list(range(2*n+3,3*n+3))
  edges=[]
  for i in range(2*n+2):edges.append((i,i+1,xL if i<n+1 else 1.))
  edges.append((n+1,2*n+3,1.))
  edges.extend((i,i+1,1.) for i in range(2*n+3,3*n+2))
  depth=np.array(list(range(2*n+3))+list(range(n+2,2*n+2)))
  self.gauge={'T':np.exp(1j*phi*depth),'H':(-1.)**depth,'D':(-1.)**depth}
  K=np.eye(self.S)
  for i,j,x in edges:K[i,j]=K[j,i]=-w*x
  self.positive_lu=splu((sp.kron(sp.csc_matrix(K),sp.eye(self.S))+sp.kron(sp.eye(self.S),sp.csc_matrix(K))).tocsc())
  self.eigs={}
  for X in 'THD':
   H=np.eye(self.S,dtype=complex)
   for i,j,x in edges:
    v=(-w*np.exp(1j*phi) if X=='T' else w)*x
    H[j,i]=v;H[i,j]=v.conjugate()
   self.eigs[X]=eigh(H)
  self.gap=min(e[0][0] for e in self.eigs.values())
 @lru_cache(None)
 def code(self,bits):
  return Chain.apply([(q,1) for a,q in enumerate(self.data) if bits>>a&1],[],self.vac)
 @lru_cache(None)
 def state(self,X,Y,bits,r,s):
  ops=[]
  for cell,pos,t in [(0,r,X),(1,s,Y)]:
   b=2*(cell*self.S+pos);a=b+1
   ops+= {'T':[(a,1),(b,0)],'H':[(b,0)],'D':[(a,1)]}[t]
  sg,st=self.code(bits);sg2,st2=Chain.apply(ops,[],st)
  return sg*sg2,st2
 def resolve(self,blocks):
  out={}
  for (X,Y,bits),B in blocks.items():
   if self.method=='positive':
    dx=self.gauge[X];dy=self.gauge[Y];src=dx.conj()[:,None]*B*dy.conj()[None,:]
    flat=src.ravel();sol=self.positive_lu.solve(flat.real)+1j*self.positive_lu.solve(flat.imag)
    out[X,Y,bits]=dx[:,None]*sol.reshape(self.S,self.S)*dy[None,:]
   else:
    ex,U=self.eigs[X];ey,W=self.eigs[Y]
    out[X,Y,bits]=U@((U.conj().T@B@W.conj())/(ex[:,None]+ey[None,:]))@W.T
  return out
 def source(self,bits):
  bR=2*self.ports[0];bL=2*(self.S+self.ports[1]);qR=bR+1;qL=bL+1
  terms=[]
  for a in [bR,qR]:
   for b in [bL,qL]:
    terms.extend([(.5,[(a,1),(b,0)]),(.5,[(b,1),(a,0)])])
  terms.extend([(.25,[(qR,1),(bR,0),(qL,1),(bL,0)])])
  sg,st=self.code(bits);out={}
  for coef,ops in terms:
   result=Chain.apply(ops,[],st)
   if result is None:continue
   sgn,st1=result
   for X,Y in [('T','T'),('D','H'),('H','D')]:
    gs,ss=self.state(X,Y,bits,*self.ports)
    if ss==st1:
     key=X,Y,bits
     if key not in out:out[key]=np.zeros((self.S,self.S),complex)
     out[key][tuple(self.ports)]+=coef*sg*sgn*gs
  return out
 def convert(self,blocks):
  out={}
  for (X,Y,bits),B in blocks.items():
   for cell in [0,1]:
    typ=[X,Y][cell]
    for a,pos in enumerate(self.ds):
     di=cell*self.n+a;c=self.data[di];b=2*(cell*self.S+pos);al=b+1
     opss=[[(c,1),(al,0)],[(al,1),(c,0)],[(b,1),(c,0)],[(c,1),(b,0)]]
     conds=[[(b,0)],[(b,0)],[(al,1)],[(al,1)]]
     for other in range(self.S):
      r,s=(pos,other) if cell==0 else (other,pos)
      if abs(B[r,s])==0:continue
      sg,st=self.state(X,Y,bits,r,s)
      for ops,co in zip(opss,conds):
       res=Chain.apply(ops,co,st)
       if res is None:continue
       sg1,st1=res;nbits=bits^(1<<di)
       occ=((st1>>b)&1,(st1>>al)&1)
       nt={(0,1):'T',(0,0):'H',(1,1):'D'}.get(occ)
       if nt is None:raise ValueError('carrier lost')
       nx,ny=(nt,Y) if cell==0 else (X,nt)
       gs,ss=self.state(nx,ny,nbits,r,s)
       assert ss==st1
       key=nx,ny,nbits
       if key not in out:out[key]=np.zeros_like(B)
       out[key][r,s]+=-sg*sg1*gs*B[r,s]
  return out
 def junction(self,blocks):
  # Q V Q within one carrier per cell. Other carrier-number sectors
  # have zero overlap with the endpoint resolvent vectors at order lambda^3.
  bR=2*self.ports[0];bL=2*(self.S+self.ports[1]);qR=bR+1;qL=bL+1
  terms=[]
  for a in [bR,qR]:
   for b in [bL,qL]:terms.extend([(.5,[(a,1),(b,0)]),(.5,[(b,1),(a,0)])])
  def dens(b,a):return [(.5,[(b,1),(b,0)]),(.5,[(a,1),(a,0)]),(.5,[(b,1),(a,0)]),(.5,[(a,1),(b,0)]),(-.5,[])]
  terms.extend((cr*cl,op+ol) for cr,op in dens(bR,qR) for cl,ol in dens(bL,qL))
  out={}
  # The projected port operator is zero when a carrier is away from
  # its port, since p and q have zero code compression.
  for (X,Y,bits),B in blocks.items():
   r,s=self.ports;val=B[r,s]
   if val==0:continue
   sg,st=self.state(X,Y,bits,r,s)
   targets={self.state(nx,ny,bits,r,s)[1]:(nx,ny,self.state(nx,ny,bits,r,s)[0]) for nx in 'THD' for ny in 'THD'}
   for coef,ops in terms:
    res=Chain.apply(ops,[],st)
    if res is None:continue
    sg1,st1=res
    if st1 not in targets:continue
    nx,ny,gs=targets[st1];key=nx,ny,bits
    if key not in out:out[key]=np.zeros_like(B)
    out[key][r,s]+=coef*sg*sg1*gs*val
  return out
 def leading(self):
  Ys=[self.convert(self.resolve(self.source(b))) for b in [0]+[1<<a for a in range(2*self.n)]]
  RYs=[self.resolve(Y) for Y in Ys]
  D=np.array([[-sum(np.vdot(A,Y[k]) for k,A in X.items() if k in Y) for Y in RYs] for X in Ys])
  assert np.max(abs(D-D.conj().T))<1e-15
  h=D[1:,1:]-D[0,0]*np.eye(2*self.n)
  n=self.n
  return h[:n,n:],h[:n,:n]+h[n:,n:]

if __name__=='__main__':
 import time
 t=__import__('time').time();r=Reduced(2,1.05,.001);T,M=r.leading()
 print('seconds',time.time()-t,'T',T,'M',M)
 ref=np.load('full_leading_n2.npz')
 # Gauge from different placement of physical dimers changes c_a signs.
 print('singular values',np.linalg.svd(T,compute_uv=False),np.linalg.svd(ref['T'],compute_uv=False))
 print('relative T',np.max(abs(T-ref['T']))/np.max(abs(T)))
 print('relative M',np.max(abs(M-ref['M']))/np.max(abs(M)))
