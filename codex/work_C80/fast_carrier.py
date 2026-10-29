"""Different fold attempt: replace the B express rail by a fast direct Laplacian."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from structured import carrier,visibility

class FastB:
 def __init__(self,n,speed):
  self.n=n;self.S=2*n+1;S=self.S
  A=carrier(n,D=10,a=4,td=1e-4,flux=0,loopflux=np.pi/2)[0]
  B=4*np.eye(S,dtype=complex);order=list(range(n))+[2*n]
  for i,j in zip(order,order[1:]):
   B[i,i]+=speed;B[j,j]+=speed;B[i,j]-=speed;B[j,i]-=speed
  H=[A,sp.csc_matrix(B)]
  self.lu={(s,t):splu((sp.kron(H[s],sp.eye(S))+sp.kron(sp.eye(S),H[t])).tocsc()) for s in range(2) for t in range(2)}
  source=np.zeros(S*S);source[(S-1)*S]=1
  self.rs={st:self.lu[st].solve(source).reshape(S,S) for st in [(0,1),(1,0)]}
 def Y(self,bits):
  n=self.n;out={};theta=-np.pi/4
  for s,t,c in [(0,1,-1j*np.exp(1j*theta)*(-1)**((bits&((1<<n)-1)).bit_count())),(1,0,np.exp(-1j*theta))]:
   R=c*self.rs[s,t]
   for cell,typ in enumerate([s,t]):
    local=(bits>>(cell*n))&((1<<n)-1)
    for j in range(n):
     bit=(local>>j)&1
     if (typ==0 and bit!=1) or (typ==1 and bit!=0):continue
     sign=(-1)**((local&((1<<j)-1)).bit_count());nb=bits^(1<<(cell*n+j))
     key=(1-typ,t,nb) if cell==0 else (s,1-typ,nb)
     if key not in out:out[key]=np.zeros((self.S,self.S),complex)
     if cell==0:out[key][j,:]+=sign*R[j,:]
     else:out[key][:,j]+=sign*R[:,j]
  return out
 def calculate(self):
  n=self.n;ys=[self.Y(v) for v in [0]+[1<<i for i in range(2*n)]];H=np.zeros((len(ys),len(ys)),complex)
  for j,Y in enumerate(ys):
   for key,v in Y.items():
    rv=self.lu[key[:2]].solve(v.ravel())
    for i,X in enumerate(ys):
     if key in X:H[i,j]-=np.vdot(X[key].ravel(),rv)
  h=H[1:,1:]-H[0,0]*np.eye(2*n)
  return visibility(n)*h[:n,n:],visibility(n)*(h[:n,:n]+h[n:,n:])

if __name__=='__main__':
 out=[]
 for n in [2,3]:
  for speed in [1,10,100,1000,10000]:
   T,M=FastB(n,speed).calculate()
   es=np.array([np.linalg.eigvalsh(M+T*np.exp(1j*k)+T.conj().T*np.exp(-1j*k)) for k in np.linspace(-np.pi,np.pi,2049)])
   lo=es.min(axis=0);hi=es.max(axis=0)
   overlap=float(min(hi)-max(lo));err=2*np.linalg.norm(T,2)*np.pi/2048
   out.append(dict(n=n,speed=speed,common_field_overlap=overlap,grid_Lipschitz_error=float(err),onsite_to_hopping_norm=float(np.linalg.norm(M,2)/np.linalg.norm(T,2)),intervals=[[float(a),float(c)] for a,c in zip(lo,hi)]))
 json.dump(out,open('fast_carrier_checks.json','w'),indent=2);print(json.dumps(out),flush=True)
