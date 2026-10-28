import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,itertools as it,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from test_fluids import matrix

def prime(q):return q>=2 and all(q%d for d in range(2,int(q**.5)+1))
def distance(x):return np.abs((np.asarray(x)+np.pi)%(2*np.pi)-np.pi)
rows=[]
for n in [1,2,4,8]:
 N=2*n;q=next(q for q in range(max(3,N),16*N+1) if prime(q))
 A=np.array([2*q*i+i*i%q for i in range(N)])
 B=int(np.ceil(np.min(np.diff(A))/2));T=2*A[-1]+4*B+2
 p=2*np.pi*(A+B+1/6)/T+1e-10*np.sqrt(np.arange(N)+2)
 bil=np.array([s*p[i]+t*p[j] for i,j in it.combinations(range(N),2) for s,t in it.product([-1,1],repeat=2)])
 d2=float(min(np.min(distance(2*p)),np.min(distance(bil))))
 d4=min((float(distance(sum(s*p[i] for i,s in zip(ii,ss)))) for ii in it.combinations(range(N),4) for ss in it.product([-1,1],repeat=4)),default=None)
 roots=np.unique((bil/(2*np.pi))%1);cuts=np.r_[.25,roots[(roots>.25)&(roots<.75)],.75];j=np.argmax(np.diff(cuts));nu=(cuts[j]+cuts[j+1])/2+1e-10*np.sqrt(7)
 dm=float(min(np.min(distance(bil+2*np.pi*nu)),np.min(distance(bil-2*np.pi*nu))))
 rows.append(dict(n=n,p=p.tolist(),nu=float(nu),delta2=d2,delta4=d4,delta_mixed=dm,conservative_mixed_bound=float(np.pi/(4*(4*(N*(N-1)//2)+1)))))
# Canonical transformation check: charge stiffness times current stiffness.
a=4.37;f=np.array([.27,.19]);v=np.array([5.1,2.3,1.7])
T=np.eye(3);T[0,0]=a;T[1:,0]=-f
A=T.T@np.diag(v)@T;Ti=np.linalg.inv(T);B=Ti@np.diag(v)@Ti.T
canon=dict(velocities_from_stiffness=np.sqrt(np.sort(np.linalg.eigvals(A@B))).tolist(),target=sorted(v.tolist()))
# A parity-even uniform p-wave pair gaps an unshifted incommensurate band.
kF=np.pi*np.sqrt(2)/4;Delta=.02;k=np.linspace(-np.pi,np.pi,200001)
E=np.sqrt((-2*np.cos(k)+2*np.cos(kF))**2+(2*Delta*np.sin(k))**2)
pair=dict(kF=float(kF),Delta=Delta,grid_size=len(k),min_gap=float(np.min(E)))
Path('additional_results.json').write_text(json.dumps(dict(margins=rows,canonical=canon,pairing=pair),indent=2))
print('MARGINS',json.dumps(rows),flush=True)
# Contiguous physical lengths, fixed charge, full graph ED.
rings=[]
for n in [1,2]:
 for L in range(14,20):
  M=round(L/2.4);r=L-2*M;Ns=[2,1][:n]
  ts=time.time();states,H,psi,phase,hops=matrix(2,M,r,Ns,True)
  e=eigsh(H,k=2,which='SA',tol=1e-10,v0=np.random.default_rng(13013).normal(size=len(psi)),return_eigenvectors=False);e.sort()
  row=dict(n=n,L=L,M=M,r=r,Ns=Ns,dim=len(states),gap=float(e[1]-e[0]),Lgap=float(L*(e[1]-e[0])),seconds=time.time()-ts)
  rings.append(row);print(json.dumps(row),flush=True)
  Path('contiguous_rings.json').write_text(json.dumps(rings,indent=2))
