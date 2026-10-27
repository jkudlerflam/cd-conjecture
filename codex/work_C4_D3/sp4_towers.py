#!/usr/bin/env python3
"""Global Sp4 quadratic Casimir of low current-momentum states."""
import json,numpy as np
from scipy import sparse,linalg
from d3 import Sector

def generators():
 J=np.zeros((4,4))
 for i in range(4):J[i,3-i]=(-1)**i
 basis=[]
 for i in range(4):
  for j in range(i+1,4):
   t=np.zeros((4,4),complex);t[i,j]=t[j,i]=.5;basis.append(t)
   t=np.zeros((4,4),complex);t[i,j]=-.5j;t[j,i]=.5j;basis.append(t)
 for k in range(1,4):
  t=np.zeros((4,4),complex)
  for i in range(k):t[i,i]=1
  t[k,k]=-k;t/=np.sqrt(2*k*(k+1));basis.append(t)
 candidates=[(t+J@t.T@J)/2 for t in basis]
 # Gram-Schmidt preserves the real vector space of Hermitian matrices.
 gs=[]
 for t in candidates:
  for g in gs:t=t-2*np.trace(g@t).real*g
  norm=np.trace(t@t).real
  if norm>1e-10:gs.append(t/np.sqrt(2*norm))
 assert len(gs)==10
 assert max(np.linalg.norm(g@J+J@g.T) for g in gs)<1e-10
 assert np.linalg.norm(sum(g@g for g in gs)-1.25*np.eye(4))<1e-10
 return gs

def casimir(sec,psivecs):
 ans=np.zeros(psivecs.shape[1]);applied=np.zeros_like(psivecs);d=4
 for g in generators():
  rows=[];cols=[];data=[]
  for site in range(sec.L):
   pi=d**site;digit=sec.basis//pi%d
   for a in range(d):
    ix=np.flatnonzero(digit==a)
    for b in np.flatnonzero(abs(g[:,a])>1e-12):
     rows.append(sec.basis[ix]+(int(b)-a)*pi);cols.append(ix);data.append(np.full(len(ix),g[b,a]))
  targets=np.concatenate(rows);ci=np.concatenate(cols);values=np.concatenate(data)
  unique,ri=np.unique(targets,return_inverse=True)
  op=sparse.coo_matrix((values,(ri,ci)),shape=(len(unique),len(sec.basis))).tocsr()
  transformed=op@psivecs
  ans+=np.sum(abs(transformed)**2,axis=0)
  applied+=op.getH()@transformed
 variance=np.sum(abs(applied-psivecs*ans)**2,axis=0)
 return ans,variance

open('sp4_towers.jsonl','w').close()
for L in [4,8,12]:
 sec=Sector('sp4',L,cartan=True,cache_blocks=False)
 # Ground is k=0 at L=8, k=pi at L=4,12 on these commensurate lengths.
 k0=0 if L%8==0 else L//2;k=(k0+1)%L;k=min(k,L-k)
 for theta in [2*np.pi/3,np.pi]:
  coeff=[np.cos(theta),np.sin(theta)];e0=sec.solve(coeff,k0,1)[0][0];sec._blocks.clear()
  e,vec,res=sec.solve(coeff,k,18);sec._blocks.clear()
  psi=np.stack([sec.expand(k,vec[:,i]) for i in range(len(e))],axis=1)
  cc,variance=casimir(sec,psi);vv=(e-e0)*L/(2*np.pi)
  row={'L':L,'theta':theta,'k0':k0,'k':k,'E0':float(e0),'residual':res,'C2':cc.tolist(),'C2_variance':variance.tolist(),'gap_times_L_over_2pi':vv.tolist()}
  for val,name in [(3,'adjoint_current_velocity'),(2,'vector_mixed_current_velocity')]:
   matches=np.flatnonzero((abs(cc-val)<1e-6)&(variance<1e-8))
   row[name]=float(vv[matches[0]]) if len(matches) else None
  with open('sp4_towers.jsonl','a') as f:f.write(json.dumps(row)+'\n')
  print(json.dumps(row),flush=True)
