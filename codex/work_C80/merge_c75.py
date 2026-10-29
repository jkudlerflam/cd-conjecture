"""Test C75 equal-charge data and express identifications on the full legal image."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
import bond18,bond15
from checks import gaussian_residual

def test(n,track):
 p=bond18.params(n,D=4,Delta=.5,g=0,td=1e-4,K=0)
 old=bond15.LET;target=(track,'A',1);removed=(track,'B',0)
 f=lambda l:target if l==removed else l
 new=list(dict.fromkeys(map(f,old)));ix={l:i for i,l in enumerate(new)};d=len(new)
 h,_,o,v=bond15.parts(p);mats=[]
 for M,arity in [(h,2),(o,1),(v,2)]:
  Q=np.zeros((d**arity,d**arity),complex);seen={}
  def idx(k):return ix[f(old[k])] if arity==1 else ix[f(old[k//15])]*d+ix[f(old[k%15])]
  for j,i in zip(*np.nonzero(M)):
   key=idx(j),idx(i)
   if key in seen:assert abs(seen[key]-M[j,i])<1e-12
   seen[key]=M[j,i]
  for key,x in seen.items():Q[key]=x
  mats.append(Q)
 h,o,v=mats
 oldws=[tuple(bond15.merge(l) for l in w) for w in bond18.words(n)]
 ws=[tuple(ix[f(l)] for l in w) for w in oldws]
 base=[x+y for x in ws for y in ws];states=list(dict.fromkeys(base));ind={w:i for i,w in enumerate(states)}
 H=sp.lil_matrix((len(states),len(states)),dtype=complex);V=H.copy();fixed=ix[f(('e','',None))]
 for col,w in enumerate(states):
  for pos,l in enumerate(w):
   for row in np.flatnonzero(o[:,l]):
    z=list(w);z[pos]=row
    if tuple(z) in ind:H[ind[tuple(z)],col]+=o[row,l]
  full=(fixed,)+w
  for pos in range(2*(n+1)):
   ci=full[pos]*d+full[pos+1]
   for ri in np.flatnonzero(h[:,ci]):
    z=list(full);z[pos]=ri//d;z[pos+1]=ri%d
    if z[0]==fixed and tuple(z[1:]) in ind:H[ind[tuple(z[1:])],col]+=h[ri,ci]
  ci=w[n]*d+w[n+1]
  for ri in np.flatnonzero(v[:,ci]):
   z=list(w);z[n]=ri//d;z[n+1]=ri%d
   if tuple(z) in ind:V[ind[tuple(z)],col]+=v[ri,ci]
 H=H.tocsc();V=V.tocsc();ref=bond18.pair(dict(p,g=.01))
 G=np.zeros((len(states),4**n),complex)
 for i,w in enumerate(base):G[ind[w]]+=ref['G'][i]
 assert np.max(abs(G.conj().T@G-np.eye(4**n)))<1e-12
 ex=np.array([sum(new[l][1] in ['A','B'] for l in w) for w in states]);sel=np.flatnonzero(ex==2)
 E=2*ref['E0'];S=(V@G)[sel];A=H[sel][:,sel]-E*sp.eye(len(sel),format='csc')
 H2=-S.conj().T@splu(A.tocsc()).solve(S)
 return dict(n=n,merge_track=track,alphabet=d,legal_image_pair_dimension=len(states),code_residual=float(np.max(abs(H@G-E*G))),g_zero_non_scalar_norm=float(np.linalg.norm(H2-H2[0,0]*np.eye(4**n),2)),g_zero_non_gaussian_relative_residual=gaussian_residual(H2,n))
if __name__=='__main__':
 out=[test(n,t) for t in ['d','x'] for n in [1,2]]
 json.dump(out,open('c75_merge_checks.json','w'),indent=2);print(json.dumps(out),flush=True)
