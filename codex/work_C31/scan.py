import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.linalg as la
from itertools import combinations
from reduced import Reduced

def roots(T,M):
 n=len(T);A=-(T+T.conj().T)/2
 e,U=la.eigh(A)
 if e[0]<=0:return None
 W=(U/np.sqrt(e))@U.conj().T
 t=W@T@W;m=W@M@W
 O=np.zeros((n,n));I=np.eye(n)
 z=la.eigvals(np.block([[-m,-t.conj().T],[I,O]]),np.block([[t,O],[O,I]]))
 p=np.sort(np.angle(z[abs(abs(z)-1)<1e-6]))
 if len(p)!=2*n:return dict(count=len(p),unit_defect=float(max(abs(abs(z)-1))))
 vel=[];labels=[];res=[]
 for k in p:
  h=m+t*np.exp(1j*k)+t.conj().T*np.exp(-1j*k)
  ev,uv=la.eigh(h);ix=np.argmin(abs(ev));v=W@uv[:,ix];v=v/la.norm(v)
  vel.append(float(np.vdot(v,(1j*T*np.exp(1j*k)-1j*T.conj().T*np.exp(-1j*k))@v).real))
  labels.append(int(ix));res.append(float(abs(ev[ix])))
 def circ(x):return abs(np.angle(np.exp(1j*x)))
 d2=min(circ(a-b) for a,b in combinations(p,2))
 pairlist=list(combinations(range(2*n),2))
 d4=min(circ(sum(p[list(a)])-sum(p[list(b)])) for a,b in combinations(pairlist,2))
 centres=[];pairs=[]
 for a in range(n):
  pp=[p[i] for i in range(2*n) if labels[i]==a]
  if len(pp)==2:
   cen=np.angle(np.exp(1j*sum(pp)))/2
   centres.append(float(cen));pairs.append(pp)
 dp=min(circ(2*(a-b)) for a,b in combinations(centres,2)) if len(centres)>1 else 0
 return dict(count=len(p),p=p.tolist(),v=vel,centres=centres,delta2=d2,delta4=d4,pair_margin=dp,residual=max(res),unit_defect=float(max(abs(abs(z)-1))))

if __name__=='__main__':
 import json,time
 out=[]
 for n in [2,3,4]:
  T0,M0=Reduced(n).leading();s=la.eigvalsh(-T0)
  print('base',n,s,'Mnorm',la.norm(M0),flush=True)
  for asym,phi in [(.05,.001),(.005,.0001),(.001,.00001),(.0001,.000001),(.00001,.0000001)]:
   T,M=Reduced(n,1+asym,phi).leading()
   for mf in [0.,.3,1.]:
    mu=mf*s[0];rr=roots(T,M+mu*np.eye(n))
    row=dict(n=n,asym=asym,phi=phi,mu=mu,smin=s[0],smax=s[-1],**rr)
    out.append(row)
    print(json.dumps(row),flush=True)
 json.dump(out,open('scan_results.json','w'),indent=2)
