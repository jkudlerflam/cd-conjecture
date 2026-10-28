import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.linalg as la
from itertools import combinations
from reduced import Reduced
from interactions import inner,vertex
import json,time

def matrices(n,third=True):
 r=Reduced(n);m=2*n;pairs=list(combinations(range(m),2));bits=[0]+[1<<i for i in range(m)]+[(1<<i)+(1<<j) for i,j in pairs]
 seq=[];ys=[]
 for b in bits:
  X=[r.resolve(r.source(b))]
  for q in range(4 if third else 2): X.append(r.resolve(r.convert(X[-1])))
  seq.append(X);ys.append(r.convert(X[1]))
 ry=[r.resolve(y) for y in ys]
 K2=-np.array([[inner(a,b) for b in ry] for a in ys])
 K3=np.zeros_like(K2)
 if third:
  for q in range(5):
   vs=[r.junction(x[4-q]) for x in seq]
   K3+=np.array([[inner(x[q],v) for v in vs] for x in seq])
 def subtract(K):
  h=K[1:1+m,1:1+m]-K[0,0]*np.eye(m)
  V=K[1+m:,1+m:].copy()-K[0,0]*np.eye(len(pairs))
  for j,(c,d) in enumerate(pairs):
   for i,(a,b) in enumerate(pairs):
    V[i,j]-=h[a,c]*(b==d)-h[a,d]*(b==c)-h[b,c]*(a==d)+h[b,d]*(a==c)
  return V,h
 V2,h2=subtract(K2);V3,h3=subtract(K3)
 return V2,V3,pairs,K2,K3

if __name__=='__main__':
 out=[]
 for n in [2,3,4]:
  t=time.time();V2,V3,pairs,K2,K3=matrices(n,True)
  T,M=Reduced(n).leading();s,U=la.eigh(-(T.real+T.real.T)/2)
  mu=.3*s[0];p=np.arccos(mu/(2*s));v=2*s*np.sin(p)
  states={(a,ch):np.concatenate([U[:,a],np.exp(1j*ch*p[a])*U[:,a]]) for a in range(n) for ch in [-1,1]}
  rows=[]
  for order,V in [(2,V2),(3,V3)]:
   g=np.array([[vertex(V,pairs,states[a,1],states[b,-1],states[a,1],states[b,-1]).real for b in range(n)] for a in range(n)])
   dK=-g/(np.pi*(v[:,None]+v[None,:]))
   vals=[]
   for a,b in combinations(range(n),2):
    l=np.eye(n)[a]-np.eye(n)[b];x=-l@dK@l
    J=vertex(V,pairs,states[a,1],states[a,-1],states[b,1],states[b,-1])
    vals.append(dict(a=a,b=b,x_per_tau2=x,J_real=J.real,J_imag=J.imag,y=abs(J)/(np.pi*np.sqrt(v[a]*v[b]))))
   rows.append(dict(lambda_order=order,g2=g.tolist(),pairs=vals))
  row=dict(n=n,p=p.tolist(),v=v.tolist(),rows=rows,seconds=time.time()-t)
  out.append(row);print(json.dumps(row),flush=True)
  np.savez('quartic_n%d.npz'%n,V2=V2,V3=V3,K2=K2,K3=K3)
 json.dump(out,open('centered_scan_results.json','w'),indent=2)
