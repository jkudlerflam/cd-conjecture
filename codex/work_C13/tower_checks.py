import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from pathlib import Path
out=[]
for common in [False,True]:
 for n in [1,2]:
  for M in [40,80,160,320]:
   b=4;r=round(.37*M);L=b*M+r;a=L/M
   Ns=np.rint(M*np.array([.27,.16][:n])).astype(int)
   theta=np.array([.31,-.57][:n]);ts=np.array([1,.73,.51][:n+1])
   beta=np.repeat((sum(Ns)-1)/2,n) if common else (Ns-1)/2
   offset=beta+(Ns-1)/2
   c=np.round(M*theta/(2*np.pi)-offset)+offset+np.arange(1,n+1)
   d=c-M*theta/(2*np.pi);D=-d;chi=np.dot(Ns,theta)/(2*np.pi)
   D0=int(round(chi-np.dot(Ns/M,D)))+1
   mpay=int(round(np.dot(Ns,c)))
   assert abs(mpay-np.dot(Ns,c))<1e-10
   j0=(D0*M-mpay)//L;ell=D0*M-mpay-L*j0
   assert 0<=ell<L
   assert (r*j0+mpay+ell)%M==0
   W=(D0-chi+np.dot(Ns/M,D))/a
   phases=np.r_[j0+ell/L,d]
   counts=np.r_[r,Ns];amp=2*ts*np.sin(np.pi*counts/M)
   exact=-sum(amp*np.cos(2*np.pi*phases/M)/np.sin(np.pi/M))
   bulk=-sum(amp)*M/np.pi;v=amp*a
   pred=bulk-np.pi*sum(v)/(6*L)+2*np.pi/L*(v[0]*W*W+np.dot(v[1:],D*D))
   out.append(dict(n=n,M=M,L=L,common_spin_structure=common,gluing_residual=int((r*j0+mpay+ell)%M),current_identity_error=float(abs(W-(j0+ell/L))),energy_error=float(exact-pred),L3_error=float((exact-pred)*L**3)))
Path('tower_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
