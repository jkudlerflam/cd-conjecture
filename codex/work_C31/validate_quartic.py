import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.sparse.linalg as sla
from model import Chain,DEFAULT
from model2 import LaneChain,segment_parts,assemble_seg
from centered_scan import matrices
import time,json
n=2;t=time.time();ch=LaneChain('switch',n,2,K=2,max_data=None,one_per_cell=True)
P=segment_parts(ch);p=dict(DEFAULT,wD=-.3,tauA=0.,tauB=0.)
H0,V=assemble_seg(P,p);H,_=assemble_seg(P,dict(p,tauA=1.,tauB=1.));C=H-H0
vac=[i for i in ch.code if all((ch.states[i]>>q)&1==0 for q in ch.data_modes)][0]
from itertools import combinations
pairs=list(combinations(range(2*n),2));bits=[0]+[1<<i for i in range(2*n)]+[(1<<i)+(1<<j) for i,j in pairs]
cols=[];signs=[]
for b in bits:
 sg,st=Chain.apply([(q,1) for a,q in enumerate(ch.data_modes) if b>>a&1],[],ch.states[vac]);cols.append(ch.index[st]);signs.append(sg)
mask=np.ones(ch.dim,bool);mask[ch.code]=False;Q=np.where(mask)[0]
B=V[Q][:,cols].toarray()*np.array(signs)[None,:];CQ=C[Q][:,Q];VQ=V[Q][:,Q];R=sla.splu(H0[Q][:,Q].tocsc())
Xs=[R.solve(B)]
for q in range(4):Xs.append(R.solve(CQ@Xs[-1]))
K2=-(CQ@Xs[1]).conj().T@Xs[2]
K3=sum(Xs[q].conj().T@(VQ@Xs[4-q]) for q in range(5))
a,b,_,ref2,ref3=matrices(n,True)
row=dict(n=n,dim=ch.dim,code_columns=len(cols),error_K2=float(np.linalg.norm(K2-ref2,2)),relative_K2=float(np.linalg.norm(K2-ref2,2)/np.linalg.norm(K2,2)),error_K3=float(np.linalg.norm(K3-ref3,2)),relative_K3=float(np.linalg.norm(K3-ref3,2)/np.linalg.norm(K3,2)),seconds=time.time()-t)
print(json.dumps(row));json.dump(row,open('quartic_validation.json','w'),indent=2)
