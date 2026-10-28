"""Local adaptation of compile model. Leading tau^2 coefficient, no scalar subtraction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import scipy.sparse.linalg as sla
from model import Chain,DEFAULT
from model2 import LaneChain,segment_parts,assemble_seg
from onebody import blocks
import time,json

def leading(n,p):
 t=time.time();ch=LaneChain('switch',n,2,K=2,max_data=3,one_per_cell=True)
 P=segment_parts(ch);p=dict(DEFAULT,**p)
 H,V=assemble_seg(P,dict(p,tauA=1.,tauB=1.))
 H0,_=assemble_seg(P,dict(p,tauA=0.,tauB=0.)); C=H-H0
 dm=ch.data_modes;vac=[i for i in ch.code if all((ch.states[i]>>q)&1==0 for q in dm)][0]
 cols=[vac];signs=[1]
 for q in dm:
  sg,new=Chain.apply([(q,1)],[],ch.states[vac]);cols.append(ch.index[new]);signs.append(sg)
 mask=np.ones(ch.dim,bool);mask[ch.code]=False; Q=np.where(mask)[0]
 B=V[Q][:,cols].toarray()*np.array(signs)[None,:]
 R=sla.splu(H0[Q][:,Q].tocsc());X=R.solve(B);Y=C[Q][:,Q]@X
 D=-Y.conj().T@R.solve(Y)
 h=D[1:,1:]-D[0,0]*np.eye(2*n)
 T,M=blocks(h,n)
 print('n dim seconds',n,ch.dim,time.time()-t,flush=True)
 return T,M
if __name__=='__main__':
 import sys
 n=int(sys.argv[1]);asym=float(sys.argv[2]) if len(sys.argv)>2 else .05;phi=float(sys.argv[3]) if len(sys.argv)>3 else .001
 T,M=leading(n,dict(wD=-.3,xL=1+asym,phiT=phi))
 print('T',T,'M',M,flush=True)
 np.savez(('full_target_n%d.npz' if len(sys.argv)>2 else 'full_leading_n%d.npz')%n,T=T,M=M)
