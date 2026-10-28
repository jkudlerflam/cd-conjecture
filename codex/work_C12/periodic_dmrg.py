"""Periodic Hamiltonian with an open MPS, retaining the closing bond explicitly."""
import numpy as np
from canonical_dmrg import CanonicalDMRG
from simple_dmrg import qutrit_model

class PeriodicDMRG(CanonicalDMRG):
    def __init__(self,L,N,h,onsite,seed=17):
        self.setting_up=True
        super().__init__(L,N,h,onsite,seed)
        self.setting_up=False
        zero=(np.zeros((1,1)),np.zeros((self.rank,1,1)),np.zeros((self.rank,1,1)),0)
        self.LE[0]=zero;self.RE[L]=zero
        self.rebuild_right()

    def grow_left(self,env,M):
        H,ops=super().grow_left(env,M)
        if len(env)<4:return H,ops
        if env[3]==0:
            first=np.einsum('asb,kst,atc->kbc',M,self.opB,M,optimize=True)
        else:
            first=np.einsum('asb,kac,csd->kbd',M,env[2],M,optimize=True)
        return H,ops,first,env[3]+1

    def grow_right(self,env,M):
        H,ops=super().grow_right(env,M)
        if len(env)<4:return H,ops
        if env[3]==0:
            last=np.einsum('asb,kst,ctb->kac',M,self.opA,M,optimize=True)
        else:
            last=np.einsum('asb,kbd,csd->kac',M,env[2],M,optimize=True)
        return H,ops,last,env[3]+1

    def closure_ops(self,i,l,d,r):
        if i==0:left=[np.kron(np.eye(l),B) for B in self.opB]
        else:left=[np.kron(B,np.eye(d)) for B in self.LE[i][2]]
        if i+1==self.L-1:right=[np.kron(A,np.eye(r)) for A in self.opA]
        else:right=[np.kron(np.eye(d),A) for A in self.RE[i+2][2]]
        return list(zip(left,right))

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=24);p.add_argument('--N',type=int,default=8);p.add_argument('--delta',type=float,default=.9);p.add_argument('--chi',type=int,default=96);p.add_argument('--sweeps',type=int,default=10);p.add_argument('--output',default='periodic_candidate');args=p.parse_args()
 h,o=qutrit_model(args.delta,mu=0)
 obj=PeriodicDMRG(args.L,args.N,h,o)
 obj.run([8,16,24,32]+[args.chi]*(args.sweeps-4),output=args.output,tol=1e-9)
