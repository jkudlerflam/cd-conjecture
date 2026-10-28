"""Charge-conserving extension of simple_dmrg, with block SVD and fixed N."""
import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import LinearOperator,eigsh
from simple_dmrg import DMRG,qutrit_model

class CanonicalDMRG(DMRG):
    def __init__(self,L,N,h,onsite,seed=17):
        super().__init__(L,h,onsite,seed)
        # Deterministic distributed charge; only a starting state.
        occ=np.diff(np.floor(np.arange(L+1)*N/L+1e-8)).astype(int)
        if np.max(occ)>=self.d:raise ValueError('density exceeds local capacity')
        self.q=[np.array([x],int) for x in np.r_[0,np.cumsum(occ)]]
        self.mps=[]
        for n in occ:
            A=np.zeros((1,self.d,1));A[0,n,0]=1;self.mps.append(A)
        self.rebuild_right()

    def optimize_bond(self,i,chi,direction,tol):
        theta=np.tensordot(self.mps[i],self.mps[i+1],axes=(2,0))
        l,d,_,r=theta.shape;shape=(l*d,d*r)
        HL=self.left_pair(self.LE[i]);HR=self.right_pair(self.RE[i+2])
        qrow=np.repeat(self.q[i],d)+np.tile(np.arange(d),l)
        qcol=np.tile(self.q[i+2],d)-np.repeat(np.arange(d),r)
        blocks=[];indices=[]
        for q in np.intersect1d(qrow,qcol):
            a=np.flatnonzero(qrow==q);b=np.flatnonzero(qcol==q)
            ix=(a[:,None]*(d*r)+b[None,:]).ravel()
            blocks.append((int(q),a,b,HL[np.ix_(a,a)],HR[np.ix_(b,b)]));indices.extend(ix)
        indices=np.asarray(indices,int)
        closure=self.closure_ops(i,l,d,r) if hasattr(self,"closure_ops") else []
        closure_blocks=[]
        for X,Y in closure:
            for oo,(_,ao,bo,_,_) in enumerate(blocks):
                for ii,(_,ai,bi,_,_) in enumerate(blocks):
                    xx=X[np.ix_(ao,ai)];yy=Y[np.ix_(bo,bi)]
                    if np.max(abs(xx))>1e-13 and np.max(abs(yy))>1e-13:
                        closure_blocks.append((oo,ii,xx,yy))
        def mv(vec):
            mat=np.zeros(shape);mat.ravel()[indices]=vec
            out=np.zeros(shape)
            subs=[mat[np.ix_(a,b)] for q,a,b,hh,rr in blocks]
            outs=[hh@sub+sub@rr.T for sub,(q,a,b,hh,rr) in zip(subs,blocks)]
            for oo,ii,xx,yy in closure_blocks:outs[oo]+=xx@subs[ii]@yy.T
            for sub,(q,a,b,hh,rr) in zip(outs,blocks):out[np.ix_(a,b)]=sub
            out+=np.tensordot(self.h4,mat.reshape(l,d,d,r),axes=((2,3),(1,2))).transpose(2,0,1,3).reshape(shape)
            return out.ravel()[indices]
        n=len(indices);v0=theta.reshape(-1)[indices]
        if n<30:
            mat=np.column_stack([mv(np.eye(n)[:,a]) for a in range(n)])
            en,vs=np.linalg.eigh(mat);e=en[0];v=vs[:,0]
        else:
            lin=LinearOperator((n,n),matvec=mv,dtype=float)
            en,vs=eigsh(lin,k=1,which='SA',v0=v0,tol=tol,maxiter=2000,ncv=min(24,n))
            e=en[0];v=vs[:,0]
        mat=np.zeros(shape);mat.ravel()[indices]=v
        pieces=[]
        for q,a,b,hh,rr in blocks:
            u,s,vh=svd(mat[np.ix_(a,b)],full_matrices=False,check_finite=False)
            for j,value in enumerate(s):
                if value>1e-14:pieces.append((value,q,a,b,u[:,j],vh[j]))
        pieces.sort(key=lambda p:-p[0]);discard=float(sum(p[0]**2 for p in pieces[chi:]));pieces=pieces[:chi]
        k=len(pieces);u=np.zeros((l*d,k));vh=np.zeros((k,d*r));s=[];qs=[]
        for j,(value,q,a,b,uu,vv) in enumerate(pieces):
            u[a,j]=uu;vh[j,b]=vv;s.append(value);qs.append(q)
        s=np.asarray(s);s/=np.linalg.norm(s);self.q[i+1]=np.asarray(qs,int)
        if direction==1:
            self.mps[i]=u.reshape(l,d,k)
            self.mps[i+1]=(s[:,None]*vh).reshape(k,d,r)
            self.LE[i+1]=self.grow_left(self.LE[i],self.mps[i])
        else:
            self.mps[i]=(u*s).reshape(l,d,k)
            self.mps[i+1]=vh.reshape(k,d,r)
            self.RE[i+1]=self.grow_right(self.RE[i+2],self.mps[i+1])
        return float(e),discard

    def run(self,*args,**kwargs):
        result=super().run(*args,**kwargs)
        output=kwargs.get('output')
        if output:np.savez(output+'_charges.npz',**{str(i):q for i,q in enumerate(self.q)})
        return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=48);p.add_argument('--N',type=int,default=16);p.add_argument('--delta',type=float,default=.9);p.add_argument('--chi',type=int,default=96);p.add_argument('--sweeps',type=int,default=8);p.add_argument('--output',default='canonical_candidate');args=p.parse_args()
    h,o=qutrit_model(args.delta,mu=0)
    obj=CanonicalDMRG(args.L,args.N,h,o)
    obj.run([8,16,24,32]+[args.chi]*(args.sweeps-4),output=args.output,tol=1e-9)
