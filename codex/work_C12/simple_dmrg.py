"""Two-site finite DMRG for real Hermitian nearest-neighbor Hamiltonians.
No external tensor-network package. NumPy/SciPy, dense renormalized blocks.
"""
import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import LinearOperator,eigsh
import time,json

def decompose(h,d):
    mat=h.reshape(d,d,d,d).transpose(0,2,1,3).reshape(d*d,d*d)
    u,s,vh=svd(mat,full_matrices=False,check_finite=False)
    keep=s>1e-12
    return (u[:,keep]*np.sqrt(s[keep])).T.reshape(-1,d,d),(vh[keep]*np.sqrt(s[keep,None])).reshape(-1,d,d)

class DMRG:
    def __init__(self,L,h,onsite,seed=17):
        self.L=L;self.d=onsite.shape[0];self.h=h;self.h4=h.reshape((self.d,)*4);self.onsite=onsite
        self.opA,self.opB=decompose(h,self.d);self.rank=len(self.opA)
        rng=np.random.default_rng(seed)
        bonds=[1]+[2]*(L-1)+[1]
        self.mps=[rng.normal(size=(bonds[i],self.d,bonds[i+1])) for i in range(L)]
        self.right_canon()
        self.LE=[None]*(L+1);self.RE=[None]*(L+1)
        self.LE[0]=(np.zeros((1,1)),np.zeros((self.rank,1,1)))
        self.RE[L]=(np.zeros((1,1)),np.zeros((self.rank,1,1)))
        self.rebuild_right()

    def right_canon(self):
        for i in range(self.L-1,0,-1):
            l,d,r=self.mps[i].shape
            q,rr=np.linalg.qr(self.mps[i].reshape(l,d*r).T,mode='reduced')
            self.mps[i]=q.T.reshape(q.shape[1],d,r)
            self.mps[i-1]=np.tensordot(self.mps[i-1],rr.T,axes=(2,0))
        self.mps[0]/=np.linalg.norm(self.mps[0])

    def left_pair(self,env):
        H,ops=env[:2];l=len(H)
        out=np.kron(H,np.eye(self.d))+np.kron(np.eye(l),self.onsite)
        for X,B in zip(ops,self.opB):out+=np.kron(X,B)
        return (out+out.T)/2

    def right_pair(self,env):
        H,ops=env[:2];r=len(H)
        out=np.kron(np.eye(self.d),H)+np.kron(self.onsite,np.eye(r))
        for A,X in zip(self.opA,ops):out+=np.kron(A,X)
        return (out+out.T)/2

    def grow_left(self,env,M):
        l,d,r=M.shape;mat=M.reshape(l*d,r)
        H=mat.T@self.left_pair(env)@mat
        acted=np.tensordot(self.opA,M,axes=(2,1)).transpose(0,2,1,3).reshape(self.rank,l*d,r)
        ops=np.matmul(mat.T,acted)
        return H,ops

    def grow_right(self,env,M):
        l,d,r=M.shape;mat=M.reshape(l,d*r)
        H=mat@self.right_pair(env)@mat.T
        acted=np.tensordot(self.opB,M,axes=(2,1)).transpose(0,2,1,3).reshape(self.rank,l,d*r)
        ops=np.matmul(mat,acted.transpose(0,2,1))
        return H,ops

    def rebuild_right(self):
        for i in range(self.L-1,-1,-1):self.RE[i]=self.grow_right(self.RE[i+1],self.mps[i])

    def optimize_bond(self,i,chi,direction,tol):
        theta=np.tensordot(self.mps[i],self.mps[i+1],axes=(2,0))
        l,d,_,r=theta.shape
        HL=self.left_pair(self.LE[i]);HR=self.right_pair(self.RE[i+2]);shape=(l*d,d*r)
        def mv(vec):
            mat=vec.reshape(shape)
            out=HL@mat+mat@HR.T
            out+=np.tensordot(self.h4,mat.reshape(l,d,d,r),axes=((2,3),(1,2))).transpose(2,0,1,3).reshape(shape)
            return out.ravel()
        n=theta.size
        if n<30:
            mat=np.column_stack([mv(np.eye(n)[:,a]) for a in range(n)])
            en,vs=np.linalg.eigh(mat);e=en[0];v=vs[:,0]
        else:
            lin=LinearOperator((n,n),matvec=mv,dtype=float)
            en,vs=eigsh(lin,k=1,which='SA',v0=theta.ravel(),tol=tol,maxiter=2000,ncv=min(24,n))
            e=en[0];v=vs[:,0]
        u,s,vh=svd(v.reshape(shape),full_matrices=False,check_finite=False,lapack_driver='gesdd')
        k=min(chi,len(s));discard=float(sum(s[k:]**2));s=s[:k];s/=np.linalg.norm(s)
        if direction==1:
            self.mps[i]=u[:,:k].reshape(l,d,k)
            self.mps[i+1]=(s[:,None]*vh[:k]).reshape(k,d,r)
            self.LE[i+1]=self.grow_left(self.LE[i],self.mps[i])
        else:
            self.mps[i]=(u[:,:k]*s).reshape(l,d,k)
            self.mps[i+1]=vh[:k].reshape(k,d,r)
            self.RE[i+1]=self.grow_right(self.RE[i+2],self.mps[i+1])
        return float(e),discard

    def diagnostics(self):
        # Called with orthogonality center on site 0.
        norm=np.ones((1,1));density=[]
        n=np.diag(np.arange(self.d));double=np.diag(np.arange(self.d)==2)
        for M in self.mps:
            weighted=np.tensordot(norm,M,axes=(1,0))
            rho=np.tensordot(M,weighted,axes=((0,2),(0,2)))
            density.append([float(np.trace(rho@n)),float(np.trace(rho@double))])
            norm=sum(M[:,s,:].T@norm@M[:,s,:] for s in range(self.d))
        prof=[]
        center=self.mps[0].copy()
        for i in range(self.L-1):
            l,d,r=center.shape
            u,s,vh=svd(center.reshape(l*d,r),full_matrices=False,check_finite=False)
            p=s*s;p=p[p>1e-14];p/=sum(p)
            prof.append(dict(cut=i+1,S=float(-sum(p*np.log(p))),S2=float(-np.log(sum(p*p)))))
            center=np.tensordot(s[:,None]*vh,self.mps[i+1],axes=(1,0))
        return dict(density=density,N=float(np.sum(np.array(density)[:,0])),entropy=prof)

    def run(self,schedule,output=None,tol=1e-8):
        history=[];start=time.time()
        for sweep,chi in enumerate(schedule):
            disc=0
            for i in range(self.L-1):e,dd=self.optimize_bond(i,chi,1,tol);disc=max(disc,dd)
            for i in range(self.L-2,-1,-1):e,dd=self.optimize_bond(i,chi,-1,tol);disc=max(disc,dd)
            row=dict(sweep=sweep+1,chi=chi,E=e,max_discard=disc,seconds=time.time()-start,**self.diagnostics())
            history.append(row)
            print(json.dumps({k:row[k] for k in ['sweep','chi','E','N','max_discard','seconds']}),flush=True)
            if output:
                with open(output+'.json','w') as f:json.dump(history,f,indent=2)
                np.savez(output+'_mps.npz',**{str(i):M for i,M in enumerate(self.mps)})
        return history

def qutrit_model(delta=0.9,mu=-1.6,V12=20,V22=20,t21=0):
    from qutrit_ed import local_h
    U=-14/9+delta
    onsite=np.diag([0,-mu,U-2*mu])
    h=local_h(U=0,V12=V12,V22=V22,t21=t21)
    return h,onsite

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=32);p.add_argument('--delta',type=float,default=.9);p.add_argument('--mu',type=float,default=-1.6);p.add_argument('--chi',type=int,default=64);p.add_argument('--output',default='dmrg_candidate');args=p.parse_args()
    h,o=qutrit_model(args.delta,args.mu)
    obj=DMRG(args.L,h,o)
    schedule=[8,16,24,32,args.chi,args.chi,args.chi,args.chi]
    obj.run(schedule,args.output)
