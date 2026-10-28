"""Complete code-connected symmetry-sector ED, no excitation truncation."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json,time,argparse,signal
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh,polar
from model import cell,parameters,combined_bond,target_data,cubic_bond


def kron(xs):
    ans=sp.csr_matrix([[1.]])
    for x in xs:ans=sp.kron(ans,x,format='csr')
    return ans


def embed(op,j,M,D,width=1):
    return kron([sp.eye(D**j),sp.csr_matrix(op),sp.eye(D**(M-j-width))])


def ed(n,M,lam,g=.2,eta=10.):
    signal.alarm(180)
    start=time.monotonic(); H,G,ports,alpha,E,q=cell(n,g)
    D,N=len(H),2**n;m=n+1
    L,R=ports
    v=sp.csr_matrix((D*D,D*D),dtype=complex)
    for r,(_,_,phi,w) in enumerate(parameters(n)):
        v+=w*(kron([R[r][0],L[r][0]])+np.exp(1j*phi)*kron([R[r][1],L[r][2]])+np.exp(-1j*phi)*kron([R[r][2],L[r][1]]))
    T,mu,W=target_data(n,g)
    Wfull=alpha**4*g*g*np.kron(np.eye(D//N),W)
    HT=sum((embed(H+lam*eta*E+lam*lam*Wfull,j,M,D) for j in range(M)),sp.csr_matrix((D**M,D**M)))
    HT+=lam*sum((embed(v,j,M,D,2) for j in range(M-1)),sp.csr_matrix(HT.shape))
    GG=kron([G]*M)
    # Bus types 0 (vacuum), 1,2. Each excited type count modulo two is conserved.
    typ=np.r_[np.zeros(m*N,dtype=int),np.ones(3*m*N,dtype=int),2*np.ones(3*m*N,dtype=int)]
    qs=np.array([0]);pa=np.array([0]);pb=np.array([0])
    for _ in range(M):
        qs=(qs[:,None]+q).reshape(-1)
        pa=(pa[:,None]^(typ==1)).reshape(-1)
        pb=(pb[:,None]^(typ==2)).reshape(-1)
    even=(pa==0)&(pb==0)
    symmetry_error=float(sp.linalg.norm(HT@sp.diags(pa.astype(float))-sp.diags(pa.astype(float))@HT)+sp.linalg.norm(HT@sp.diags(pb.astype(float))-sp.diags(pb.astype(float))@HT))
    qlog=np.array([i.bit_count() for i in range(N**M)])
    K=np.zeros((N**M,N**M),complex);res=0.;overlap=1.;low=[];nxt=[];sizes=[]
    for charge in range(n*M+1):
        inds=np.flatnonzero(even&(qs==charge));logs=np.flatnonzero(qlog==charge);k=len(logs)
        B=HT[inds][:,inds].tocsr();sizes.append(len(inds))
        if len(inds)<200:
            ev,U=eigh(B.toarray(),subset_by_index=(0,k))
        else:
            rng=np.random.default_rng(100+n*10+M+charge)
            ev,U=eigsh(B,k=k+1,which='SA',tol=1e-12,ncv=min(len(inds),max(4*(k+1)+10,50)),v0=rng.normal(size=len(inds)))
            U,_=np.linalg.qr(U);ev,rot=eigh(U.conj().T@(B@U));U=U@rot
        res=max(res,float(np.linalg.norm(B@U-U*ev,axis=0).max()))
        low.extend(ev[:k].tolist());nxt.append(float(ev[k]))
        C=GG[inds][:,logs].conj().T@U[:,:k];Q,S=polar(C)
        overlap=min(overlap,float(np.linalg.eigvalsh(S)[0]))
        K[np.ix_(logs,logs)]=(Q*ev[:k])@Q.conj().T
    bond=alpha**4*combined_bond(n,g)
    K2=sum((embed(bond,j,M,N,2).toarray() for j in range(M-1)),np.zeros_like(K))
    K2+=alpha**4*g*g*sum((embed(W,j,M,N).toarray() for j in range(M)),np.zeros_like(K))
    K3=sum((embed(alpha**4*cubic_bond(n,g,eta),j,M,N,2).toarray() for j in range(M-1)),np.zeros_like(K))
    diff=K-lam*lam*K2;tr=diff-np.trace(diff)*np.eye(len(K))/len(K)
    row=dict(n=n,cells=M,lambda_=lam,g=g,eta=eta,full_legal_dimension=D**M,even_type_dimension=int(even.sum()),charge_sizes=sizes,code_levels=N**M,full_error=float(np.linalg.norm(diff,2)),traceless_error=float(np.linalg.norm(tr,2)),error_over_lambda3=float(np.linalg.norm(diff,2)/lam**3),cubic_norm=float(np.linalg.norm(K3,2)),after_cubic_error=float(np.linalg.norm(diff-lam**3*K3,2)),eigen_residual=res,overlap_min=overlap,next_gap=min(nxt)-max(low),bandwidth=max(low)-min(low),hermiticity=float(sp.linalg.norm(HT-HT.conj().T)),charge_error=float(sp.linalg.norm(HT@sp.diags(qs.astype(float))-sp.diags(qs.astype(float))@HT)),type_parity_error=symmetry_error,elapsed_seconds=time.monotonic()-start)
    signal.alarm(0)
    print(json.dumps(row),flush=True)
    return row

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--cells',type=int,required=True);p.add_argument('--lam',type=float,required=True);p.add_argument('--eta',type=float,default=10.)
    a=p.parse_args();r=ed(a.n,a.cells,a.lam,eta=a.eta)
    with open(f'ed_n{a.n}_M{a.cells}_l{a.lam}_eta{a.eta}.json','w') as f:json.dump(r,f,indent=2)
