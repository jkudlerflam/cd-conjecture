"""Full legal-sector ED of the C28 frustration-free bus with local densities."""
import json
import time
import argparse
import numpy as np
from scipy import sparse as sp
from scipy.linalg import eigh,polar
from scipy.sparse.linalg import eigsh
from bus_base import h2,data_ops
from ff_bus import cell,field,charges,PHASES
from explore import coefficients,normal_order,roots,scattering

def kron(xs):
    out=sp.csr_matrix([[1.]])
    for x in xs:out=sp.kron(out,x,format='csr')
    return out

def embed(op,j,M,D,width=1):return kron([sp.eye(D**j),sp.csr_matrix(op),sp.eye(D**(M-j-width))])

def parameters(n):
    coeff=coefficients(n,phases=PHASES,ratio=1.)
    _,h,_,_=normal_order(coeff[2],2*n)
    M=h[:n,:n]+h[n:,n:];T=h[:n,n:]
    if n==1:mu=.37*abs(T[0,0])
    else:
        b,c=T[0,1].real,T[1,0].real
        mu=np.sqrt(b*b+c*c+2*b*c*np.cos(1.4))
    values=np.full(n,mu)-M.diagonal().real
    return values,coeff,h

def ed(n,M,lam,g=.2):
    start=time.monotonic();H,G,ports,alpha=cell(n,g=g);D,N=len(H),2**n
    values,_,_=parameters(n);Fp,Fl=field(n,values)
    scale=alpha**4*g*g
    L,R=ports
    v=np.kron(R[0],L[0])+1j*np.kron(R[1],L[2])-1j*np.kron(R[2],L[1])
    HT=sum((embed(H+lam**2*scale*Fp,j,M,D) for j in range(M)),sp.csr_matrix((D**M,D**M)))
    V=sum((embed(v,j,M,D,2) for j in range(M-1)),sp.csr_matrix(HT.shape));HT+=lam*V
    GG=kron([sp.csr_matrix(G)]*M)
    q=charges(n);qs=np.array([0])
    for _ in range(M):qs=(qs[:,None]+q[None,:]).ravel()
    qlog=np.array([i.bit_count() for i in range(N**M)])
    K=np.zeros((N**M,N**M),complex);lows=[];nexts=[];res=0.;overlap=1.;sector_dims=[]
    for charge in range(n*M+1):
        inds=np.flatnonzero(qs==charge);logs=np.flatnonzero(qlog==charge);k=len(logs)
        B=HT[inds][:,inds].tocsr();sector_dims.append(len(inds))
        if len(inds)<200:e,U=eigh(B.toarray(),subset_by_index=(0,k))
        else:
            rng=np.random.default_rng(2800+n*100+M*10+charge)
            e,U=eigsh(B,k=k+1,which='SA',tol=1e-12,ncv=min(len(inds),max(4*(k+1)+10,50)),v0=rng.normal(size=len(inds)))
            U,_=np.linalg.qr(U);e,W=eigh(U.conj().T@(B@U));U=U@W
        res=max(res,float(np.linalg.norm(B@U-U*e,axis=0).max()));nexts.append(float(e[k]));lows.extend(e[:k])
        C=GG[inds][:,logs].conj().T@U[:,:k];Q,S=polar(C)
        overlap=min(overlap,float(np.linalg.eigvalsh(S)[0]));K[np.ix_(logs,logs)]=(Q*e[:k])@Q.conj().T
    Kb,_=h2(n,g=g,phases=PHASES)
    pred=alpha**4*sum((embed(Kb,j,M,N,2).toarray() for j in range(M-1)),np.zeros_like(K))
    pred+=scale*sum((embed(Fl,j,M,N).toarray() for j in range(M)),np.zeros_like(K))
    err=K-lam**2*pred;tc=err-np.trace(err)/len(err)*np.eye(len(err))
    row=dict(n=n,cells=M,lambda_=lam,g=g,dimension=D**M,cell_dimension=D,sector_dimensions=sector_dims,
             full_error=float(np.linalg.norm(err,2)),traceless_error=float(np.linalg.norm(tc,2)),
             eigenpair_residual=res,overlap_min=overlap,next_gap=min(nexts)-max(lows),bandwidth=float(max(lows)-min(lows)),
             physical_field_coefficients=(scale*values).tolist(),elapsed_seconds=time.monotonic()-start,
             charge_commutator=float(sp.linalg.norm(HT@sp.diags(qs.astype(float))-sp.diags(qs.astype(float))@HT)))
    np.savez(f'effective_n{n}_M{M}_lam{lam}.npz',exact=K,second=pred)
    with open(f'ed_n{n}_M{M}_lam{lam}.json','w') as f:json.dump(row,f,indent=2)
    print(json.dumps(row),flush=True)

def algebra():
    out=[]
    for n in (1,2):
        values,K,h=parameters(n);_,h4,U4,high=normal_order(K[4],2*n)
        _,_,U2,high2=normal_order(K[2],2*n);rr=roots(h,values)
        H,G,ports,a=cell(n);F,Fl=field(n,values);q=charges(n)
        row=dict(n=n,phases=PHASES,fields=values.tolist(),h2_real=h.real.tolist(),h2_imag=h.imag.tolist(),
                 h4_real=h4.real.tolist(),h4_imag=h4.imag.tolist(),quartic4_real=U4.real.tolist(),quartic4_imag=U4.imag.tolist(),
                 g2_quartic_norm=float(np.linalg.norm(U2)),g4_degree_above_four=high,
                 g_odd_norms=[float(np.linalg.norm(K[r])) for r in (1,3)],
                 cell_gap=float(eigh(H,eigvals_only=True)[2**n]),code_error=float(np.linalg.norm(H@G)),
                 field_invariance=float(np.linalg.norm(F@G-G@Fl)),charge_error=float(np.linalg.norm(q[:,None]*H-H*q[None,:])),
                 roots=[dict(k=r[0],v=r[1]) for r in rr])
        if n==2:
            R=[r for r in rr if r[1]>0];L=sorted([r for r in rr if r[1]<0],key=lambda r:-r[0])
            fwd=np.array([[scattering(U4,r,l,r,l).real for l in L] for r in R])
            pair=scattering(U4,R[1],L[1],R[0],L[0]);v=np.array([r[1] for r in R])
            x=-(fwd[0,0]/(2*np.pi*v[0])+fwd[1,1]/(2*np.pi*v[1])-(fwd[0,1]+fwd[1,0])/(np.pi*sum(v)))
            y=abs(pair)/np.pi*np.sqrt(1/(2*v[0]*v[1])+2/sum(v)**2)
            row.update(forward_RL=fwd.tolist(),pair_transfer=[pair.real,pair.imag],x_over_g2=x,y_over_g2=y,
                       bkt_margin_over_g2=-x-y,backscatter_mismatch=float(np.angle(np.exp(1j*(R[0][0]+L[1][0]-R[1][0]-L[0][0])))))
        np.savez(f'coefficients_n{n}.npz',K0=K[0],K2=K[2],K4=K[4],fields=values)
        out.append(row);print(json.dumps(row),flush=True)
    with open('final_algebra.json','w') as f:json.dump(out,f,indent=2)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--algebra',action='store_true');p.add_argument('--n',type=int);p.add_argument('--cells',type=int);p.add_argument('--lam',type=float,default=.08)
    args=p.parse_args()
    if args.algebra:algebra()
    else:ed(args.n,args.cells,args.lam)
