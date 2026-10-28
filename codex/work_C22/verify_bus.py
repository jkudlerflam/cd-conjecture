"""Full ED of legal, 17-letter walker-bus cells and analytic diagnostics."""
import json
import math
import time
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh, polar
from scipy.integrate import quad
from spatial_bus import walker_cell, h2, h2_g2, quadratic_fit, data_ops, kron_all


PHASES=(.41,.31,-.07)


def skron(xs):
    out=sp.csr_matrix([[1.]])
    for x in xs:
        out=sp.kron(out,x,format='csr')
    return out


def embed(op,j,M,D,width=1):
    return skron([sp.eye(D**j),sp.csr_matrix(op),sp.eye(D**(M-j-width))])


def local_charges(n):
    m=n+1
    q=np.array([i.bit_count() for i in range(2**n)])
    return np.concatenate([np.tile(q,m),np.concatenate([q+a for _ in range(m) for a in (0,1,-1)])])


def ed(n,M,lam,g=.2):
    start=time.monotonic()
    H,G,ports,vis=walker_cell(n,g=g,phases=PHASES)
    D,N=len(H),2**n
    L,R=ports
    v=np.kron(R[0],L[0])+1j*np.kron(R[1],L[2])-1j*np.kron(R[2],L[1])
    HT=sum((embed(H,j,M,D) for j in range(M)),sp.csr_matrix((D**M,D**M)))
    V=sum((embed(v,j,M,D,2) for j in range(M-1)),sp.csr_matrix(HT.shape))
    HT=HT+lam*V
    GG=skron([sp.csr_matrix(G)]*M)
    # Elementwise charge labels respect the tensor ordering of the Hamiltonian.
    q=local_charges(n)
    qs=np.array([0])
    for _ in range(M):
        qs=(qs[:,None]+q[None,:]).reshape(-1)
    qlog=np.array([i.bit_count() for i in range(N**M)])
    K=np.zeros((N**M,N**M),complex)
    all_low=[]; next_levels=[]; res=0.; min_overlap=1.
    for charge in range(n*M+1):
        inds=np.flatnonzero(qs==charge)
        logs=np.flatnonzero(qlog==charge)
        k=len(logs)
        B=HT[inds][:,inds].tocsr()
        if len(inds)<200:
            e,U=eigh(B.toarray(),subset_by_index=(0,k))
        else:
            rng=np.random.default_rng(100+n*10+M+charge)
            e,U=eigsh(B,k=k+1,which='SA',tol=2e-12,
                      ncv=min(len(inds),max(4*(k+1)+10,50)),
                      v0=rng.normal(size=len(inds)))
            ix=np.argsort(e); e=e[ix];U=U[:,ix]
            # Reorthogonalize near-degenerate eigenvectors before Ritz.
            U,_=np.linalg.qr(U)
            e,rot=eigh(U.conj().T@(B@U));U=U@rot
        res=max(res,float(np.linalg.norm(B@U-U*e,axis=0).max()))
        next_levels.append(float(e[k]))
        all_low.extend(e[:k].tolist())
        C=GG[inds][:,logs].conj().T@U[:,:k]
        Q,S=polar(C)
        min_overlap=min(min_overlap,float(np.linalg.eigvalsh(S)[0]))
        K[np.ix_(logs,logs)]=(Q*e[:k])@Q.conj().T
    Kbond,_=h2(n,g=g,phases=PHASES)
    Kpred=sum((embed(Kbond,j,M,N,2).toarray() for j in range(M-1)),np.zeros_like(K))
    Kpred*=lam**2*vis**4
    diff=K-Kpred
    trdiff=diff-np.trace(diff)*np.eye(len(diff))/len(diff)
    row=dict(n=n,cells=M,dimension=D**M,cell_dimension=D,lambda_=lam,g=g,
             elapsed_seconds=time.monotonic()-start,eigen_residual=res,
             code_overlap_min=min_overlap,next_gap=min(next_levels)-max(all_low),
             full_error=float(np.linalg.norm(diff,2)),traceless_error=float(np.linalg.norm(trdiff,2)),
             bandwidth=float(max(all_low)-min(all_low)),visibility=vis,
             phases=PHASES,commutator_charge=float(sp.linalg.norm(HT@sp.diags(qs.astype(float))-sp.diags(qs.astype(float))@HT)),
             hermiticity=float(sp.linalg.norm(HT-HT.conj().T)))
    print(json.dumps(row),flush=True)
    return row


def analytic_checks():
    rows=[]
    for n in (1,2,3):
        H,G,ports,a=walker_cell(n,g=.1,phases=PHASES)
        ev=eigh(H,eigvals_only=True)
        comp=max(np.linalg.norm(G.conj().T@F@G) for p in ports for F in p)
        q=local_charges(n)
        K=h2_g2(n,phases=PHASES)
        fit,one,_=quadratic_fit(K,2*n)
        T=one[:n,n:];M=one[:n,:n]+one[n:,n:]
        bands=np.array([np.linalg.eigvalsh(M+np.exp(1j*k)*T+np.exp(-1j*k)*T.conj().T)
                        for k in np.linspace(-np.pi,np.pi,20001)])
        # independent heat-kernel integral formula
        m=n+1
        A=4*np.eye(m)-np.diag(np.ones(m-1),1)-np.diag(np.ones(m-1),-1)
        e,U=eigh(A)
        def f(s,i,tau):
            # analytic time convolution, avoiding nested quadratures
            weights=U[s]*U[i]
            ea=e[:,None];eb=e[None,:]
            dif=ea-eb
            out=np.empty_like(dif)
            same=abs(dif)<1e-10
            out[same]=np.broadcast_to(tau*np.exp(-ea*tau),dif.shape)[same]
            out[~same]=((np.exp(-eb*tau)-np.exp(-ea*tau))/np.where(same,1,dif))[~same]
            return np.sum(weights[:,None]*weights[None,:]*out)
        kap=PHASES[0]-(PHASES[1]+PHASES[2])/2
        chi=(PHASES[2]-PHASES[1])/2
        Tf=np.zeros((n,n),complex)
        for i in range(n):
            for j in range(n):
                F=quad(lambda tau:f(n,i,tau)*f(0,j,tau),0,np.inf,epsabs=1e-15,epsrel=1e-9)[0]
                Tf[i,j]=2j*1j*F*np.exp(1j*chi*(n-i+j))*np.sin(kap*(n-i-j))
        Kreal=h2_g2(n)
        fr,_,_=quadratic_fit(Kreal,2*n)
        row=dict(n=n,cell_dimension=len(H),code_dimension=int(sum(abs(ev)<1e-10)),
                 code_residual=float(np.linalg.norm(H@G)),gap=float(ev[2**n]),
                 first_order=float(comp),parity_error=float(np.linalg.norm(H@np.diag((-1.)**q)-np.diag((-1.)**q)@H)),
                 g2_quadratic_residual=fit['residual'],g2_field_norm=float(np.linalg.norm(M)),
                 no_flux_quadratic_norm=fr['quadratic_norm'],heat_formula_error=float(np.linalg.norm(T-Tf)),
                 visibility=a,T_real=T.real.tolist(),T_imag=T.imag.tolist(),
                 band_min=bands.min(axis=0).tolist(),band_max=bands.max(axis=0).tolist(),
                 smallest_abs_band=float(np.abs(bands).min()),
                 pole_count_g0=m)
        print(json.dumps(row),flush=True);rows.append(row)
    return rows


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--algebra',action='store_true')
    p.add_argument('--n',type=int);p.add_argument('--cells',type=int);p.add_argument('--lam',type=float,default=.01)
    a=p.parse_args()
    if a.algebra:
        result=analytic_checks();name='analytic_results.json'
    else:
        result=ed(a.n,a.cells,a.lam);name=f'ed_n{a.n}_M{a.cells}_lam{a.lam}.json'
    with open(name,'w') as f:json.dump(result,f,indent=2)
