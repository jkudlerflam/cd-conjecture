"""C124 direct-sum, basis and tensor-product checks."""
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp
from checks import ROOT, histories_dense, embed_bond
import xx_kernels as xx


def main():
    out = dict(sectors=[], rotated=[], spectator=[])
    x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex)
    z=np.diag([1.,-1.])
    ks={(0,0):-.3*x-.2*z,(0,1):.1*z+.07*y+.08*np.eye(2),
        (1,0):np.array([[.17]]),(1,1):np.array([[-.09]])}
    h=np.zeros((9,9),complex)
    for a,b in product(range(3),repeat=2):
        alpha=int(a==2); beta=int(b==2)
        for c in ([0,1] if alpha==0 else [2]):
            h[3*c+b,3*a+b]=ks[alpha,beta][c if alpha==0 else 0,a if alpha==0 else 0]
    A=np.kron(h,np.eye(3)); B=np.kron(np.eye(3),h)
    assert np.linalg.norm(A@B-B@A)<1e-12
    H={n:sum((embed_bond(h,n,3,i,(i+1)%n) for i in range(n)),
             np.zeros((3**n,3**n),complex)) for n in [2,4]}
    U=np.exp(2j*np.pi*np.arange(3)[:,None]*np.arange(3)[None,:]/3)/np.sqrt(3)
    hu=np.kron(U,U)@h@np.kron(U,U).conj().T
    Hu={n:sum((embed_bond(hu,n,3,i,(i+1)%n) for i in range(n)),
              np.zeros((3**n,3**n),complex)) for n in [2,4]}
    for beta in [.2,2.,6.]:
        for N in [1,2,3]:
            hist=list(product(range(3),repeat=N))
            zz={ab:float(np.exp(logsumexp(-beta*eigh(k,eigvals_only=True))))
                for ab,k in ks.items()}
            fs={ab:(histories_dense(k,2,1,beta,N,[0]) if ab[0]==0
                    else np.array([np.exp(-beta*k[0,0])])) for ab,k in ks.items()}
            BB=np.zeros((3**N,3**N))
            for si,s in enumerate(hist):
                al=int(s[0]==2)
                if any(int(v==2)!=al for v in s): continue
                ind=sum(v*2**(N-1-j) for j,v in enumerate(s)) if al==0 else 0
                for ti,t in enumerate(hist):
                    be=int(t[0]==2)
                    if any(int(v==2)!=be for v in t): continue
                    BB[si,ti]=sum(fs[al,ga][ind]*zz[ga,be] for ga in [0,1])
            q=histories_dense(H[2],3,2,beta,N,[0])
            qd=histories_dense(H[4],3,4,beta,N,[0,2])
            pred=[]
            for labels in product(range(3),repeat=2*N):
                s=sum(labels[2*j]*3**(N-1-j) for j in range(N))
                t=sum(labels[2*j+1]*3**(N-1-j) for j in range(N))
                pred.append(BB[s,t]*BB[t,s])
            err=max(np.max(abs(q-np.diag(BB)))/max(q),np.max(abs(qd-pred))/max(qd))
            G=float(np.log(q@q)-np.log(qd.sum()))
            R=float(2*logsumexp(-beta*eigh(H[2],eigvals_only=True))
                    -logsumexp(-beta*eigh(H[4],eigvals_only=True)))
            T=np.array([[zz[a,b] for b in [0,1]] for a in [0,1]])
            rT=float(2*np.log(np.trace(T@T))-np.log(np.trace(np.linalg.matrix_power(T,4))))
            assert err<1e-10 and G<1e-10 and R<=np.log(2)+1e-10 and abs(R-rT)<1e-10
            out['sectors'].append(dict(beta=beta,N=N,G=G,R=R,error=float(err)))
            qu=histories_dense(Hu[2],3,2,beta,N,[0])
            qdu=histories_dense(Hu[4],3,4,beta,N,[0,2])
            Gu=float(np.log(qu@qu)-np.log(qdu.sum()))
            assert Gu<=np.log(2)+2*N*np.log(3)+1e-10
            out['rotated'].append(dict(beta=beta,N=N,G=Gu))
    # Critical d=4 chain: stoquastic XX times a slow on-site qubit.
    # Compute the exact product G using the free-fermion and on-site formulas.
    for L in [4,16,128,1024,8192]:
        eps=L**-2
        for a in [.25,1.,4.]:
            tau=2*a*L/(1+eps)
            logz=lambda t:float(np.logaddexp(0.,-eps*t))
            Cslow=logz(2*tau)+(2*L-2)*logz(tau)
            for N in [1,2,4]:
                gx=xx.xx_case(L,a/(1+eps),N)['G']
                gs=logz(2*N*tau)-2*logz(N*tau)
                out['spectator'].append(dict(L=L,a=a,N=N,G=gx+gs,
                    spectator_G=gs,spectator_boundary=Cslow,
                    # C_XX >= G_XX at one insertion.
                    boundary_lower=Cslow+xx.xx_case(L,a/(1+eps),1)['G']))
    (ROOT/'extra_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('extra checks', {k:len(v) for k,v in out.items()})
    print('sector separator error',max(r['error'] for r in out['sectors']))


if __name__=='__main__': main()
