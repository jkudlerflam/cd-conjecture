"""Ring ED of the blockaded second-order model with fermion boundary condition, momentum distribution n(k).
Fermion c's defined by JW from site 0; a hop crossing the bond (L-1, 0) gets the sign (-1)^(N-1) * bc
(bc = +1 periodic, -1 antiperiodic fermions). Units t2 = -1."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from eff_ed import basis_fixed

def ring_basis(L,N,blockade=True):
    b=basis_fixed(L,N,blockade)
    if blockade: b=b[~(((b>>(L-1))&1).astype(bool) & ((b&1).astype(bool)))]
    return b

def popcount(a):
    a=a.copy(); c=np.zeros_like(a)
    while np.any(a): c+=a&1; a>>=1
    return c

def build_ring(L,N,w,t2=-1.0,V2=-2.0,V3=-2/3,bc=1,blockade=True):
    b=ring_basis(L,N,blockade); n=len(b); bit=lambda j:(b>>(j%L))&1
    diag=np.zeros(n)
    for j in range(L):
        diag+=V2*(bit(j)&bit(j+2))+V3*(bit(j)&bit(j+3))
    rows=[np.arange(n)];cols=[np.arange(n)];vals=[diag.astype(complex)]
    for dist,amp in ((1,w),(2,t2)):
        for j in range(L):
            a=j; c=(j+dist)%L
            ok=(bit(a)==1)&(bit(c)==0)
            if blockade:
                ok&=(bit(c+1)==0)
                if dist==2: ok&=(bit(j+1)==0)
            src=np.nonzero(ok)[0]
            if len(src)==0: continue
            ts=(b[src]&~(1<<a))|(1<<c)
            tgt=np.searchsorted(b,ts)
            sgn=np.ones(len(src))
            if c<a:  # crossed the boundary
                sgn*= bc*(-1)**(N-1)
            if dist==2 and not blockade:
                sgn*=1-2*bit(j+1)[src]
            rows+=[tgt,src];cols+=[src,tgt];vals+=[amp*sgn,np.conj(amp)*sgn]
    H=sp.csr_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n))
    return b,H

def obdm(b,psi,L,N,bc=1):
    """<c+_i c_j> on the ring (JW ordering from site 0)."""
    idx={int(s):k for k,s in enumerate(b)}
    G=np.zeros((L,L),complex)
    for i in range(L):
        G[i,i]=np.sum(np.abs(psi)**2*((b>>i)&1))
    for j in range(L):
        occ=np.nonzero((b>>j)&1)[0]
        for i in range(L):
            if i==j: continue
            sel=occ[((b[occ]>>i)&1)==0]
            if len(sel)==0: continue
            ts=(b[sel]&~(1<<j))|(1<<i)
            pos=np.searchsorted(b,ts); pos=np.minimum(pos,len(b)-1)
            good=b[pos]==ts
            lo,hi=min(i,j),max(i,j)
            mask=((1<<hi)-1)&~((1<<(lo+1))-1)
            between=popcount(b[sel]&mask)
            sg=(-1.0)**between
            G[i,j]=np.sum(np.conj(psi[pos[good]])*psi[sel[good]]*sg[good])
    return G

def nk(G,L,bc=1,N=None):
    sh=0.0 if (bc*(-1)**(N-1))==1 else 0.5
    ks=2*np.pi*(np.arange(L)+sh)/L
    x=np.arange(L)
    out=[]
    for k in ks:
        ph=np.exp(1j*k*(x[:,None]-x[None,:]))
        out.append((ph*G).sum().real/L)
    return ks,np.array(out)

def low(L,N,w,k=4,**kw):
    b,H=build_ring(L,N,w,**kw)
    if H.shape[0]<600:
        e,v=np.linalg.eigh(H.toarray()); return b,e[:k],v[:,:k]
    e,v=spla.eigsh(H,k=k,which='SA',tol=1e-10); o=np.argsort(e); return b,e[o],v[:,o]
