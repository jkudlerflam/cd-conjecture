"""Ring ED at the symmetric phase phi = pi/2: low levels, total momentum, valley imbalance of each level."""
import numpy as np, scipy.sparse.linalg as spla, sys
from ring_ed import build_ring, obdm, nk
def translate(b,L):
    return ((b<<1)|(b>>(L-1)))&((1<<L)-1)
def run(L,N,phi,wa=1.0,k=12):
    w=wa*np.exp(1j*phi)
    b,H=build_ring(L,N,w,bc=1)
    e,v=spla.eigsh(H,k=k,which='SA',tol=1e-10); o=np.argsort(e); e=e[o]; v=v[:,o]
    tb=translate(b,L); perm=np.searchsorted(b,tb)
    # fermion sign for translation: particle crossing from site L-1 to 0 picks (-1)^(N-1) for periodic fermions
    sg=np.where((b>>(L-1))&1, (-1.0)**(N-1), 1.0)
    out=[]
    i=0
    while i<len(e):
        j=i
        while j+1<len(e) and abs(e[j+1]-e[i])<1e-8: j+=1
        V=v[:,i:j+1]
        # diagonalize translation in the degenerate block
        TV=np.zeros_like(V); TV[perm]=V*sg[:,None]
        Tm=V.conj().T@TV
        ev,U=np.linalg.eig(Tm); Vd=V@U
        for m in range(Vd.shape[1]):
            psi=Vd[:,m]/np.linalg.norm(Vd[:,m])
            G=obdm(b,psi,L,N); ks,n=nk(G,L,1,N)
            imb=n[np.cos(ks)>0].sum()-n[np.cos(ks)<0].sum()
            K=np.angle(ev[m])/np.pi
            out.append((e[i+m] if False else e[i], K, imb))
        i=j+1
    return out
if __name__=="__main__":
    for (L,N) in ((24,3),(28,4),(32,4)):
        for phi in (np.pi/2,):
            res=run(L,N,phi)
            e0=res[0][0]
            print(f"L={L} N={N} phi={phi:.4f}")
            for (E,K,imb) in res: print(f"   E-E0={E-e0:.5f}  K/pi={K:+.3f}  valley imbalance={imb:+.2f}")
