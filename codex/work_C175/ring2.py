"""Independent brute-force check: two hard-core bosons on a ring of L sites, full two-particle basis,
model (3.2) with blockade (no two adjacent), NN hop w, NNN hop -1 (only over an empty site), V2=-2, V3=-2/3.
Returns the lowest levels relative to 2 eps_min."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, sys
from scat import valleys
def ring2(L, w, t2=-1.0, V2=-2.0, V3=-2/3, k=6):
    conf=[(x,y) for x in range(L) for y in range(x+2,L) if (y-x)%L>=2 and (x-y)%L>=2]
    idx={c:i for i,c in enumerate(conf)}
    def key(a,b): return (min(a,b),max(a,b))
    rows=[];cols=[];vals=[]
    for i,(x,y) in enumerate(conf):
        d=min((y-x)%L,(x-y)%L)
        dv={2:V2,3:V3}.get((y-x)%L,0)+{2:V2,3:V3}.get((x-y)%L,0)
        rows.append(i);cols.append(i);vals.append(dv)
        for p,q in ((x,y),(y,x)):
            for s,amp in ((1,w),(-1,np.conj(w)),(2,t2),(-2,t2)):
                np_=(p+s)%L
                if np_==q: continue
                if abs(s)==2 and ((p+s//2)%L)==q: continue   # no passing
                # blockade of new position
                if (np_-q)%L in (1,L-1): continue
                j=idx[key(np_,q)]
                rows.append(j);cols.append(i);vals.append(amp)
    H=sp.csr_matrix((vals,(rows,cols)),shape=(len(conf),)*2)
    assert abs(H-H.getH()).max()<1e-12
    e=np.sort(spla.eigsh(H,k=k,which='SA',tol=1e-12)[0])
    return e
if __name__=="__main__":
    for wa in (0.5,1.0,2.0):
        w=1j*wa; em=valleys(w)[0][1]
        for L in (120,240):
            e=ring2(L,w)
            print(f"|w|={wa} L={L}: lowest E - 2eps_min =", np.round(e[:4]-2*em,5))
