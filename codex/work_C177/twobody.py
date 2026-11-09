"""Two-body problem of the blockaded model on a long ring: bound states below the two-particle continuum.
Relative coordinate r in [2, L-2] (blockade), total momentum K. Units t2 = -1, V2 = -2, V3 = -2/3."""
import numpy as np
def eps(k,w,t2=-1.0): return (w*np.exp(-1j*k)+np.conj(w)*np.exp(1j*k)).real+2*t2*np.cos(2*k)
def two_body(K,w,L=400,t2=-1.0,V2=-2.0,V3=-2/3):
    rs=np.arange(2,L-1); n=len(rs); H=np.zeros((n,n),complex)
    amp={1:w,-1:np.conj(w),2:t2,-2:t2}
    for i,r in enumerate(rs):
        H[i,i]+= V2*((r==2)+(r==L-2)) + V3*((r==3)+(r==L-3))
        for d,a in amp.items():
            # particle 2 moves by d: r -> r+d
            r2=r+d
            if 2<=r2<=L-2: H[r2-2,i]+=a
            # particle 1 moves by d: x -> x+d, r -> r-d, Bloch phase e^{iKd}
            r1=r-d
            if 2<=r1<=L-2: H[r1-2,i]+=a*np.exp(1j*K*d)
    return np.linalg.eigvalsh(H)[0]
if __name__=="__main__":
    ks=np.linspace(-np.pi,np.pi,2001)
    for wa in (0.0,0.5,1.0,2.0):
      for phi in (1.3,):
        w=wa*np.exp(1j*phi); e1=eps(ks,w); emin=e1.min()
        worst=0; rows=[]
        for K in np.linspace(-np.pi,np.pi,49):
            cont=(eps(ks,w)+eps(K-ks,w)).min()
            eb=two_body(K,w)
            rows.append((K,eb-cont)); worst=min(worst,eb-cont)
        Kb=[r[0] for r in rows if r[1]<-1e-3]
        print(f"|w|={wa} phi={phi}: 2*eps_min={2*emin:.4f}, min over K of (E_2body - continuum) = {worst:.4f}; bound at K in [{min(Kb) if Kb else None}, {max(Kb) if Kb else None}] ({len(Kb)}/49);  min E_2 - 2eps_min = {min(two_body(K,w) for K in np.linspace(-np.pi,np.pi,49))-2*emin:.4f}")
