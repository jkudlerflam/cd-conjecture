"""C175: energy-dependent two-body S-matrix at phi = pi/2 and the derived 'scattering-length matrices'
A(q) = (i/2q) log(-S(q)) for the unlike pair (K = k_a + k_b, channels by valley of the right particle)
and a_aa(q) for the like pair (K = 2 k_a).  Flux-normalized open channels; evanescent modes kept if decaying."""
import numpy as np, scipy.linalg as sl
from scat import valleys, eps, rel_coeffs, boundary_rows, free_roots
def vel(p,K,w,h=1e-6):
    E=lambda x: eps(K-x,w)+eps(x,w)
    return (E(p+h)-E(p-h))/(2*h)
def smat(K,E,w):
    z=free_roots(K,E,w)
    opn=[x for x in z if abs(abs(x)-1)<1e-7]; dec=[x for x in z if abs(x)<1-1e-7]
    ps=np.angle(opn); vs=np.array([vel(p,K,w) for p in ps])
    inc=[i for i in range(len(ps)) if vs[i]<0]; out=[i for i in range(len(ps)) if vs[i]>0]
    modes=[lambda r,x=x: x**r for x in opn]+[lambda r,x=x: x**r for x in dec]
    M=boundary_rows(modes,K,E,w)
    # unknowns: incoming amps (given), outgoing amps, decaying amps
    no=len(out); nd=len(dec)
    cols_o=[out[i] for i in range(no)]; cols_d=list(range(len(opn),len(opn)+nd))
    Mu=M[:,cols_o+cols_d]; S=np.zeros((no,len(inc)),complex)
    for j,ii in enumerate(inc):
        sol=np.linalg.lstsq(Mu,-M[:,ii],rcond=None)[0]
        S[:,j]=sol[:no]*np.sqrt(abs(vs[out]))/np.sqrt(abs(vs[ii]))
    return S,ps[inc],ps[out]
if __name__=="__main__":
    qs=np.linspace(0.01,0.6,60)
    for wa in (0.5,1.0):
        w=1j*wa; (k0,e0),(k1,e1)=valleys(w)[:2]; ka,kb=k0,k1
        th_aa=[]; th_ab=[]; U=0
        for q in qs:
            S1,_,_=smat(2*ka,eps(ka+q,w)+eps(ka-q,w),w); th_aa.append(np.angle(-S1[0,0]))
            S,pin,pout=smat(ka+kb,2*eps(kb+q,w),w)
            U=max(U,np.linalg.norm(S.conj().T@S-np.eye(S.shape[1])))
            th_ab.append(np.sort(np.angle(np.linalg.eigvals(-S))))
        th_aa=np.unwrap(th_aa); th_ab=np.unwrap(np.array(th_ab),axis=0)
        print(f"# |w|={wa}: eigenphases theta of -S (S = -e^{{i theta}}, theta = -2 q a); a = -theta/(2q); max unitarity err {U:.1e}")
        print("#   q    a_aa(q)   a_1(q)   a_2(q)   [unlike-pair eigen-channels, branch continuous from q->0]")
        for j,q in enumerate(qs):
            if j%6==0: print(f"  {q:5.2f} {-th_aa[j]/(2*q):8.3f} {-th_ab[j,0]/(2*q):8.3f} {-th_ab[j,1]/(2*q):8.3f}")
