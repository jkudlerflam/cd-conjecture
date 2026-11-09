"""Coordinator calibration: C180 criterion A_d - |A_o| < a_s < A_d + |A_o| (c = 2), with V2, V3 switched on separately.
(V2, V3) = (0, 0) at |w| = 1 is the V -> infinity limit of Mishra-Greschner-Santos (arXiv:1411.7270), expected two-component.
Also: lowest two-body level over all K minus 2 eps_min (negative = bound pair below the gas threshold)."""
import numpy as np, scipy.linalg as sl
from scat import A_matrix, a_like, valleys, eps, rel_coeffs
def E0(w,K,V2,V3,Lr=1500):
    c=rel_coeffs(K,w); n=Lr; pot=np.zeros(n); pot[0]=V2; pot[1]=V3
    ab=np.zeros((3,n),complex); ab[0,:]=c[0]+pot; ab[1,:-1]=c[-1]; ab[2,:-2]=c[-2]
    return sl.eig_banded(ab,lower=True,select='i',select_range=(0,0),eigvals_only=True)[0]
print("# |w|  V2    V3     a_s      A_d     |A_o|   window(A_d-|A_o|, A_d+|A_o|)  Delta   verdict(Tonks)   min_K E2-2eps_min")
for wa in (0.5,1.0,2.0):
    w=1j*wa; em=valleys(w)[0][1]
    for V2,V3 in ((0,0),(-2,0),(0,-2/3),(-2,-2/3)):
        A,ka,kb,K=A_matrix(w,V2=V2,V3=V3); AH=(A+A.conj().T)/2
        a_s=a_like(w,0,V2=V2,V3=V3)[0].real; Ad=AH[0,0].real; Ao=abs(AH[0,1])
        D=(a_s-Ad)/Ao
        v='c=2' if abs(D)<1 else ('vector chiral (one valley)' if D<=-1 else 'valley crystal')
        Lr=3000 if wa<1 else 1500
        eb=min(E0(w,KK,V2,V3,Lr) for KK in np.linspace(-np.pi,np.pi,49))-2*em
        print(f"{wa:4.1f} {V2:5.1f} {V3:6.3f} {a_s:8.4f} {Ad:8.4f} {Ao:8.4f}   ({Ad-Ao:8.4f}, {Ad+Ao:8.4f})   {D:7.3f}  {v:28s} {eb:9.5f}")
