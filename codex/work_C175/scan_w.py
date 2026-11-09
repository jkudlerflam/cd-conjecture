"""Scan |w| at phi = pi/2: threshold data (A-matrix Hermitian part, a_aa), Delta, and the bound pair at K = pi."""
import numpy as np, scipy.linalg as sl
from scat import *
def EK(w,K,Lr=1500):
    c=rel_coeffs(K,w); n=Lr
    # banded solver: lowest eigenvalue of the relative problem on r=2..Lr+1
    ab=np.zeros((3,n),complex)
    ab[0,:]=c[0]+np.array([Vpot(r) for r in range(2,n+2)])
    ab[1,:-1]=c[-1]; ab[2,:-2]=c[-2]
    return sl.eig_banded(ab,lower=True,select='i',select_range=(0,0),eigvals_only=True)[0]
print("# |w|  k_a  eps''  a_aa  A11  |A12|  arg(A12)  a_minus a_plus  Delta  E_b(K=pi)  1/(m a_plus^2)")
for wa in [0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,3.5,3.9]:
    w=1j*wa
    A,ka,kb,K=A_matrix(w); AH=(A+A.conj().T)/2
    aa,_=a_like(w,0)
    (k0,e0)=valleys(w)[0]
    epp=8-wa**2/2
    ev=np.linalg.eigvalsh(AH)
    Lr=max(1500,int(40*ev[1]))
    Eb=2*e0-EK(w,np.pi,Lr)
    D=(aa.real-AH[0,0].real)/abs(AH[0,1])
    print(f"{wa:5.2f} {k0:7.4f} {epp:6.3f} {aa.real:8.4f} {AH[0,0].real:8.4f} {abs(AH[0,1]):8.4f} {np.angle(AH[0,1]):7.4f} {ev[0]:8.4f} {ev[1]:9.4f} {D:8.4f} {Eb:10.6f} {epp/ev[1]**2:10.6f}")
