"""Unlike-pair binding off the symmetric line: lowest two-body level at K = k_a + k_b against eps_a + eps_b and 2 eps_a."""
import numpy as np, scipy.linalg as sl
from scat import valleys, eps, rel_coeffs, Vpot
def E0(w,K,Lr=1500):
    c=rel_coeffs(K,w); n=Lr
    ab=np.zeros((3,n),complex); ab[0,:]=c[0]+np.array([Vpot(r) for r in range(2,n+2)])
    ab[1,:-1]=c[-1]; ab[2,:-2]=c[-2]
    return sl.eig_banded(ab,lower=True,select='i',select_range=(0,0),eigvals_only=True)[0]
print("# |w| phi  delta=eps_b-eps_a  E0(K=ka+kb)-(eps_a+eps_b)  E0-2eps_a  min_K E0(K)-2eps_a")
for wa in (1.0,2.0):
  for phi in (np.pi/2,1.568,1.565,1.56,1.55,1.54):
    w=wa*np.exp(1j*phi); (ka,ea),(kb,eb)=valleys(w)[:2]
    e=E0(w,ka+kb)
    Ks=np.linspace(-np.pi,np.pi,73); emin=min(E0(w,K,800) for K in Ks)
    print(f"{wa:4.1f} {phi:.4f} {eb-ea:9.5f} {e-ea-eb:10.5f} {e-2*ea:10.5f} {emin-2*ea:10.5f}")
