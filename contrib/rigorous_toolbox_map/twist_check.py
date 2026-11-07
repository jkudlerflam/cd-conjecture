import numpy as np
from itertools import product
def sector_H(L,n,Delta,Phi,J=1.0):
    states=[s for s in range(2**L) if bin(s).count('1')==n]
    idx={s:i for i,s in enumerate(states)}
    H=np.zeros((len(states),)*2,dtype=complex)
    for i,s in enumerate(states):
        for x in range(L):
            y=(x+1)%L
            bx=(s>>x)&1; by=(s>>y)&1
            H[i,i]+=J*Delta*(0.25 if bx==by else -0.25)
            if bx!=by:
                t=s^(1<<x)^(1<<y)
                ph=np.exp(1j*Phi*(1 if bx==1 else -1)) if y==0 else 1.0
                H[idx[t],i]+=0.5*J*ph
    H=(H+H.conj().T)/2  # symmetrize (already Hermitian up to convention)
    return np.linalg.eigvalsh(H)[0]
worst=0
for L in [4,5,6,7,8]:
  for n in range(0,L//2+1):
    for Delta in [-0.9,0,0.5,1,2]:
      Ph=np.linspace(0,2*np.pi,121)
      E=np.array([sector_H(L,n,Delta,p) for p in Ph])
      r=(E.max()-E.min())/(np.pi**2/(4*L))
      worst=max(worst,r)
print("max ratio (E(Phi)-minE)/(pi^2 J/(4L)) =",worst)
