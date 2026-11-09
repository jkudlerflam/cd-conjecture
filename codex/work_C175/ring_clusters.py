"""Independent ring check of the cluster energies (N = 2, 3, 4) at |w| = 1 and 2, phi = pi/2, large rings."""
import numpy as np
from ring_k import levels
from scat import valleys
for wa,cases in ((1.0,((2,60),(3,45),(3,60),(4,48),(4,60))),(2.0,((2,60),(3,45),(3,60),(4,40),(4,52)))):
    w=1j*wa; em=valleys(w)[0][1]
    for N,L in cases:
        E=min((levels(L,N,w,m,k=1)[0],m) for m in range(L))
        print(f"|w|={wa} N={N} L={L}: min_K E0 - N eps_min = {E[0]-N*em:.6f} at m={E[1]} (K={2*np.pi*E[1]/L:.4f})",flush=True)
