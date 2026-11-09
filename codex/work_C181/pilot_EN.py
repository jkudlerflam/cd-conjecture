"""C181 design pilot: E(N) on a short open chain, to locate the density range (not a physics test)."""
import numpy as np, sys
from fibdmrg import ground
wa=float(sys.argv[1]); phi=float(sys.argv[2]); L=int(sys.argv[3]); Ns=[int(x) for x in sys.argv[4].split(',')]
w=wa*np.exp(1j*phi)
for N in Ns:
    dm,h=ground(L,N,w,[16,32,48,64,64,96,96],ch_kw=dict(gauge=1))
    d=np.array(dm.measure()['dens']) if False else None
    print(L,N,h[-1][1],h[-1][2],flush=True)
