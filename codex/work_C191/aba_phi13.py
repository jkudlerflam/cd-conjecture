"""C191 check: one-valley ABA at phi = 1.3 (only the lower valley occupied) against C179's ring K_rho."""
import numpy as np
from scat191 import Model
from aba_valley import ABA
m = Model(w=np.exp(1.3j))
vs = m.valleys()
print("valleys", [(k, m.eps(k)) for k in vs])
ka = vs[0]
aba = ABA(m, ka, 1.0, npts=141)
for rho in (2/15, 1/6, 0.2):
    for L in (240, 480):
        N = int(round(rho*L))
        K, u = aba.KU(L, N)
        print(f"phi=1.3 rho={rho:.4f} L={L}: ABA K={K:.4f} u={u:.4f}")
