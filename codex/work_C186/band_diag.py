"""C186: dressed one-particle band of the mixed model on an N-capped ring, all grid momenta (diagnostic)."""
import numpy as np, sys
from ring_mixed import Ring, model
L, Nmax, u, v, phi, mu = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
ring = Ring(L, Nmax)
m = model(u, phi, v, mu)
for i in range(L):
    q = 2*np.pi*i/L
    e, Nx = ring.lowest(m, q, k=3)
    eps = 2*u*np.cos(q-phi) - 2*np.cos(2*q) - mu
    print(f"q={q:.3f} bare={eps:+.4f} levels " + "  ".join(f"{a:+.4f}(N={b:.2f})" for a, b in zip(e, Nx)), flush=True)
