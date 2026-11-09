"""C186: the N-capped ring (Ring(L, Nmax)) must equal the full ring when Nmax >= max N, and the full ring restricted to N <= Nmax."""
import numpy as np
from ring_mixed import Ring, model
L = 16
m = model(1.2, 1.5, 0.3, -2.3)
full = Ring(L); cap = Ring(L, 8); cap3 = Ring(L, 3)
for kk in (0, 3, 8):
    K = 2*np.pi*kk/L
    ef = np.linalg.eigvalsh(full.H(m, K)[0].toarray())
    ec = np.linalg.eigvalsh(cap.H(m, K)[0].toarray())
    Hf, Nf = full.H(m, K)
    Hf = Hf.toarray()[np.ix_(Nf <= 3, Nf <= 3)]
    e3 = np.linalg.eigvalsh(cap3.H(m, K)[0].toarray())
    print(kk, len(ef), len(ec), np.abs(ef-ec).max(), "| cap3 vs restricted full:", len(e3), np.abs(np.sort(np.linalg.eigvalsh(Hf)) - e3).max())
