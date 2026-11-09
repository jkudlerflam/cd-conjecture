"""C186: particle-number content of the dressed single particle in the mixed second-order model (truncation check)."""
import numpy as np, scipy.sparse.linalg as spla, sys
from mixed_fewbody import Model, bare_valleys
for (u, v) in [(1.0, 0.25), (1.0, 0.5), (2.0, 1.0), (2.0, 2.0), (0.5, 0.0625)]:
    m = Model(u, np.pi/2, v)
    ka, kb = bare_valleys(m)[:2]
    m.mu = float(m.eps(ka))    # bare band bottom
    for k, lab in ((ka, 'a'), (kb, 'b')):
        for Nmax in (1, 2, 3, 4, 5):
            R1 = 40 if Nmax <= 4 else 26
            H, basis = m.build(k, Nmax, R1)
            if H.shape[0] < 400:
                e, V = np.linalg.eigh(H.toarray()); e0, x = e[0], V[:, 0]
            else:
                e, V = spla.eigsh(H, k=1, which='SA', tol=1e-11); e0, x = e[0], V[:, 0]
            Ns = np.array([len(c) for c in basis])
            wts = [float(np.sum(abs(x[Ns == N])**2)) for N in range(1, Nmax+1)]
            print(f"u={u} v={v} r={m.r:.3f} valley {lab} k={k:.4f} Nmax={Nmax} R1={R1} e-eps={e0:.6f}  N-weights " +
                  " ".join(f"{w_:.4f}" for w_ in wts), flush=True)
