"""C176 validation, part 2: correlators Cnn and Cpp against ED, and the c = 2 control at even N."""
import numpy as np
from fibdmrg import ground
from validate import ed
from eff_ed import build
cases = [(20, 4, 1.0, np.pi/2, {}, {}, 1), (18, 4, 1.0, 1.45, {}, {}, 0),
         (16, 4, 0.0, 0.0, dict(V2=0.0, V3=0.0, blockade=False), dict(V2=0.0, V3=0.0, blockade=False), 0)]
for (L, N, wa, phi, ckw, ekw, g) in cases:
    w = wa*np.exp(1j*phi)
    b, e0, psi, e1 = ed(L, N, w, **ekw)
    from eff_ed import entropies
    Sed = entropies(b, psi, L)
    dm, hist = ground(L, N, w, chis=[16, 32, 64, 128], ch_kw=dict(gauge=g, **ckw), tol=1e-12)
    i0 = L//2 - 2
    m = dm.measure(gwin=(0, L), refs=[i0])
    nj = lambda j: ((b >> j) & 1).astype(float)
    p2 = np.abs(psi)**2
    dens = np.array([np.sum(p2*nj(j)) for j in range(L)])
    cnn = np.array([np.sum(p2*nj(i0)*nj(j)) - dens[i0]*dens[j] for j in range(i0+1, L)])
    # P_j = c+_{j+1}c_j + h.c. as a sparse operator in the ED basis (adjacent hop, no sign)
    import scipy.sparse as sp
    def Pop(j):
        rows, cols = [], []
        for (a, c) in ((j, j+1), (j+1, j)):
            ok = (((b >> a) & 1) == 1) & (((b >> c) & 1) == 0)
            src = np.nonzero(ok)[0]
            ts = (b[src] & ~(1 << a)) | (1 << c)
            pos = np.searchsorted(b, ts); pos = np.minimum(pos, len(b)-1); good = b[pos] == ts
            rows += list(pos[good]); cols += list(src[good])
        return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(b), len(b)))
    Pexp = np.array([np.real(np.vdot(psi, Pop(j) @ psi)) for j in range(L-1)])
    Pi = Pop(i0) @ psi
    cpp = np.array([np.real(np.vdot(Pop(j) @ psi, Pi)) - Pexp[i0]*Pexp[j] for j in range(i0+2, L-1)])
    print(f"L={L} N={N} phi={phi:.3f} {ckw}: dE={hist[-1][1]-e0:+.1e} gap={e1-e0:.2e} max|dS|={np.max(np.abs(dm.S-Sed)):.1e} "
          f"max|dCnn|={np.max(np.abs(np.array(m['Cnn'][str(i0)])-cnn)):.1e} max|dCpp|={np.max(np.abs(np.array(m['Cpp'][str(i0)])-cpp)):.1e} "
          f"max|Cpp|={np.max(np.abs(cpp)):.2e}")
    print("   per-j |dCpp|:", np.round(np.abs(np.array(m["Cpp"][str(i0)])-cpp), 12)[:6])
