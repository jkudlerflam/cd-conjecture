"""C185: three-body bound states of the cusp model with V2 = -2 + v, V3 = -2/3, phi = pi/2, via C175's
nbody_box.build (unmodified). Threshold for v above the pair-binding boundary is 3 eps_min (no bound pair).
Prints min over K of E3 - 3 eps_min at box R, and the value at the minimizing K for a larger R."""
import sys, numpy as np, scipy.sparse.linalg as spla
sys.path.insert(0, "../work_C175")
from nbody_box import build
from scat import valleys
def low(N, K, w, R, V2, V3):
    H = build(N, K, w, R, V2=V2, V3=V3)
    return spla.eigsh(H, k=1, which='SA', return_eigenvectors=False)[0]
for u, R1, R2 in ((1.0, 90, 140), (2.0, 70, 110)):
    w = 1j*u; em = valleys(w)[0][1]
    for v in [float(x) for x in sys.argv[1].split(',')]:
        V2, V3 = -2+v, -2/3
        Ks = np.linspace(-np.pi, np.pi, 33)
        vals = [low(3, K, w, R1, V2, V3) - 3*em for K in Ks]
        i = int(np.argmin(vals)); K = Ks[i]
        e2 = low(3, K, w, R2, V2, V3) - 3*em
        print(f"u={u} v={v}: min_K(E3-3eps) R={R1}: {vals[i]:+.5f} at K={K:+.3f}; R={R2}: {e2:+.5f}", flush=True)
