import numpy as np
from ringk import reps, build_k
from ringv import build
for (L, N, w, Phi) in ((14, 3, 1j, 0.0), (15, 4, 0.7 + 0.4j, 0.37), (16, 4, 1j, 1.1)):
    b, H = build(L, N, w, Phi); er = np.sort(np.linalg.eigvalsh(H.toarray()))
    data = reps(L, N); ek = np.sort(np.concatenate([np.linalg.eigvalsh(build_k(L, N, 2*np.pi*m/L, w, Phi=Phi, data=data)[0].toarray()) for m in range(L)]))
    print(f"L={L} N={N} w={w} Phi={Phi}: max |diff| = {np.abs(er-ek).max():.2e}")
