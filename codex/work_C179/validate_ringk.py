"""Validate ringk.py against C173's real-space ring_ed.py (bosonic seam sign) at small L."""
import sys, numpy as np
sys.path.insert(0, '../work_C173')
from ring_ed import build_ring
from ringk import reps, build_k
for (L, N, w) in ((14, 3, 1j), (15, 4, 0.7 + 0.4j), (16, 4, 1j)):
    b, H = build_ring(L, N, w, bc=(-1) ** (N - 1))
    e_real = np.sort(np.linalg.eigvalsh(H.toarray()))
    data = reps(L, N); ek = []
    herm = 0
    for m in range(L):
        Hk, R, P = build_k(L, N, 2 * np.pi * m / L, w, data=data)
        A = Hk.toarray(); herm = max(herm, np.abs(A - A.conj().T).max())
        ek += list(np.linalg.eigvalsh(A))
    ek = np.sort(ek)
    print(f"L={L} N={N} w={w}: dim {len(e_real)} vs {len(ek)}, max |diff| = {np.abs(e_real-ek).max():.2e}, max non-hermiticity {herm:.1e}")
