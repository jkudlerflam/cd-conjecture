"""Calibration of the K_rho estimator of manybody.py on hard rods (blockade, NN hop only, no t2/V2/V3), where the
exact Luttinger parameter in boson convention is K = (1 - n)^2 (exclusion length 1)."""
import numpy as np, scipy.sparse.linalg as spla
from ringv import build
def E(L, N, Phi=0.0):
    b, H = build(L, N, 1.0, Phi, t2=0.0, V2=0.0, V3=0.0)
    return spla.eigsh(H, k=1, which='SA', tol=1e-12)[0][0] if H.shape[0] > 800 else np.linalg.eigvalsh(H.toarray())[0]
for (L, N) in ((24, 4), (36, 6)):
    dP = 0.05
    vK = np.pi * L * (E(L, N, dP) + E(L, N, -dP) - 2 * E(L, N)) / dP ** 2
    vk = (L / np.pi) * (E(L, N + 2) + E(L, N - 2) - 2 * E(L, N)) / 4
    n = N / L
    print(f"hard rods L={L} N={N}: K_rho estimate = {np.sqrt(vK/vk):.3f}, exact (1-n)^2 = {(1-n)**2:.3f}")
