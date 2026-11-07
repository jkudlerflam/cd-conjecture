"""Rank of the frustration-free projector at KTM points: eigenvalues of the gauge-optimized bond h + o.(s1 - s2)
above its minimum; for rank one, the Bravyi-Gosset T matrix of the excited eigenvector psi."""
import numpy as np
from ed import bond, gauge_lmin, PAULI, I2, spectrum, all_levels
def T_of(psi):
    c = np.conj(psi); return np.array([[c[1], c[3]], [-c[0], -c[2]]])
for (jx, jy, jz) in [(0.6, 0.2, 0.3), (0.4, -0.2, 0.3), (0.3, -0.25, 0.3), (0.5, 0.5, 0.3), (0.6, 0.0, 0.3)]:
    b = 2 * np.sqrt((jx + jz) * (jy + jz))
    h = bond(np.diag([jx, jy, jz]), (0, 0, b))
    lmin, o = gauge_lmin(h)
    G = [np.kron(P, I2) - np.kron(I2, P) for P in PAULI]
    hg = h + sum(o[a] * G[a] for a in range(3))
    w, v = np.linalg.eigh(hg)
    rank = int(np.sum(w - w[0] > 1e-6))
    msg = f"jx={jx} jy={jy} jz={jz}: bond eigenvalues {np.round(w - w[0], 5)}, projector rank {rank}"
    if rank == 1:
        mu = np.linalg.eigvals(T_of(v[:, -1]))
        msg += f", |mu| = {np.round(np.abs(mu), 5)}, arg(mu2/mu1)/pi = {np.angle(mu[1]/mu[0])/np.pi:.4f}"
    print(msg)
# larger rings for the ambiguous first case
jx, jy, jz = 0.6, 0.2, 0.3; b = 2 * np.sqrt((jx + jz) * (jy + jz))
h = bond(np.diag([jx, jy, jz]), (0, 0, b))
for L in [13, 14, 15, 16]:
    w, _ = all_levels(spectrum(h, L, nev=6)); w = np.sort(w)
    print(f"  (0.6,0.2,0.3) L={L}: lowest excitations {np.round(w[1:6]-w[0], 5)}, L*gap3 = {L*(w[2]-w[0]):.3f}", flush=True)
