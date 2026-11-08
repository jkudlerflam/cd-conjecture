"""Open-chain entanglement profile S(l) of the free (V = 0) parity-class chain from the BdG ground state, and the same
central-window c fits as dmrg_pt.py. Fermion units: H = sum [t e^{i phi} c+_{j+1} c_j + h.c. - Delta (c_j c_{j+1} + h.c.) - mu n_j].
Usage: python3 free_open.py phi s mu L1 L2 ...   (t = 1, Delta = s sin phi)."""
import sys
import numpy as np
from dmrg_pt import cfit
def profile(phi, Dl, mu, L):
    w = np.exp(1j * phi)
    A = np.zeros((L, L), complex); B = np.zeros((L, L), complex)
    for j in range(L - 1):
        A[j + 1, j] += w; A[j, j + 1] += np.conj(w)
        B[j, j + 1] += -Dl; B[j + 1, j] += Dl           # H = c+ A c + (1/2)(c B c + h.c.) convention
    A -= mu * np.eye(L)
    H = np.block([[A, B], [-B.conj(), -A.conj()]])     # BdG in (c, c+) basis, H = (1/2) Psi+ H Psi
    e, U = np.linalg.eigh(H)
    V = U[:, e < 0]                                      # occupied quasiparticle modes of the BdG ground state
    G = V @ V.conj().T                                   # <Psi Psi+> structure
    S = []
    for l in range(1, L):
        idx = list(range(l)) + list(range(L, L + l))
        g = np.linalg.eigvalsh(G[np.ix_(idx, idx)])
        g = np.clip(g, 1e-14, 1 - 1e-14)
        S.append(float(-0.5 * np.sum(g * np.log(g) + (1 - g) * np.log(1 - g))))
    return np.array(S)
if __name__ == "__main__":
    phi, s, mu = (float(x) for x in sys.argv[1:4]); Dl = s * np.sin(phi)
    for L in (int(x) for x in sys.argv[4:]):
        S = profile(phi, Dl, mu, L)
        print(f"phi={phi} s={s} mu={mu} L={L}: c25 = {np.round(cfit(S, L, .25, .75), 3)}, c12 = {np.round(cfit(S, L, .125, .875), 3)}, Smid = {S[L//2-1]:.4f}")
