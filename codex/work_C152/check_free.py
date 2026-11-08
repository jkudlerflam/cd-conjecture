"""Check of the fermion parametrization: at V = 0 the open-chain BdG single-particle gap equals the ED gap, and the
infinite-chain band minimum of E_+ (free_band) is approached. Also ring ED gap vs min over allowed momenta."""
import numpy as np
from model import *
def bdg_open(t, phi, Dl, mu, L):
    w = t * np.exp(1j * phi)
    A = np.zeros((L, L), complex); B = np.zeros((L, L), complex)
    for j in range(L - 1):
        A[j + 1, j] += w; A[j, j + 1] += np.conj(w)
        B[j, j + 1] += -Dl; B[j + 1, j] -= -Dl
    for j in range(L): A[j, j] -= mu if 0 < j < L - 1 else mu / 2   # bond-split field: edge sites carry half
    # H = c+ A c + 1/2 (c B c + h.c.) with B antisym: BdG [[A, B^dag],[B, -A^T]] on (c, c+)
    Hb = np.block([[A, B.conj().T], [B, -A.T]])
    e = np.linalg.eigvalsh(Hb)
    return np.sort(e[e > -1e-12])[:3]
for (t, phi, Dl, mu) in [(1, 0.4, 0.5, 0.7), (1, 1.2, 0.3, 0.2), (1, 0.0, 0.6, 2.5)]:
    p = spin_from_fermion(t, phi, Dl, 0.0, mu)
    L = 12
    lev = []
    for par in (0, 1):
        H = ham(bond(*p), L, par, periodic=False)
        lev += list(np.linalg.eigvalsh(H.toarray())[:3])
    lev = np.sort(lev)
    k = np.linspace(-np.pi, np.pi, 4001)
    print((t, phi, Dl, mu), 'ED open gaps', np.round(lev[1:3] - lev[0], 5), 'BdG modes', np.round(bdg_open(t, phi, Dl, mu, L), 5), 'band min', round(free_band(k, t, phi, Dl, mu).min(), 4))
