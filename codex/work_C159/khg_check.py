"""C159: the Kitaev-Heisenberg-Gamma chain (alternating x, y bonds) as a uniform B_2 bond, and its stratum.
Spectra of the alternating chain and of the uniform chain J' = J0 Q are compared on 6-site rings (dense)."""
import numpy as np
import tk
from est import PAULI


def ring_H(Jfun, L):
    N = 2 ** L
    H = np.zeros((N, N), complex)
    for j in range(L):
        J = Jfun(j)
        for a in range(3):
            for b in range(3):
                if J[a, b] == 0:
                    continue
                ops = [np.eye(2)] * L
                ops = list(ops); ops[j] = PAULI[a]; ops[(j + 1) % L] = PAULI[b]
                M = ops[0]
                for o in ops[1:]:
                    M = np.kron(M, o)
                H += J[a, b] * M
    return H


Q = np.zeros((3, 3)); Q[2, 0] = 1; Q[1, 2] = 1; Q[0, 1] = 1
for (K, G, Jh) in [(-0.6, 0.5, 0.3), (0.4, -0.7, 0.25), (-1.0, 0.3, -0.4)]:
    Jx = np.array([[K, 0, 0], [0, 0, G], [0, G, 0]]) + Jh * np.eye(3)
    Jy = np.array([[0, 0, G], [0, K, 0], [G, 0, 0]]) + Jh * np.eye(3)
    best = None
    for J0 in (-np.array([[K + Jh, 0, 0], [0, G, Jh], [0, Jh, G]]), -np.array([[K + Jh, 0, 0], [0, Jh, G], [0, G, Jh]])):
        Ju = J0 @ Q
        eo = np.linalg.eigvalsh(ring_H(lambda j: Jx if j % 2 == 0 else Jy, 6))
        eu = np.linalg.eigvalsh(ring_H(lambda j: Ju, 6))
        d = np.max(np.abs(eo - eu))
        if best is None or d < best[0]:
            best = (d, Ju)
    d, Ju = best
    syms = tk.symmetry_check(Ju, np.zeros(3))
    par = tk.spiral_parent(Ju)
    sv = np.linalg.svd(Ju, compute_uv=False)
    print(f"K={K} G={G} J={Jh}: spectral diff {d:.1e}; J'={np.round(Ju,3).tolist()}; singular values {np.round(sv,3).tolist()};"
          f" uniform symmetries {[(s[1]) for s in syms]}; spiral theta/pi={par['theta']/np.pi:.3f} m={np.round(par['m'],3).tolist()}")
