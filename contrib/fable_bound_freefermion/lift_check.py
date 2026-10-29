"""Numerical check of the lifted J-orthogonality proof of  sum |v| <= 4 ||B||_1  (Theorem 4).
For random chains: build x_a = (u_a, z_a u_a)/sqrt2, J = [[0, iB], [-iB^dag, 0]];
verify x_b^dag J x_a = (v_a/2) delta_ab and sum_a |x_a^dag J x_a| <= tr|J| = 2 ||B||_1."""
import numpy as np
from ff_bound import fermi_points, random_nc, random_mj, quad_coeffs, trace_norm

rng = np.random.default_rng(9)
for kind in ("nc", "mj"):
    worst_off, worst_ratio = 0.0, 0.0
    for n in (1, 2, 3, 4):
        for _ in range(300):
            A0, B = random_nc(n, rng) if kind == "nc" else random_mj(n, rng)
            Bc = B if kind == "nc" else 1j * B            # complex hopping matrix of h(k) = A0 + z Bc + z^* Bc^dag
            J = np.block([[np.zeros((n, n)), 1j * Bc], [-1j * Bc.conj().T, np.zeros((n, n))]])
            xs = []
            for k, U in fermi_points(kind, A0, B):
                z = np.exp(1j * k)
                dh = 1j * (Bc * z - Bc.conj().T / z)
                V = U.conj().T @ dh @ U
                w, R = np.linalg.eigh((V + V.conj().T) / 2)
                Ud = U @ R                                 # kernel basis diagonalizing the velocity matrix
                for i in range(Ud.shape[1]):
                    xs.append(np.concatenate([Ud[:, i], z * Ud[:, i]]) / np.sqrt(2))
            if not xs:
                continue
            X = np.array(xs).T
            Gram = X.conj().T @ J @ X
            off = np.abs(Gram - np.diag(np.diag(Gram))).max()
            worst_off = max(worst_off, off)
            worst_ratio = max(worst_ratio, np.abs(np.diag(Gram)).sum() / (2 * trace_norm(Bc)))
    print(f"{kind}: max off-diagonal |x_b^dag J x_a| = {worst_off:.2e};  max sum|d_a| / tr|J| = {worst_ratio:.5f} (must be <= 1)")
