"""C160 (Claude Opus agent, 2026-11-08). Check of the refined Lemma 3 (PROOF section 3.1).

For a uniform polarized reference along n (local frame u, v, n) of the qubit bond (J, B), the
linear-spin-wave Hamiltonian is H = sum_k [A_k a_k^+ a_k + (g cos k a_k a_-k + h.c.)] with
  A_k = A0 + 2 Re(t) cos k - 2 Im(t) sin k,  A0 = -4 J'_nn - 2 B'_n,
  t = J'_uu + J'_vv + i (J'_vu - J'_uv),     g = J'_uu - J'_vv - i (J'_uv + J'_vu),
and omega(k) = -2 Im(t) sin k + sqrt(Abar_k^2 - 4|g|^2 cos^2 k), Abar_k = A0 + 2 Re(t) cos k.
Stationarity of n fixes B'_u, B'_v = -(J'_un + J'_nu), -(J'_vn + J'_nv).

Claims tested:
 (i)  a_s = 0 (J'_uv = J'_vu, no DM along n) and energetic stability (Abar_k >= 2|g||cos k|) imply that
      every minimum of omega(k) lies at k = 0 or pi (and no soft pair +-k0 occurs).
 (ii) the BdG spectrum of a ring (direct diagonalization of the 2L x 2L bosonic problem) agrees with the
      formula, so the coefficients are right.
"""
import numpy as np

rng = np.random.default_rng(1)


def coeffs(Jp, Bn):
    A0 = -4 * Jp[2, 2] - 2 * Bn
    t = Jp[0, 0] + Jp[1, 1] + 1j * (Jp[1, 0] - Jp[0, 1])
    g = Jp[0, 0] - Jp[1, 1] - 1j * (Jp[0, 1] + Jp[1, 0])
    return A0, t, g


def omega(k, A0, t, g):
    Ab = A0 + 2 * t.real * np.cos(k)
    disc = Ab ** 2 - 4 * abs(g) ** 2 * np.cos(k) ** 2
    return -2 * t.imag * np.sin(k) + np.sqrt(np.maximum(disc, 0)), Ab, disc


def bdg_ring(A0, t, g, L):
    """Direct bosonic BdG on a ring, H = sum_j [A0 a+a + (t a_j^+ a_j+1 + h.c.) + (g a_j a_j+1 + h.c.)]."""
    Amat = np.zeros((L, L), complex); Bmat = np.zeros((L, L), complex)
    for j in range(L):
        Amat[j, j] += A0
        Amat[j, (j + 1) % L] += t; Amat[(j + 1) % L, j] += np.conj(t)
        Bmat[j, (j + 1) % L] += g; Bmat[(j + 1) % L, j] += g     # g a_j a_j+1 written symmetrically
    # H = 1/2 (a^+, a) [[A, B^*],[B, A^*]] (a, a^+)^T + const, with the pair term (1/2) a B a + h.c.
    M = np.block([[Amat, np.conj(Bmat)], [Bmat, np.conj(Amat)]])
    eta = np.diag(np.r_[np.ones(L), -np.ones(L)])
    ev = np.linalg.eigvals(eta @ M)
    return np.sort(ev.real[ev.real > -1e-9])


# (ii) coefficient check
Jp = rng.normal(size=(3, 3)); Bn = rng.normal()
A0, t, g = coeffs(Jp, Bn); A0 += 4 * abs(t) + 4 * abs(g) + 1.0   # shift into the stable region
L = 24
ks = 2 * np.pi * np.arange(L) / L
w, Ab, disc = omega(ks, A0, t, g)
if np.all(disc > 0) and np.all(Ab > 0):
    print('BdG ring check, max |formula - direct| =', np.max(np.abs(np.sort(w) - bdg_ring(A0, t, g, L)[:L])))

# (i) random a_s = 0 references
kk = np.linspace(-np.pi, np.pi, 2001)
nstable = 0; bad = 0; soft_interior = 0
for trial in range(30000):
    Jp = rng.normal(size=(3, 3)); Jp[0, 1] = Jp[1, 0]      # no DM along n
    Bn = 3 * rng.normal()
    A0, t, g = coeffs(Jp, Bn)
    Ab = A0 + 2 * t.real * np.cos(kk)
    if np.any(Ab < 2 * abs(g) * np.abs(np.cos(kk)) - 1e-12):
        continue                                              # not energetically stable
    nstable += 1
    w = np.sqrt(np.maximum(Ab ** 2 - 4 * abs(g) ** 2 * np.cos(kk) ** 2, 0))
    i = np.argmin(w)
    kmin = abs(kk[i])
    if min(kmin, abs(np.pi - kmin)) > 2e-3 and w[i] < min(w[0], w[1000]) - 1e-9:
        bad += 1
    if w[i] < 1e-3 and min(kmin, abs(np.pi - kmin)) > 0.05:
        soft_interior += 1
print(f'a_s = 0: {nstable} energetically stable references, interior global minima: {bad}, interior soft points: {soft_interior}')

# near h_F: J = -1 + eps M, a_s = 0, minimum location and gap scale
for eps in (0.3, 0.1, 0.03):
    worst = 0.0
    for trial in range(3000):
        M = rng.normal(size=(3, 3)); M[0, 1] = M[1, 0]
        Jp = -np.eye(3) + eps * M
        Bn = eps * rng.normal()
        A0, t, g = coeffs(Jp, Bn)
        Ab = A0 + 2 * t.real * np.cos(kk)
        if np.any(Ab < 2 * abs(g) * np.abs(np.cos(kk)) - 1e-12):
            continue
        w = np.sqrt(np.maximum(Ab ** 2 - 4 * abs(g) ** 2 * np.cos(kk) ** 2, 0))
        worst = max(worst, abs(kk[np.argmin(w)]))
    print(f'near h_F eps={eps}: largest |k_min| over stable a_s=0 references = {worst:.4f}')
