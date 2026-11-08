"""Checks: (1) ED ground energy at Jz = 0 against the free BdG ground energy on rings (both parity sectors with the
correct JW boundary condition); (2) timing at L = 16, 20."""
import numpy as np, time
from model import *

def free_E0_sector(Jx, Jy, D, hz, L, par):
    # fermion ring: parity even (number of down spins even) -> antiperiodic fermions, odd -> periodic.
    # Build the BdG matrix in real space with boundary sign s, and take the ground state in the given parity sector.
    s = -1.0 if par == 0 else 1.0
    w = Jx + Jy + 2j * D; dl = Jx - Jy
    A = np.zeros((L, L), complex); B = np.zeros((L, L), complex)
    for j in range(L):
        k = (j + 1) % L; sg = s if k == 0 else 1.0
        A[k, j] += w * sg; A[j, k] += np.conj(w) * sg   # w c+_{j+1} c_j + h.c.
        # -(dl)(c_j c_k + h.c.)  ->  (1/2) sum c_a B_ab c_b form: B antisymmetric
        B[j, k] += -dl * sg; B[k, j] -= -dl * sg
        A[j, j] += -2 * hz
    # H = sum A_ab c+_a c_b + sum_{a<b}... write H = 1/2 Psi+ Hbdg Psi + const with Psi = (c, c+)
    # pairing term P = sum_{j} B[j,k] c_j c_k (one per bond, written as 1/2 sum_ab B_ab c_a c_b) + h.c.
    Hb = np.block([[A, B.conj().T * 0 + (-B.conj())], [ -B * 0 + B, -A.conj()]])
    # careful construction: H = sum A_ab c+a cb + 1/2 sum (B_ab c_a c_b + h.c.)
    Hb = np.block([[A, -B.conj()], [B, -A.conj()]])  # placeholder, verified numerically below
    e = np.linalg.eigvalsh(Hb)
    ep = np.sort(e)[L:]
    E_vac = -0.5 * np.sum(ep) + 0.5 * np.trace(A).real
    return E_vac, ep

L = 10
rng = np.random.default_rng(1)
for trial in range(4):
    Jx, Jy, D, hz = rng.normal(size=4)
    h = bond(Jx, Jy, D, 0.0, hz)
    for par in (0, 1):
        H = ham(h, L, par)
        w = np.linalg.eigvalsh(H.toarray())[:1]
        # spin constant: ZZ term zero; hz Z = hz (1 - 2n): constant L*hz
        Ef, ep = free_E0_sector(Jx, Jy, D, hz, L, par)
        print(trial, par, round(w[0], 8), round(Ef + L * hz, 8))
for L in (16, 20):
    t = time.time(); r = ring_data(bond(1.2, 0.8, 0.3, 0.7, 0.4), L); print(L, r['E0'], r['S'], round(time.time() - t, 1), 's')
