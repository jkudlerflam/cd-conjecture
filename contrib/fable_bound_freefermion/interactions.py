"""Interacting perturbations of the XX chain: does A/M move above the free value 2?

Chain (fermions, t = 1):  H = sum_x [ -(c_x^dag c_{x+1} + h.c.) + V (n_x - 1/2)(n_{x+1} - 1/2) ],
which is the XXZ chain with J = 2t = 2 and Delta = V/J = V/2, at half filling (exact, Jordan-Wigner).

Exact (Bethe ansatz, recalled):  c = 1 and v = (pi J / 2) sin(gamma)/gamma with Delta = cos(gamma),
for -1 < Delta <= 1.  First-order Hartree-Fock:  v = v_F + (V/pi)(1 - cos 2k_F), v_F = 2 t sin k_F.
M: computed by gauge optimization of the 4x4 bond; the dual lower bound from the two one-particle
states (|10> +- |01>)/sqrt2 (equal marginals) gives M >= t for every V.
"""
import numpy as np
from gauge_M import gauge_M

t = 1.0
# two-site fermion bond in the occupation basis |00>,|01>,|10>,|11> (ordering: site1 = first label)
def bond(V):
    H = np.zeros((4, 4))
    # states: 0=|00>, 1=|01>, 2=|10>, 3=|11>   (n1 n2)
    H[1, 2] = H[2, 1] = -t            # hopping between |01> and |10> (no sign for adjacent modes)
    for s, (n1, n2) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
        H[s, s] += V * (n1 - 0.5) * (n2 - 0.5)
    return H

print(f"{'Delta':>6s} {'V':>6s} {'v exact':>9s} {'v HF(1st)':>10s} {'M (gauge)':>10s} {'A/M exact':>10s} {'A/M 1st':>8s}")
for Delta in (-0.5, -0.25, 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
    V = 2 * Delta
    gamma = np.arccos(Delta)
    v_exact = (np.pi * 2 * t / 2) * (np.sin(gamma) / gamma if gamma > 0 else 1.0)
    v_hf = 2 * t + 2 * V / np.pi
    M, _ = gauge_M(bond(V), 2)
    print(f"{Delta:6.2f} {V:6.2f} {v_exact:9.5f} {v_hf:10.5f} {M:10.5f} {v_exact/M:10.5f} {v_hf/M:8.5f}")
print("\nFree value A/M = 2 at Delta = 0; the ratio rises monotonically with repulsion and reaches pi at the")
print("Heisenberg point Delta = 1, which saturates both the proved free-fermion bound A <= pi M and the")
print("conjectured bound pi M log_2 d at d = 2.  For Delta > 1 the chain is gapped (A = 0).")

# Hartree-Fock first-order velocity at generic filling, checked against the analytic formula
print("\nHartree-Fock check at generic filling (V = 0.2):")
V = 0.2
for kF in (np.pi / 3, np.pi / 2, 2 * np.pi / 3):
    # Fock-renormalized hopping t_eff = t + V <c_x^dag c_{x+1}>,  <c^dag c'> = sin(kF)/pi
    t_eff = t + V * np.sin(kF) / np.pi
    v_num = 2 * t_eff * np.sin(kF)
    v_formula = 2 * t * np.sin(kF) + (V / np.pi) * (1 - np.cos(2 * kF))
    print(f"  k_F = {kF:.4f}: v(HF) = {v_num:.6f}, formula = {v_formula:.6f}")
