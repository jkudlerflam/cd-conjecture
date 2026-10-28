# Route A at finite coupling: chains of identical second-order exchange cells

## Claim

Each site is one identical 48-dimensional C8 cell (`codex/work_C8/two_port_N3.npz`, C8 §4.2, cell gap 1). The chain is translation invariant, nearest neighbour and Hermitian on d = 48,

H = Σ_j H0_j + λ Σ_j Σ_a F_{R,a}(j) F_{L,a}(j+1),

with F the Paulis on the right and left qubit ports. At second order it is the SU(3) Uimin–Lai–Sutherland (ULS) chain b Σ P_{j,j+1}, b = 27λ²/8192, with IR SU(3)_1 and c = 2. We studied λ = 0.02 to 0.2.

1. The logical chain at finite λ is not SU(3) symmetric. The breaking starts at order λ³ (relative size about λ). Its two-body part is SO(3) symmetric to a few percent (the logical qutrit is a spin 1), with dipole exchange above quadrupole exchange, which puts the chain at the bilinear–biquadratic angle θ = π/4 − (0.06 to 0.08)πλ, on the Haldane side of ULS. Each cell also feels a uniform logical field with three distinct eigenvalues of order 0.4λb to 0.75λb, and the three-cell terms are 3λ²b to 4.5λ²b.
2. Up to correlation lengths of about 30 cells, the chain is indistinguishable from c = 2 for λ ≤ 0.1, by two estimators calibrated on the Heisenberg and ULS chains. z = 1 holds, and there is no dimerization.
3. SU(3)_1 is not the IR fixed point of this family. At λ = 0.2 the finite-size spectrum drifts away from SU(3)_1 with L, and amplifying the measured breaking by 8 or 16 gives a uniform gapped state. Whether the physical couplings end gapped or in a critical window is beyond the reachable scales. This fits the separatrix picture of SU(N)_1.
4. Question (3). With an incommensurate free-fermion target (XX chain in a field, c = 1) the chain is gapped at every coupling tested, ξ ≈ 1.25/λ cells, because the cell generates a strictly relevant pairing γ ≈ 0.4λ. The IR has c = 0.

## Status

Numerical. An RG integration of the measured couplings would give the IR at λ ≤ 0.1, which is far beyond reach numerically.

## How it was established and checked

- Author(s): Claude (Opus 5.5) agent, 2026-10-27, for Jonah Kudler-Flam.
- AI assistance: all code and runs by Claude. No Codex.
- Checked by: internal cross-checks only (third-order perturbation theory, calibrations). Not independently rerun.

## Argument or method

**Exact logical couplings** (`sw_cluster.py`, `sw_table.py`, `pt3_check.py`). For open clusters of 2 and 3 cells we take the 3^n lowest eigenvectors of the complete Hamiltonian and form the des Cloizeaux (direct-rotation Schrieffer–Wolff) effective Hamiltonian on the code. It is linked-cluster additive, so h2 = H_eff(2) and h3 = H_eff(3) − h2⊗1 − 1⊗h2 are the bond and three-cell terms of the logical chain. The band stays isolated (next level 0.74 or higher). b/λ² → 27/8192, and the λ³ term agrees with third-order perturbation theory to 5 × 10⁻⁶ relative [checked].

Couplings in units of b (g_D, g_Q are the dipole and quadrupole exchange, "other" the remaining SO(3)-breaking two-body part, and the field collects the one-cell parts of h2 and h3):

| λ | b/λ² | (g_D − g_Q)/b | θ/π | field eigenvalues / b | other / b | ‖h3‖/b |
| --: | --: | --: | --: | :-- | --: | --: |
| 0.02 | 0.003370 | 0.016 | 0.2488 | −0.007, −0.001, 0.008 | 0.001 | 0.001 |
| 0.05 | 0.003482 | 0.041 | 0.2467 | −0.019, −0.002, 0.021 | 0.003 | 0.008 |
| 0.1 | 0.003670 | 0.087 | 0.2430 | −0.044, −0.005, 0.048 | 0.006 | 0.036 |
| 0.2 | 0.004040 | 0.199 | 0.2337 | −0.121, −0.030, 0.152 | 0.012 | 0.181 |

The three-cell term is mostly one-cell. The cell's QR completion commutes with no SO(3) action, yet the λ³ exchange comes out SO(3) symmetric.

**Estimators** (`effchain.py`, `analyze.py`, reusing `contrib/d2_followup_220_1282/vumps.py`). (i) Finite-entanglement scaling (FES), VUMPS with a six-site cell, c = 6 dS/d ln ξ over the last four χ. (ii) Periodic ED (L = 6 to 15) with momenta, the Casimir fit E0/L = e∞ − πcv/(6L²) on consecutive sizes with v from the k = 2π/L level, and the ratio x(2π/3)/x(2π/L), which is 2/3 in SU(3)_1. Calibration: FES gives 1.00 for Heisenberg (ξ = 623) and 2.02 for ULS (ξ = 30); ED gives 1.007 for Heisenberg (L = 18, 20) and 2.09 for ULS (L = 9, 12).

**Results.** Rows marked h3 include the three-cell term, the others use h2 only. ED c is the fit on L = 9 and 12.

| chain | FES c (χ ≤, ξ) | ED c | x-ratio, L = 9, 12 (15) | v, L = 12 |
| :-- | :-- | --: | :-- | --: |
| ULS | 2.02 (128, 30) | 2.09 | 0.652, 0.642 | 2.08 |
| λ = 0.02 | 2.12 (96, 23) | 2.11 | 0.649, 0.638 | 2.07 |
| λ = 0.05 | 2.05 (96, 24) | 2.14 | 0.644, 0.633 (0.628) | 2.04 |
| λ = 0.1 | 2.02 (128, 31) | 2.18 | 0.635, 0.623 | 2.00 |
| λ = 0.2 | 1.97 (128, 32) | 2.28 | 0.612, 0.597 (0.588) | 1.91 |
| λ = 0.02, h3 | | 2.11 | 0.649, 0.638 | 2.07 |
| λ = 0.1, h3 | 1.98 (64, 16) | 2.20 | 0.633, 0.619 | 1.98 |
| λ = 0.2, h3 | 1.95 (64, 16) | 2.40 | 0.600, 0.575 | 1.79 |

Over the same χ range (24 to 64) ULS gives 2.03. Up to λ = 0.1 both estimators match the ULS calibration. At λ = 0.2 the Casimir estimate is 0.2 to 0.3 above ULS and the x-ratio falls steadily with L, whereas in ULS it sits just below 2/3 (logarithmic correction). Since FES still gives about 2, we read this as a flow away from the fixed point rather than a different central charge.

Bond energies in the six-site cell show a period-3 modulation of 1 to 2 × 10⁻² at χ = 128 (ULS 3 × 10⁻³), decreasing with χ, and no period-2 modulation. Velocities are L independent to 3 percent (z = 1). The BLBQ chain at the angle of λ = 0.2 (θ = 0.2337π) also looks like c = 2 here (FES 1.94, ξ = 32), while θ = 0.2π is gapped (ξ = 9.8), so the Haldane length of the SO(3)-symmetric part alone exceeds the reachable scale.

**Amplified breaking.** We replaced the λ = 0.1 logical chain (with h3) by (P − 1/3) + κ[h_eff − (P − 1/3)]. At κ = 8 and 16 the state is gapped and translation invariant (ξ = 5.9 and 4.9, S saturated to 10⁻³, no modulation in a six-site cell). At κ = 2 and 4 it still looks critical at the χ reached (local FES slopes 1.9 and 1.8 at ξ = 16 and 13).

**Free-fermion target** (`part3_ff.py`). `cell.build_twoport(K, 2)` builds an identical two-qubit-port cell (D = 32, C8 geometry with G1's syndrome trick, energies 1 and 3 only; C9's recipe with m + 2 energies is ill conditioned, its Cauchy matrix reaching eigenvalues 10⁻¹⁷) with second-order bond τ²K, K = XX + YY + 0.6(Z1 + 1Z) [checked to 2 × 10⁻¹⁵], so k_F = 0.295π. At finite λ the logical chain acquires the pairing γ = (J_XX − J_YY)/(J_XX + J_YY) = 0.0072, 0.019, 0.040, 0.089 and transverse fields up to 0.07. FES (two-site cell, χ ≤ 128):

| λ | 0.02 | 0.05 | 0.1 | 0.2 |
| :-- | --: | --: | --: | --: |
| ξ (cells) | 57 | 25 | 12.7 | 6.5 |
| S (saturated) | 0.672 | 0.518 | 0.401 | 0.289 |
| v/(2Δ), pairing alone | 69 | 27 | 12.6 | 5.6 |

The last row uses the XY-chain gap Δ = min_k 2[(J cos k − h)² + J²γ² sin²k]^{1/2} (J = J_XX + J_YY, h the field per site) and agrees with ξ within 20 percent. The U(1)-symmetric part alone stays critical (FES 0.98 to 0.99, as for the bare target).

**Not done.** (i) Direct VUMPS on the 48-dimensional sites (`bigsite.py`, kept for reference). The logical scale is 10⁻⁴ to 10⁻⁶ of the cell gap and the preconditioned Davidson solves did not converge within budget, so the big-site statements rest on the exact cluster couplings, whose expansion converges like λ². A four-cell exact check of the omitted connected term (`sw4_check.py`, job 17513021, output `data/sw4_λ.json`) was still running when this was written. (ii) Open-chain entropy fits were not converged in χ for SU(3) (and give 0.92 for Heisenberg at L = 192), so they are not used.

## Reproducing

numpy and scipy only. `cell.py` checks the cell. `sw_cluster.py --n3 --lams …` gives the couplings (minutes per λ), `sw_table.py` the table, `pt3_check.py` the third-order check, `sw4_check.py λ` the four-cell check (32 GB, hours). `effchain.py MODEL {fes,ed,obc}` runs the estimators (models `heis`, `uls`, `blbq:T`, `eff:λ`, `eff3:λ`, `amp:λ:κ`, `file:…`; FES to χ = 128 takes hours at d = 3), and `analyze.py` prints every estimate (`analyze_stdout.txt`). `part3_ff.py` builds the free-fermion cell. Slurm scripts `sw3.slurm`, `sw4.slurm`, `eff.slurm` with `tasks_*.txt`. Typhon jobs 17512603, 17512641 (couplings), 17512642, 17512643, 17513016, 17513017, 17513079 (estimators), 17513021 (four-cell check). The amplified scan, the h3 runs and the free-fermion FES ran on the Typhon login node at low priority.

## Relation to the rest of the repository

This tests REPORT §9.1–9.3 (G1, C8 §4, C9 §2) at finite coupling, alongside route B (`contrib/encoded_chain_numerics/`) and the incommensurate multi-flavour target (`contrib/routeA_incommensurate/`). The second-order mechanism survives finite coupling quantitatively (b/λ² is 2 percent above 27/8192 at λ = 0.02 and 23 percent above at λ = 0.2), but the cell deforms the target by O(λ) in directions it does not control. SU(3)_1 is not protected against them, and a single-flavour free fermion is gapped by them. This supports the switch to incommensurate multi-flavour targets with an exact logical fermion parity (REPORT §9.9). It does not change "Where things stand".
