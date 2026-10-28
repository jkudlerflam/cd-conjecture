# Route A at finite coupling, with the incommensurate two-flavour target

## Claim

Class: route A of STATUS step 3, a chain of identical big cells (one site each, D = 64), H = Σ_j H0_j + λ Σ_j Σ_{a=X,Y,Z} σ^R_a(j) σ^L_a(j+1), cell gap Δ = 1. It is translation invariant, nearest neighbour in the big sites, Hermitian, and has an exact on-site Z₂ acting on the code as logical fermion parity. First order vanishes, and the second-order bond is exactly two complex fermion flavours with distinct hopping phases and fillings (the `../ir_control/` target), c = 2.

1. Such a cell exists at D = 64. C9's recipe has to be changed to respect the grading [checked to 10⁻¹⁵].
2. For λ/Δ ≤ 0.2 the logical Hamiltonian, computed from two- and three-cell diagonalization of the full cells (exact through O(λ⁴)), is the target plus parity-even corrections of relative size ≲ 0.5λ, including pairing, with Fermi points moved by < 0.06 and no resonance (margin ≥ 0.54) [numerical].
3. On rings of three full cells, this logical Hamiltonian reproduces the spacings of all 64 low levels to 10⁻⁴ J (λ = 0.05) and 10⁻³ J (λ = 0.1), against 0.4 J and 0.85 J for the second-order chain [numerical].
4. On rings up to L = 9 the logical chain at λ = 0.1 and 0.2 tracks the finite-size spectrum of the free target (same quantum numbers, gaps within about 30%), with L·gap of order one. The symmetric target instead acquires L-independent splittings (0.03 J and 0.13 to 0.15 J at λ = 0.2), masses from resonant pairing [numerical].
5. A calibrated large-scale c for the incommensurate chain was not obtained (below).

Verdict: route A's incommensurate target appears to survive finite coupling at n = 2. The cells generate only small, non-resonant, parity-even corrections, so the infrared statement reduces to `../ir_control/` with measured couplings. The flow to c = 2 itself is supported only up to L = 9.

## Status

Numerical. A calibrated c from open chains at L ≥ 64 (array 17514542, pending at writing) or from U(1)-symmetric MPS would upgrade item 4.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27/28. No Codex/GPT.
- New code, reusing `../d2_followup_220_1282/{vumps.py, fdmrg.py}` unchanged.
- Checks: the second-order bond from the spectrum of H0 alone equals the target to 1.5 × 10⁻¹⁵, and so does the three-cell ring operator (4 × 10⁻¹³ J). Estimators were run on XX, two decoupled XX chains and the free targets. Not rerun by anyone else.

## The cell (`cell.py`)

Qubit ports L, R, each in a Bell pair with a partner, and a logical ℂ⁴ (two modes, intra-cell Jordan–Wigner), so D = 64 and G|i⟩ = |Φ⟩_{L,pL}|Φ⟩_{R,pR}|i⟩. Every port Pauli compresses to zero. U = Z_L Z_R Z_pL Z_pR P_log commutes with H0 and the bond and equals P on the code, with X, Y odd and Z even.

The grading forces a change. With U exact, M_e[a,b] = G†σ_a P_e σ_b G is parity even when a and b have equal parity, so the hopping (a†P) ⊗ a can only come from the mixed blocks [X,Z] and [Y,Z]. C9's block-diagonal recipe puts it in one syndrome block and breaks the grading. Here M_e = p_e 1 + ρ z_e 𝒜, 𝒜 = |X⟩⟨Z| ⊗ e^{−iθ₁}a₁†P + |Y⟩⟨Z| ⊗ e^{−iθ₂}a₂†P + h.c. (right port, with t_α a_α on the left), with excited energies 1 and 5, z ⊥ (1,1) normalized in C_ef = 1/(e_e + e_f), and chemical potentials in the Z–Z block (they enter linearly in ρ). Distinct odd channels per flavour forbid cross-flavour hopping. A graded Naimark dilation with random partner vectors realizes the weights, so higher orders are generic apart from the Z₂. H0 has levels 0 (×4), 1 (×24), 2.5 (×12), 5 (×24), and

h₂ = λ²[−0.7 + J₀ Σ_α (−e^{−iθ_α} c†_{j,α}c_{j+1,α} + h.c. − μ_α n_{j,α})],  J₀ = 0.0222,

with (θ, k_F) = (−2.617, 0.962), (0.298, 0.742), Fermi points 1.040, 2.704, 4.628, 5.839 and δ₂ = 0.596. Both Fermi-sea centres are displaced. The logical scale is 3% of the leakage scale, the price of positive weights.

## Effective bond at finite λ (`ed_check.py`, `analyze_heff.py`)

Two cells (dimension 4096) and three cells (262144, block Davidson) give the des Cloizeaux Hamiltonian of the low band, h₂(λ) and the connected h₃(λ). The traceless residual ‖h₂(λ) − λ²h₂‖ is 2.34λ J (power 3.00, λ = 0.0025 to 0.3), and ‖h₃‖ ≈ 1.4λ² to 1.8λ² J (0.014 J at λ = 0.1, 0.072 J at 0.2). Four-cell terms start at λ⁵. In global fermions, in units of J = J₀λ²:

| λ | t₁, t₂ | pairing D₁₁, D₂₂ | on-site D₁₂ | t₁₂ | degree ≥ 4 | Fermi points | δ₂ |
|---|---|---|---|---|---|---|---|
| 0.02 | 1.018, 1.019 | 0.006, 0.007 | 0.005 | 0.0006 | 0.006 | 1.044 2.700 4.628 5.838 | 0.598 |
| 0.05 | 1.046, 1.047 | 0.015, 0.019 | 0.014 | 0.0015 | 0.016 | 1.051 2.693 4.629 5.836 | 0.603 |
| 0.1 | 1.094, 1.096 | 0.034, 0.042 | 0.032 | 0.003 | 0.034 | 1.063 2.680 4.631 5.831 | 0.589 |
| 0.2 | 1.194, 1.200 | 0.084, 0.106 | 0.079 | 0.007 | 0.084 | 1.099 2.647 4.640 5.817 | 0.544 |
| 0.3 | 1.30, 1.32 | 0.158, 0.199 | 0.150 | 0.015 | 0.157 | 1.149 2.600 4.660 5.796 | 0.475 |

Phases move by < 0.04 rad. Leakage adds pairing and interactions of the same order as the hopping renormalization, none at a resonant momentum. The same-flavour pairing (0.34λ, 0.42λ) is the term that gapped the undisplaced target of `../routeA_numerics/`. Here it pairs momenta summing to 2θ_α, far from zero.

## Full chain on rings (`bigchain.py ring`, `analyze_ring.py`)

Rings of three cells, the whole 64-level logical band in both parity sectors (residual < 10⁻¹⁰). Compared with h₂(λ) + h₃(λ) on the same ring, level spacings agree to 1.0 × 10⁻⁴ J at λ = 0.05 and 9.3 × 10⁻⁴ J at 0.1, with a uniform shift of 0.03 J and 0.14 J (loops around a three-cell ring, absent from the infinite chain). The second-order chain misses the spacings by 0.40 J and 0.85 J. So leakage only renormalizes the logical couplings, and the cluster chain captures this. The λ = 0.2 ring and an L = 4 ring (job 17514543) were still running.

## Logical-chain numerics (`chains.py`, `analyze.py`, `freeobc.py`)

Calibration. VUMPS finite-entanglement scaling gives c = 0.986 for XX and 2.007 for the free symmetric target (χ ≤ 192), but fails at the incommensurate fillings. On the free target the iMPS is non-injective for χ ≥ 8 (ξ diverges), and two decoupled XX chains jump between filling sectors (c = 1.99 or 2.78 depending on the χ window). Open-chain DMRG without U(1) at χ = 32 misses the exact free entropy by 0.5 with fractional particle numbers, while the exact free profile fits c = 1.98 at L = 64. Neither estimator is calibrated here yet.

Rings (Lanczos, L = 4 to 9). Scaled lowest gaps L(E₁ − E₀), in units of the bare hopping:

| L | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|
| free target | 0.79 | 0.36 | 1.70 | 0.97 | 0.74 | 1.52 |
| logical chain, λ = 0.1 | 0.88 | 0.39 | 1.85 | 0.89 | 0.71 | 1.81 |
| logical chain, λ = 0.2 | 0.91 | 0.60 | 2.13 | 0.71 | 0.60 | 2.16 |

The momenta and parity sectors of the ground state and first excitation agree with the free target at every L, and the next three gaps track it to within about 30%. The irregular L dependence comes from the incommensurate filling, which rules out Casimir fits at these sizes. Nothing grows with L, consistent with z = 1 and no gap.

## Symmetric contrast

The SO(4) Majorana point (θ = π/2, half filling, both flavours) from the same construction has the same pairings (D₁₁ ≈ 0.39λ, D₂₂ ≈ 0.44λ), inter-flavour hopping 0.19λ and on-site pairing D₁₂ ≈ 0.27λ, but now all Fermi points sit at 0 and π and the on-site pairing is resonant. On rings at λ = 0.2 the free chain's fourfold ground state at odd L splits as L·ΔE = (0.17, 0.65, 0.74), (0.23, 0.89, 1.04), (0.28, 1.14, 1.34) at L = 5, 7, 9, i.e. by amounts constant in L to within 10%, 0.03 J and 0.13 to 0.15 J (the latter ≈ 2D₁₂). Read as masses, they give correlation lengths of roughly 80 and 17 sites, beyond VUMPS up to χ = 48, where the symmetric chain still tracks the free one (ξ = 10.9 against 11.1). So what goes wrong is a gap at long distances, invisible to entanglement scaling at these χ, as `../ir_control/` anticipated for symmetric points.

## Why VUMPS on the full chain failed

`bigchain.py fes` (sparse MPO, Davidson preconditioned by the cell and environment energies, start dressed by the exact first-order Schrieffer–Wolff MPO of `dressing.py`) does not converge. Near the optimum the lowest eigenvector of H_AC has overlap below 1% with the current A_C (its lowest levels are split by about J/2), so VUMPS jumps between logical states, and under-relaxed and preconditioned-gradient variants stall (`dbg_hac.py`). Leakage Schmidt directions also carry far more energy than the logical tail, so a faithful state needs several times the logical χ (an estimate).

## Three flavours

`cell.py --n3` builds n = 3 at D = 152, exact to 9 × 10⁻¹⁵, but positivity limits J₀ to 3.7 × 10⁻⁴, so finite-coupling runs were not attempted. Its trial momenta (δ₂ = 0.12) should give way to a Sidon-type choice (`../routeA_assembly/`).

## Reproducing

numpy and scipy. `cell.py`, `ed_check.py --target INC|SYM` (writes `data/heff_T.*`, log in `logs/`), `analyze_heff.py`, `chains.py MODEL fes|obc|ed` (MODEL: xx, xx2, free:T, eff:T:λ, eff3:T:λ), `freeobc.py`, `bigchain.py ring|fes`, `analyze_ring.py`, `analyze.py`. Typhon jobs 17514062 (calibrations), 17514542 (logical chains), 17514543 (L = 4 ring), plus login-node runs, all writing to `/scratch/lustre/jkudlerflam/routeA_incommensurate/`.

## Relation to the rest of the repository

This addresses STATUS help item 2 and gap G4 of `../routeA_assembly/`. A route-A cell can carry the exact parity grading at D = 64 for n = 2, leakage only renormalizes the couplings, and each flavour's Fermi-sea centre must be displaced, since the cells generate pairing of order 0.4λ J. The open numerical gap is a calibrated large-scale c at incommensurate filling.
