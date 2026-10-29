# Logical vacancy model at n = 2: open-chain DMRG and large-χ FES

Follow-up to `../vacancy_mobile_numerics/` (section "Late results and calibrated c (2026-10-28, evening)"), which left the n = 2 member of the logical vacancy model LV_n unresolved. Target c = 3 (charge plus two flavours), control c = 2 (one flavour gapped).

## Claim

[numerical, open-chain U(1)³ DMRG, L = 40 to 100, χ ≤ 1024 (1536 at L = 60), interval-entropy profiles calibrated on the exact three-mode free reference] LV_2 at generic fillings has c = 2.9 ± 0.1, and the control with a Néel-gapped second flavour has c = 1.85 ± 0.15 at L = 75 and 105 (1.69 at L = 45). At every size the model sits about one unit above the control, so the third mode is resolved. The model stays about 0.1 below 3 between L = 40 and 100, which is the same size of deficit the control shows against 2 at L ≥ 75. I read both deficits as finite-size corrections to the fit [heuristic]. They are not a DMRG artefact, since the matched-velocity free references run through the same code and fit converge to their exact profiles (to 0.012 in c at χ = 1024 for three chains at L = 40, and to 0.000 for two chains at L = 75).

[numerical, VUMPS FES, χ = 16 to 160 for the model] Finite-entanglement scaling now gives the model 2.8 to 3.2 (calibrated, depending on the χ window and the variant, fit errors 0.05 to 0.2) and the control 1.95 ± 0.14 on its final branch (χ = 96 to 192). This is consistent with the DMRG result, but the large-χ part of the FES (256 to 512) was not reached for the model when this was written, and on its own the FES does not settle n = 2.

## Status

The interval estimator settles n = 2 at the physics standard (c = 3, not 2), with a residual 0.1 deficit at L ≤ 100 that I have not removed. Not done: χ = 256 to 512 FES for the model (the four VUMPS jobs were still running, the model at χ = 192), and open chains beyond L = 100 at converged χ. Nobody else has rerun anything.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28/29. No Codex call (the work is repository-bound numerics and code).
- `qdmrg.py` is a new block-sparse two-site DMRG with three conserved charges (cell number M and the two flavour numbers). `check_ed.py` compares it with exact diagonalization of open LV_2 chains (L = 8 and 9, including a control with the Ising term on). Energies agree to 2 × 10⁻¹¹ and every S(ℓ) to 6 × 10⁻⁹ (`data/check_ed.json`).
- The same code, run on the three decoupled free-fermion chains with matched velocities, converges to the exact profiles as χ grows (table 3). This calibrates the truncation bias, which turns out to be the main systematic.
- A direct fit of the exact three-chain profile (no subtraction) gives 2.98 to 3.01 for L ≥ 80 at lmin = 8, but 2.4 to 3.0 at L = 40 depending on the variant, which is why the subtraction estimator below is the main one.

## Method

**Model.** LV_2 of `../vacancy_mobile_numerics/logical_vac.py` (bond `model_bond(2, t, J, 0, τ, 0, vz)`), t = 1, J = 0.5, τ = (1, 0.5), on open chains. A site is a vacancy or a cell with two Jordan–Wigner flavours (d = 5). Cells hop through vacancies with amplitude t, flavours hop with Jτ_f between adjacent cells, and the bond carries the intra-cell strings, so the open-chain MPS needs no extra signs. In a fixed sector μ_v is a constant. The flavour fields are set to zero, since the sector fixes the fillings (with fields on they act as an adjacency-weighted chemical potential, which is how the FES sets its fillings).

**Fillings.** ρ = M/L = 3/5, ν₀ = N₀/M = 1/3, ν₁ = N₁/M = 1/4, at L = 40, 60, 80, 100 (M = 0.6L, N₀ = 0.2L, N₁ = 0.15L). No flavour is half filled, ν₀ + ν₁ ≠ 1 (no inter-flavour umklapp), and the lowest commensurate harmonics are of high order [heuristic, bosonization]. The control adds an Ising term vz = 3 (XXZ anisotropy Δ = vz/(2τ₁) = 3) on the second flavour, which needs half filling. It uses odd M with N₁ = (M + 1)/2, so the squeezed Néel state is unique, at L = 45, 75, 105 (M = 27, 45, 63, N₀ = M/3).

**DMRG** (`qdmrg.py`). Every bond carries the charge vector accumulated to its left, encoded as one integer. The bond h is split by charge transfer into Σ_k A_k ⊗ B_k (12 terms for the model). The two-site effective Hamiltonian is written as a part that conserves the right bond charge (block Hamiltonian, coupling to site 1, and h on the two sites) plus a part that conserves the left bond charge. Each part is assembled as one dense matrix per bond sector, so a Lanczos step costs two dense products per sector. Runs start from an evenly spread product state and grow χ through 64, 128, 256, 384, 512, 768, 1024, with 4 to 6 sweeps per χ. The second pass adds a random admixture of relative norm 10⁻³ and 10⁻⁴ in the first two sweeps at every χ and checkpoints the MPS to scratch.

**Estimator.** The main one is the difference estimator c = m + c_Δ, with c_Δ the slope of the profile difference ΔS(ℓ) = S_model(ℓ) − S_ref(ℓ). Here S_ref is the exact profile of m independent open free-fermion chains of length L with particle numbers (M, N₀, N₁) (m = 3), or (M, N₀) for the control (m = 2), which reproduces the physical Friedel wavevectors. The slope comes from a linear least-squares fit of ΔS over lmin ≤ ℓ ≤ L − lmin to

(c_Δ/6) X + s₀ + Σ_q [a_q cos(qℓ) + b_q sin(qℓ)] e^{−X}, with X = ln[(2(L + 1)/π) sin(πℓ/(L + 1))],

where q runs over 2πM/(L + 1), 2πN_f/(L + 1) and, in half the variants, their sums and differences. The six variants (lmin = 6, 8, 12, with and without combination wavevectors) give the quoted range, and the central value is their median. Subtracting the reference removes most of the open-chain finite-size structure, which a direct fit of S(ℓ) does not (the exact three-chain reference itself gives 2.4 to 3.0 from a direct fit at L = 40). Also recorded are the ratio calibration m c_model/c_ref of direct fits (it agrees with the difference estimator to 0.01 once χ is converged), the same difference taken against the matched-velocity reference computed by DMRG at the same χ, and half-chain pairs 6[S(L₂/2) − S(L₁/2)]/ln[(L₂ + 1)/(L₁ + 1)] calibrated by the exact reference.

**Truncation.** For these c = 3 states the energy-optimal MPS at fixed χ has much less entanglement than the truncated exact state (at L = 80 and χ = 256 the DMRG reference has S(L/2) 0.35 below exact, while keeping the 256 largest exact Schmidt values costs only 0.01). So χ convergence is checked run by run, and the value at the largest χ is extrapolated linearly in the discarded weight from the last two χ. The same extrapolation on the DMRG references underestimates the exact value by 0.000 (L = 40, χ = 768 to 1024, discarded weight 4 × 10⁻⁷), 0.07 (L = 60, 3 × 10⁻⁶) and 0.15 (L = 80, 5 × 10⁻⁶), so the extrapolated model values below are, if anything, low by a few hundredths.

**Matched velocities.** The reference chains have hoppings (1, 0.549, 0.290), which reproduce the first-order velocities of the three modes, v_c = 2 sin πρ and v_f = 2Jατ_f sin(πν_f)/ρ with α = ρ − sin²(πρ)/(π²ρ) = 0.447 [estimate, first order in J/t, as in the previous README]. The exact reference does not depend on the hoppings, but the finite-χ DMRG reference does.

**FES** (`fes2.py`). One-site VUMPS of `../d2_followup_220_1282/vumps.py` (dense, no symmetry) on LV_2 with flavour fields (0.313, 0.3) and μ_v = 0.3, which gives ρ = 0.631, ν₀ = 0.455, ν₁ = 0.409 (generic, no flavour half filled, ν₀ + ν₁ ≠ 1). The control has fields (0.313, 0) and the same Ising term, and lands at ρ = 0.765, ν₀ = 0.45, ν₁ = 1/2. References are decoupled XX chains with velocities (1, 0.40, 0.19), and (1, 0.40) for the control, at fillings (0.61, 0.275, 0.244), the first-order estimates for the model at ρ ≈ 0.61. c comes from S = (c/6) ln ξ with ξ from the leading transfer-matrix eigenvalue (m = 1), or with ln ξ replaced by the mean of ln ξ_j over the m leading distinct transfer-matrix scales (m = 2, 3, the "one ξ per species" variant). Calibration is m_modes c_model/c_ref over the same χ window.

## Results

**Table 1. LV_2 (target 3), difference estimator against the exact reference.**

| L | χ | discarded weight | c at χ_max | variant range | c at previous χ | extrapolated |
|---|---|---|---|---|---|---|
| 40 | 1024 | 1.6e-7 | 2.859 | 2.45 to 2.94 | 2.849 (768) | 2.864 |
| 60 | 1536 | 1.9e-7 | 2.890 | 2.87 to 2.91 | 2.840 (1024) | 2.909 |
| 80 | 1024 | 1.4e-6 | 2.759 | 2.75 to 2.78 | 2.627 (768) | 2.893 |
| 100 | 1024 | 2.0e-6 | 2.635 | 2.60 to 2.67 | 2.447 (768) | 2.866 |

At L = 40 the variant range is wide because lmin = 12 leaves few points. Against the DMRG reference at the same χ the model gives 2.94 to 3.16 already at χ = 64 to 512 (L = 60: 3.14, 3.14, 3.09, 3.05, 2.99 at χ = 128 to 768), approaching the converged values from above.

**Table 2. Control, Néel second flavour (target 2).**

| L | χ | discarded weight | c at χ_max | variant range | c at previous χ | extrapolated |
|---|---|---|---|---|---|---|
| 45 | 1024 | 2.5e-8 | 1.692 | 1.57 to 1.75 | 1.690 (768) | 1.693 |
| 75 | 1024 | 9.5e-8 | 1.832 | 1.76 to 1.88 | 1.823 (768) | 1.837 |
| 105 | 1024 | 2.2e-7 | 1.871 | 1.82 to 1.91 | 1.846 (768) | 1.888 |

Calibrated half-chain pairs give 2.00 (L = 45 → 75) and 1.97 (75 → 105) for the control, against 2.81 (40 → 60) for the model, the only model pair with both ends converged.

**Table 3. DMRG on the matched-velocity references: c_Δ against the exact profile, and the half-chain entropy error.**

| reference | χ = 128 | 256 | 384 | 512 | 768 | 1024 |
|---|---|---|---|---|---|---|
| three chains, L = 40 | −1.28 (−0.26) | −0.47 (−0.09) | −0.20 (−0.04) | −0.10 (−0.02) | −0.03 (−0.005) | −0.012 (−0.002) |
| three chains, L = 60 | −2.07 (−0.37) | −1.20 (−0.19) | −0.73 (−0.11) | −0.46 (−0.07) | −0.20 (−0.03) | |
| three chains, L = 80 | −2.37 (−0.59) | −1.60 (−0.35) | −1.11 (−0.23) | −0.79 (−0.15) | −0.42 (−0.08) | |
| two chains, L = 75 | −0.22 (−0.04) | −0.03 (−0.005) | −0.007 | −0.002 | 0.000 | |

The three-mode reference needs several times the χ of the two-mode one for the same accuracy, and the model behaves the same way (tables 1 and 2), which is why the earlier χ ≤ 96 work saw about two modes.

**Table 4. FES (VUMPS), calibrated c.** Model to χ = 160, its reference to χ = 128, control to χ = 192, its reference to χ = 160.

| window | estimator | LV_2 (raw, reference, calibrated) | control (raw, reference, calibrated) |
|---|---|---|---|
| χ ≥ 16 | one ξ | 2.80 ± 0.12, 2.94, 2.85 | 1.75, 1.90, 1.83 |
| χ ≥ 32 | one ξ | 3.19 ± 0.15, 3.00, 3.18 | 1.75, 1.98, 1.77 |
| χ ≥ 48 | one ξ | 2.99 ± 0.13, 2.89, 3.11 | 1.73, 2.01, 1.72 |
| χ ≥ 64 | one ξ | 3.05 ± 0.22, 3.24, 2.82 | 1.69, 1.87, 1.81 |
| χ ≥ 32 | two scales | 3.07 ± 0.07, 3.02, 3.05 | 2.30, 1.98, 2.32 |
| χ ≥ 32 | three scales | 2.90 ± 0.06, 3.10, 2.81 | 2.43, 1.96, 2.48 |
| χ = 96 to 192 | one ξ | not reached | 1.96 ± 0.10, 2.01, 1.95 |

LV_2 reached ξ = 20.3 at χ = 160 (S = 2.18). The control changes branch between χ = 64 and 96 (ξ jumps from 12.7 to 51, and ν₀ moves from 1/2 to 0.45), so windows that straddle the jump mix two states and only the χ ≥ 96 row is clean for it. The multi-scale variants are not meaningful for the control, whose subleading transfer-matrix scales include the gapped flavour. [numerical, negative] For the model the FES windows scatter between 2.8 and 3.2, and the χ = 256 to 512 runs had not arrived.

## Interpretation

[numerical] Three things point to c = 3. First, the model is converged or extrapolated to 2.86 to 2.91 at four sizes, while the control is 1.69 to 1.89, so the gap between them is about 1.0 to 1.2 at every comparable L. Second, the matched-velocity calibration at fixed χ sits at 2.8 to 3.16 from χ = 64 on, so the third mode is present even where the uncalibrated profile still looks like c ≈ 1. Third, the control's half-chain pairs return 2.0 once both sizes are converged.

[heuristic] The leftover deficit of about 0.1 seems to be a finite-size effect of the interacting charge mode, whose Friedel terms decay with a Luttinger exponent different from the free reference's, so subtracting the reference leaves a small slope. The control's deficit shrinks from 0.31 to 0.11 between L = 45 and 105, which fits this reading. The model's does not visibly shrink between L = 40 and 100, so this part is not established.

## Reproducing

numpy and scipy only, `python3 -B` from this folder.

- `check_ed.py` (a minute) writes `data/check_ed.json`.
- `run_obc.py NAME --model lv|kxx --L L --M M --N N0 N1 [--vz 0 3] [--vels ...] --chis ... [--state FILE]` writes `data/obc/NAME.json` after every χ. Task lists are `tasks_obc.txt` (first pass: `lv2_L40` and `ctrl_L45`, `ctrl_L75`, `ctrl_L105` come from it, with mixing only at the first χ and 6 sweeps per χ) and `tasks_obc2.txt` (second pass: `lv2_L60`, `L80`, `L100` and the references, with mixing at every χ, 4 sweeps per χ, and checkpoints under `state/` on scratch, not in the repository). The runs took 4 to 6 hours each to χ = 1024 on the Typhon login node at nice 19. The first apollo array (17522348, 17522349, 17522358) was killed by the partition after 26 minutes, which is why most of the work ran on the login node.
- `fes2.py NAME ...` (`tasks_fes.txt`, Typhon array 17522449 on the typhon partition, qos long) writes `data/fes/NAME.json`. These jobs were still running, as were `lv2_L80`, `lv2_L100`, `ref3_L60` and `ref3_L80` (toward χ = 1536) on the login node. Their output lands in the campaign scratch area on Typhon (`vacancy_n2_dmrg/`), and rerunning `analyze.py` on it updates every table.
- `analyze.py` prints every table above and writes `data/summary.json`. `fitlib.py` holds the profile fit, `models.py` the charges, product states, kxx bonds and exact references.

No MPS is stored in the repository.

## Relation to the rest of the repository

This answers item (i) of "What remains" in `../vacancy_mobile_numerics/README.md` for the interval estimator and partly for the FES. It does not touch item (ii) (the bond's own T with on-site μ = iA, U(1) only) or the first-order step from LV_n back to the real lift, which stays heuristic. The truncation result of table 3 also explains the earlier n = 2 FES deficit. For a fixed χ, energy-optimal MPS give up the slowest mode first.
