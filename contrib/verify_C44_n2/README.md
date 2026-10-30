# The n = 2 member of C44: projected chain and c = 3

## Claim

`contrib/verify_C44/` measured c = 2 for the second-order projected n = 1 chain of the fixed-d₀ construction (C39 core, C44 mobile wrapper). Here the same is done for n = 2, where C44 predicts c = n + 1 = 3. At a scheduled point with a weak length interaction (V = −0.16 j_l) and comparable velocities, the projected n = 2 chain gives c = 3.00 ± 0.02 from open-chain interval entanglement and c = 3.0 ± 0.1 from finite-entanglement scaling, each calibrated on the exactly solvable reference of three free modes with the same velocities. Freezing the length mode gives c = 2.00 from both estimators, and multiplying the length-payload coupling by 1000 does not move either result.

## Status

Numerical (open chains of L = 32 and 48 cells at bond dimensions 512 to 1536, infinite chains at bond dimensions 120 to 288), with a heuristic caveat on C39's flavour margin (last subsection of Method). As at n = 1, this is the central charge of the projected model. It tests the Gaussian stability of the retained theory (part (b) of the referee's input S), not the elimination of the auxiliary sectors.

Both estimators are strongly truncation-limited before calibration. Raw values lie between 1.7 and 2.7 for a true c = 3, and the calibrated numbers rest on the exact reference absorbing that bias. The evidence that it does is that the calibrated profile values stay within 0.04 of 3 while the raw values move by 0.9 as the bond dimension triples, and that the reference itself approaches its exact profile (deviation 0.21, 0.073, 0.032 at χ = 512, 1024, 1536). The errors quoted in the Claim are the spread over fit windows, calibration forms, bond dimensions and sizes (profile, χ ≥ 1024), and over subsets of bond dimensions (FES). They are not statistical errors.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-30.
- AI assistance: all code and runs by Claude, partly on Typhon. No Codex, since every step needs the repository's rule text and the verify_C44 code.
- Checked by: internal cross-checks only. The bond is evaluated at 34 and 50 digits and against double precision, the block-sparse DMRG against exact diagonalization, and every truncated profile of the references against its exact free-fermion profile. Nothing from C39, C44 or C45 was rerun, and the verify_C44 rebuild is used as is.

## Method and results

### The projected n = 2 chain

Each cell is one site of dimension 8, the length bit l (lengths n + 3 and n + 4) times the two payload modes z = z₀ + 2z₁ (flavour f is data bit f, in the Jordan-Wigner order the rebuilt bond carries). The bond is

    h = Σ_{l,l'} Π_l Π_l' ⊗ κ_(l+1)(l'+1) + V n_l n_l' − j_l (|10⟩⟨01| + h.c.)_length,

with κ_pq the traceless part of the rebuilt s/δ → 0 bond K_pq (`../verify_C44/precise.py:kbar_bonds`), V = c₂₂ − c₂₁ − c₁₂ + c₁₁ the length density interaction, and the linear length terms absorbed into the length chemical potential. This is the n = 1 model of `verify_C44/effmodel.py` with a 4-dimensional payload (`model2.py`). Dropped, as there, are the O(s/δ) remainder, the λ⁴, λ²u_mov and u_mov² terms, readout leakage, the ring gluing and the physical coordinate.

Precision. At n = 2 the payload hopping sits 10 to 13 orders below the scalar part of K_pq at C44's parameters, so `mpkbar.py` re-evaluates `kbar_bonds` in mpmath from the exact binary values of the rebuilt core matrices and subtracts the scalar parts before rounding. At the scheduled point the 34- and 50-digit bonds agree to 5e-37 (absolute, against a scalar part of 0.05), and the double-precision bond agrees to 8e-7 of the payload velocity. Payload number violation is below 2e-26 and the one-body onsite part below 4e-18 of the payload velocity, as the exact symmetries require.

Structure of the payload. A Hilbert-Schmidt (Majorana) decomposition (`ferm.py`) splits κ₁₁ into a nearest-neighbour hopping matrix T (2 × 2 in flavour space), no onsite one-body term, and an interaction. T is close to rank one, with singular values 14.5 and 0.35 in units of v_p, the mean payload Fermi velocity. The two bands therefore have widths of about 29 v_p but cross zero only where the small singular value matters. The four Fermi points sit in two pockets of width 1.536, one per band (centred at 0 and π after the gauge below), with Fermi velocities 0.857 and 1.144 v_p (the two chiralities of each pocket differ, since T is complex). The fast scale produces no gapless mode, but it produces a large entanglement offset (the exact half-chain entropy of the free payload at L = 32 is 2.09, against 1.0 from the logarithm alone) and makes the problem stiff.

Schedule. At the default core of `core46.py` (w = 0.8, a = 4, η = 0.6, v = 0.65, D = 0.2, g = 0.4) the same structure has a fast-to-slow ratio of about 1700 and an interaction 36 times the payload velocity. A scan (`scan.py`) of core parameters inside the frustration-free, gapped regime (the core gap stays 0.509 and the code dimension 4, `checks.json`) gives a much milder hierarchy at w = 1.8, a = 3.5 (C39's diagnostic value), η = 2, v = 0.3 and D = 1, with q = 0.8, K = 10 and U = 1 unchanged. The conversion g = 0.22 sets the largest payload interaction element to 0.175 v_p (it scales as g²), and the vacancy ratio ρ = 0.008 sets V = −0.085 v_p (it scales as ρ⁴, checked over a factor 10 in ρ). With length density ν = 3/8 and j_l chosen so that the free length velocity 2 j_l sin πν equals v_p (j_l = 0.541 v_p), V/j_l = −0.157, well inside |V| < 2 j_l. The traceless length-payload terms are 2.6e-4 v_p. An earlier choice ν = 0.3 was dropped because 2πν came within 0.02 of a payload wavevector.

The uniform payload hopping phase is removed by the exact gauge c_j → e^{ijφ} c_j (φ = −0.649), which centres the pockets at 0 and π.

| Quantity (units v_p = 3.19e-11 λ²/s) | Value |
| :-- | :-- |
| Payload Fermi velocities | 0.857, 1.144 (ratio 1.34) |
| Length velocity | 1 |
| Fast hopping, singular values of T | 14.5, 0.35 |
| Largest payload interaction element | 0.175 |
| V, and V/j_l | −0.085, −0.157 |
| Traceless length-payload terms, M1 and M2 | 2.6e-4, 0.26 |

### Models and estimators

M0 is the exactly solvable reference, the quadratic payload part (with its two unequal velocities) plus a free length XX chain with the same j_l, so three decoupled free modes with c = 3. M1 is the scheduled point. M2 is M1 with every traceless length-payload term multiplied by 1000 (to 0.26 v_p). C0 is M0 with the length frozen (every cell short, c = 2 exactly) and C1 is M1 with the length frozen, the control expected to give c = 2.

Profile. The open chain needs bond dimensions near 10³ at c = 3 (the exact M0 half-chain entropy at L = 32 is 2.92), so `bsdmrg.py` is a two-site DMRG with both U(1) charges (payload number, length number) in block-sparse form, at N_p = L and N_l = 3L/8. At L = 4 it reproduces exact diagonalization of M0, M1, M2 and C1 in three charge sectors to 5e-11 (`edcheck.json`). The fit is

    S(ℓ) = (c/6) ln x + s₀ + Σ_k [A_k cos kℓ + B_k sin kℓ]/x,  x = (2L/π) sin(πℓ/L),

with Friedel wavevectors k from all differences of the four payload Fermi points (1.536, 1.606, π) and 2πN_l/L, over ℓ₀ ≤ ℓ ≤ L − ℓ₀ with ℓ₀ = 3, 4, 6. On the exact M0 and C0 profiles it gives 2.83 and 1.90 at L = 32 (ℓ₀ = 4) and 2.92 and 1.95 at L = 128, a slow finite-size bias from the shallow pockets. Two calibrations are reported, at identical L and χ. The ratio form is c = 3 f(M)/f(M0). The difference form, c = 3[1 + (f(M) − f(M0))/f_exact(M0)], removes the truncation deficit of the reference additively and the fit bias multiplicatively. For C1 the reference is C0 and the target 2.

FES. One-site VUMPS (`../verify_C44/dmrg.py`) and S = (c/6) ln ξ + b, with ξ from the leading transfer-matrix gap. From random states VUMPS lands in a metastable state with payload density 0.755 (both pockets emptied), which is variationally favoured at bond dimension 32 or less. `run_fes.py` therefore starts from the product of a half-filled payload state (C0, accepted only with density 1 to 1e-2) and a free length state, χ = 3χ_p, and the full model is then optimized from there. Only points with gauge error below 1e-6 enter the fits. This excludes M0 and M2 at χ = 96 and M1 at χ = 144 (errors 2e-3 to 1e-2, also after a restart with a different split of χ), and M1 at χ = 288 did not finish.

### Results

Profile estimator. Entries give ratio / difference calibration at ℓ₀ = 4, with the range over ℓ₀ = 3, 4, 6 and both forms in brackets.

| Model | L = 32, χ = 512 | L = 32, χ = 1024 | L = 32, χ = 1536 | L = 48, χ = 1024 |
| :-- | :-- | :-- | :-- | :-- |
| M0 raw (exact 2.83; 2.86 at L = 48) | 1.77, max deviation 0.21 | 2.43, deviation 0.073 | 2.65, deviation 0.032 | 2.02, deviation 0.19 |
| M1 | 2.98 / 2.99 [2.98, 3.01] | 3.00 / 3.00 [3.00, 3.01] | 3.00 / 3.00 [3.00, 3.01] | 2.99 / 2.99 [2.99, 3.00] |
| M2 | 2.96 / 2.98 [2.96, 3.00] | 2.99 / 3.00 [2.99, 3.01] | 3.00 / 3.00 [3.00, 3.01] | 2.99 / 2.99 [2.99, 3.00] |

| Control | L = 32, χ = 512 | L = 48, χ = 512 | L = 64, χ = 512 |
| :-- | :-- | :-- | :-- |
| C0 raw (exact 1.90, 1.92, 1.93) | 1.88, deviation 0.004 | 1.83, deviation 0.019 | 1.75, deviation 0.045 |
| C1 | 2.003 [2.003, 2.004] | 2.002 [2.002, 2.003] | 2.002 [2.002, 2.003] |

Payload densities are 1.000 in the bulk of every run and length densities 0.378 to 0.380 (N_l/L = 0.375), and the truncated weights at the largest χ are 3.4e-6 (L = 32) and 1.2e-5 (L = 48).

FES estimator. Calibrated values use the bond dimensions common to the model and M0. The range is over all subsets of at least three of them.

| Model | χ | raw c | calibrated c |
| :-- | :-- | :-- | :-- |
| M0 (c = 3) | 120, 144, 192, 240, 288 | 2.44 | |
| M1 | 120, 192, 240 | 2.46 (M0 on the same set 2.46) | 3.00 (two-point values 2.86 to 3.07) |
| M2 | 120, 144, 192, 240, 288 | 2.43 | 2.99 [2.90, 3.09] |
| C0 (c = 2) | 32, 48, 64, 96, 128, 160 | 2.04 | |
| C1 | 32, 48, 64, 96, 128, 160 | 2.04 | 2.001 [1.999, 2.003] |

The correlation lengths only reach ξ = 6 to 11.5 sites, set by the payload (its four-fold degenerate transfer-matrix eigenvalue). The raw M0 slope is still rising across this window (its local value between χ = 240 and 288 is 3.2), so the raw numbers are pre-asymptotic and the FES calibration is the weaker of the two. It is also less independent than it looks. M1 and M0 differ by weak terms, and a calibration that divides M1 by M0 tests that these terms do not change the scaling, not the scaling itself.

What the results test. The calibrated profile values for M1 and M2 lie between 2.99 and 3.01 at both sizes and every χ ≥ 1024, so neither the payload interaction (0.175 v_p), nor the attractive length interaction, nor a length-payload coupling 1000 times the scheduled one (0.26 v_p) changes the count of gapless modes at these sizes. The frozen-length control returns 2.00 from both estimators, so the third unit of c is carried by the length fluid. The payload interaction barely moves the entanglement (the C1 and C0 slopes agree to 2e-3, relative), which reflects that at g = 0.22 the interacting payload sits close to its free reference.

### Flavour margin (heuristic)

C39 removes the residual half-filling coincidences of the payload pair sums with a weak density readout, μ q^{2j+2}/S n_j. It is not included here. At this point the two pockets have equal widths exactly at μ = 0 (h(k + π) = −h(k)), so the two flavours share 2k_F = 1.536 and a pair-transfer process is resonant. Adding the readout with μ up to 2 v_p changes the width difference only to below 4e-5 (`checks.json`), since C39's slopes s_j are tiny at these parameters. A marginal process at coupling 0.175 v_p would, heuristically, open a gap of order exp(−πv_p/g_eff), far below anything a chain of 48 cells resolves. These sizes therefore cannot tell the scheduled point from one with C39's margin, and the c = 3 found here holds only on scales below that crossover.

## Reproducing

numpy, scipy and mpmath, `python3 -B` from this folder. The scripts add `../verify_C44` to the path.

1. `mpkbar.py 0.008 0.22 1.8 3.5 2.0 0.3 1.0` writes `kbar_n2_rho0.008_g0.22_w1.8_a3.5_eta2.0_v0.3_D1.0.npz` (about a minute, output in `log_mpkbar.txt`). `checks.py` writes `checks.json`, `edcheck.py` writes `edcheck.json` (seconds), and `scan.py` reproduces the parameter scan.
2. `run_profile.py MODEL L CHI` writes `results*/MODEL_L{L}_chi{CHI}.json`. At L = 32 and χ = 1024 a run takes 10 to 30 hours on 8 cores, and `jobs_prof.sh A|B|C` packs the runs into Typhon allocations.
3. `run_fes.py MODEL CHI [CHI_L]` writes `results*/MODEL_vumps_chi{CHI}*.json` (one to twenty hours each), with `jobs_fes.sh` and `jobs_fes_redo.sh` the Typhon packing.
4. `fit2.py` reads `results/` and `results_typhon/` and writes `c_results.json`.

## Relation to the rest of the repository

This extends `verify_C44` item (i) from n = 1 to n = 2 at the level stated. It does not test S(a), the microscopic chain, the readout margin of C39 or the lacunary schedule at large n.
