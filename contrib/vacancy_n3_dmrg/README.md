# Logical vacancy model at n = 3: open-chain DMRG

Follow-up to `../vacancy_n2_dmrg/`, which measured c = 2.9 ± 0.1 for the n = 2 member of the logical vacancy model (target 3). This folder takes the same code and estimator to n = 3, where a site is a vacancy or a cell with three Jordan–Wigner flavours (local dimension 9). Target c = 4 (charge plus three flavours), control c = 3 (one flavour gapped).

## Claim

[numerical, open-chain U(1)⁴ DMRG, χ ≤ 1024, interval-entropy profile calibrated on the exact four-chain free reference] At L = 30, the only size converged in χ within this run, LV_3 at generic fillings with three distinct flavour hoppings has c = 3.98 ± 0.05. The value has stopped moving with χ (3.973 at χ = 768, 3.979 at χ = 1024, discarded weight 5.5 × 10⁻⁶) and the extrapolation in discarded weight gives 3.99. At L = 40 the run reached χ = 768 and gives 3.82, extrapolated 3.88, still rising with χ. Calibrating instead against the matched-velocity four-chain reference run through the same DMRG at the same χ gives 4.00 at L = 30 (χ = 256 to 512) and 3.95 to 3.99 at L = 40 (χ = 256 to 384).

[numerical] The control, with the third flavour Néel-gapped, has c = 2.89 at L = 35 and χ = 768 (variants 2.81 to 3.06, extrapolated 2.93), with its three-chain reference converged to 3.00 at the same size. So the model sits about one unit above the control, and the fourth mode is resolved at L = 30 to 35.

[heuristic] Larger L (45 to 60) was not converged in χ in the time available, and the four-mode state needs much more χ per unit length than the three-mode one, so nothing here tests whether the L = 30 value survives at larger sizes. The control's deficit of about 0.1 at L = 35 has the same sign and size as the n = 2 control's at comparable L, which I read as a finite-size effect of the gapped flavour's boundary region.

## Status

n = 3 is settled at L = 30 (model) and L = 35 (control) at the physics standard, and the trend at L = 40 is consistent with 4. Not done: convergence at L ≥ 40 for the model and L ≥ 45 for the control, and the four-chain DMRG references beyond χ = 512. Runs were still going when this was written (see Reproducing), so the tables can be extended by rerunning `analyze.py` on their output. Nobody else has rerun anything.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-29, in about two hours. No Codex call (repository-bound numerics with code reuse).
- `qdmrg.py`, `models.py`, `fitlib.py` are unchanged copies of the n = 2 files (the DMRG takes any number of abelian charges, here M, N₀, N₁, N₂). `run_obc.py` differs only in its defaults and docstring, and `analyze.py` in its model and reference names, the dropped FES part, and one extra estimator (below).
- `check_ed.py` compares the DMRG with exact diagonalization of open LV_3 chains (L = 7, and L = 8 including the control with the Ising term on). Energies agree to 2 × 10⁻¹¹ and every S(ℓ) to 1.2 × 10⁻⁹ (`data/check_ed.json`).
- The matched-velocity references (four decoupled XX chains for the model, d = 16, and three for the control, d = 8) are run through the same code, and their convergence to the exact profile (table 3) calibrates the truncation bias of the fit at each L and χ.
- The direct fit of the exact four-chain profile, with no subtraction, ranges from 0.5 to 3.0 at L = 30 and from 0.05 to 5.4 at L = 40 across the variants, which is why the subtraction estimator below is the only usable one at these sizes.

## Method

**Model.** LV_3 of `../vacancy_mobile_numerics/logical_vac.py` (bond `model_bond(3, t, J, 0, τ, 0, vz)`), t = 1, J = 0.5, τ = (1, 0.5, 0.75), on open chains, as in the n = 2 README with one more flavour. Cells hop through vacancies with amplitude t, flavours hop with Jτ_f between adjacent cells, and the bond carries the intra-cell strings. Fields are zero since the sector fixes the fillings.

**Fillings.** ρ = M/L = 3/5 and ν_f = N_f/M close to (1/3, 1/4, 1/6), at L = 30, 40, 50, 60 with (M, N₀, N₁, N₂) = (18, 6, 4, 3), (24, 8, 6, 4), (30, 10, 8, 5), (36, 12, 9, 6). No flavour is half filled and no pair sums to 1. The control adds an Ising term vz = 4.5 on the third flavour (XXZ anisotropy Δ = vz/(2τ₂) = 3, as for the n = 2 control), at half filling with odd M and N₂ = (M + 1)/2, at L = 35, 45, 55 with (M, N₀, N₁, N₂) = (21, 7, 5, 11), (27, 9, 7, 14), (33, 11, 8, 17).

**DMRG.** As for n = 2. χ grows through 64, 128, 256, 384, 512, 768, 1024 (1536 queued), with three sweeps per χ, the first carrying a random admixture of relative norm 10⁻³, and a checkpoint after every χ. Three sweeps per χ instead of four to six was the price of the time limit; at L = 30 the energy changes by less than 10⁻⁶ between the second and third sweep at every χ ≥ 384.

**Estimator.** The difference estimator of the n = 2 README: c = m + c_Δ, with c_Δ the slope of ΔS(ℓ) = S_model(ℓ) − S_ref(ℓ), S_ref the exact profile of m open free-fermion chains with particle numbers (M, N₀, N₁, N₂) (m = 4) or (M, N₀, N₁) for the control (m = 3), fitted over lmin ≤ ℓ ≤ L − lmin to (c_Δ/6) X + s₀ + Friedel terms with X the chord log. The fit wavevectors are 2πM/(L + 1) and 2πN_f/(L + 1) for all three flavours (the gapped one included), and in half the variants their pairwise sums and differences. The variants are lmin = 6, 8, 12 with and without the combination wavevectors.

New here is that I quote as the main value the median of the three variants that include the combination wavevectors. The reason is the control. Its half-filled gapped flavour produces strong oscillations at the sums and differences of its wavevector with the others, and without them the fit fails (lmin = 12 gives 0.83 at L = 35 and χ = 768 while the combination variants give 2.81 to 3.06). For the model at L = 30 the two choices agree (3.98 against the six-variant median 3.975), and for the converged references the combination variants are the closer ones to the exact answer (table 3). The six-variant median and range are listed as well. The extrapolation is linear in the discarded weight from the last two χ, as for n = 2.

**Matched velocities.** The reference chains have hoppings (1, 0.549, 0.290, 0.452) for the model and (1, 0.549, 0.290) for the control, which give the first-order velocities v_c = 2 sin πρ and v_f = 2Jατ_f sin(πν_f)/ρ with α = 0.447 [estimate, first order in J/t, as in the n = 2 README]. The exact reference does not depend on them, but the finite-χ DMRG reference does.

## Results

**Table 1. LV_3 (target 4), difference estimator against the exact four-chain reference.**

| L | χ | discarded weight | c (combination variants) | six-variant median and range | c at previous χ | extrapolated | vs DMRG reference at the same χ |
|---|---|---|---|---|---|---|---|
| 30 | 1024 | 5.5e-6 | 3.979 | 3.975 (3.68 to 4.06) | 3.973 (768) | 3.986 | 4.00 (χ = 256, 384, 512) |
| 40 | 768 | 1.9e-5 | 3.821 | 3.41 (1.80 to 4.09) | 3.754 (512) | 3.880 | 4.05, 3.99 (χ = 256, 384) |
| 50 | 384 | 7.1e-5 | 2.757 | 2.06 (1.07 to 3.07) | 2.558 (256) | 3.02 | 4.53 (χ = 256) |
| 60 | 768 | 2.6e-5 | 2.096 | 2.09 (1.92 to 2.39) | 1.505 (512) | 2.68 | not run at matching χ |

At L = 30 the three combination variants at χ = 1024 are 4.006, 3.965, 3.979, and the χ change from 768 is +0.006, so I quote 3.98 ± 0.05 [numerical], the error covering the variant spread, the χ change and the extrapolation. Only L = 30 is converged. L = 50 and 60 are far from converged and are listed only to show where the runs stood.

**Table 2. Control, Néel third flavour (target 3), against the exact three-chain reference.**

| L | χ | discarded weight | c (combination variants) | six-variant median and range | c at previous χ | extrapolated | vs DMRG reference at the same χ |
|---|---|---|---|---|---|---|---|
| 35 | 768 | 8.3e-6 | 2.887 | 2.45 (0.83 to 3.06) | 2.835 (512) | 2.927 | 2.89 (χ = 768) |
| 45 | 768 | 1.2e-5 | 2.526 | 2.11 (1.59 to 2.83) | 2.262 (512) | 2.760 | 2.58 (χ = 768) |
| 55 | 512 | 2.7e-5 | 1.611 | 1.61 (1.51 to 1.65) | 1.313 (384) | 2.04 | 1.97 (χ = 384) |

The control converges more slowly in χ than its own three-chain reference, presumably because the gapped flavour multiplies every Schmidt value of the gapless part by its own short-range spectrum, so fewer of the χ states go to the gapless modes. The same-χ calibration therefore undercorrects the control at L ≥ 45, and only L = 35 is usable. There the combination variants give 2.81, 2.89, 3.06, the change from χ = 512 is +0.05, and I quote 2.9 ± 0.1 [numerical].

**Table 3. DMRG on the matched-velocity references: c_Δ against the exact profile (six-variant median) and the half-chain entropy error.**

| reference | χ = 128 | 256 | 384 | 512 | 768 |
|---|---|---|---|---|---|
| four chains, L = 30 | −0.17 (−0.52) | −0.18 (−0.27) | −0.10 (−0.17) | −0.10 (−0.11) | |
| four chains, L = 40 | −1.19 (−0.84) | −1.61 (−0.57) | −1.18 (−0.38) | | |
| four chains, L = 50 | −3.29 (−1.00) | −2.50 (−0.57) | | | |
| three chains, L = 35 | −0.89 (−0.14) | −0.32 (−0.05) | −0.12 (−0.02) | −0.05 (−0.008) | −0.014 (−0.002) |
| three chains, L = 45 | −1.57 (−0.23) | −0.64 (−0.10) | −0.29 (−0.05) | −0.15 (−0.02) | −0.05 (−0.009) |
| three chains, L = 55 | −2.02 (−0.37) | −1.07 (−0.18) | −0.62 (−0.10) | | |

With the combination variants the four-chain reference at L = 30 and χ = 512 reads 3.931 (extrapolated 3.996), and the three-chain reference at L = 35 and χ = 768 reads 2.998. At fixed χ and L the model at L = 30 is closer to its exact answer than the four-chain reference is (3.93 against 3.90 at χ = 512), which is why the same-χ calibration lands at 4.00 already from χ = 256.

## Interpretation

[numerical] Three observations point to c = 4 at n = 3. The model at L = 30 is converged at 3.98 with χ = 1024. The same-χ calibration against the matched four-chain reference gives 3.95 to 4.05 at L = 30 and 40 from χ = 256 on, before either is converged. And the control at L = 35 reads 2.89, about one unit lower, with its reference converged to 3.00 at the same size. Together with n = 1 (c = 2, `../vacancy_mobile_numerics/`) and n = 2 (c = 2.9 ± 0.1, `../vacancy_n2_dmrg/`), this fits c = n + 1.

[heuristic] The four-mode state needs much more χ than the three-mode one. At L = 40 and χ = 384 the four-chain reference still has its half-chain entropy 0.38 below exact, against 0.05 for three chains at L = 45, so χ of a few thousand would be needed for L = 60. The L = 30 value carries an unknown finite-size correction, which at n = 2 was about −0.1 at L = 40 to 100. The fact that L = 30 here lands closer to the target than n = 2 did at L = 40 may be partly luck of the fillings at this size, and I would not read anything into it.

## Reproducing

numpy and scipy only, `python3 -B` from this folder (`run_obc.py` imports `logical_vac.py` from `../vacancy_mobile_numerics/`).

- `check_ed.py` (a minute) writes `data/check_ed.json`.
- `tasks_obc.txt` lists all fourteen runs (`run_obc.py NAME --model lv|kxx --L L --M M --N N0 N1 N2 [--vz 0 0 4.5] [--vels ...] --chis ... --state FILE`), each writing `data/obc/NAME.json` after every χ and a checkpoint under `state/` (on scratch, not in the repository). `launch_login.sh TASKFILE T` starts the lines in the background at nice 19 with T threads, and `node.slurm TASKFILE` runs them in parallel on one Typhon node, resuming from the checkpoints.
- How the data in `data/obc/` were produced. `lv3_L30`, `lv3_L40`, `ref4_L30`, `ref4_L40` ran on a 10-core laptop with two threads each (L = 40 at χ = 768 took about 40 minutes). The other ten started on the Typhon login node at nice 19, which was overloaded by other jobs (load above 100 on 48 cores), and were moved after about 50 minutes to a single node of the typhon partition (job 17526254, 96 cores, 3 hour limit), resuming from their checkpoints. When this was written that job was still running toward χ = 1024 and beyond, in the campaign scratch area on Typhon under `vacancy_n3_dmrg/slurm/`, and the four laptop runs had been stopped and queued there from their checkpoints (job 17526438, pending on the per-user job limit; lv3_L40 at χ = 1024 had not finished a full set of sweeps and restarts that χ). Copying the newer JSON files into `data/obc/` and rerunning `analyze.py` extends every table.
- `analyze.py` prints every table above and writes `data/summary.json`. The half-chain pairs it also prints are meaningless here, since the fillings are rounded differently at each L.

No MPS is stored in the repository.

## Relation to the rest of the repository

This extends `../vacancy_n2_dmrg/` from n = 2 to n = 3 with the same model, code and estimator, and adds one more point to the c = n + 1 count for the logical vacancy model. It does not touch the first-order step from LV_n back to the real lift, which stays heuristic, nor any FES at n = 3.
