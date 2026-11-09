This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C176 result: a blockade-exact DMRG and converged open chains of the Fibonacci-cusp model

Derivations and the code design are in PROOF.md, the log in PROGRESS.md. Code: fibdmrg.py (DMRG), validate.py and validate2.py (checks against C173's ED, copied unmodified as eff_ed.py and ring_ed.py), run.py (production driver), analyze.py (estimators).

## Pre-registered criteria (written before any production output was looked at)

Registered 2026-11-09, after validation runs at L ≤ 24 and timing and density scans at L = 48 only.

Definitions.
- S(l) is the von Neumann entropy at cut l (left block = sites 0..l−1) of the final right-to-left sweep.
- Window average S̄_W(L) is the mean of S(l) over the cuts l ∈ [L/2 − W, L/2 + W]. W = 1 (three middle cuts, C173's estimator), W = 6 and W = 12 are all reported.
- Two-size estimate c_two(L₁ → L₂) = 6 [S̄_W(L₂) − S̄_W(L₁)] / ln(L₂/L₁), at fixed density n = N/L.
- Profile estimate c_prof(L): S(l) smoothed with a running mean of 12 cuts, fitted to (c/6) ln[(2L/π) sin(πl/L)] + a over L/6 ≤ l ≤ 5L/6.
- Truncation is the largest discarded weight in the final sweep.
- Valley polarization: P_j = c†_{j+1}c_j + h.c. in the c frame (the vector chirality of the J1–J2 picture, odd under k → π − k). C_P(r) is its connected correlation from the four sites next to the centre, averaged over r ∈ [L/4 − 4, L/4 + 4]. P_pol = 2 n sin(k₀) with k₀ the d-frame valley position is the value of ⟨P⟩ for a fully valley-polarized gas (about 2n).
- Two-valley momentum distribution: the fermion n(k) (c frame) from the one-body density matrix on the central half of the chain. W_a and W_b are its weights on cos k > 0 and cos k < 0.

First-stage criterion (C173 PROOF 5.4, as the swarm spec states it), at the main point φ = π/2, |w| = 1, n = 3/16. All five must hold.
1. Truncation below 1e-7 at L = 128.
2. c_two(64 → 128) > 1.2 for each of W = 1, 6, 12.
3. c_prof(128) > 1.2.
4. min(W_a, W_b)/max(W_a, W_b) ≥ 0.8 at L = 128 (two occupied valleys).
5. No long-range valley polarization. C_P(L/4) at L = 128 is below 0.05 P_pol², and it is smaller at L = 128 than at L = 64.

Calibration requirement. The estimator counts only if the controls read correctly at the same sizes: c_two(64 → 128) of the c = 1 control within 0.15 of 1, and of the c = 2 control within 0.25 of 2. If not, the pass or fail is reported as uncalibrated.

Registered predictions for the large sizes, stated now so the reading at L = 192 and 256 can be scored.
- If the point is a c = 2 phase: c_two(128 → 256) ≥ 1.6 for W = 6, and c_two(128 → 256) not smaller than c_two(64 → 128) by more than 0.3. C_P(r) decays at least as fast as 1/r.
- If it is a c = 1 phase with a pseudospin gap below the scale reached by C173 (an antiferromagnetic or polarized pseudospin): c_two falls monotonically with L, and c_two(128 → 256) < 1.4.
- If it is valley-polarized (vector chiral, c = 1): C_P(L/4) saturates at a nonzero fraction of P_pol² and does not fall from 128 to 256.

Density convention (decided before production, see PROGRESS). In the exact model μ = −2.2 is only 0.05 above the band bottom (−2.25), and the L = 48 scan gives N* = 6, n = 1/8, at μ = −2.2. C173's n = 0.17 at μ = −2.2 came from the penalty model (U = 40), whose virtual terms shift the density at this dilute point. So the runs are done at fixed even N, with n = 3/16 (closest to C173's reading density with an even integer N at L = 64 to 256) as the main point, n = 1/8 (the exact model's μ = −2.2) and n = 1/4 (about the exact model's μ = −2.0) as the two other densities. μ = −2.4 is below the band bottom and gives the empty chain in the exact model.

### Addendum, registered 2026-11-09 before job 17604009 started (coordinator note on C175's bound pair)

Correction to the density note above. The L = 48 statement "n = 1/8 at μ = −2.2" is about N/L of an open chain, whose edges are depleted. The timing run (local/t128.log, L = 128, N = 24) has a mid-chain density of 0.204 against N/L = 0.1875, and C175's ring ED gives n = 0.167 at μ = −2.2. So the penalty in C173 is not shown to shift the density, and the claim above is withdrawn. A series at n = 5/32 (L = 64, 128, 192, 256, bulk density expected near 0.17) was added. The registered first-stage criterion stays on n = 3/16. The same five tests are applied to n = 5/32 as a secondary reading, with the same thresholds.

Parity gap. Δ_P(L) = E(N+1) + E(N−1) − 2E(N) on open chains, N even, for n = 3/16, 5/32 and 1/8, with all three energies taken at the same χ (the largest χ common to the three runs, 512 for L ≤ 128 and 768 for L ≥ 192).
- Two-component liquid (c = 2): L Δ_P tends to a constant (finite-size spacing of the pseudospin sector).
- Pair liquid with a pseudospin gap Δ_s (c = 1): Δ_P tends to a nonzero constant of order Δ_s, so L Δ_P grows linearly once L ≫ v_s/Δ_s.
- Registered reading. If L Δ_P increases by more than 50 percent from the smallest to the largest L of a series and the increase steepens, that series is called pair-gapped at these sizes. If L Δ_P is flat or falling within 20 percent, the parity gap is gapless at these sizes. A gap of 1e-4 (C175's estimate at n ≈ 0.17) would not be visible below L ≈ 1e4, so a flat L Δ_P at n = 5/32 and 3/16 does not exclude it.

Pair against single-particle correlations. D(r) = ⟨(c_i c_{i+2})† c_j c_{j+2}⟩ (and separation 3), r = j − i, from the centre, against |G_f(r)| = |⟨c†_i c_j⟩|, both on the same chain.
- Pair liquid: |G_f(r)| decays exponentially (faster than any power) while |D(r)| decays as a power, with D(r) carrying the sign pattern of momentum π, so |D(r)/G_f(r)| grows with r.
- Two-component liquid: both decay as powers, so the log-log slopes of both are bounded and |D/G_f| does not grow exponentially.
- Registered reading at L = 256: compute the effective exponents of |G_f| and |D| (envelope fits over r ∈ [8, 64]) and the ratio of |G_f| envelopes between r = 16 and r = 64. A |G_f| envelope that falls by more than a factor 100 between r = 16 and 64 while |D| falls by less than a factor 10 is read as pairing.

## Status at write-up

The production runs are in Typhon job 17604009 (pending at write-up, scheduler estimate 2026-11-09 12:30). No production output has been read. What exists now is the code, its validation, and a local preliminary set at L = 48 and 64 (χ ≤ 256), which cannot test the registered criterion (it needs L = 128) and is reported only as orientation.

## Verdict

A blockade-exact, N-conserving DMRG for C173's model (3.2) now exists and is validated. The registered first-stage test of the candidate (C173 PROOF 5.4, at L = 128) is queued, not yet decided. At L ≤ 64 the converged readings (truncation 1e-8 to 4e-8 at χ = 256) are what C173 saw: two-size estimates near 2 on the main point at the bulk density of C173's reading, controls that read 1 and 2, no valley polarization, and a parity gap that scales like 1/L. These sizes cannot separate c = 2 from C175's pair liquid with a small pseudospin gap.

## Key claims

1. [proved] Label scheme (PROOF 1.1, 1.2). A bond label (N, n_{i−1}) used on both sides of the centre gives a non-orthonormal right basis, so its Schmidt values are wrong. Left-physical labels on left-canonical tensors and right-physical labels on right-canonical tensors give orthonormal bases on both sides. In that scheme the blockade and N are structural zeros, the effective Hamiltonian is P H P on the variational space, and the entanglement from the full mixed-label centre matrix is exact for the represented state.
2. [numerical, validate.py, validate2.py, validate3.py] Against C173's ED (unmodified), open chains L = 12 to 24, ten parameter sets including both controls and both gauges, the code agrees to machine precision. Energies agree to 3e-12, S(l) to 6e-11, densities to 5e-11, the full fermion one-body matrix to 3e-9, the ⟨n n⟩, ⟨P P⟩ (|i − j| ≥ 3) and pair correlators to 1e-9, and ⟨n_j n_{j+1}⟩ = 0 exactly.
3. [numerical, ringcheck.out, ringcheck_open.out] The open-chain bulk energy density at n = 1/6 is −0.37327. Ring ED (C173's ring_ed.py) gives −0.3722 to −0.3749 at L = 24 to 36, oscillating with the valley shell filling.
4. [numerical, scan_L48_phipi2.out; correction in the addendum] In the exact model the band bottom is −9/4. μ(N) at L = 48 shows an even-odd alternation of about 0.05. Open-chain edges are depleted, with n_mid ≈ N/L + δ/L and δ ≈ 2, so N/L is not the bulk density. My first note, that μ = −2.2 means n = 1/8 in the exact model, confused N/L with the bulk density and is withdrawn. μ = −2.4 is the empty chain [derived]. Runs therefore fix even N.
5. [numerical, preliminary, local/pre/summary.txt, χ = 256, two sizes only] The mid-chain entanglement (S̄ over W = 6 cuts each side) and the two-size estimates for 48 → 64 are:
   - c = 1 control: 0.98 (0.77 → 0.82).
   - c = 2 control: 1.92 (1.54 → 1.64).
   - φ = π/2, n = 1/8 (bulk density 0.164 → 0.154): 1.79 at W = 6 and 2.05 at W = 12, but 0.34 at W = 1. At this density the Friedel period is about 16 sites, so a three-cut window is meaningless here.
   - At a fixed bulk density, interpolating in ln n_mid across the φ = π/2 series (an unregistered estimator, PROOF 3), the reading is 2.26 (W = 6) and 2.67 (W = 12) at n_bulk = 0.17, and 2.08 and 2.71 at 0.19.

   A ratio L₂/L₁ = 4/3 makes all of these noisy at the ±0.3 level. They agree with C173 and say nothing yet about L ≥ 128.
6. [numerical, preliminary] Valley polarization. C_P(r), the connected correlation of the vector chirality P_j = c†_{j+1}c_j + h.c., is negative and decays at L = 64. For n = 3/16 it is −3.4e-2, −6.3e-3 and −1.7e-3 at r = 4, 8 and 16, against a long-range-order scale P_pol² ≈ 0.19. No vector-chiral order is visible. The negative sign (unlike valleys favoured at short distance) is what a valley-singlet pair liquid or an easy-plane pseudospin would give, and it is the opposite of the vector-chiral side. The valley weights are equal (W_a/W_b = 0.998 to 0.999), but at φ = π/2 the symmetry forces that for any unbroken state, so test 4 is informative only off φ = π/2.
7. [numerical, preliminary, local/pre/summary_parity_density.txt] The parity gap Δ_P = E(N+1) + E(N−1) − 2E(N), converged in χ to 1e-6, is:
   - n = 1/8: 0.0563 at L = 48 and 0.0375 at L = 64, so L Δ_P = 2.70 and 2.40, falling;
   - n = 1/6: 0.0581 at L = 48 (L Δ_P = 2.79).

   By the registered reading this is gapless at these sizes. C175's pair binding E_b = 0.0195 is the two-body scale. A many-body pseudospin gap below about 0.02 would not show at L ≤ 64, and C175's 1e-4 estimate at n ≈ 0.17 would not show below L ~ 1e4.

## The Typhon job (17604009)

One node, 48 cores, 320 GB, 3 days. It runs 72 tasks, 24 at a time with 2 BLAS threads each. Code and configs are in typhon/ (tasks.txt, cfg/, c176.slurm), on scratch at /scratch/lustre/jkudlerflam/c176/.
- φ = π/2, |w| = 1, fixed even N:
  - n = 3/16 at L = 64, 96, 128, 192, 256 (the registered main point);
  - n = 5/32 at L = 64, 128, 192, 256 (bulk density nearest C173's 0.17);
  - n = 1/6 at L = 48, 96, 192 (C173's sizes);
  - n = 1/8 at L = 48 to 256;
  - n = 1/4 at L = 64, 128, 256.
- φ = 1.54, 1.52, 1.50 at n = 3/16, L = 64, 128, 256 (complex arithmetic, χ ≤ 768).
- c = 1 control (hard rods, NN hop only) and c = 2 control (two XX legs) at n = 1/4, L = 48 to 256.
- Parity gap: N ± 1 for n = 3/16, 5/32 and 1/8 at every size, χ ≤ 512 (L ≤ 128) or 768.
- χ schedule: 16 sweeps at χ ≤ 64, then 128, 192, 256, 384, 512, 768, 1024 (two sweeps each, capped by L). Measurements at every sweep with χ ≥ 512: density, ⟨n_j n_{j+1}⟩, the one-body matrices on the central half (boson and fermion, hence n(k) and valley weights), C_nn and C_P from four central reference sites, and pair correlations D(r) with separations 2 and 3.
- None of this duplicates C173's job 17603887 (penalty model, grand canonical, χ ≤ 192) or C175's (C173's penalty code at fixed μ).

Expected outputs: /scratch/lustre/jkudlerflam/c176/out/<tag>.jsonl (one line per sweep), logs/<tag>.log, logs/c176_17604009.out, summary.txt (analyze.py: the registered criterion, the predictions at 128 → 256, the parity table, pair against single-particle correlations, the secondary n = 5/32 reading and the density-corrected estimator). Checkpoints in out/ckpt/ are large and need not be copied.

Collection (run from codex/work_C176/):

    mkdir -p typhon_out/out typhon_out/logs
    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c176/out/*.jsonl' typhon_out/out/
    scp -r -o LogLevel=ERROR typhon:/scratch/lustre/jkudlerflam/c176/logs typhon:/scratch/lustre/jkudlerflam/c176/summary.txt typhon_out/
    python3 analyze.py typhon_out/out > typhon_out/summary_local.txt

## What is left open

1. The registered first-stage test at L = 128 and the predictions at 128 → 256. All of it is in job 17604009.
2. Whether a pseudospin gap of C175's size (about 1e-4 at n ≈ 0.17) exists. Open chains at L ≤ 256 cannot see it. The job's discriminating content is at n = 1/8, where C175's dilute pairing should be strongest: the parity gap L Δ_P against L, and the decay of |G_f| against |D| at L = 256.
3. FES at χ ≥ 256 (C173's second stage) is not part of this section. The block-sparse code here is open-chain only.
