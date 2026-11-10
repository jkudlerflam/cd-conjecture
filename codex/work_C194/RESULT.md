This run was done by a Claude Opus 5.5 agent (2026-11-10), not by Codex (Astra out of credits).

# C194: follow the pseudospin gap in the effective model, and locate the boundary

Spec: `specs/C193_C194_sign_tests.md`, section C194. Model: C173 PROOF eq. (3.2) at bare Schrieffer–Wolff couplings (t₂ = −1, V₂ = −2, V₃ = −2/3), symmetric line φ = π/2, |w| = 1, blockade exact, C176's conventions. Derivations in PROOF.md, log in PROGRESS.md.

## Registration (copied from the spec before any production run, 2026-11-10)

- (G0) Calibration. On the two-leg XX control Δ_pb extrapolates to zero within the fit error and ξ_1/ξ_2 is stable over the last χ doubling. On the hard-rod control the estimators read their known values.
- (G1) Positive control. At one or more densities the gap is seen directly, meaning Δ_pb(L → ∞) is nonzero by more than three fit errors, and on the infinite chain ξ_1 grows by less than 1.15 over the last χ doubling while ξ_2 grows by more than 1.5.
- (G2) Candidate on the gapless side at a measured distance. G0 and G1 hold, both estimates of n_c lie below 0.11, no density at or above 0.11 shows a saturating ξ_1 or a nonzero Δ_pb, and K_s(n) rises with n from the boundary to 3/16.
- (G3) Gap risk stands. A saturating ξ_1 or a nonzero Δ_pb is seen at any n ≥ 0.11, or either estimate of n_c lies above 0.11, or K_s at n = 1/8 or 3/16 drifts below 1.
- (G4) C191's flow is checked separately. Its separatrix prediction passes if both estimates of n_c lie in 0.05 to 0.11, and fails otherwise. Say which.
- Anything between G2 and G3 is reported as undecided, with the numbers.

### Operational definitions fixed with the registration (before any production output)

These make the registered words computable. Fit ranges and windows may still be fixed on the controls (the spec allows it), and are then frozen in the section "Frozen after the controls" below.

1. Δ_pb(n, L) = E(N + 1) + E(N − 1) − 2E(N) at even N, all three energies at the same final χ, on open chains. The bulk density of an open chain is not N/L (the edges hold fewer particles), so N is chosen at each L as the even integer nearest to n·L, and the realized N/L is reported.
2. Extrapolation. Δ_pb(L) = Δ_∞ + a/L + b/L², least squares over the available sizes (PROOF §2 derives the form). The fit error is the standard error of Δ_∞ from the fit, enlarged by the spread of Δ_∞ over the two sub-fits (drop the smallest L; drop the b term). "Nonzero" means Δ_∞ > 3 × that error.
3. Infinite chain. ξ_1 is the charge-1 correlation length with the fermion string (C189's T1f), ξ_2 the charge-2 one (T2), ξ_0 the neutral one, all from C189's `run_inf.py` with its convergence rule. The "last χ doubling" is between the largest χ and half of it. "Saturating ξ_1" is the G1 condition on ξ_1 (growth below 1.15) with ξ_2 growing by more than 1.5. Every run reports its iteration count and whether it stopped on tolerance or on the cap.
4. K_s(n) on the gapless side. C192's estimators as frozen in §9.161 and §9.162: the boson one-body exponent η_b at the valley momentum and the inter-valley pair exponent at momentum π, both with the open-chain conformal forms (`offline.py`, largest window, spread over the four windows as the error), and K_s = 1/(4η_b − η_pair) on the symmetric line. "Drift" is K_s(256) − K_s(128).
5. n_c from the gapped side. Fit Δ_∞(n) and 1/ξ_1(n) to the KT form A exp(−B/√(n_c − n)) (PROOF §3), with n_c, A, B free, error from the fit covariance and from dropping one point. From the gapless side: fit K_s(n) linearly (and quadratically as a check) near the boundary and find where it crosses the separatrix value K_s* = 1 + |y_1|/2, with y_1 from C191's scattering data, and also where it crosses 1. Both crossings are reported.
6. "Rises with n" in G2 means K_s at every measured gapless density is not smaller than at the next lower one by more than its error.

## Frozen after the controls

Written 07:25, before any candidate output was opened (`controls_exact.py`, `ctl_check.py`, `out_ctl/`).

- Sizes. At each density the sizes are those in 48 to 256 at which n·L is an even integer (`tasks_main.json`), so N/L = n exactly at every size: 1/32 and 1/16 and 3/32 at L = 64, 128, 192, 256; 1/24 and 1/12 at 48, 96, 144, 192, 240; 1/20 at 80, 120, 160, 200, 240; 3/40 at 80, 160, 240; 1/10 at 60, 120, 180, 240. The densities 1/8, 5/32, 3/16 are taken from C176's converged runs (L = 48 or 64 to 256, χ to 768), not rerun.
- Sweeps. χ = 16 (6 sweeps), 32 (4), 64 (4), 96 (3), 128 (4), plus 192 (3) at L ≥ 160 and n ≥ 1/16. Convergence is reported per run as the change between the last two sweeps and between the last χ = 128 and χ = 192 sweeps.
- Δ_pb fit. As registered (definition 2), all sizes. On exact energies of both controls at n = 1/32 to 1/4 it gives |Δ_∞| ≤ 8e-5 with errors 7e-5 to 2.3e-3 (0.0 to 0.1 errors from zero). The fit error at low density is dominated by 1/L² curvature and is the honest resolution: about 2e-3 (two-leg XX) and 6e-4 (hard rods) at n = 1/32.
- DMRG against exact. Hard rods: energies to 3e-10 at L = 64, 128, Δ_pb to 3e-11. Two-leg XX: to 1e-12 at low density and 4e-5 at L = 128, n = 1/8, χ = 128 (Δ_pb off by 1.5e-5, L·Δ_pb off by 0.002). A c = 2 state at L = 128 and χ = 128 is thus accurate to about 1e-5 in Δ_pb, far below the fit errors.
- K_s. K_s = 1/(4η_b − η_pair) with `offline.py` (C192): η_b the boson exponent at the valley momentum, η_pair the mean of the d = 2 and d = 3 inter-valley pair exponents at momentum π. Central value window 4 (r from 24 to 0.7W for η_b, 24 to 0.4L for the pair), error the spread over windows 2 to 4. On the two-leg XX control (C176's c2 runs) the same fits read η_b = 0.498, 0.500 and the inter-leg pair (d = 3) 0.988, 0.998 at L = 128, 256 (exact 1/2 and 1). On that control the d = 2 pair is the same-leg pair (2.0), so the d = 2/d = 3 average applies only to the candidate, where both are inter-valley (C192).
- Infinite chain. C189's `run_inf.py` unchanged, with max_it raised to 800 and χ = 16 to 192 locally. The two-leg XX control cannot run in this code (legs conserve number separately, C189), so the G0 ratio check uses C189's free two-valley run (ξ_1/ξ_2 changes by 1.005 over 384 → 768) and its hard-rod run, plus a hard-rod rerun here with the higher cap if cores allow.

## Results (local machine only, 2026-11-10, stopped 11:41)

Every number below comes from this machine (`out_main/`, `inf/out/`) or from C176's stored runs (marked). Nothing from Typhon or from the desktop copy had finished when this was written. Output tables: `gap_table.out`, `gap_chi.out`, `inf_table.out`, `ks_table.out`, `nc_fit.out`, `nc_fit_minN.out`, `nc_ks.out`, `predict194.json`.

### Step 1. Predictions and the dilute check

C191's flow on the grid (`predict194.py`): ξ_s = 100 to 140 sites at n = 1/32, 230 to 300 at 1/24, 700 to 800 at 1/20, 2e4 to 6e4 at 1/16, beyond 1e7 at 0.07, gapless from 0.075 to 0.08 on (variant-dependent). The dilute limit (PROOF §1) is E_b − 2E_3 = 0.0183 ≤ Δ_pb(n → 0) ≤ E_b = 0.0195, pair size 21 sites.

### Step 2. Direct gap on open chains (registered fit, all sizes)

| n | sizes | Δ_∞ | error | Δ_∞/error | N > 2 sizes only (unregistered) |
|---|---|---|---|---|---|
| 1/48 | 96, 192, 288 | 0.0129 | 0.0028 | 4.6 | (two sizes) |
| 1/32 | 64 to 256 | 0.0085 | 0.0029 | 2.9 | 0.0094 ± 0.0021 (4.4) |
| 1/24 | 48 to 240 | 0.0049 | 0.0034 | 1.5 | 0.0056 ± 0.0021 (2.6) |
| 1/20 | 80 to 240 | 0.0028 | 0.0008 | 3.3 | same |
| 1/16 | 64 to 256 | 0.0013 | 0.0013 | 1.0 | same |
| 1/8 (C176) | 48 to 256 | 0.0020 | 0.0056 | 0.4 | |
| 5/32 (C176) | 64 to 256 | −0.0005 | 0.0007 | −0.8 | |
| 3/16 (C176) | 64 to 256 | 0.0031 | 0.0069 | 0.4 | |

- L·Δ_pb rises with L at low density (1/32: 2.48, 2.90, 3.42 at L = 128, 192, 256; 1/48: 2.45, 3.34, 4.46 at 96, 192, 288), is nearly flat at 1/16 (2.19, 2.22, 2.28) and falls at 1/8 (2.28 to 2.18). The registered error is dominated by the sub-fit without the 1/L² term, which the N = 2 sizes (pair size 21 comparable to L/N) pull hard.
- Dilute check. Δ_∞ = 0.0129 (1/48) and 0.0085 to 0.0094 (1/32), linear in n, extrapolate to about 0.020 at n = 0, against 0.0183 to 0.0195. Consistent.
- Cross-check with the infinite chain. At 1/32, Δ_s = Δ_∞/2 ≈ 0.0045 and C191's u_s ≈ 0.24 give ξ_s ≈ u_s/Δ_s ≈ 53 sites, against the saturating ξ_1 ≈ 57 to 61 below.
- Convergence. Last-sweep energy changes are below 4e-6 everywhere. The χ dependence is not negligible at the largest L: L·Δ_pb falls with χ (1/20, L = 240: 2.658, 2.534, 2.517 at χ = 64, 96, 128; 1/32, L = 256: 3.633, 3.434, 3.418), so the χ = 128 values carry an upward bias growing with L, of order 0.01 to 0.02 in L·Δ_pb (a few 1e-5 in Δ_pb) at the largest sizes. That is small against the fit errors, but it acts in the direction of a fake gap. Unverified at χ = 384 (those states are in the Typhon and desktop payloads).

### Step 3. Infinite chain (C189's code, all stages stopped on tolerance, none on the cap)

| n | χ range | ξ_1/ξ_2 against χ | last doubling: ξ_1, ξ_2 growth | c_fit (four largest χ) |
|---|---|---|---|---|
| 1/32 | 16 to 128 (192 snapshot) | 0.70, 0.64, 0.61, 0.35, 0.29, 0.21, 0.15 (χ = 16 to 128); 0.106 at 192 | 64 → 128: 1.19, 2.25; 96 → 192 snapshot: 1.12, 2.17 | 1.08 |
| 1/20 | 16 to 192 | 0.83, 0.77, 0.69, 0.58, 0.47 (χ = 48 to 192) | 96 → 192: 1.45, 2.14 | 1.29 |
| 1/16 | 16 to 384 | 1.14, 1.13, 0.97, 0.93, 0.84 (χ = 96 to 384) | 192 → 384: 1.58, 1.80 | 1.72 (1.52 with ξ_2) |
| 3/40 | 16 to 96 | 1.47, 1.45, 1.51 (χ = 48 to 96) | 48 → 96: 1.76, 1.72 | 2.01 |
| 1/12 | 16 to 256 | 1.61, 1.67, 1.72, 1.64, 1.68 (χ = 64 to 256) | 128 → 256: 1.82, 1.86 | 2.01 |

- The 1/32 χ = 192 entry is a snapshot (`snap_inf.py`) of the checkpoint at iteration 392, 43 iterations into the stage, not converged by run_inf's rule (two-iteration entropy change about 2e-4 against the tolerance 1e-5, energy per site steady to 3e-10). ξ_1 = 60.5, ξ_2 = 571. Unverified until the stage converges.
- Largest ξ reached: 1/32 ξ_2 = 375 (571 snapshot), 1/20 273, 1/16 338, 3/40 114, 1/12 266. For comparison C189 has 1/8 and 3/16 at ξ ≈ 510 to 590 with ξ_1/ξ_2 stable to 1% and c_fit 2.0.
- Reading. A gap signature (ξ_1 slowing, ξ_2 growing, the ratio falling) is clear at 1/32 and 1/20, starts at 1/16 only beyond ξ ≈ 180, and is absent at 3/40 (short ξ) and 1/12 (to ξ ≈ 270). So the infinite chain brackets the boundary between 1/16 and roughly 1/12, with the usual KT caveat that 3/40 and 1/12 may still be gapped beyond the ξ reached. The ratio drift is an unregistered diagnostic. The registered 1/ξ_1 KT fit cannot be done, since only 1/32 saturates.
- Against C191. Its ξ_s is too long at 1/32 (100 to 140 predicted, ξ_1 ≈ 60 measured) and much too long at 1/16 (2e4 to 6e4 predicted, a drift already visible at ξ ≈ 200 to 340). Its separatrix location is not contradicted.

### Step 4. K_s on the gapless side (frozen estimator)

| n | L | K_s (window 4) | error (spread windows 2 to 4) | windows 1 to 4 | source |
|---|---|---|---|---|---|
| 1/12 | 144 | 0.89 | 0.25 | 1.66, 1.39, 1.20, 0.89 | here, χ = 128, N = 12 |
| 3/32 | 128 | 0.95 | 0.25 | 1.73, 1.44, 1.24, 0.95 | here, χ = 128, N = 12 |
| 1/10 | 120 | 0.99 | 0.23 | 1.78, 1.45, 1.27, 0.99 | here, χ = 128, N = 12 |
| 1/8 | 128, 256 | 1.286, 1.265 | 0.05, 0.03 | | C176 |
| 5/32 | 128, 256 | 1.408, 1.383 | 0.02, 0.02 | | C176 |
| 3/16 | 128, 256 | 1.500, 1.462 | 0.02, 0.02 | | C176 |

- Drift K_s(256) − K_s(128): −0.021 (1/8), −0.025 (5/32), −0.037 (3/16). Small, downward, nowhere near 1.
- Below 1/8 only L = 120 to 144 exist (the L = 192 to 256 measured states did not run), with N = 12. There K_s falls steeply from the short to the long window, which is either the running of K_s toward the boundary or a failure of the estimator at 12 particles. The estimator was calibrated at N ≥ 16, so these three values are unvalidated.
- K_s rises with n from 1/12 to 3/16 (0.89, 0.95, 0.99, 1.29, 1.41, 1.50 at L ≈ 128), within errors.

### Step 5. The boundary

- Gapped side, registered KT fit of Δ_∞(n): n_c = 0.138, with the χ² + 1 interval 0.066 to 0.200, the whole grid. The fit is not constrained (two densities above three errors, and the leave-one-out fits are undefined). The N > 2 variant gives 0.179, interval 0.056 to 0.200. The 1/ξ_1 fit is not possible.
- Gapless side (K_s at L ≈ 128, weighted fits, Monte Carlo errors): crossing K_s = 1 at 0.050 ± 0.020 (linear) or 0.090 ± 0.016 (quadratic), crossing 1 + |y_1|/2 at 0.076 ± 0.012 (linear) or 0.098 ± 0.013 (quadratic). By PROOF §3 (ii) the first pair are lower estimates and the second upper ones, so n_c ≈ 0.05 to 0.10 from this side, with the fit form as the dominant error and the low-n inputs unvalidated.
- Infinite chain (unregistered ratio diagnostic): 1/16 < n_c ≲ 1/12, i.e. 0.0625 to about 0.083.
- C191: 0.07 to 0.085.
- Distance of the candidates. n = 1/8 sits about 0.04 to 0.06 above the bracket from the infinite chain and the K_s side, with K_s = 1.27 to 1.29 and KT margin 2(K_s − 1) − |y_1| ≈ 0.47. n = 3/16 sits about 0.1 above it, with K_s = 1.46 to 1.50.

### Registered outcome

- G0. Partly met. Δ_pb on both controls extrapolates to zero within the fit error (exact energies and DMRG to 1e-5). ξ_1/ξ_2 stability rests on C189's free two-valley run (1.005), since the two-leg XX control cannot run in the infinite code. The hard-rod infinite rerun with the higher cap was not done locally (it is in the Typhon payload), so the c = 1 side of the infinite-chain calibration is still C189's 0.847 reading.
- G1. Met only across densities, provisionally. The finite chain shows Δ_∞ above three errors at 1/48 (4.6) and 1/20 (3.3), and 1/32 at 2.9 (4.4 without the N = 2 sizes). The infinite chain meets the ξ_1 < 1.15, ξ_2 > 1.5 condition only at 1/32 and only in the unconverged χ = 192 snapshot (1.12, 2.17). No single density meets both parts on converged data.
- G2. Not met, since the gapped-side estimate of n_c is not below 0.11 (it is unconstrained).
- G3. Formally triggered by the same estimate: its central value, 0.138, lies above 0.11. This does not discriminate, because its interval is the whole grid. Every other G3 condition fails. No density at or above 0.11 shows a saturating ξ_1 or a nonzero Δ_pb (C176 and C189 data), and K_s at 1/8 and 3/16 drifts by −0.02 and −0.04 from 1.29 and 1.50.
- G4. Fails formally (the gapped-side estimate is not inside 0.05 to 0.11). The gapless-side estimates (0.05 to 0.10) and the infinite-chain bracket (0.0625 to 0.083) do lie inside it.
- Reading, stated separately from the letters. The test now has a positive control. A gap is seen at n ≤ 1/20 on both the finite and the infinite chain, with the right dilute limit and a consistent ξ_s. Every density at or above 1/12 reads gapless to the ξ reached. The boundary sits between 1/16 and about 1/12, below the candidates by about 0.04 at n = 1/8. The registered gapped-side KT fit cannot place n_c from the gaps visible at L ≤ 256, which is the KT problem the brief anticipated. Taken as registered, the outcome is G3 by the letter, and undecided between G2 and G3 in substance.

### What failed or is missing

- No converged infinite-chain stage meets G1. The 1/32 χ = 192 stage was stopped unconverged.
- The KT fit from the gapped side is unconstrained, and 1/ξ_1 has one saturating density.
- K_s below 1/8 comes from three L ≈ 128, N = 12 runs at χ = 128 with large window spread, and there is no L = 256 drift there.
- The local finite queue did not reach 3/40, 1/12, 1/10, 3/32 Δ_pb at the larger sizes, or the 1/64 sizes beyond L = 128. Of the high-density L ≥ 144 states only 1/12 at L = 144 ran.
- Local runs were stopped at 11:41. The iDMRG checkpoints in `inf/out/ckpt/` resume with `python3 -E -s run_inf.py inf/cfg/<tag>.json inf/out`. The finite pool resumes with `python3 -E -s run194.py tasks_main3.json out_main 3` (finished states are skipped, unfinished ones restart).

## Typhon

- Job 17611350 (c194), submitted 07:19, package in `typhon/` (one run script `typhon/tasks.txt` driven by `typhon/c194.slurm`), scratch `/scratch/lustre/jkudlerflam/c19x/c194/`. At 09:29 it was pending on priority with a scheduler start estimate of 14:43 EST today (it read 2026-11-14 at 07:20).
- Overlap step, not started. The coordinator suggested (09:28) starting the same payload now as an `srun --overlap` step inside the running c189_s316 allocation 17606915. I did not. Both my brief and the spec's common rules say not to start steps inside other allocations, and a relayed request is not a change to those rules. The job stays queued. If the coordinator decides to place it, the single command is, run on the login node and on its own line:

  `nohup setsid srun --jobid=<donor> --overlap --nodes=1 --ntasks=1 --cpus-per-task=45 bash -c 'cd /scratch/lustre/jkudlerflam/c19x/c194 && mkdir -p logs out_inf out_fin && cat typhon/tasks.txt | xargs -P 10 -I CMD sh -c CMD' > /scratch/lustre/jkudlerflam/c19x/c194/logs/overlap.log 2>&1 &`

  Caveats for whoever places it. The iDMRG runs checkpoint every 10 iterations and resume (C189's `run_inf.py`), but the finite-chain pool (`run194.py`) restarts an unfinished state from scratch, so a step that is cancelled loses its finite runs in progress. The step has no self-cancel when 17611350 starts, so it would need the TYPHON_OVERLAP.md poll wrapper, otherwise two copies of the payload would write to the same files. Peak memory has not been measured here (unverified; the χ = 512 iDMRG at q = 32 to 40 is the largest item).
- At 11:41 job 17611350 was still pending, with a start estimate of 16:32 EST today.
- Collection (the multi-source form tested at 11:43 on the two folders that exist already; the brace form does not work with this scp), run from `codex/work_C194/`:

  `scp -o LogLevel=ERROR -r "typhon:/scratch/lustre/jkudlerflam/c19x/c194/out_inf" "typhon:/scratch/lustre/jkudlerflam/c19x/c194/out_fin" "typhon:/scratch/lustre/jkudlerflam/c19x/c194/logs" "typhon:/scratch/lustre/jkudlerflam/c19x/c194/summary_inf.txt" typhon_out/`

  then `python3 -E -s inf_table.py inf/out typhon_out/out_inf`, `python3 -E -s gap_table.py out_main typhon_out/out_fin` and `python3 -E -s ks_table.py` (which reads `typhon_out/out_fin` itself). The desktop copy is read the same way, with `../work_C194_desktop/out_inf` and `../work_C194_desktop/out_fin` as extra directories; numbers from it must be labelled by machine. Where a state exists in two places, gap_table.py keeps the one at larger χ; report the energy difference.
