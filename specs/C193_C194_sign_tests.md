# C193 and C194: the sign of the pseudospin coupling, on the qubit chain and across the gap boundary (2026-11-10)

Coordinator: Claude Fable 5. Jonah approved both on 2026-11-10. Context: REPORT §9.139, §9.160 to §9.163, §9.174 to §9.180, STATUS item 1, `specs/C189_C191_decisive_tests.md`.

## Where things stand

The full qubit chain reads c ≈ 2 from L = 128 to 256 at six tuned bonds (C190, §9.178), and C173's effective model reads c = 2.0 on the infinite chain out to ξ ≈ 500 to 600 sites (C189, §9.179). Larger sizes cannot exclude the one alternative left, a pseudospin gap far below the finite-size spacing. What can decide it is the sign of the pseudospin backscattering flow.

- C191 (§9.160) predicts K_s = 1.34 to 1.67 at n = 1/8 to 3/16 on the symmetric line, with the separatrix at n ≈ 0.07 to 0.085 and a gap below it.
- C192 (§9.161) measures K_s ≈ 1.2 to 1.5 in the effective model from the boson one-body exponent η_b and the inter-valley pair exponent 1/K_c, through 1/K_s = 4η_b − 1/K_c. The exponents move by at most 0.004 from L = 128 to 256.

Two things are missing. First, K_s has only been measured in the effective model, not on the qubit chain. Second, the gap has never been seen where it is predicted to exist, so the tests have no positive control and the boundary has never been located. A measured K_s > 1 at finite size is also not sufficient by itself, since a strong enough backscattering still flows to a gap. The drift of K_s with scale bounds that coupling, and neither test has used it yet.

## Common rules

- Work only in your own folder, `codex/work_C193/` or `codex/work_C194/`. Everything else is read-only. Keep RESULT.md (main file, with the registration), PROOF.md (derivations, every step) and PROGRESS.md (log).
- Copy the registration below into RESULT.md before you open any production output that decides it. You may fix fit ranges and windows on the controls first, and then they are frozen. Report the registered outcome by its letter, also when you think it does not discriminate, and say why separately.
- Physics standard (RG, perturbation theory, calibrated numerics). Show every step and run consistency checks, with limits and controls of known exponents.
- Local machine: at most 8 cores in total. Don't pip install anything. numpy and scipy are available. Use `python3 -E -s` when a script imports sibling modules.
- Typhon: `ssh -o LogLevel=ERROR -o ConnectTimeout=20 typhon`, retried with waits of about 30 s inside the command (the bastion rate-limits, and a local hook blocks `sleep N; command` chains). SLURM only, never compute on the login node. Partition typhon, CPU. Job names must start with c193 or c194. Scratch `/scratch/lustre/jkudlerflam/c19x/<yours>/`. Submitting is always OK. Never cancel, hold, release, pause or signal any job or step you did not submit. Do not start steps inside other allocations. The CPU queue is days deep, so plan on local runs for the first answer.
- A hook blocks `grep -r`, `find .` and `ls -R`. Use explicit paths. The local shell is zsh, so quote globs in scp.
- Don't commit or push. No probabilities in any file. No absolute local paths in files you write. No Fable subagents. Astra (gpt_solve.sh) is out of credits, don't call it.
- Prose: commas, not em-dashes, few colons, no italics, American spelling.
- Final report: what you established with numbers, the registered outcome, what failed, what is left, and any Typhon job ids with tested collection commands.

## C193: K_s on the qubit chain, from the stored C190 states

Goal: measure η_b, 1/K_c and K_s on the nearest-neighbour qubit chain itself, at all six tuned bonds, and bound the backscattering from the drift with L. No new DMRG runs are needed.

Data.
- GPU runs: `codex/work_C190/gpu/typhon_out/outG/<key>_L<L>_P<0|1>_<init>.json`, keys p153c, sym40_n316, rat_n316, sym40_n18, p153_n316, sym80_n316, at L = 64, 128, 256. Use only files with status done.
- CPU runs and both controls (single XX chain and two-leg XX): `codex/work_C190/typhon_out/outA/` and `outB/`.
- Each file has `meas` with `Gb_re`, `Gb_im`, `Gf_re`, `Gf_im` (one-body matrices on a 48-site window given by `gwin`), `CP` (pair function from site `CP_i0`), `dens`, `Z`, `P_bond`, `e_bond`. Read the measurement code in `codex/work_C190/z2dmrg.py` and `run190.py` to see which operators these are, and `codex/work_C190/PROOF.md` and RESULT.md for the map from qubit operators to the effective model.
- Effective-model references: `codex/work_C192/CHECKPOINT.md`, `onebody.py`, `pairconf.py`, `sq2.py`, and C191's frame (§9.160, `codex/work_C191/PROOF.md`).

Steps.
1. Derive in PROOF.md which effective operators the stored qubit correlators represent at leading order in Δ/V (boson, Jordan–Wigner fermion, inter-valley pair), which momenta they carry on and off the symmetric line, and the open-chain conformal forms used for the fits.
2. Calibrate on the two controls with the same 48-site window. Exact values: two-leg XX has η_b = 1/2 and pair exponent 1, the single XX chain has its own known exponents (state them). Freeze the fit ranges here.
3. Apply to the six keys, both parity sectors and the ground branch, at L = 128 and 256 (L = 64 as a check). Report η_b, the pair exponent, 1/K_c and K_s = 1/(4η_b − 1/K_c), with the spread over fit ranges as the error.
4. Compare with the effective model at matched density: η_b = 0.326 at n = 3/16 and 0.351 at n = 1/8, pair exponent 0.61 to 0.64 (§9.161).
5. Drift. Write the Kosterlitz–Thouless flow for (y_∥, y_1) with y_∥ ∝ K_s − 1 in PROOF.md, with the normalization fixed. From K_s at L = 128 and 256, bound |y_1| through dy_∥/dl = −y_1², and state whether y_∥ > |y_1| holds (flow to the gapless line) with the measured numbers and their errors. Do the same for the effective model from C192's numbers as a cross-check.

Registration (copy into RESULT.md before opening any candidate correlator).
- (Q0) Calibration. With the frozen ranges the two-leg control reads η_b within 0.03 of 1/2 and the pair exponent within 0.08 of 1 at L = 128 and 256. If not, everything below is uncalibrated.
- (Q1) Tracks the effective model. At L = 256 on the ground branch, every n = 3/16 key has η_b within 0.05 of 0.326 and the n = 1/8 key within 0.05 of 0.351, and every pair exponent lies in 0.50 to 0.75.
- (Q2) Gapless sign. K_s ≥ 1.1 at L = 128 and 256 for every key in both sectors, K_s(256) − K_s(128) ≥ −0.05, and the drift bound gives y_∥ > |y_1| at every key.
- (Q3) Gap sign. Any key has K_s(256) < 1, or K_s falls by more than 0.1 from 128 to 256, or the drift bound allows |y_1| ≥ y_∥.
- Anything between Q2 and Q3 is reported as undecided, with the numbers.
- If the 48-site window is too short for a stable η_b on the controls, say so, stop the candidate analysis for that estimator, and give the exact command that would recompute the correlators from the Typhon checkpoints under `/scratch/lustre/jkudlerflam/c18x/c190g/` (read only) in a SLURM job named c193.

## C194: follow the gap in the effective model, and locate the boundary

Goal: see the pseudospin gap where C191 predicts one, locate the boundary in density on the symmetric line at |w| = 1, and show on which side of it and how far from it the candidate densities n = 1/8 and 3/16 sit.

Tools and inputs.
- The model is C173 PROOF eq. (3.2), with C176's blockade-exact conventions (`codex/work_C176/fibdmrg.py`, `analyze.py`).
- Infinite chain: `codex/work_C189/idmrg.py`, `ivumps.py`, `run_inf.py`, `summarize_inf.py` (U(1) at fixed density, correlation lengths by charge sector). Validated in `validate_inf.py`.
- K_s estimators: `codex/work_C192/onebody.py`, `pairconf.py`, `sq2.py`, and `codex/work_C189/kest*.py`.
- C191's flow and its predictions: `codex/work_C191/RESULT.md`, `ktflow.py`, `predict.json`. C175's two-body data: `codex/work_C175/RESULT.md`.

Steps.
1. From C191 and C175, tabulate the predicted pseudospin gap and correlation length ξ_s against n below the separatrix, and the two-body binding energy and pair size. Choose the densities where the predicted ξ_s is within reach (at most a few hundred sites). In the dilute limit the gap should tend to half the two-body binding energy. Derive that limit in PROOF.md and use it as a check.
2. Direct gap on open chains. With the finite-chain code, measure the pair-breaking gap Δ_pb(n, L) = E(N + 1) + E(N − 1) − 2E(N) at even N, at L = 64, 128, 192, 256, for a density grid from about 1/32 to 3/16 that includes 1/16, 1/12, 1/10, 1/8, 5/32 and 3/16. Separate the gap from the compressibility term, which scales as 1/L (derive the form, and check it on the two-leg XX control where the gap is zero and on the hard-rod control). Converge each run and report residuals.
3. Infinite chain. At the same densities, run the χ sequence and report the correlation lengths in the charge-1 sector (ξ_1), the pair sector (ξ_2) and the neutral sector, with c_fit. A pseudospin gap shows as ξ_1 saturating while ξ_2 keeps growing, with the slope tending to 1/6.
4. K_s(n) on the gapless side, with the same estimators and frozen ranges as §9.161, at L = 128 and 256, and its drift.
5. Locate the boundary n_c two ways. From the gapped side, fit Δ_pb(n) and 1/ξ_1(n) to the Kosterlitz–Thouless form. From the gapless side, find where K_s(n) extrapolates to its critical value. Give both with errors and compare with C191's 0.07 to 0.085.
6. If no reachable density shows a gap, build a positive control by adding an inter-valley attraction to (3.2) until the gap is visible, show that every estimator detects it, and follow it back toward zero added attraction at n = 3/16. Say plainly that the density scan then has no positive control.
7. Prepare a Typhon package (one run script, one sbatch file named c194) for the points that need χ above what the local machine reaches, submit it, and put the single command an overlap step would need in RESULT.md. The coordinator decides where to place it.

Registration (copy into RESULT.md before any production run).
- (G0) Calibration. On the two-leg XX control Δ_pb extrapolates to zero within the fit error and ξ_1/ξ_2 is stable over the last χ doubling. On the hard-rod control the estimators read their known values.
- (G1) Positive control. At one or more densities the gap is seen directly, meaning Δ_pb(L → ∞) is nonzero by more than three fit errors, and on the infinite chain ξ_1 grows by less than 1.15 over the last χ doubling while ξ_2 grows by more than 1.5.
- (G2) Candidate on the gapless side at a measured distance. G0 and G1 hold, both estimates of n_c lie below 0.11, no density at or above 0.11 shows a saturating ξ_1 or a nonzero Δ_pb, and K_s(n) rises with n from the boundary to 3/16.
- (G3) Gap risk stands. A saturating ξ_1 or a nonzero Δ_pb is seen at any n ≥ 0.11, or either estimate of n_c lies above 0.11, or K_s at n = 1/8 or 3/16 drifts below 1.
- (G4) C191's flow is checked separately. Its separatrix prediction passes if both estimates of n_c lie in 0.05 to 0.11, and fails otherwise. Say which.
- Anything between G2 and G3 is reported as undecided, with the numbers.
