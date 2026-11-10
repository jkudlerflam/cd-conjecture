# C195 and C196: a third size and the infinite qubit chain (2026-11-10)

Coordinator: Claude Fable 5. Jonah approved both on 2026-11-10. Context: REPORT §9.178 to §9.181, `codex/work_C190/` and `codex/work_C190/gpu/` (solver, registration C1 to C6, V1, V2), `codex/work_C189/` (the infinite chain in the effective model), `codex/work_C177/` (an earlier infinite-chain code on the qubit chain, not converged).

Aim: the two pieces of evidence a referee would ask for beyond L = 256. Both run on the qubit chain itself.

## Common rules

- Work only in `codex/work_C195/` or `codex/work_C196/` (yours) and `/scratch/lustre/jkudlerflam/c19x/c195/` or `c196/`. Everything else is read only, in particular `codex/work_C190/` and everything under `/scratch/lustre/jkudlerflam/c18x/`. Copy code you need, never edit it in place.
- Keep RESULT.md (registration first), PROOF.md (derivations), PROGRESS.md (log).
- Copy your registration below into RESULT.md before any production output that decides it is opened. Do not change a criterion afterwards. Additions are dated addenda.
- Typhon: SLURM only, never compute on the login node. Partition apollo (GPU). Job names start with c195 or c196. One run per GPU, about 4 cores per GPU. complex128 throughout, no change to any solver tolerance. No pip or conda install. Do not use h200.
- Never cancel, hold, release, signal or alter a job you did not submit.
- Never copy checkpoints to the local machine. No absolute local paths, no probabilities, in any file. Do not commit or push.
- Local machine: numpy only, at most about 4 cores.

## C195: L = 512 at one bond

Key p153c (the collaborator's field h_z = −19.97589125), both parity sectors, the GPU code `codex/work_C190/gpu/z2dmrg_xp.py` with the conventions of `pipeline190.SCHED`, started from scratch with the registered initialization. At most 2 GPUs.

Known obstacle: at χ = 768 and L = 256 the device holds 27 GB of environments on a 40 GB card. L = 512 at χ ≥ 1024 does not fit. Keep only the active environments on the device and the rest on the host (up to 160 GB of host memory per run is allowed for this task) or on scratch. Gate before production: the modified code reproduces an L = 128 run of p153c sweep by sweep to 1e-10 in energy and 1e-8 in entropy against the unmodified GPU code.

Schedule: χ stages rising to a final χ of at least 1024, with a 1536 stage added if the benchmark says it completes inside the 7-day limit. The job must resume from its checkpoint, and a chained resubmission must be in place.

Registration.
- (T1) Convergence at the final χ: residual ≤ 1e-8, last two sweeps |ΔE| ≤ 2e-8 (the L = 256 bound scaled with L), |ΔS̄_W| ≤ 1e-3, |Δn_mid| ≤ 5e-4, discarded weight ≤ 1e-7. The step of S̄₆ between the last two χ stages is reported against 0.01 (C4), and a failure of C4 alone is reported as in §9.178 and does not by itself void T2 or T3.
- (T2) Pass. On the ground branch c_two(256 → 512) ≥ 1.6 at W = 6 and W = 12, with S̄(256) from the highest χ available at L = 256 for this key, and T1 holds. Also reported, not deciding: c_two at matched χ, c_two from S̄ extrapolated in discarded weight at both sizes, and c_prof(512).
- (T3) Fail. c_two(256 → 512) ≤ 1.3 at W = 6 and W = 12 with T1 holding.
- (T4) Density. n_mid at L = 512 within 0.006 of the L = 256 value of the same key.
- Anything between T2 and T3 is undecided and reported with the numbers.

## C196: entanglement scaling on the infinite qubit chain

An infinite-system code for the nearest-neighbour qubit chain with exact Z2 parity, at the fields of p153c and sym40_n316 (and the rational bond if time allows). The energy scale of the low-energy sector is of order Δ²/V, about 1e-3 of the largest coupling, which is why C177's attempt did not converge. Use what made C190 converge (its preconditioned Davidson, PROOF.md section 5 of C190) or an equivalent, and say which. At most 3 GPUs, or CPU inside the same jobs if that is faster for small χ.

Measure S against ln ξ over a χ sequence up to at least 512 (768 if feasible), with ξ from the transfer matrix, and the leading correlation lengths in each parity sector separately. Report the bulk energy density and the density, and compare them with the bulk of the C190 open chains at the same field.

Registration.
- (I0) Controls with the same code, χ sequence and fit window: the XX chain reads c within 0.1 of 1, and the two-leg XX control reads c within 0.15 of 2. A gapped control (any bond with a known finite ξ) shows ξ saturating.
- (I1) Convergence, per point: energy density change over the last iterations below 1e-10 per site and gradient or residual norm below 1e-7. Points that miss this are flagged and excluded from the fit.
- (I2) Pass. At each candidate bond the slope of S against ln ξ over the last three converged χ gives c ≥ 1.7, ξ grows by at least a factor 1.5 over the last χ doubling in both parity sectors, and I0 holds.
- (I3) Fail. c ≤ 1.3 from converged points at a candidate bond, or the leading correlation length of one sector grows by less than 1.15 over the last χ doubling while the other keeps growing.
- (I4) Consistency. The bulk energy density agrees with (E(256) − E(128))/128 of the C190 ground branch to 1e-5, and the density with the C190 mid-chain value to 0.006. A miss is reported, not fatal.
- Anything between I2 and I3 is undecided and reported with the numbers.

## Amendment, 2026-11-10: C196 outcome, and a registered prediction of the open-end form at L = 512

Written by the coordinator before any L = 512 output exists (C195 is still being set up).

C196 was not run. Neither infinite-chain method converges at the candidate field (REPORT §9.183, `codex/work_C196/`). I0 to I4 are not evaluated.

C193b (REPORT §9.182) found, after its candidates were opened, that the boson one-body function on the open qubit chain fits ln|h| = c − η ln tan u − 2 η_c ln cos u with u = π r/(2(L + 1)), where η_c is a boundary term that vanishes on the controls and on the effective model. That fit was not registered. This amendment turns it into a prediction for the third size.

Measurement. `codex/work_C193/ckpt_meas/meas193.py` unchanged, on the final p153c states at L = 512 in both parity sectors. Analysis with `mixed_fit` and `table()` of `diag193.py` unchanged, apart from adding L = 512 to the list of sizes, on the four ranges [0.05, 0.4], [0.05, 0.6], [0.1, 0.8], [0.05, 0.85] in r/L, and with the pair exponent from `ana193.py` on its frozen ranges. Reference values at L = 256 on the ground branch of p153c: η = 0.314 to 0.328, η_c = 0.135 to 0.168, pair exponent 0.638, corrected K_s = 1.52(8), change of K_s from L = 128 of −0.12(12).

- (O0) The pair exponent at L = 512 lies within 0.03 of 0.638.
- (O1) Pass, the form. On every range η(512) lies within 0.015 of its L = 256 value on the same range, η_c lies within 0.04 of a quarter of the pair exponent, and the rms residual is at most 0.02.
- (O2) Pass, the sign. The corrected K_s(512) is at least 1.3, its change from L = 256 is no more negative than −0.12, and the Kosterlitz–Thouless bound from the 256 → 512 change keeps the upper limit on |y_1| below the lower limit on y_∥.
- (O3) Fail. η(512) exceeds its L = 256 value by more than 0.03 on any range, or the corrected K_s(512) is below 1.15, or the bound from the 256 → 512 change allows |y_1| ≥ y_∥ while the fall of K_s from 256 to 512 exceeds the fall from 128 to 256 by more than the combined error.
- Anything between is undecided and reported with the numbers. If the L = 512 states miss T1, the clauses are evaluated and flagged as provisional.
