This run was done by a Claude Opus 5.5 agent (C196b, 2026-11-10), not by Codex. Astra was not used (out of credits).

# C196b result: second attempt at entanglement scaling on the infinite qubit chain

Spec: `specs/C196b_infinite_chain_second_attempt.md`. Method and derivations in PROOF.md, the log in PROGRESS.md.

## Registration

### I0 to I4, copied verbatim from `specs/C195_C196_third_size_and_infinite_chain.md` (section C196) before any production output exists

- (I0) Controls with the same code, χ sequence and fit window: the XX chain reads c within 0.1 of 1, and the two-leg XX control reads c within 0.15 of 2. A gapped control (any bond with a known finite ξ) shows ξ saturating.
- (I1) Convergence, per point: energy density change over the last iterations below 1e-10 per site and gradient or residual norm below 1e-7. Points that miss this are flagged and excluded from the fit.
- (I2) Pass. At each candidate bond the slope of S against ln ξ over the last three converged χ gives c ≥ 1.7, ξ grows by at least a factor 1.5 over the last χ doubling in both parity sectors, and I0 holds.
- (I3) Fail. c ≤ 1.3 from converged points at a candidate bond, or the leading correlation length of one sector grows by less than 1.15 over the last χ doubling while the other keeps growing.
- (I4) Consistency. The bulk energy density agrees with (E(256) − E(128))/128 of the C190 ground branch to 1e-5, and the density with the C190 mid-chain value to 0.006. A miss is reported, not fatal.
- Anything between I2 and I3 is undecided and reported with the numbers.

### Additions of the C196b spec (verbatim in substance), and the details fixed here before any production output

- χ sequence 64, 96, 128, 192, 256, 384, 512, and 768 if feasible, each started from the previous one.
- The sector clause of I2 and I3 is meaningful only for a state with a fixed number of particles per cell. If the converged state mixes particle numbers at leading order, this is said, c is reported from S against ln ξ, and the sector clause is marked not evaluable.
- I4's density clause is weak here because the cell pins the density. It is reported with the number per cell and its variance.
- If no route reaches I1 at a candidate bond, stop and report with numbers. I1 is not relaxed.

Estimator details (fixed here, used identically for the controls).
- S̄. Bond entropy averaged over the q inequivalent bonds of one cell (bonds q, …, 2q − 1 of the window, from the singular values of the final left pass of an iteration), then averaged over the last two iterations of the χ point (C189's convention for a state that may alternate between iterations).
- ξ. From the two-cell transfer map X ↦ M†XM with M = A[0..q−1] Λ_t B[q..2q−1] Λ_{t−1}^{-1} of the last iteration, resolved by r = q_bra ⊕ q_ket: ξ_even from the second eigenvalue of r = 0, ξ_odd from the leading eigenvalue of r = 1 (the sector of a single-spin operator), ξ_r = −2q/ln(|λ_r|/λ₀) in sites. ξ = max(ξ_even, ξ_odd) is the abscissa of the fit; the fits with ξ_even and ξ_odd alone are also reported.
- c. 6 × least-squares slope of S̄ against ln ξ over the last three converged χ (I2). Local slopes between consecutive χ are reported.
- I1 as implemented in run196b.py. Energy per site e_t = (window eigenvalue)/(2q) + e_off, the energy of the 2q sites added at iteration t. Converged when, after at least min_it iterations, |e_t − e_{t−2}| < 1e-10 (twice the change of the two-iteration mean of e, so a state that alternates between iterations, C189's case, is covered), |S̄_t − S̄_{t−2}| < 1e-6, and every two-site Davidson solve of the last two iterations has residual below 1e-7 (its target is 1e-9). The energy of a point is the mean of its last two iterations. (Amended 2026-11-10 before any production output: the first version also required |e_t − e_{t−1}| < 1e-10, which a converged period-two state of the two-leg control cannot meet, its e alternating by 1e-6 with |e_t − e_{t−2}| = 4e-12.) The residual of I1 is that local residual, since in iDMRG the local eigenproblem is the whole update. The residual of the McCulloch prediction and its fidelity with the optimized window are logged as diagnostics, not used as criteria (on the XX control they stay at 3e-3 and 1 − 4e-7 at a point that meets I1, because the prediction is built from truncated cells).
- Number per cell. ⟨N_cell⟩ and Var N_cell of the left cell of the last iteration (closed with Λ_t²). "Stable" means |⟨N_cell⟩ − p| < 0.02 for the mean of the last two iterations and a change of that mean below 1e-3 over the last ten iterations.
- Reference for I4 and G3 at sym40_n316. The spec quotes (E(256) − E(128))/128 = −10.0011861. From the C190 GPU files (ground branch P = 1 at both sizes, E_MPO = −13.1890724885 and −16.4089104718, e_off = −9.976047251 per site) I get −10.0012022. Both are reported. I could not reproduce −10.0011861 from any pair of C190 runs at this key.

## Gates

### G1, the q-site-cell code against ED: passed [numerical, local]

`python3 -E -s validate196b.py ed` (output in local/G1_ed.out). Untruncated growth from a random open window, q = 2, 3, 4, chains up to L = 16, with three extra complete window sweeps per iteration (needed only in the first growth steps, PROOF 2). The accumulated window eigenvalues against sparse-Lanczos ED of the same open chain (bulk MPO, P = 0): worst deviation 2.7e-11 in the first run and 1.3e-11 in the saved rerun (the random start differs), over sym40_n316 and p153c (all deviations of order 1e-11, set by the Lanczos and Davidson tolerances) and the TFIM at h = 1.05 (below 2e-13). Tolerance 1e-10. Pass. G1 was not run on a separate route-2 code path beyond the XX check below.

### G2, controls (I0): passed at χ ≤ 96 on the local machine, not run on the production χ sequence [numerical]

`python3 -E -s summarize196b.py local/out c1loc c2loc t_tfg` (saved in local/summary_controls.txt). Cell q = 8, p = 2 (n = 1/4), μ = −√2, N-definite start from fibdmrg (PROOF 1), production estimator.
- XX chain (c1loc), χ = 16 to 64, every point meets I1. c = 0.987 over χ = 32, 48, 64 (ξ = max; 0.989 with ξ_even, 0.987 with ξ_odd), local slopes 1.017, 0.996, 0.961, 1.038. ⟨N_cell⟩ = 2.0000 to 2.0003, Var N_cell = 0.398 to 0.404 (0.47 for a 16-site cell, the logarithmic growth of a U(1)-symmetric liquid). Within 0.1 of 1. Pass.
- Two-leg XX (c2loc), χ = 16 to 96, mixer in the first 40 iterations, every point meets I1. c = 2.148 over χ = 48, 64, 96 with ξ = max (2.019 with ξ_even), local slopes 1.84, 1.85, 2.39 (ξ = max). Within 0.15 of 2. Pass, with the χ = 96 local slope high.
- Gapped control, TFIM h = 1.2, q = 8: I1 in 10 iterations at every χ, ξ_odd = 5.238 and ξ_even = 2.436 from χ = 24 to 64 (unchanged to 1e-11), against the exact asymptotic 1/ln 1.2 = 5.48 and 2.74. Saturated. Pass.
- The two-leg control on Typhon (job 17612141, c_xx2) did not meet I1 at χ = 32 within 400 iterations (|e_t − e_{t−2}| = 2.0e-10). It is still running. The controls were not run on χ = 64 to 768 as the registration requires for the final I0, because no candidate passed G3 (below).

### G3, sym40_n316 at χ = 64 and 128: FAILED at χ = 32 and χ = 64 [numerical]

Route 1 (q = 32, p = 6, effective-model start, registered field h_z = −19.976047251). No point reaches I1.
- χ = 32 (local, 50 and 31 iterations in two runs): |e_t − e_{t−2}| scatters between 1e-9 and 1e-7 with no downward trend, |S̄_t − S̄_{t−2}| about 1e-3, prediction fidelity 1 − 3e-5 to 1 − 5e-5. e = −10.0011944.
- χ = 64 (local run G3loc continued from the χ = 32 state for 84 iterations, and Typhon run g3_s316 from the effective start, 18 iterations at collection, the same numbers to 1e-9): |e_t − e_{t−2}| between 5e-10 and 3e-8, |S̄_t − S̄_{t−2}| between 2.5e-5 and 1.6e-3, fidelity 1 − 5e-6 to 1 − 1.2e-5, every local residual below 1e-9. e = −10.0011998 (two-iteration mean), S̄ = 1.389.
- The number per cell is not pinned. At χ = 64 the two-iteration mean of ⟨N_cell⟩ is 6.124 (stable to 1e-3 over 40 iterations, but 0.12 above p, so the registered "stable" rule fails), which is the C190 open-chain bulk density 0.1917 × 32 = 6.13. ⟨N_cell⟩ alternates by ±0.065 between iterations. At χ = 32 a field scan Δh_z = 0, −2e-5, −6e-5, −1.5e-4 moves ⟨N_cell⟩ continuously through 6.00, 5.89, 5.63, 5.05, so the commensurate cell does not lock the density, and at the field where ⟨N_cell⟩ = 6.00 the scatter is the same. The allowed field retune therefore cannot help, and I did not record a retuned field.
- Energy against I4/G3's reference. e(χ = 64) = −10.0011998 is 2.4e-6 above my reading of the C190 ground branch (−10.0012022) and 1.37e-5 below the spec's value (−10.0011861). Not evaluated as a pass because I1 fails.
- sym40_n18 (q = 16, p = 2, h_z = −19.976333694) behaves the same: χ = 32 scatter 5e-8 in e, χ = 64 (Typhon g3_s18, 19 iterations) |e_t − e_{t−2}| 8e-9 to 2e-7, ⟨N_cell⟩ = 2.12 (C190 bulk 0.1324 × 16 = 2.12).
- A 64-site cell (q = 64, p = 12) at χ = 32 removes the alternation (⟨N_cell⟩ = 11.98 to 12.00) but not the scatter (3e-8). One extra full window sweep per iteration was prepared as a variant but not run (the first Typhon job died on a code mismatch, see below).

Why it fails [numerical, with a heuristic reading]. The iDMRG state at the candidate breaks the approximate U(1) of the real particles. The anomalous correlator |⟨b_i b_{i+r}⟩| (b = |0⟩⟨1|, local/u1diag.py) stays as large as the normal one out to r = 56: 0.080 against 0.080 at χ = 32 and 0.081 against 0.084 at r = 48, χ = 64 (local/u1diag_s316_chi64.out). On the XX control with the same code it is 1e-6 against 0.09 (local/u1diag_xx_chi64.out). The number variance of a cell grows linearly with the cell, 1.60 at 32 sites and 2.77 at 64 sites (XX: 0.40 at 8 and 0.47 at 16, logarithmic). So the finite-χ state is a pair condensate (a cat of two U(1)-broken components under the exact parity), its phase is pinned only by the weak third-order pair term, and the iteration wanders in that soft direction. This is the case the spec anticipated ("a cat of U(1)-broken components"), and it sets in at χ = 32 from a number-definite start within the first iterations (Var N_cell = 1.08 at the first iteration, 1.6 after five).

### Route 2 (VUMPS with a multi-site cell): FAILED at χ = 32

vumps196b.py, started from the 2q-site closing cell of the iDMRG state (PROOF 4). On the XX control (16-site cell, χ = 32) it converges (ε from 6e-5 to 4e-7 in 15 steps, energy 1e-7 below the iDMRG value). At sym40_n18 (32-site cell, χ = 32) the first A_C solves end at residual 5e-6 after 300 Davidson iterations, and within four steps the state collapses to the dressed vacuum (N per 32 sites 4.04 → 0.06), after which ε stays at 4e-3 to 9e-3 for 40 steps. This is the first attempt's VUMPS failure, now from a liquid start.

### Route 3 (start cut from a C190 open chain): not run

The spec makes route 3 conditional on the embedding being the problem. It is not: the embedded start has N = 12 exactly, the energy reaches −10.0011944 (χ = 32) within ten iterations, and the failure is the U(1) breaking that sets in during the iteration whatever the start (the χ = 64 runs from the effective start and from the χ = 32 state reach the same state). A start from the C190 chain would meet the same iteration. I record this as a judgement, not a test.

### G4: not reached

No production was started. Timing for reference: one iteration of the 64-site window at χ = 64 takes 40 to 50 s on 4 host cores (Typhon and local alike).

## Outcome

Every route fails at the gates. I1 is not reached at either candidate bond at χ = 32 or 64, so I2, I3 and I4 are not evaluated and no production exists. G1 passes, G2 passes at χ ≤ 96 with the production estimator, G3 fails, G4 is not reached.

## Typhon

- Job 17611816 (c196b_g3, gpu-medium, apollo02) started 10:30 and ended after 96 s: the uploaded run196b.py passed a warm-up argument that the uploaded idmrg196b.py (an older copy) did not accept. My error, no output.
- Job 17612141 (c196b_g3, gpu-medium, 3 A100 on apollo02, 1-day limit), submitted and started 11:06 EST with the corrected code: g3_s316 (χ = 64, 96, 128, at most 300 iterations each), g3_s18 (same at sym40_n18), c_xx2 (two-leg control, χ = 32 to 256). It will not change the verdict unless the scatter falls below 1e-10 at χ = 96 or 128. It can be cancelled without loss.

Layout on /scratch/lustre/jkudlerflam/c19x/c196b/: code at the top level, cfg/<tag>.json, out_g3/<tag>.jsonl (one record per finished χ), out_g3/<tag>.iter.jsonl (one line per iteration), out_g3/ckpt/ (do not copy), logs/<tag>_<jobid>.log, logs/c196b_<jobid>.out. Slurm script typhon/c196b_run3.slurm. Collection (tested, run from codex/work_C196/b/):

    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c19x/c196b/out_g3/*.jsonl' typhon_out/out_g3/
    python3 -E -s local/show.py typhon_out/out_g3/g3_s316.iter.jsonl | tail
    python3 -E -s summarize196b.py typhon_out/out_g3 c_xx2

## Files

idmrg196b.py (q-site-cell iDMRG), run196b.py (driver), validate196b.py (G1), summarize196b.py, vumps196b.py (route 2), local/ (configs, outputs, show.py, u1diag.py, vtest.py, G1_ed.out, summary_controls.txt, the U(1) diagnostics), typhon/ (slurm script, configs). Copied unchanged: z2dmrg.py, z2dmrg_xp.py, models.py, idmrg196.py, vumps196.py, fibdmrg.py.

## What a further attempt would need

(i) A way to keep the finite-χ state U(1)-symmetric in the real particles when H conserves them only up to virtual pairs. No exact label exists (PROOF 5), so this would have to be a symmetric ansatz or an explicit symmetrization, not a label. (ii) Alternatively, accept the U(1)-broken state and read c from it as the spec allows, but that needs a convergence criterion other than I1, which is registered and not relaxed here.
