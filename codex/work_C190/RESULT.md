This run was done by a Claude Opus 5.5 agent (C190, 2026-11-09), not by Codex. Astra was not used (out of credits).

# C190 result: a converged full-qubit test at the candidate bonds

Derivations (the map to the effective model, density matching, the parity-gap bound, the solver) are in PROOF.md, the log in PROGRESS.md.

## Registration (written before any production output)

Written after the validation and the local pilots listed in PROGRESS.md (L ≤ 64, collaborator's field only), and before the Typhon job started. Nothing below was changed after production output existed.

### Runs

Bonds (PROOF 1): p153 (the collaborator's φ = 1.53 bond, unrounded), rat (their rational simplification), sym40 (φ = π/2, V/Δ = 40), sym80 (φ = π/2, V/Δ = 80, optional). Field on every site, open chains, both parity sectors, L = 64, 128, 256.

Fields.
- p153c: the collaborator's own field h_z = −19.97589125, the direct reproduction test (μ_eff = −2.143).
- Tuned keys p153_n316, rat_n316, sym40_n316, sym40_n18, sym80_n316. Three L = 128 even-sector runs (χ ≤ 192) at μ_eff = μ* − 0.03, μ*, μ* + 0.03, with μ* = −2.165 (n = 3/16) or −2.215 (n = 1/8) from C176's E(N ± 1). A least-squares line n_mid(μ_eff) through the three final middle-half densities gives the μ_eff that hits the target n_mid = 0.1957 (3/16) or 0.1326 (1/8), clipped to ±0.05 outside the grid, and h_z follows from PROOF eq. (2.3). Then this h_z is held fixed for all L and both parities (grand canonical).

Initializations. Random Z₂-symmetric MPS (rand:1) everywhere; in addition a product state at the target density (prod:N) at L = 64 and 128 for p153c, at L = 128 for p153_n316, and at L = 128 and 256 for sym40_n316. Where two initializations exist, the lower final energy is used.

Bond dimension schedules (pipeline190.SCHED; mixer amplitude α in the early stages):
- L = 64: χ = 32, 64, 128 (mixer), 192, 256, then final χ = 384.
- L = 128: χ = 32, 64, 128, 192 (mixer), 256, 384, then final χ = 512.
- L = 256: χ = 32, 64, 128, 192, 256 (mixer), 384, 512, then final χ = 768.
Two sweeps per stage. The final stage runs at least 3 sweeps (2 at L = 256) and at most 8 (6), and stops early when ΔE and ΔS̄₆ between consecutive final sweeps fall below (1e-9, 2e-4), (1e-9, 3e-4), (2e-9, 5e-4) at L = 64, 128, 256. Local solves target a residual of 1e-9, with at most 300 Davidson iterations.

Controls at the same L, the same χ schedules and the same code: the XX chain (c = 1) and two decoupled XX legs coupled only by the next-nearest-neighbour hop (c = 2), each at density 0.1957 (per leg for the two-leg control), with μ midway between open-chain levels so that the grand canonical ground state is unique. Exact free-fermion entanglement for the same chains is computed in the analysis, which gives the DMRG bias of the controls directly.

### Estimators (identical to C176)

- S̄_W(L): mean entanglement over cuts l with |l − L/2| ≤ W, for W = 1, 6, 12.
- c_two(L₁ → L₂; W) = 6 [S̄_W(L₂) − S̄_W(L₁)] / ln(L₂/L₁).
- c_prof(L): C176's analyze.c_prof (12-cut moving average, fit window [L/6, 5L/6]).
- n_mid: mean of ⟨(1 − Z_j)/2⟩ over j ∈ [L/4, 3L/4).
- Ground branch: at each L, the parity sector with the lower energy. All four parity pairings are also reported.
- Parity: Δ_GC(L) = E_odd − E_even at the same field, reported as L·Δ_GC/|t₂|.
- Valley weights W_a, W_b: from n(k) of the fermion one-body matrix on the central 48 sites (32 at L = 64), weights with cos k > 0 and < 0, as in C176.
- Density-corrected secondary estimator (added before any production output, after the local pilot showed that the two parity sectors differ in n_mid by about 1/L, 0.218 against 0.202 at L = 64): at each L, S̄_W is interpolated linearly in n_mid between the two parity sectors (which differ by one particle) to a common density n*, the mean n_mid of the ground-branch runs at L = 256 (or at the largest L available). This gives c_dc(L₁ → L₂; W). If the ground-branch c_two and c_dc differ by more than 0.3 at W = 6, the verdict is reported as estimator-dependent.
- Effective-model comparison: ΔS̄₆(L) = S̄₆^full(L) − S̄₆^eff(L, n_mid), with S̄₆^eff interpolated linearly in ln n_mid across C176's φ = π/2 families at the same L (sym18, sym532, sym16, sym316, sym14, table c176_ref.json).

### Convergence targets (per run)

- (C1) Local residual: every two-site solve of the last final-χ sweep reaches ‖Hθ − Eθ‖ ≤ 1e-8 (Δ = 1 units, 8e-7|t₂|). The maximum is reported.
- (C2) Stationarity: between the last two final-χ sweeps, |ΔE| ≤ 1e-8, |ΔS̄_W| ≤ 1e-3 for W = 1, 6, 12, and |Δn_mid| ≤ 5e-4.
- (C3) Truncation: maximum discarded weight of the last sweep ≤ 1e-7.
- (C4) χ convergence: |S̄₆(final χ) − S̄₆(previous stage)| ≤ 0.01. This is a bound of about 0.09 on c_two(128 → 256) from χ alone.
- (C5) Initializations: where two exist and their energies agree within 1e-7, their S̄₆ must agree within 0.01. Otherwise the higher-energy state is reported as metastable and not used.
- (C6) Density (tuned keys): n_mid at L = 256 within ±0.006 of the target.

A run failing C1 to C4 is kept and flagged, and any reading that depends on it is labelled unconverged.

### Verdict criteria

- (V1) The full chain tracks the effective model for a key if the ground-branch c_two(64 → 128) and c_two(128 → 256) are both ≥ 1.6 at W = 6 and W = 12, c_prof(256) ≥ 1.6, all runs used pass C1 to C4, and the two-leg control reads 2 ± 0.25 at the same size pairs.
- (V2) The full chain falls toward 1 if the ground-branch c_two(128 → 256) ≤ 1.3 at W = 6 and 12, or c_two(64 → 128) ≤ 1.3 at W = 6 and 12, with the same convergence conditions.
- Anything else is reported as unresolved (a crossover or an unconverged reading).
- (V3) The collaborator's downturn (c from 64 to 128 of 0.74 to 0.93, three central cuts) reproduces if the converged p153c ground-branch c_two(64 → 128, W = 1) ≤ 1.2. It does not reproduce if that value is ≥ 1.6.
- (V4) The full chain is faithful to the effective model in the size dependence if ΔS̄₆(256) − ΔS̄₆(64) lies within ±0.05 (equal slopes to about 0.2 in c over 64 → 256). A constant offset is allowed (C177 found about −0.2 at V/Δ = 8).
- (V5) Parity: if |L·Δ_GC|/|t₂| ≤ 1.5 at all three sizes with no growth, no pair gap is resolved (PROOF 4: an unpaired liquid is bounded by L·D_P/2 ≈ 1.2). If L·Δ_GC/|t₂| > 0 at all three sizes and grows by more than a factor 2 from 64 to 256, that is a pair-gap signature.
- (V6) Two valleys if W_min/W_max ≥ 0.8 at L = 128 and 256.
- Controls calibrate a size pair if XX reads 1 ± 0.15 and two-leg XX 2 ± 0.25 at W = 6. A size pair that fails calibration is not used for V1 or V2.

## Typhon jobs

Submitted 2026-11-09, both pending at submission (Priority). Partition typhon, one node each, 48 cores as 8 slots of 6 BLAS threads, 240 GB, 3-day limit. Directory /scratch/lustre/jkudlerflam/c18x/c190/, scripts typhon/c190a.slurm and typhon/c190b.slurm (copies in this folder). Each job first reruns validate.py (logs/validate_A.log, logs/validate_B.log), then pipeline190.py, then analyze190.py.

- Job 17606939 (c190a, node A, output outA/): p153c (collaborator's field) at L = 64, 128, 256, both parities, rand:1 plus prod at L ≤ 128; tuning and production of p153_n316.
- Job 17606940 (c190b, node B, output outB/): tuning and production of sym40_n316, sym40_n18, sym80_n316, rat_n316, and the XX and two-leg XX controls.

Overlap steps (coordinator). The queue was memory-bound (start estimate 2026-11-10 06:29 or later), so the coordinator launched the same two bodies as steps inside two existing allocations with `srun --jobid=<donor> --overlap` (typhon/overlap_run.sh): node A (8 slots of 6 threads) in 17604455 on node18, node B (6 slots of 6 threads) in 17602270 on node26. They write the same outA/ and outB/, with logs logs/validate_A_ov.log, logs/pipeline_A_ov.log, logs/overlap_A.log and the same for B. Each step cancels itself when its queued job (17606939 or 17606940) starts, and the queued job then resumes every run from its checkpoint. A donor can also end first (17602270 lives at most about 22 h), so killed-and-resumed runs are expected. analyze190.py uses only runs with status done unless --partial is given. Do not cancel, resubmit or edit 17606939 and 17606940, and do not start a second pipeline on outA or outB.

The xx2 row of validate.py (coordinator's question). On both nodes validate.py 10 prints a worst-case xx2 row with E = 3.3e-1, S = 8.3e-1. That row is the maximum over both parity sectors and both initial states, and the mismatch is entirely in the odd sector: xx2 P=0 agrees with ED to 3e-15 in energy and 1e-15 in every observable from both the random and the product state, while xx2 P=1 is off. This is expected and does not affect the control. At L = 10 the two legs have 5 sites each, so the odd sector's ground state is exactly degenerate (the extra particle on either leg). The random start finds the right energy (1.7e-15) but a different state in the degenerate pair, hence S = 0.17. The two leg numbers are separately conserved by this Hamiltonian, so the product start (one particle on one leg) cannot reach the other leg-number sector at all, hence E = 0.33. The production two-leg control runs only in the even sector, from random states. Equal leg lengths and equal μ give N_A = N_B, so the grand canonical ground state is even and nondegenerate there (pipeline190.free_control_mu: N = 12, 26, 50 at L = 64, 128, 256), and analyze190.py compares every control run with the exact free-fermion entropies at all cuts. Any failure of this kind would show there as a large DMRG − exact difference.

Outputs per run: out*/<tag>.json (configuration, per-sweep history with energies, residuals, discarded weights, S̄ and density profiles, and at the end the exact entropy profile, density, bond observables, one-body matrices and C_P), out*/tune/<key>.json (the tuned field), out*/logs/<tag>.log, out*/summary.txt, logs/pipeline_*.log. Checkpoints out*/ckpt/*.npz are large (up to about 5 GB each at L = 256) and should not be copied. A run that is cut by the time limit resumes from its checkpoint if the same script is resubmitted.

To collect, from inside codex/work_C190 (no brace expansion, no checkpoints). The remote globs must be quoted, since zsh otherwise stops with "no matches found". These commands were tested on the partial output on 2026-11-09 (typhon_out/ holds that snapshot, 0.8 MB).

    mkdir -p typhon_out/outA typhon_out/outB typhon_out/logs
    scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190/outA/*.json' typhon_out/outA/
    scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190/outB/*.json' typhon_out/outB/
    scp -r typhon:/scratch/lustre/jkudlerflam/c18x/c190/outA/tune typhon_out/outA/
    scp -r typhon:/scratch/lustre/jkudlerflam/c18x/c190/outB/tune typhon_out/outB/
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/outA/summary.txt typhon_out/outA/
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/outB/summary.txt typhon_out/outB/
    scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190/logs/*.log' typhon_out/logs/

The two summary.txt files exist only after a pipeline has finished. Then

    python3 analyze190.py typhon_out/outA > typhon_out/summaryA.txt
    python3 analyze190.py typhon_out/outB > typhon_out/summaryB.txt

(add --partial to include runs still in progress at their last completed sweep). Entropies are recorded only on sweeps with the mixer off, so with --partial a run appears only once it has reached χ = 192 (L = 64), χ = 256 (L = 128) or χ = 384 (L = 256). Before that only the "field tuning" block is printed. The summary prints, per run, the convergence flags C1 to C4 and the χ-stage history, and per key the parity table, c_two for all size pairs (ground branch and all parity pairings), c_prof, valley weights, ΔS̄₆ against C176 at equal n_mid, C5, C6 and V4. For the controls it adds the exact free-fermion entropies and the DMRG bias.

## Results (status at the end of the agent session, 2026-11-09)

One production run is done (the XX control at L = 256, R4). No candidate-bond run has reached its registered final χ, so the registered criteria V1 to V6 are not evaluated. R4 gives an interim reading at χ = 256 for the collaborator's bond and field at L = 64 and 128.

### R0. Validation and solver [numerical; PROOF 5, PROGRESS 3 to 6]

- validate.py at L = 10 passes on both Typhon nodes (logs/validate_A_ov.log, logs/validate_B_ov.log) as it does locally. Energies agree with ED to 7e-12, entropies to 6e-10 and the one-body matrices to 3.5e-10 for p153, rat and sym40 in both sectors. The xx2 row is the expected odd-sector degeneracy explained above.
- The (s₁, s₂)-block Kronecker preconditioner (PROOF eq. 5.1) is what makes the full chain tractable. On Typhon every local solve of every sweep so far meets its residual target (1e-8 or 1e-7 in the growth stages) with 70 to 90 matvecs on average. The diagonal preconditioner (C177, C183, the collaborator) stalls at 1e-7 to 1e-6.

### R1. Local pilot at the collaborator's field, L = 64 [numerical, validation only, PROGRESS 7]

At χ = 128, 192, 256 the even sector gives S̄₁ = 1.5993, 1.6125, 1.6142, with E = −650.0414814402, −650.0414849170, −650.0414852592 and n_mid = 0.2183. The collaborator's χ = 128 value is 1.60106. So even at L = 64 the entropy moves by 0.015 between χ = 128 and 256 while the discarded weight at χ = 128 is only 9e-9. The even sector is lower, with L·Δ_GC/|t₂| = +0.25. Against the effective model at equal n_mid, S̄₆ is lower by 0.23, close to C177's offset at V/Δ = 8.

### R2. Field tuning [numerical, out*/tune/]

- sym40_n316: at L = 128 (χ = 192), μ_eff = −2.195, −2.165, −2.135 give n_mid = 0.1561, 0.2038, 0.2375. The fitted slope is dn_mid/dμ_eff = 1.357, against about 1.25 for the effective model. The target n_mid = 0.1957 needs μ_eff = −2.1676, that is h_z = −19.976047251. At equal μ_eff the full chain is denser than the effective model by about 0.009, a chemical-potential shift of about 0.007|t₂|, of the size expected from the next even order.
- All five keys are tuned (out*/tune/<key>.json, copies in typhon_out/). Each line gives the three L = 128, P = 0 points as (μ_eff: n_mid, N_tot), then the fit.

  | key | points | slope dn_mid/dμ_eff | tuned μ_eff | tuned h_z |
  | :-- | :-- | :-- | :-- | :-- |
  | p153_n316 | −2.195: 0.1553, 22.06; −2.165: 0.1884, 26.06; −2.135: 0.2210, 30.05 | 1.096 | −2.1582 | −19.975988793 |
  | rat_n316 | −2.195: 0.1554, 22.07; −2.165: 0.1884, 26.07; −2.135: 0.2211, 30.05 | 1.094 | −2.1582 | −19.975988986 |
  | sym40_n316 | −2.195: 0.1561, 22.10; −2.165: 0.2038, 27.92; −2.135: 0.2375, 32.04 | 1.357 | −2.1676 | −19.976047251 |
  | sym40_n18 | −2.245: 0.0800, 13.09; −2.215: 0.1380, 19.91; −2.185: 0.1724, 24.07 | 1.539 | −2.2134 | −19.976333694 |
  | sym80_n316 | −2.195: 0.1554, 22.03; −2.165: 0.2045, 28.01; −2.135: 0.2372, 32.01 | 1.363 | −2.1674 | −39.988023246 |

- Limit of the tuning. At fixed L and parity the particle number is close to an integer of that parity (26.06, 30.05, and so on), so n_mid(μ_eff) is a staircase with steps of two particles, about 0.016 in n_mid at L = 128. A linear fit through three such points fixes the field to about half a step, not better. The steps halve at L = 256, where C6 is evaluated, and the density-corrected estimator covers what is left. Some points are not near an integer (13.09 and 19.91 in the even sector), which means the state mixes N and N ± 2 there. Only parity is exact in the full chain.
- The rat and p153 bonds give the same tuning points to 1e-4 in n_mid, as expected for a rational simplification of the same bond.

### R3. Production snapshot, 2026-11-09 21:02 EST [numerical, typhon_out/, logs on scratch]

- Both overlap pipelines started at 20:23. Node A runs the three p153_n316 tuning runs and five p153c runs (L = 128: P0 rand1, P0 prod28, P1 rand1, and L = 256: P0 rand1, P1 rand1). The L = 64 runs start when the tuning slots free up. Node B has finished the sym40_n316 tuning and is running the rat_n316 and sym40_n18 tuning. No production run is done.
- Every local solve so far meets its residual target, with 70 to 120 matvecs on average (largest single solve 232, cap 300).
- Measured sweep times with 6 threads per run. Complex bonds (p153), L = 128: 320 to 680 s at χ = 128 and 720 to 1120 s at χ = 192. Complex, L = 256: 570 to 860 s at χ = 128. Real (sym40 in the staggered gauge), L = 128: 225 to 280 s at χ = 192.
- The donor nodes have 96 cores each and the donor jobs use few of their 48 (sstat gives the node A step about 46 cores of CPU time). So the overlap steps already run at close to dedicated speed, and the queued jobs (start estimate 2026-11-10 05:48 when this was written, see R4 for the later ones) should not be expected to be much faster.
- Wall-time extrapolation (an estimate, from the measured 95 GFLOP/s matvec at D = 512 in PROGRESS 6 and the times above, uncertain by a factor of about 2): a complex L = 128 run needs about 1 h per sweep at χ = 512 and 4 to 13 h in total, and a complex L = 256 run needs 3 to 6 h per sweep at χ = 768 and 10 to 40 h in total. Real runs are about 4 times cheaper. (R4 revises this upward with later sweep times. The χ = 768 figure here is too low.) Checkpoints are written once per sweep, so each hand-over from an overlap step to its queued job (and each time-limit kill) costs every run its sweep in progress.

### R4. Interim reading at χ = 256 at the collaborator's bond and field, snapshot 2026-11-09 22:00 EST [numerical, interim; typhon_out/, local/out/]

Status. One production run is done (the XX control at L = 256). No p153c run has reached its registered final χ (384 at L = 64, 512 at L = 128, 768 at L = 256), so V1 to V6 are not evaluated, and everything below is a reading at χ = 256. It comes from `python3 analyze190.py typhon_out/outA --partial` and from the per-stage records in the json files, read with interim190.py (which also pairs the local L = 64 pilot with the Typhon L = 128 runs at equal χ, `python3 interim190.py local/out/p153c_L64_P0.json local/out/p153c_L64_P1.json typhon_out/outA/p153c_L128_P0_rand1.json typhon_out/outA/p153c_L128_P0_prod28.json typhon_out/outA/p153c_L128_P1_rand1.json`).

XX control. xx1_L256_P0_rand1 is done at χ = 768 and meets C1 to C4 (residual 1e-9, discarded weight 1.5e-13, ΔE = 1.7e-13 and ΔS̄ = 1.7e-11 between the last two sweeps). It gives S̄₁, S̄₆, S̄₁₂ = 1.1186, 1.1206, 1.1207, each 3.6e-8 below the exact free-fermion value, and c_prof = 0.993.

p153c (φ = 1.53, h_z = −19.97589125), last mixer-off sweep at χ = 256. Every local solve met its residual target (1e-8 on Typhon, 1e-9 in the local pilot).

| L | sector, init | E (Δ = 1) | S̄₁ | S̄₆ | S̄₁₂ | c_prof | n_mid | N_tot | discarded |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 64 | even, rand:1 (local pilot) | −650.0414852455 | 1.6140 | 1.6085 | 1.5929 | 2.57 | 0.2183 | 16.014 | 3.3e-10 |
| 64 | odd, rand:1 (local pilot) | −650.0414362271 | 1.5937 | 1.5881 | 1.5723 | 2.50 | 0.2018 | 15.025 | 2.8e-10 |
| 128 | even, rand:1 | −1290.1134121766 | 1.8476 | 1.8386 | 1.8314 | 2.09 | 0.2211 | 30.045 | 3.5e-9 |
| 128 | even, prod28 | −1290.1134121773 | 1.8476 | 1.8386 | 1.8314 | 2.09 | 0.2211 | 30.045 | 3.5e-9 |
| 128 | odd, rand:1 | −1290.1134075210 | 1.8403 | 1.8316 | 1.8247 | 2.10 | 0.2130 | 29.052 | 3.4e-9 |

- Reproducibility. The Typhon L = 64 even run at χ = 192 agrees with the local pilot at χ = 192 (E = −650.0414849190 against −650.0414849170, S̄₁ = 1.6125 in both).
- Initializations (C5, read at χ = 256). The random start and the product start with 28 particles agree to 7e-10 in energy and to four digits in S̄ and c_prof. Both end at N_tot = 30.05.
- Ground branch and parity. The even sector is lower at both sizes, with E_odd − E_even = 4.90e-5 at L = 64 and 4.66e-6 at L = 128, that is L·Δ_GC/|t₂| = +0.251 and +0.048. Both are below the unpaired-liquid bound of about 1.2 and the value decreases with L, so there is no pair-gap signature in this size pair. V5 needs L = 256.

c_two(64 → 128), with χ = 256 at both sizes.

| pairing (64 → 128) | W = 1 | W = 6 | W = 12 | Δn_mid |
| :-- | :-- | :-- | :-- | :-- |
| even → even (ground branch) | 2.022 | 1.991 | 2.064 | +0.0028 |
| even → odd | 1.959 | 1.931 | 2.006 | −0.0053 |
| odd → even | 2.198 | 2.168 | 2.243 | +0.0193 |
| odd → odd | 2.135 | 2.108 | 2.184 | +0.0112 |

- The density-corrected estimator at n* = 0.2211 gives c_dc = 1.99, 1.96, 2.03 (W = 1, 6, 12), within 0.03 of the ground-branch c_two.
- The collaborator's numbers for the same quantity (W = 1) are 0.79 on their ground branch (even at L = 64, odd at L = 128 by their energies) and 0.74 to 0.93 over the four pairings. At χ = 256 the downturn does not reproduce. The value 2.02 is on the "does not reproduce" side of V3 (≥ 1.6), but V3 is stated for the registered final χ and stays unevaluated until those stages finish.

Where the two computations differ. The energies below are E − E_off in my convention (the collaborator's shifted energy plus 10, checked at both sizes).

| L, sector | collaborator, χ = 128 | this work, χ = 256 | collaborator − this work |
| :-- | :-- | :-- | :-- |
| 64, even | −11.584442504 | −11.584445246 | 2.74e-6 |
| 64, odd | −11.584393385 | −11.584396227 | 2.84e-6 |
| 128, even | −13.199287936 | −13.199332177 | 4.42e-5 |
| 128, odd | −13.199296931 | −13.199327521 | 3.06e-5 |

- At L = 128 the collaborator's states lie 3e-5 to 4e-5 above mine. That is 7 to 9 times the parity splitting, and their branch order is reversed (odd lower by 9.0e-6, against even lower by 4.7e-6 here).
- My own mixer-on sweep at χ = 128 has E − E_off = −13.1992881, next to their even-sector value. So their L = 128 states appear to be what χ = 128 gives, not the product of a failed solver. Their entropies there are S̄₁ = 1.686 (even) and 1.692 (odd), below my χ = 256 values by 0.16 and 0.15.
- χ dependence of S̄₁ for the even N = 30 state at L = 128. It is 1.686 at χ = 128 (collaborator), 1.817 at χ = 192 (the tuning run tune_p153_n316_2, at a field higher by 4.75e-5, with the same N_tot = 30.05, discarded weight 9e-9, stationary to 2e-4 over its last three sweeps) and 1.848 at χ = 256. The steps are 0.131 and 0.031. At L = 64 the same sequence is 1.599, 1.613, 1.614, with steps 0.013 and 0.0015. So a χ = 128 truncation costs about ten times more entropy at L = 128 than at L = 64. Relative to χ = 256 that shifts c_two(64 → 128) down by 6 × 0.15/ln 2 ≈ 1.3, which is the size of the disagreement. The χ = 128 point is the collaborator's and the χ = 192 point is from a neighbouring field, so this sequence is an indirect comparison. The χ = 384 and 512 stages of the same run give the direct one.
- What is still open on my side. The step from χ = 192 to 256 at L = 128 is 0.031, so S̄(128) at χ = 256 is not converged at the level of C4 (0.01). The steps shrink by a factor of about 4 per stage (9 at L = 64), so I expect about 0.01 more from the χ = 384 and 512 stages. An increase of 0.01 in S̄(128) raises c_two(64 → 128) by 0.09, away from the collaborator's value.

Against the effective model (C176, interpolated in ln n_mid at each L). ΔS̄₆ = −0.232 at L = 64 (n_mid = 0.2183) and −0.246 at L = 128 (n_mid = 0.2211). The effective model's own c_two between these two points is 2.11 at W = 6, against 1.99 for the full chain at χ = 256. The difference of −0.014 in ΔS̄₆ is inside V4's ±0.05 window, and the remaining χ deficit at L = 128 has the same sign and about the same size. The L = 128 reference point is interpolated between N = 24 and N = 32 (n_mid 0.204 and 0.268), which adds an error of about 0.005 to 0.01 in S̄. So over 64 → 128 the full chain at this bond follows the effective model's c ≈ 2 to within about 0.1 to 0.2, with a constant offset of about −0.24 in S̄₆ (C177 found about −0.2 at V/Δ = 8). The decisive pair for V1 and V2 is 128 → 256, which is not available yet.

Tuning runs as χ = 192 data at L = 128 (not used for any verdict). The fifteen tuning runs (even sector, final χ = 192, discarded weight 6e-10 to 1e-8) give S̄₆ between 1.45 and 1.82 and c_prof between 1.6 and 2.1. Against the effective model at equal n_mid, ΔS̄₆ runs from −0.20 at n_mid ≈ 0.155 to −0.28 (p153, rat), −0.29 (sym40) and −0.32 (sym80) at n_mid ≈ 0.22 to 0.24. Part of that growth with density is the χ = 192 deficit, since the p153 point with N = 30 moves from −0.276 at χ = 192 to −0.246 at χ = 256. The rat and p153 bonds agree to 3e-3 in S̄ at every point. One tuning point (sym40_n316 at μ_eff = −2.165, N_tot = 27.92) was not stationary, with S̄₆ still moving by 0.015 in its last sweep. Its n_mid = 0.204 agrees with the V/Δ = 80 point at the same μ_eff (0.2045), so the tuned field is not affected beyond the half-step limit of R2.

L = 256 so far. The two p153c runs are in the χ = 192 mixer-on stage, with N_tot = 57.16 (even) and 57.15 (odd), n_mid = 0.2155 and 0.2154, and E − E_off = −16.42889214 and −16.42890298. The even sector has no integer close to 57.16 with the right parity, so at this field it sits between N = 56 and N = 58 or has not settled, and the odd sector is lower by 1.1e-5 at this stage. The ground branch at L = 256 may therefore be odd. This is unconverged and is recorded only so that the final parity table is read with it in mind.

Timing, revised (this replaces the estimate in R3).

- Measured sweep times for complex bonds with 6 threads: 1170 to 1650 s at L = 128, χ = 256, and 1790 to 2690 s at L = 256, χ = 192. The observed scaling lies between χ² and χ³.
- A complex L = 128 run then needs 0.9 to 1.3 h per sweep at χ = 384 and 1.6 to 3.1 h at χ = 512, so 7 to 27 h more for its remaining 2 + (3 to 8) sweeps.
- A complex L = 256 run needs 7 to 40 h per sweep at χ = 768, and 1.5 to 6 days for its whole schedule. R3's 3 to 6 h per sweep was too low. The registered χ = 768 at L = 256 for the complex bonds (p153c, p153_n316, rat_n316) may not be reached within one 3-day job. The real bonds (sym40, sym80) and the controls are 3 to 4 times cheaper.

Allocation lifetimes at 22:05 EST (squeue).

- Donor 17602270 (node B, node26) has 20 h 45 min left, so the node B step ends by 2026-11-10 18:50. Donor 17604455 (node A, node18) has 40.5 h left, so the node A step ends by 2026-11-11 14:37.
- The start estimate of the queued jobs 17606939 and 17606940 moved during the evening from 2026-11-10 05:48 to 2026-11-10 16:09 and then to 2026-11-14 14:52 (reason Priority for both). If the last estimate holds, node B stops for almost four days when its donor ends and node A for three. I did not touch the jobs or start anything else. Whether to move the payload to another allocation is the coordinator's decision.

### Coordinator note (2026-11-09, after the agent's interim stop)

- R1's χ = 256 pilot values (S̄₁ = 1.6142, E = −650.0414852592) are from the first sweep at that χ. The last sweep in local/out/p153c_L64_P0.json gives S̄₁ = 1.6140 and E = −650.0414852455. The lead is unchanged.
- The interim state is written up in REPORT §9.164. No production run had finished at that point.
