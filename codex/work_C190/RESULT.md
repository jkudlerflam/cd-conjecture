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

To collect, from inside codex/work_C190 (no brace expansion, no checkpoints):

    mkdir -p typhon_out/outA typhon_out/outB typhon_out/logs
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/outA/*.json typhon_out/outA/
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/outB/*.json typhon_out/outB/
    scp -r typhon:/scratch/lustre/jkudlerflam/c18x/c190/outA/tune typhon_out/outA/
    scp -r typhon:/scratch/lustre/jkudlerflam/c18x/c190/outB/tune typhon_out/outB/
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/outA/summary.txt typhon_out/outA/
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/outB/summary.txt typhon_out/outB/
    scp typhon:/scratch/lustre/jkudlerflam/c18x/c190/logs/*.log typhon_out/logs/

Then

    python3 analyze190.py typhon_out/outA > typhon_out/summaryA.txt
    python3 analyze190.py typhon_out/outB > typhon_out/summaryB.txt

(add --partial to include runs still in progress at their last completed sweep). The summary prints, per run, the convergence flags C1 to C4 and the χ-stage history, and per key the parity table, c_two for all size pairs (ground branch and all parity pairings), c_prof, valley weights, ΔS̄₆ against C176 at equal n_mid, C5, C6 and V4. For the controls it adds the exact free-fermion entropies and the DMRG bias.

## Results (status at the end of the agent session, 2026-11-09)

No production result exists yet. The registered criteria V1 to V6 are not evaluated.

### R0. Validation and solver [numerical; PROOF 5, PROGRESS 3 to 6]

- validate.py at L = 10 passes on both Typhon nodes (logs/validate_A_ov.log, logs/validate_B_ov.log) as it does locally. Energies agree with ED to 7e-12, entropies to 6e-10 and the one-body matrices to 3.5e-10 for p153, rat and sym40 in both sectors. The xx2 row is the expected odd-sector degeneracy explained above.
- The (s₁, s₂)-block Kronecker preconditioner (PROOF eq. 5.1) is what makes the full chain tractable. On Typhon every local solve of every sweep so far meets its residual target (1e-8 or 1e-7 in the growth stages) with 70 to 90 matvecs on average. The diagonal preconditioner (C177, C183, the collaborator) stalls at 1e-7 to 1e-6.

### R1. Local pilot at the collaborator's field, L = 64 [numerical, validation only, PROGRESS 7]

At χ = 128, 192, 256 the even sector gives S̄₁ = 1.5993, 1.6125, 1.6142, with E = −650.0414814402, −650.0414849170, −650.0414852592 and n_mid = 0.2183. The collaborator's χ = 128 value is 1.60106. So even at L = 64 the entropy moves by 0.015 between χ = 128 and 256 while the discarded weight at χ = 128 is only 9e-9. The even sector is lower, with L·Δ_GC/|t₂| = +0.25. Against the effective model at equal n_mid, S̄₆ is lower by 0.23, close to C177's offset at V/Δ = 8.

### R2. Field tuning [numerical, outB/tune/sym40_n316.json]

- sym40_n316: at L = 128 (χ = 192), μ_eff = −2.195, −2.165, −2.135 give n_mid = 0.1561, 0.2038, 0.2375. The fitted slope is dn_mid/dμ_eff = 1.357, against about 1.25 for the effective model. The target n_mid = 0.1957 needs μ_eff = −2.1676, that is h_z = −19.976047251. At equal μ_eff the full chain is denser than the effective model by about 0.009, a chemical-potential shift of about 0.007|t₂|, of the size expected from the next even order.
- The other four keys (p153_n316, rat_n316, sym40_n18, sym80_n316) were still tuning when the session ended. Their fields will be in outA/tune/ and outB/tune/.

### Coordinator note (2026-11-09, after the agent's interim stop)

- R1's χ = 256 pilot values (S̄₁ = 1.6142, E = −650.0414852592) are from the first sweep at that χ. The last sweep in local/out/p153c_L64_P0.json gives S̄₁ = 1.6140 and E = −650.0414852455. The lead is unchanged.
- The interim state is written up in REPORT §9.164. No production run had finished at that point.
