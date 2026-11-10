This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C175 result: the valley-pseudospin sign of the Fibonacci-cusp gas, analytically

Derivations are in PROOF.md, the log in PROGRESS.md. Scripts: scat.py (exact lattice threshold scattering), analytic_A.py and analytic_res.py (closed forms, sympy), scan_w.py, calib.py, phi_scan.py, smatrix.py, ring2.py (brute-force two-particle ring ED), nbody_box.py (exact few-body bound states), ring_k.py, ring_eos.py and ring_clusters.py (momentum-resolved ring ED). block_dmrg.py and dmrg_d.py are C173's, copied unmodified for the Typhon job. Outputs: scan_w.out, calib.out, phi_scan.out, smatrix.out, ring2.out, nbody3.out, nbody3b.out, nbody45.out, nbody4_w1_R110.out, ring_clusters.out, ring_eos_w1.jsonl, ring_eos_w1b.jsonl, local/bd_pilot*.jsonl.

## Pre-registered predictions and pass criteria (written before computing, unchanged)

P1. Dilute-limit criterion. In the impenetrable dilute gas of model (3.2) at φ = π/2 the valley pseudospin is an XXZ chain on the squeezed lattice with J_z = 2c (a_aa − A_11), J_⊥ = 2c |A_12|, c = π² n³/(3m), where a_aa is the like-valley threshold scattering length and A is the 2 × 2 threshold scattering-length matrix of the unlike pair (channels ab, ba). The dilute gas is two-component (c = 2) if and only if Δ ≡ (a_aa − A_11)/|A_12| lies in (−1, 1], and in addition no two-body bound state lies below 2ε_min at any total momentum. Δ ≤ −1 means one valley (immiscible, c = 1), Δ > 1 a pseudospin Néel state (valley alternation, c = 1).

P2. Consistency checks the two-body code must pass before its output is used: (i) A is Hermitian to 1e-8 (current conservation); (ii) for pure nearest-neighbour hard rods (no NNN hop, no V2, V3) the scattering length is a = 1; (iii) the threshold data reproduce two-particle ring levels at L ≥ 200 to relative accuracy O(1/L).

P3. Ring ED cross-check (N = 3 to 6). The effective theory (Tonks charge sector on the squeezed ring, twisted pseudospin XXZ ring with the computed A, plus the valley-odd cubic dispersion term) must reproduce the lowest 2^N-manifold of ED levels, with the error in pseudospin splittings falling with L at fixed N. Pass: relative error in the splittings below 20% at the largest L for N = 3 and 4, and decreasing with L.

No expectation is registered for the sign of Δ itself.

Report against them: P1's second condition fails everywhere at φ = π/2 (a bound pair exists), so by the registered criterion the dilute gas is not two-component. P2(i) failed as stated, for a kinematic reason that is understood (below), and the Hermitian part is used. P2(ii) passed. P2(iii) passed for the bound-state level. P3 could not be tested as registered, because its effective theory presupposes no bound state; the ring ED was used instead to test the pairing picture (PROOF 7).

## Verdict

At φ = π/2 the dilute two-valley gas of (3.2) is not a two-component liquid. The question "easy-plane or not" is pre-empted, because one particle from each valley binds. The bound pair sits at total momentum k_a + k_b = π, it exists for every 0 < |w| < 4, and in the dilute limit the gas is a liquid of valley-singlet pairs with a pseudospin gap of order E_b, so c = 1. Among scattering states the pseudospin is ferromagnetic anyway (Δ_σ < −1, the vector-chiral side), so neither dilute instability is easy-plane. The predicted c = 2 window at small n is empty, at φ = π/2 and around it. The candidate density n ≈ 0.17 lies outside the dilute regime (n a_+ ≈ 3.5 at |w| = 1). There the ring ED excludes the vector-chiral state and is consistent with pairing, but it cannot separate a small pseudospin gap from c = 2 at L ≤ 36. A Yang–Gaudin-type estimate puts the pair gap at about 1e-4 |t₂| (correlation length of order 1e4 sites), which would make C173's c ≈ 2 reading at L ≤ 96 a crossover, but that estimate is heuristic. This is a strong negative indication for the candidate, not a closure at n = 0.17.

## Key claims

1. [derived] Dilute mapping (PROOF 1). The impenetrable two-valley gas maps to the XXZ chain J_z = 2c(a_s − A_d), J_⊥ = 2c|A_o|, c = π²n³/(3m), so c = 2 iff a_− < a_s < a_+ with a_± = A_d ± |A_o|. This agrees with C180's criterion. The lattice threshold matrix A has a potential-independent anti-Hermitian part diag(∓iγ), γ = ε‴/(6ε″), from the cubic relative dispersion. It is valley-odd and drops out, and the Hermitian part is the physical one [derived, checked numerically].
2. [proved, algebra, matched to the numerics to all printed digits] Closed form at φ = π/2 (PROOF 2.1). With u = |w|, A_d = 5/2 + 8/u², |A_o| = √(u⁴ + 160u² + 256)/(2u²), a_∓ = [5u² + 16 ∓ √(u⁴ + 160u² + 256)]/(2u²). At |w| = 1, a_− = 0.289, a_+ = 20.71.
3. [proved, algebra] An exact zero-energy resonance (PROOF 2.2). The unlike pair has a scattering-length pole when D = 2V₂V₃ + 4V₂ + V₃u² + 4V₃ + 8 = 0. With V₃ = 0 this is V₂ = 2t₂, the Schrieffer-Wolff value, for every |w|. The generated V₃ = (2/3)t₂ detunes it to D = −2u²/3, on the bound side. Higher Schrieffer-Wolff orders move D by relative O((Δ/V)²) and cannot undo this [derived].
4. [numerical, two independent codes; derived for the whole interval] The bound pair (PROOF 2.3). E_b = 2ε_min − E₂ = 1.2e-4, 1.7e-3, 0.0195, 0.064, 0.1225, 0.195 at |w| = 0.25, 0.5, 1, 1.5, 2, 3 (scan_w.out). It is confirmed by ring2.py (brute-force ring ED, −0.01953 at L = 240, |w| = 1). At small |w|, E_b ≈ ε″/a_+² ≈ |w|⁴/32, the shallow bound state of the large a_+. Because D ≠ 0 for 0 < u < 4, it cannot leave through threshold, so it persists on the whole interval. The pair size is a_+ ≈ 21 sites at |w| = 1 and 8 at |w| = 2.
5. [numerical] Coordinator's calibration (PROOF 2.4, calib.out). V₂ = V₃ = 0 is two-component (Δ_σ = 0.25 at |w| = 1, 0.125 and 0.50 at |w| = 0.5 and 2), reproducing the Mishra–Greschner–Santos limit (arXiv:1411.7270, as cited by C180). V₂ alone sits on the exact resonance (Δ_σ = 1 boundary). V₃ alone stays two-component (Δ_σ = 0.38 at |w| = 1). Together they bind the pair and give a_s = −3.17, far below a_− = 0.29 (Δ_σ = −1.34). So C180's Born direction is right for the combination, but the bound state dominates.
6. [numerical] Sign of J_z among scattering states. Δ_σ = −1.07, −1.16, −1.34, −1.54, −1.80, −2.89 at |w| = 0.25, 0.5, 1, 1.5, 2, 3. So J_z < 0 and |J_z| > J_⊥ (ferromagnetic, vector chiral) for all |w| < 3.6. The like pair has a_s < 0 there, with a like-pair resonance near |w| ≈ 3.63.
7. [numerical] Clusters (PROOF 4). Atoms bind to dimers: the trimer lies 6.0e-4 (|w| = 1, converged) and 0.018 (|w| = 2, converged) below dimer + atom. At |w| = 2 dimers bind too: the tetramer lies at least 0.009 below two dimers, and the pentamer below every split. The |w| = 1 tetramer is not resolved by R = 110. [heuristic] At |w| = 2 the dilute gas is unstable to clustering, a pair liquid with a first-order onset. A local DMRG pilot at |w| = 2, μ = −2.9 gives n = 0.355 with the full V's, against 0.066 with V₂ = V₃ = 0 (local/bd_pilot.jsonl).
8. [derived at leading order in n, with numerical inputs] Dilute phase diagram and the c = 2 window (PROOF 5). Pairs win while the valley offset δ < E_b (|φ − π/2| ≲ 0.005 at |w| = 1, ≲ 0.017 at |w| = 2, phi_scan.out), and a single valley wins beyond. The dimer + excess-atom mixture in between is not two-component: trimers bind, and mass-imbalanced impenetrable composites are immiscible as n → 0 (C164's interface-cost result, an input here). Predicted c = 2 window in (φ, μ) at small n: empty. This replaces C173's wedge |δ| ≲ n³.
9. [numerical, momentum-resolved ring ED, L = 12 to 36] At n ≥ 0.13 (PROOF 6a, 6b). The vector-chiral state is excluded: even-N ground states sit at K = Nπ/2 (the balanced pair sector), with the polarized sectors 0.46 to 0.66 higher. The grand-canonical density at μ = −2.2 is n = 0.167 at L = 24, reproducing C173. The valley-flip energy at n = 1/6 is 0.210, 0.157, 0.117 at L = 12, 24, 36 (L·E = 2.5, 3.8, 4.2, growing), and the odd-even energy is 0.013 to 0.027. Both are consistent with a pair gap, but these sizes cannot separate a gap from c = 2.
10. [heuristic] At n = 0.17, |w| = 1 (PROOF 6c). On the BCS side (n a_+ = 3.5), a Yang–Gaudin template gives Δ_s ≈ 1e-4 |t₂| and ξ_s of order 1e4 sites. A strongly attractive like-valley channel pushes the other way (toward a gapless or ferromagnetic pseudospin), so a gapless window at intermediate density is not excluded by anything derived here.

## Answers to the spec's items

- Two-body problem, exact. Done in closed form for the unlike pair (claims 2 to 4) and numerically for the like pair (scan_w.out). Zero-energy resonances: the unlike pair at w → 0 and exactly at V₃ = 0 (claim 3), the like pair at |w| ≈ 3.63. Near-threshold structure at φ = π/2: a true bound state, not a resonance.
- J_⊥ and J_z to leading order in n. Formally J_z = 2c(a_s − A_d) < 0 and |J_z| > J_⊥ (claim 6), but the mapping does not apply because the pair binds (claim 4). Three-body input: the trimer is bound (claim 7).
- Map against |w| and n. The pairing holds for all 0 < |w| < 4 at φ = π/2, with crossover density n_* ≈ 1/a_+ ≈ u²/16 at small u (0.05 at |w| = 1, 0.12 at |w| = 2). Clustering is strong at |w| = 2 and weak at |w| = 1.
- c = 2 window in (φ, μ) at small n. Empty (claim 8).
- Ring ED N = 3 to 6. Done (claim 9). It confirms pairing and excludes vector chiral, but it is not decisive at n = 0.17.
- Carry to n ≈ 0.17. Not by a controlled route. The finite-q Tonks extrapolation (smatrix.out) predicts vector chiral and is contradicted by the ED (PROOF 6d). What remains is claim 10 and the Typhon test below.

## What is left open

1. The phase at |w| = 1, n ≈ 0.17. Pairing with a gap too small to see at L ≤ 192 is the most natural reading of everything here, but a gapless two-component window at intermediate density is not excluded. Decisive tests (PROOF 8): the parity gap E(N+1) + E(N−1) − 2E(N) on open chains against L, pair correlations at momentum π against single-particle correlations, and c readings at lower density (n ≲ 0.05) or at larger |w|, where the pair gap is large.
2. Whether the |w| = 1 dilute gas clusters beyond trimers (the tetramer is unresolved). The |w| = 2 gas clusters.
3. The interpretation of the exact V₃ = 0 resonance as a remnant of the legs' SU(2) is heuristic.
4. The full qubit chain. The binding is a robust order-one feature of the second-order model (claim 3), so it should be present at V/Δ = 40. C177 can test it directly by looking for a two-magnon bound state at momentum π below twice the one-magnon minimum.

## Typhon job (pre-registration written before submission)

Effective model (3.2) at φ = π/2 with C173's block_dmrg.py (unmodified, penalty U = 40 across blocks), open chains L = 96 and 192 (Lb = 32, 64), χ up to 128. The two-size estimator is c = 6[S̄(192) − S̄(96)]/ln 2 with S̄ the mean of the three middle block cuts, as in C173. Tasks (typhon/tasks.txt, configs typhon/cfg/):
- W1_225, W1_224: |w| = 1, μ = −2.25 and −2.24 (pilot at Lb = 12: n = 0.091 at μ = −2.24), below the candidate density, where n a_+ ≈ 1 to 2.
- W15_255: |w| = 1.5, μ = −2.55 (a_+ = 11.6, E_b = 0.064).
- W2_290: |w| = 2, μ = −2.9 (pilot n = 0.355, a_+ = 8.3, E_b = 0.12).
- Controls: W1M_215 and W2M_250 (V₂ = V₃ = 0, the MGS limit, two-component by PROOF 2.4), and W1R_220 (V₃ = 0, the exact resonance, no prediction registered).

Registered prediction of the pairing picture. W1_225, W1_224, W15_255 and W2_290 read c(96 → 192) below 1.4, and below the C173 reading at n = 0.17 for the same |w|. The V₂ = V₃ = 0 controls read above 1.6. Report only runs whose final truncation is below 1e-6; if a pairing-picture point reads at or above 1.7 with truncation below 1e-7, the pairing picture fails at these sizes and densities.

Job 17604038 (submitted 2026-11-09, partition typhon, one 48-core node, 24 h limit), scratch /scratch/lustre/jkudlerflam/c175/. Expected outputs: bd_W1_225_{32,64}.jsonl, bd_W1_224_{32,64}.jsonl, bd_W15_255_{32,64}.jsonl, bd_W2_290_{32,64}.jsonl, bd_W1M_215_{32,64}.jsonl, bd_W2M_250_{32,64}.jsonl, bd_W1R_220_{32,64}.jsonl, logs/. Collect with

    mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c175/{*.jsonl,logs} typhon_out/

Each jsonl record has E, n_mid, S (all block cuts), Smid and the χ history (truncation in hist).

## Typhon outcome (coordinator, 2026-11-10)

The job ended at its time limit with eleven of fourteen runs finished. The L = 192 runs of W1_225, W1_224 and W1M_215 hit the 1400 minute cap and left empty files. Readings and the evaluation against the registration above are in REPORT §9.185, from `summarize_typhon.py` (output in `typhon_out/summary.txt`). In short, the test is inconclusive: W2_290 reads 1.98 on a profile that is not of Calabrese–Cardy form, the V₂ = V₃ = 0 control reads 1.19 on an unconverged run, and W15_255 is above the truncation limit.
