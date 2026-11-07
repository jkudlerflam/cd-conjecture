# C132 state

Updated 2026-11-07. Task-scoped memory because writes outside work_C132 are forbidden.

[Proved] Claims 1 to 3 of cd_free_route survive. The dual-M norm on marginal differences equals trace norm. Periodization works for any CPTP recovery fixing the prescribed marginal. PROOF.md gives the support-reduced SDP, its dual, and admissibility of rotated and twirled Petz maps.

[Proved] The source entropy-continuity argument does not prove a trace-norm no-go. Its log L and crossover dependence prevent the alleged contradiction.

[Proved] For the d = 2 ferromagnetic Heisenberg full ground-space mixture, Petz is SDP-optimal, delta = 1/[2l(l+1)]. With l = L-1, Q = L/(L-1). K = 2 includes the L = 2 product convention. Applies also to even-ring XXZ Delta = -1 via staggered rotation. Actual energy bound on this frustration-free class has K = 0.

[Proved] For the quadratic class, momentum-box smoothing plus parity projection gives ||h_cross|| <= M, hence variant (a) with K = pi^2/6. Smallest dimension 2; extends to d = 2^m. This does not prove P prime on that class.

[Numerical] Original quick checks all pass. Small recovery SDPs show no material gain for XX, Heisenberg or SU(3); small gains for XXZ; appreciable gains for TFIM. L = 8 central SCS objectives are not feasible enough to use as optima. Their dual diagnostics and Petz upper bounds are retained. Full large-L coefficient remains unresolved.

[Numerical] Search includes 16 seeded generic M = 1 qubit bonds, XXZ near -1 and TFIM g = 1+c/L. Random bond 13 grows from L = 4 to 6 then decreases through 12. No divergent route constant found. Explicit bond and tables are in RESULT.md. Initial sparse FM endpoint mixtures were incomplete at L = 10,12 and were excluded, then replaced by exact staggered Dicke projectors.

[Assumed, open] Remaining lemma is the uniform upper bound on the optimal prescribed-marginal recovery defect p_L <= K(2) M/L^2, or its all-d version. General interacting qubit case is not closed.

Files: RESULT.md is the task report; PROOF.md contains derivations; PROGRESS.md records work; reproduction.txt and verification.txt are passing checks. JSONL files retain solver status and residuals. random_bonds.npy is task-generated reproducibility data. long_checks.slurm, sdp_large.slurm and falsify_large.slurm are prepared but not submitted, since Typhon is unavailable.
