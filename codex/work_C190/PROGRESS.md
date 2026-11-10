This run was done by a Claude Opus 5.5 agent (C190, 2026-11-09), not by Codex.

# C190 progress log

1. Read the spec (specs/C189_C191_decisive_tests.md, section C190), REPORT 9.147, 9.152, 9.155, 9.158, C177 RESULT and dmrg2.py, C173 PROOF, the collaborator's contrib/d2_dmrg_priority_20261109/RESULT.md and SIMPLE_PARITY_CANDIDATE.md, and contrib/d2_tenpy_20261109/RESULT.md. TeNPy is not installed locally, so a new solver was written on numpy (no installs).

2. z2dmrg.py: two-site DMRG with exact Z₂ parity blocks (no penalty term), block-sparse environments and matvec, White mixer, checkpoints, exact final entropies, density and bond observables, one-body matrices, C_P, energy variance. models.py: the four bonds, the second-order map, the staggered real gauge for the symmetric line, the XX and two-leg XX controls.

3. validate.py, L = 10, both sectors: all candidate bonds agree with ED to 3e-12 in energy and 5e-10 in entropy (details in PROOF 5). The two-leg control's odd sector is degenerate and is not used.

4. Benchmark against the collaborator, rational bond, L = 32, χ = 128: E_even = −330.0058807140 (theirs −330.0058806527), E_odd = −330.0059531605 (theirs −330.0059530758), parity splitting −7.24e-5 (same), S̄₁ (odd) 1.35352 (theirs 1.35324). But the diagonal-preconditioned Davidson needed 200 to 400 matvecs per solve with residuals stuck at 1e-7 to 1e-6, the same problem as the collaborator's.

5. Local-solver diagnosis (diag_local.py, diag_sigma.py, on the centre two-site problem of that state, D = 128): the lowest local gaps are 5.5e-4, 1.1e-3, 1.4e-3 (0.044, 0.089, 0.11 |t₂|), the diagonal spans 186. Diagonal preconditioning reaches 1e-10 only after 1200 to 1500 matvecs from a converged start, independent of the floor (1e-4 to 1). A Kronecker-sum preconditioner without the central ZZ was worse (it puts doubly occupied central bonds at low energy). The (s₁, s₂)-block Kronecker preconditioner of PROOF eq. (5.1) reaches 1e-9 in 95 matvecs from a perturbed start and 150 from a random one. With it, every local solve of a sweep reaches its target with 60 to 70 matvecs on average. Validation repeated with it (same accuracy).

6. Matvec profiling: about 95 GFLOP/s on 2 local threads at D = 512 complex after replacing the W₁W₂ tensordot by a batched 20 × 20 matmul.

7. Local pilot (validation only, not part of the registered criterion): p153 at the collaborator's field, L = 64, both sectors, rand:1, χ up to 256 (local/out/, local/*.log). Last sweep of each χ stage (energies in the physical convention, Δ = 1):

   | χ | E_even | S̄₁, S̄₆, S̄₁₂ (even) | E_odd | S̄₁, S̄₆, S̄₁₂ (odd) |
   | :-- | :-- | :-- | :-- | :-- |
   | 128 | −650.0414814402 | 1.5993, 1.5944, 1.5804 | −650.0414332127 | 1.5811, 1.5761, 1.5617 |
   | 192 | −650.0414849170 | 1.6125, 1.6070, 1.5916 | −650.0414359528 | 1.5924, 1.5868, 1.5712 |
   | 256 | −650.0414852592 | 1.6142, 1.6087, 1.5930 | −650.0414362376 | 1.5938, 1.5882, 1.5724 |

   All local solves reach their targets (1e-8 at χ ≤ 128, 1e-9 at 192 and 256), with 50 to 70 matvecs on average, and the discarded weight is 2e-10 at χ = 256. The even sector is lower (L·Δ_GC/|t₂| = +0.25), n_mid = 0.2183 (even, N_tot = 16.014) and 0.2018 (odd, N_tot = 15.025). The collaborator's χ = 128 value, S̄₁ = 1.60106 (even), sits next to this χ = 128 state, but S̄₁ still rises by 0.015 between χ = 128 and 256 while the discarded weight at χ = 128 is only 9e-9. So at these bonds a small discarded weight does not certify the entanglement, and the χ = 128 states at L = 128 behind the collaborator's downturn are likely to underestimate S̄(128) by more than this. That is the first thing the Typhon job tests. Against the effective model at equal n_mid (C176 interpolated), S̄₆ is lower by 0.23 at L = 64, close to C177's offset at V/Δ = 8. PROOF 2 gives one reason: the open ends of the qubit chain carry a chemical potential lower by V/2.

   The two parity sectors differ in n_mid by about 1/L, which at dS̄/dn ≈ 2.4 moves S̄ by about 0.04. That is why the density-corrected secondary estimator was added to the registration before any production output.

8. Registration written in RESULT.md, then the Typhon jobs 17606939 (c190a) and 17606940 (c190b) submitted.

9. Overlap steps and first data (2026-11-09, 20:23 to 22:05 EST). The coordinator launched both pipeline bodies as overlap steps inside two donor allocations (RESULT, Typhon jobs). I answered the xx2 validate question (odd-sector degeneracy at L = 10, RESULT), tested the collection commands on the partial output (the remote globs need quotes under zsh) and took two snapshots into typhon_out/. All five field tunings finished (RESULT R2). The tuning points show that n_mid(μ_eff) is a staircase in steps of two particles at fixed parity, so the linear fit fixes the field to about half a step. The XX control at L = 256 finished and agrees with the exact entropies to 4e-8. The first candidate-bond entropies appeared when the p153c L = 128 runs reached the mixer-off χ = 256 stage. With the local L = 64 pilot at the same χ they give a ground-branch c_two(64 → 128) of 2.02, 1.99, 2.06 at W = 1, 6, 12, two initializations agreeing to 7e-10 in energy, and energies 3e-5 to 4e-5 below the collaborator's χ = 128 states at L = 128 (RESULT R4). The collaborator's downturn is therefore not reproduced at χ = 256, and the entropy sequence in χ (1.686, 1.817, 1.848 at χ = 128, 192, 256) points to their χ = 128 truncation as the cause. This is an interim reading. The registered final stages (χ = 512 at L = 128, 768 at L = 256) were still running, and the sweep times measured so far put the complex L = 256 runs at 7 to 40 h per χ = 768 sweep.
