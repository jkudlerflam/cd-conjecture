# cd_crossed_fock progress log

Round on C(d), free-boson majorization of the crossed spectrum (companion to C127 and C128). Opus 5.5 agent, 2026-11-06.

## Log

- 2026-11-06. Read `specs/_Cd_target.md`, `specs/C127_crossed_majorization.md`, `specs/C128_worldline_majorization.md`, `specs/_bound_background.md`, `specs/_reconnection_lemma.md`, REPORT.md sections 9.76 to 9.88, and `codex/work_C125/` (RESULT, PROOF, CHECKS, checks.py, audit_checks.py). Folder created.
- 2026-11-06. Built and validated the checkerboard crossed transfer (cflib.py): crossed traces equal the direct Trotter traces to 1e-14 for XX, XXZ 0.5, Heisenberg, XXZ 2, TFIM, XYZ and staggered-field XXZ; realigned gates PSD.
- XX exact (m = infinity) majorant scan done (xx_exact_majorant.py). Finite-m full spectra m = 2..7 (all models), m = 8 (XX). Zero-T deficits D(L) = L(L e0 - E_L) for L <= 18 (zeroT.py): sup at L = 2 for every model, D decreasing in L.
- Single-insertion analysis (insertion.py, analysis_insertion.py, two_insertion.py): Fock model over single-insertion levels dominates Z~ in every case so far; level-by-level counting fails (Heisenberg dimer singlet, SU(2) selection rule).
- Scope update from coordinator: d = 2 positive class has priority. Added TFIM, XY on Ising line, XYZ, staggered-field XXZ (models.py).
- Local machine under heavy memory pressure (swap full from other processes): m = 8, 9 spectra and m = 12, 13 Lanczos moved to Typhon (jobs 17577755, 17577756, queued). eps0 Richardson extrapolation validated on XX against exact p_c (1e-6 at beta = 8 with m <= 11).
- XX majorization theorem with explicit constants (xx_theorem.py), checked against the exact crossed product.
- m = infinity small-ring tables (infinite_m.py) and minimal majorant constants (minimal_constants.py): S = 2, V'/M = 1.2 to 6.9 at K = log 2 over the eight d = 2 chains, binding at L = 2, beta 3 to 8.
- Proved: XX majorization theorem (PROOF 4), conditional Theorem 5.1 ((H1)+(H2)+(H3) imply the target), SU(2) selection rule for single insertions (column representation checked, run_checks I3), Proposition 6.1 (level-by-level Fock domination fails on Heisenberg, with two stated numerical facts).
- Numerical: trace-level Fock domination holds everywhere tested; degraded level-by-level domination holds with lambda >= 0.871 (Heisenberg singlet, m = 8, beta = 8 worst case); two-insertion energies far from additive.
- Typhon jobs 17577755, 17577756, 17577782 never started (queue estimate 2026-11-10); 17577756 cancelled when superseded, the other two cancelled unrun. Heisenberg m = 8 insertion analysis run locally instead.
- RESULT.md, PROOF.md, run_checks.py written; run_checks.py passes all 12 checks in about 90 s.
