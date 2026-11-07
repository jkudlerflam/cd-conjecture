# cd_free_route progress log

Opus 5.5 agent, 2026-11-06. Free choice of route for the C(d) round (companion to C126, C127, C128, cd_crossed_fock).

## Log

- Read specs/_Cd_target.md, specs/_bound_background.md, specs/_reconnection_lemma.md, REPORT.md 9.68 to 9.88, codex/work_C125/RESULT.md, specs C126 to C128, and the RESULT.md files of work_C126, work_C127, work_C128 (to avoid duplicating them).
- Route chosen: variant (a) through periodization of ring marginals (see RESULT.md once written). First step is an exact kinematic reformulation, second a universal periodization channel (iterated Petz recovery), tested by exact diagonalization.
- Proved (PROOF.md section 1): variant (a) with uniform K(d) is equivalent to the kinematic inclusion Omega_L in Omega_inf + (K/L^2) B_* of nearest-neighbour marginal sets (ring versus Z), by support functions.
- Proved and checked: momentum box smoothing of a quasi-free ring state gives E_L >= -csc(pi/L) for XX, equality on every even L (ED L = 2 to 12).
- Route: universal periodization by one-step Petz recovery. With A one site, B = l - 1 sites, C one site of the ring marginal rho_l, the Z state built by iterating R_{B->BC} has every nearest-neighbour marginal equal to the (A, B_1) marginal of (id_A (x) R)(rho_AB). Key inequality: L^2 (e_trial - E_L/L)/M <= K(d), which with e_0 <= e_trial gives variant (a).
- Numerical: bounded and decreasing in L at the best l (near L/2) for XX (2.12 at L = 14), XXZ 0.5 (2.77), Heisenberg (3.36), XXZ 2 (2.20), SU(3) swap on L = 3, 6, 9 (9.0, 7.49, 5.39), FM Heisenberg with the full ground multiplet (about 1.1 at l = L - 1), TFIM near criticality (massive modes exponentially suppressed: 0.019 at g = 1.6, L = 14), C106 two-band long-crossover chain (2.0 to 3.3 uniformly in epsilon down to 0.03, L <= 8, l <= 4).
- Numerical: on Heisenberg the trace-norm recovery error is about 4.8/L^2 and the CMI I(A:C|B) about 3.36/L^2 (CFT value (c/3) pi^2/L^2 = 3.29/L^2). So the trace-norm bound already gives 1/L^2 there, but the CMI scales with c, not with c v.
- Running: coupled fast/slow ladder (d = 4) to test velocity weighting.
- Ladder: bounded uniformly in eps (1.37 at eps = 0.01, L = 8). Large L: XX 2.109 and Heisenberg 3.333 at L = 20, approaching pi A/(3M). SU(3) L = 12: 4.93.
- Proved obstruction to closing the key lemma through trace norm (PROOF section 5), and velocity weighting on decoupled stacks (section 6).
- Wrote RESULT.md, PROOF.md (sections 1 to 8), run_checks.py (all 12 checks pass, --long in about one minute). Status: route reduced to Lemma P, open. No commits made.
