# C151 progress log (Claude Opus agent, 2026-11-07, not Codex)

Route: the Destri-de Vega (DdV) contour representation of root sums. Other routes (discrete comparison principle, L versus 2L) are run by other agents and are not touched here.

## Reading (done)

specs/_strip_context.md, specs/C151_ddv.md, specs/_Cd_target.md, specs/_bound_background.md; C138 RESULT and PROOF (Theorems A and B), C142 RESULT and PROOF Part II (Theorem II.6), C145 RESULT and PROOF Parts A and D, C146 RESULT and PROOF Parts 3 and 4, C147 REFEREE (verdict table), C136 REFEREE (head). Solvers bethe_ferro.py and continuum.py copied from work_C146.

## Plan

1. Exact contour identity for root sums on [-Q_L, Q_L] at every filling (argument principle on a rectangle whose vertical sides sit at the counting boundary, where 1 + (-1)^n e^(2 pi i theta) = 2).
2. Positivity: Re sigma_L > 0 in the strip |Im x| < min(eta, pi/4). This gives the DdV hypotheses (no complex zeros) and, through Harnack's inequality, a pointwise lower bound on 2 pi Im theta(x + i eps) by eta L sigma_L(x) at the same real point.
3. Energy: choose a contour that is lifted over the dense set {lambda >= M_0} and stays on the real axis in the sparse tail. Bound the four pieces (tail by the sawtooth, lifted segment, two vertical connectors or Fermi-edge pieces, mismatch).
4. Identify what is still assumed, compare with C146's (G) and (G_c), give margins.
5. Numerics: identity at several heights, free point, CFT finite-size term, Harnack bound, shell margins, and a DdV NLIE solve at half filling against Bethe roots.

## Log

- First pass of the analysis (in head, written up in PROOF.md): the contour identity and the Harnack bound close the core and the Fermi edge without (D1), (W1) or any variation bound. The shell still needs a pointwise lower bound on sigma_L, but only a logarithmic one, eta L sigma_L >= 0.52(1 + delta) log lambda. A bootstrap for that lower bound meets the same geometric obstruction as C145 D.1 (first-order charge of cells adjacent to the hypothesis region), see PROOF Part 5.
- Identity checks (check_half.py, check_finite.py): the contour formula reproduces F - L e_bulk at half filling (free point and small eta) and at finite filling to quadrature accuracy, independently of the height; Fermi-edge pieces give the CFT term when lambda(Q) is large.
- Shell and Harnack checks (check_shell.py, 46 cases to L = 16384, 1 + Delta to 1e-6): Harnack ratio >= 1.53, min sigma/r on the shell 0.99.
- NLIE (nlie.py): Picard iteration diverges for eta <= 0.1; Newton-Krylov with continuation in eta converges and matches Bethe. eta = 0.01 converges at L = 16 (ratio 1.25) and L = 64 (ratio 1.1); L = 256 not run.
- Constants (constants_C151.py): K = 76.72 (C = 386.1) under C146's (G_c); K(c) for c = 1/2, 1/4, 1/10; unconditional K_0(eta) about 31.4 log(1/eta) + 140.
- Coordinator message (C150 reduces the strip to (H_I), implied by any positive comparability on the shell): checked whether the DdV representation can give a positive lower bound on sigma/r on the shell. It cannot by itself (PROOF Part 5: Hilbert-sum log at the bootstrap, geometric obstruction at the tail interface, phase regularity needs spacing regularity). Recorded that (S_c) implies (H_I) with C_I <= 6.07/c and that Theorem 2.1 closes with any c > 0, so the two routes need the same missing input.
- PROOF.md and RESULT.md written.
