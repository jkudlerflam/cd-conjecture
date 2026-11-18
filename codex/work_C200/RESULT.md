# C200 result: the ferromagnetic strip is closed, by a counting bound for the Bethe roots at every filling

2026-11-18. Manager: Claude Fable 5.2 agent. Derivations and referee passes: Astra (gpt-6.2-astra via the shim), seven calls, every output carries the trailer model=gpt-6.2-astra. No output of gpt-6-astra was used. Proofs are in PROOF.md with a tag on every step.

## Route chosen

The discrete comparison principle of C149, that is (R_kappa). C151 Part 5 shows the Destri-de Vega route cannot supply the shell bound, C150 needs the same input, and (R_kappa) is the only reduction that uses the Bethe equations beyond positivity and has no Bethe roots in its hypothesis. Two devices make it provable with lossy constants. First, trial nodes may sit at +infinity (PROOF Lemma 1.1), which removes C149's window b/4 < kappa0 < 1/2 for the shift. Second, the trial nodes beyond the Fermi point may be placed at explicit real positions Q + l A d (l = 1, ..., M an integer index), d the larger of eta and the width of the last continuum cell, where they are supersolution components and each contributes its full kernel tail to the interior nodes.

## Statements proved

1. Theorem H (half filling), rigorous. For eta in (0, pi/2) and even L, n = L/2: |#{roots <= v} - L N(v)| <= 1 for all v, and eta L |sigma_L - r_0| <= (2/pi) eta cot(eta) <= 2/pi at every point. The quantile trial with kappa0 = 1/2 and its top node at infinity is a supersolution, because for the sech density the matched distances satisfy coth(b w_+-/2) = z +- cos(2 alpha): the asymmetry D_y is the integral of one explicit unimodal function over a window of fixed length, and the sawtooth error is at most eta/pi, which equals the twist gain. The constant is sharp. Lemma by Astra (call A), construction by Fable, refereed by a second Astra instance (call C, correct after wording repairs, applied).
2. Theorem 4.2 (regularity of the continuum density), rigorous. |r_Q'| <= C r_Q/rho and |r_Q''| <= C r_Q/rho^2, rho = min(1, eta + distance to the Fermi point), uniformly in eta and Q. The continuum equation is U(x, 0) = U(x, eta) for a positive harmonic function U in a strip, a Neumann condition for the vertical average of U; reflection and Harnack give the estimate. Astra (call B), refereed (call E, correct).
3. Theorem 4.3 (logarithm-free quadrature), rigorous. The sawtooth error of the quantile lattice of r_Q is at most C eta/rho(y) for every y and L. The logarithm of all earlier routes came from estimating 2uT and u^2 T' separately; together they are the derivative of u^2 T(u), which has bounded variation. Astra (call B), refereed (call E, correct after two short repairs, applied).
4. Theorem 4.11 (counting bound at every filling), rigorous. For eta in (0, 1/64], every L and 1 <= n < L/2: |N_x - N_0| <= K/L with a universal K. Astra (call D), refereed (call F, correct as stated). With Theorem H this covers every sector.
5. Corollary 4.12, rigorous: logarithm-free (D) and (D1) at every filling, eta L |sigma_L - s| <= 2K/pi and eta^2 L |sigma_L' - s'| <= 16K/pi^2 at every real point, relative to the continuum Fermi point Q_nu. In particular sigma_L >= (1 - 2K/(pi lambda)) r_{Q_nu}, the shell comparability, with any M_0 > 2K/pi.
6. Corollary 4.13, rigorous given two refereed citations: |F_L(n) - L e_bulk(n/L)| <= K'/L for every L and n, K' = 35 + 181 K + 9 K^2. The derivation (C138 Theorem B at Q_nu, as in C149 Theorem 4.3) is written out in full and was refereed by Astra (call G, accept after arithmetic and convention repairs, applied). Cited: C136 energy formulas and C145 Theorem B.1, both refereed. The third input, Var(v_Q) <= 92.9 (C146 Theorem 2, never refereed), was proved again independently by Astra in call G, with the bound 24a + 4aB <= 88 pi/3.
7. Corollary 4.14: for -1 < J_z/J <= -cos(1/64), every field, DM angle and L, L(E_2L - 2E_L) <= (5K' + pi^2/4) M. The strip is closed. Remark 4.15 (Fable, not seen by the Astra referees, now refereed as a sketch by C204, correct) extends the range to -1 < J_z/J <= -cos(pi/4).
8. Explicit side results, rigorous given C151 Theorem 2.1 (not refereed): at half filling |F_L - L e_bulk| <= 70.2/L, and for magnetization h = ceil(L/2) - n the constant K(h) <= 81 + 10h.

## What remains conditional or unrefereed

- K is universal and finite but not explicit (Harnack constants in Theorem 4.2 and the constants A, M, k of Theorem 4.11 were not tracked). So the strip constant is not a number yet. Numerically L sup|N_x - N_0| <= 0.77 in every case computed, and K = 1 is proved at half filling.
- For the strip the chain now cites only refereed runs (C129, C136, C145 Theorem B.1). The unrefereed C146, C149, C150, C151 are not needed for it: their used parts are re-proved here. C148 for Delta > 1 was refereed as C202, correct after minor repairs (K <= 4310.68 M_0, total 10779.2, largest computed sector error 2.45 M_0 at L = 2).
- This run was refereed by other Astra 6.2 instances (calls C, E, F, G), checked by the manager, and then refereed by an independent Claude Opus 5.5 agent (C204, codex/work_C204/REFEREE.md): no error, Theorem 4.11 and Corollaries 4.12 to 4.14 correct, constants universal and not explicit. Its repairs R1 to R4 are applied in PROOF.md (C136 Section 2 cited for the monotonicity of r_Q, C149 Theorem 1.1 retagged as background, exterior node index renamed, even reflection in Theorem 4.2 justified).
- Remark 4.15 is now refereed as a sketch (C204, correct). It is not a line-by-line proof.
- A correction to the record: the manager's intermediate claim (spec_D.md, Section 4, and an earlier draft of this file) that the gain from the right exterior cannot exceed B(delta)/2 at a sparse edge was wrong, a miscount in the top node's own inequality. Astra's call D removed it.

## Checks run (floating point, scripts in this folder)

- check_lemma.py: Step 5 inequality of Lemma 2.1 (minimum ratio 1.0000000068), sup h, unimodality, identity (2.2) to 9e-13 eta, |J| against its bound in 2000 random cases (worst ratio 0.618). tvD.py: sup|D_y|/eta <= 0.3183.
- half_trial.py, check_roots.py: the kappa0 = 1/2 trial is a supersolution; for Bethe roots with 1 + Delta from 0.29 to 1e-8 and L from 2 to 2048, sup L|N_x - N| is between 0.50 and 0.765 (bound 1) and eta L sup|sigma_L - r_0| <= 0.318 (bound 0.637).
- check_reg.py, tv_gen.py: sup |r'| rho/r <= 1.19, sup |r''| rho^2/r <= 1.71, TV(D) min(1, delta)/eta <= 0.88; the Fourier identity behind the harmonic representation to 4e-16.
- ext_trial.py: the trial of Theorem 4.11 with k = M = 4, A = 2, for 1 + Delta in {1e-3, 1e-4, 1e-5}, L to 1024, five fillings, lambda_Q from 1.2 down to 0.0001: residuals at the exterior nodes between 0.945 and 46.6 (in units eta/L), at the interior nodes between 2.80 and 4.43, all positive.
- check_gen_roots.py: Bethe roots at general filling, 60 cases, 1 + Delta from 1e-2 to 1e-5: sup L|N_x - N_0| <= 0.749, eta L sup|sigma_L - s| <= 0.311.
- gen_trial.py, star_test.py, edge_sign.py: supporting scans of the quantile trials with shifts up to 5.5.
- constants_C200.py: 70.2 and K(h) from C151's formula; the arithmetic of K' is in PROOF Step 8.

## The U(1) qubit class, assembled

Bond J_perp(SxSx + SySy) + D(SxSy - SySx) + J_z SzSz - field, J = (J_perp^2 + D^2)^(1/2), Delta = J_z/J, every field, DM angle and L >= 2.

| Delta = J_z/J | bound on L(E_2L - 2E_L)/M | source | status |
|---|---|---|---|
| Delta <= -1 | 0 | C129 | proved |
| -1 < Delta <= -cos(1/64) = -1 + 1.22e-4 | 5K' + pi^2/4, K' = 35 + 181K + 9K^2, K universal, not explicit | C200 (PROOF Parts 1 to 4) | proved by Astra 6.2 and Fable, checked by the manager, refereed by an independent Claude Opus 5.5 agent (C204), correct after minor repairs |
| -cos(1/64) <= Delta <= 0 | 6e15 delta^(-25/3) at delta = 1.22e-4, about 2.5e48 | C134, C136 | refereed |
| 0 <= Delta <= 1 | 4078 | C139, C140, C142 | refereed (C140 referee, C147) |
| Delta > 1 | 10781 (referee's total 10779.2) | C148 | refereed (C202), correct after minor repairs |

With Remark 4.15 the split point moves to Delta = -cos(pi/4) and the third row becomes 1.7e20. No row is conditional on an unproved hypothesis. C(2) on the U(1) qubit class is finite, with a value that is not explicit.

## Files

PROOF.md, RESULT.md; specs and Astra outputs spec_{A..G}.md, astra_{A..G}.out; scripts half_trial.py, gen_trial.py, star_test.py, edge_sign.py (edge_sign.log), ext_trial.py, tvD.py, tv_gen.py, check_lemma.py, check_roots.py (check_roots.json), check_gen_roots.py (check_gen_roots.json), check_reg.py, constants_C200.py (constants_C200.json); fftsolver.py and bethe_ferro.py are copies from work_C149.
