# C140: The two Bethe-root estimates that close the U(1) spin-1/2 class

Your folder is `../work_C140/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md`, then all of `codex/work_C138/` and `codex/work_C139/` (RESULT.md, PROOF.md: these are your starting points), and `codex/work_C134/`, `codex/work_C136/` (RESULT.md, PROOF.md, REFEREE.md; both refereed correct), and `codex/work_C129/RESULT.md`. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md).

## Setting

For the XXZ ring H = J sum (Sx Sx + Sy Sy + Delta Sz Sz) at minimizing flux, the whole U(1) spin-1/2 class (all fields, DM angles and lengths) with J_z/J in [-1, 1] now reduces to two estimates on the ground-state Bethe roots.

(D), from C138, at the ferromagnetic end (Delta = -cos eta near -1): the pointwise density estimate |sigma_L - r_{Q_L}| <= C_D/(eta L) wherever L eta r_{Q_L} >= 4 C_D, together with C138's four continuum bounds (its PROOF section 4.2). Given these, K_* <= 125 and the attractive class closes with C <= 628. Numerically eta L sup|sigma_L - r| lies in [0.23, 0.32] for L = 16 to 8192 and 1 + Delta from 1e-2 to 1e-6.

(M), from C139, on the repulsive side 0 <= Delta <= 1: the counting density sigma_L(x) = a_1(x) - (1/L) sum_k K(x - x_k) is nonnegative on the Fermi interval [-Q_L, Q_L] at the ground-state roots. Given (M), the ordered-root hypothesis (I) holds at both parities and L(E_{2L} - 2E_L) <= (2215 + pi^2/4) M on J_z/J in [0, 1], Heisenberg included. Numerically (M) holds for Delta up to 1, L up to 4096 and all fillings, thinnest at the outermost root.

Both are pointwise comparisons between the discrete counting density and its continuum limit, for the same Bethe equations on the two sides of Delta = 0. C139's hole kernel (I + K)^(-1) K, positive with norm below 1/2 for Delta >= 0, and C138's observation that the pointwise error is the uniform quantity while the cumulative error is not, are the natural tools.

## Tasks

1. Prove (M) for all 0 <= Delta <= 1, all L and n, with any explicit margin. Candidate routes: a discrete maximum principle for the linearized counting equation; a pointwise comparison of sigma_L with the continuum density r_Q using the proved bounds r_0 <= r_Q <= (1 + 4/pi) r_0 and the hole-kernel series; an induction on n or a continuity argument in Delta from Delta = 0, where sigma_L is explicit. Treat the outermost roots separately if needed.
2. Prove (D) and C138's continuum bounds near Delta = -1, uniformly in eta and density. If a single pointwise discrete-continuum lemma covers both (D) and (M), state and prove it once.
3. Assemble the complete theorem: L(E_{2L} - 2E_L) <= C M on the whole U(1) spin-1/2 class, J_z/J in [-1, 1], every field, DM angle and L >= 2, with one explicit C, listing exactly which earlier results (C129, C134, C136, C138, C139) each region uses.
4. Check every new inequality numerically with the existing Bethe solvers (Delta from -1 + 1e-6 to 1, all fillings, L up to a few thousand), and ED at L <= 12.
