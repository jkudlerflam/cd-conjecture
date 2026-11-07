# C139: The repulsive U(1) spin-1/2 class, 0 < Delta <= 1

Your folder is `../work_C139/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md`, then all of `codex/work_C134/` and `codex/work_C136/` (RESULT.md, PROOF.md and REFEREE.md; both proofs were found correct by independent referees), `codex/work_C129/` (RESULT.md, PROOF.md: the reductions for field, DM and twist), `codex/work_C138/RESULT.md` if present (the ferromagnetic crossover, a parallel run), and `contrib/rigorous_toolbox_map/TOOLBOX.md` sections 4 and 7. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md).

## Setting

C134 and C136 proved, for the XXZ ring H = J sum (Sx Sx + Sy Sy + Delta Sz Sz) at minimizing flux and Delta = cos gamma in [-1 + delta, 0], that every sector satisfies |F_L(n) - J L e_bulk(Delta, n/L)| <= K_all(delta) J/L, which with the C129 reductions gives L (E_{2L} - 2E_L) <= C(delta) M on the whole attractive U(1) class. Their real-axis argument uses that the scattering kernel is negative. For 0 < Delta < 1 (0 < gamma < pi/2) the kernel is positive, with L1 norm 1 - 2 gamma/pi < 1, so the Neumann series still converges, and the norm tends to one only as Delta -> 1 (Heisenberg), where marginally irrelevant logarithmic corrections appear. Exact anchors: Delta = 0 (free fermions, C129) and Delta = 1/2 (supersymmetric twisted identity, C134, verified numerically). Yang–Yang convexity identifies the ground state for all -1 < Delta < 1.

## Goal

Prove L (E_{2L} - 2E_L) <= C M on the U(1) class with J_z/J in [0, 1 - delta], all fields, DM angles and L, with explicit C(delta), and then decide the Heisenberg endpoint.

## Tasks

1. List every place where C134 and C136 use the sign of the kernel (monotonicity of the counting function, positivity in the quadrature, comparison and clipping steps, the resolvent bound) and, for each, give a replacement valid for positive kernels with norm q < 1, or show it fails. A two-sided Neumann estimate with |K|, or the Yang–Yang convexity (which gives a monotone counting function for all Delta in (-1, 1)), are natural candidates.
2. Prove the compact theorem on Delta in [0, 1 - delta] at all fillings with explicit constants, uniform in density, and assemble it with the C129 reductions into a statement on the repulsive U(1) class. Check it against the exact anchors Delta = 0 and Delta = 1/2.
3. Heisenberg. Only an upper bound on the doubling defect is needed, and the data decrease (4, 2.79, 2.62, 2.56 at L = 2, 4, 6, 8). Decide whether a bound uniform as Delta -> 1 can be obtained (rescaled rapidities, as C138 does near the ferromagnet, or a direct argument at Delta = 1 with continuity), and state precisely what fails if not. Note that a common bound for Delta < 1 extends to Delta = 1 by finite-matrix continuity.
4. Numerics: with the C134/C136 Bethe solvers, check every intermediate inequality on Delta in {0.1, 0.5, 0.9, 0.99, 0.999}, densities 1/L to 1/2, L up to a few thousand, with ED cross-checks at L <= 12, and report L (F_{2L} - 2F_L)/M at half filling as Delta -> 1.
