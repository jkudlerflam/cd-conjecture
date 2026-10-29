## Background for the runs on the velocity-weighted bound (C98 to C101, and the Fable agents)

Setting. A one-site translation-invariant nearest-neighbour chain H = sum_i h_{i,i+1} on (C^d)^L, with h a Hermitian d^2 x d^2 matrix. Its low-energy theory is a finite set of unitary CFT sectors a with Virasoro central charges c_a and velocities v_a (lattice units, hbar = 1; if left and right speeds differ, use the average). The gauge-invariant bond scale is M = min over h -> h + alpha 1 + o (x) 1 - 1 (x) o of (lambda_max(h) - lambda_min(h))/2.

The proposed bound (consistent with every chain checked so far):

    A = sum_a c_a v_a  <=  pi M log_2 d.

A is exactly the Casimir coefficient, E_0(L) = e_0 L - (pi/(6L)) A + o(1/L) on a ring, and also the CFT entanglement growth rate after a boundary-state quench with equal extrapolation length ell, dS/dt = pi A/(12 ell). Checks: the Heisenberg chain h = S.S has M = 1/2, v = pi/2, c = 1, so A = pi/2 = pi M log_2 2 (saturated); the SU(N) swap chain h = P has M = 1, v = 2 pi/N, c = N - 1, A = 2 pi (N-1)/N <= pi log_2 N. The cleaner Casimir bound A <= 2 pi M is false (SU(N) chains with rank-k symmetric representations, for example SU(8), k = 2 exceeds it by about 5 percent): see `../../contrib/refined_conjecture/`, `../../contrib/casimir_survey/`, `../../contrib/lit_casimir_d2/`, `../work_C14/`, `../work_C35/`. The best rigorous result so far is C69 (`../work_C69/RESULT.md`): with the small-incremental-entangling bound (rate at most 18 M log d per cut) and a boundary-state preparation of length ell_0 realized to relative accuracy alpha, A <= [216/(pi alpha)] ell_0 M log d; the length ell_0 cannot be removed by that route. Also proved there: a crossover bound c <= O(xi_cross log d) from the entropy density.

Important constraint from this campaign: c itself is not bounded by any f(d) (constructions with c = n + 1 at fixed d exist, conditionally; see `../../STATUS.md`), because each added mode can be made exponentially slow. So a proof must use the velocities, and any argument that would bound sum_a c_a alone is wrong. The constructions do satisfy A <= pi M log_2 d.

Rules: tag claims proved, numerical, heuristic, assumed. No probabilities, odds or percent-confidence in any file. Python (numpy, scipy) only, inside your folder; no network commands (a browsing tool may be used for primary literature); never modify files outside your folder. Plain prose without em-dashes. Keep PROGRESS.md and write RESULT.md and PROOF.md.
