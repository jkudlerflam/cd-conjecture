# C132: Lemma P, the energy defect of a periodized ring ground state

Your folder is `../work_C132/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md` first, then all of `contrib/cd_free_route/` (RESULT.md, PROOF.md, run_checks.py), `specs/_bound_background.md` with all updates, and REPORT.md sections 9.71, 9.76 to 9.79 and 9.88. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; no network commands except a browsing tool for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## Setting

This is variant (a) of `_Cd_target.md`, E_L >= e_0 L - K(d) M/L, which gives A <= (6/pi) K(d) M. A bound at a single fixed d, d = 2 in particular, is already a result Jonah wants, so state every result at the smallest d for which you can prove it, and say whether it extends.

`contrib/cd_free_route/` (a Claude agent, 2026-11-06) proved the following, which you should re-verify before building on it.

1. Variant (a) at fixed d is equivalent to a statement about quantum marginals. Every nearest-neighbour marginal of a translation-invariant L-ring state is within dual-M distance K/L^2 of the nearest-neighbour marginal of a translation-invariant state on Z (PROOF section 1).
2. A periodization channel. For a translation-invariant ring state rho with l-site marginal rho_l (sites A = first site, B = the other l - 1) and the Petz map R: B -> BC with R(rho_{l-1}) = rho_l, iterating R produces a translation-invariant state on Z whose every nearest-neighbour marginal is sigma_2*, the (A, B_1) marginal of (id_A (x) R)(rho_l). Hence e_0 <= Tr[h sigma_2*] (PROOF section 3).
3. With delta_L = Tr[h sigma_2*] - E_L/L, delta_L = Tr[h~_{A B_1} ((id (x) R)(rho_AB) - rho_ABC)] for every gauge representative h~ (PROOF section 4).
4. Obstruction. A trace-norm bound on the recovery error cannot give K(d), because it is bounded below by I(A:C|B), which in CFT is (c/3) pi^2/L^2 and carries the central charge c rather than A, and c is unbounded at fixed d. Any proof must use that rho is the ground state of the same h (PROOF section 5).
5. Numerics. L^2 delta_L/M is bounded and decreasing on all tests (Heisenberg 4.05 to 3.33 for L = 6 to 20, SU(3) swap 9.00 to 4.93, long-crossover chains uniformly bounded), and on critical chains it tends to pi A/(3M), twice the best possible constant pi A/(6M).

Lemma P. There are K(d) and l(L) such that for every bond h on C^d (x) C^d and every L >= 2, the translation-invariant ground state (or the equal mixture over the ground space) of H_L satisfies delta_L <= K(d) M(h)/L^2.

## A freedom the agent did not use

Check this first. The proof of Lemma 3.1 in PROOF section 3 appears to use only that R is completely positive and trace preserving from B to BC with R(rho_{l-1}) = rho_l, not that R is the Petz map. If so, then

    e_0 <= min over such R, and over l, of Tr[h ((id_A (x) R)(rho_l))_{A B_1}],

the minimum is a semidefinite program in the Choi matrix of R at fixed rho, and R may depend on h. The rotated and twirled recoveries of Junge, Renner, Sutter, Wilde and Winter (arXiv:1509.07127) are admissible choices with an explicit formula. Prove or refute that the construction extends, and if it does, restate Lemma P with the optimal R as Lemma P'.

## Tasks

1. Verify claims 1 to 3 and the extension above. Reproduce the numerics of claim 5 with `contrib/cd_free_route/run_checks.py`.
2. Numerically, at small L, solve the SDP of Lemma P' (cvxpy or a hand-written solver, whichever is available in the sandbox without installation; if neither is, write the SDP as a SLURM-ready script). Compare the optimal L^2 delta_L/M with Petz on XX, Heisenberg, XXZ at Delta = 0.5 and 2, the SU(3) swap chain, TFIM at g = 1 and g = 1.6, and the C106 two-band chain. Decide whether the optimal R removes the factor 2 on critical chains (L^2 delta_L -> pi A/(6M)) and whether its constants stay bounded as L grows. The optimizing R is evidence of what a proof should construct, so describe its structure.
3. Look for an energy-weighted proof of Lemma P or P', at d = 2 first. Candidate mechanisms, each to be proved or obstructed with a precise reason. (a) Ground-state stationarity. Since rho minimizes Tr[H_L rho], first and second order variations of rho by local channels near the A-B boundary give inequalities involving h itself, which may bound the energy change produced by R. (b) Modular flow. The Petz map acts near the A-B boundary as an imaginary modular flow by K_BC - K_B, which in CFT is a weighted energy density with weight of order x^2/L^2 and inverse velocity (PROOF section 8, heuristic). Decide whether any lattice statement of this kind can be proved without locality of block modular Hamiltonians, which is not available at criticality. (c) An h-dependent R built from the ring dynamics, for example from the ground state of a longer ring or an open chain, for which the defect can be computed in terms of energies alone. (d) Any rigorous lattice tool from `contrib/rigorous_toolbox_map/TOOLBOX.md` that applies.
4. Falsification. On d = 2 bonds, search for growth of L^2 delta_L/M (Petz and, where feasible, optimal R) with L. Include the ferromagnetic end of XXZ, Delta -> -1+ with the ground-space mixture, which is the d = 2 long-crossover face of lesson 3; TFIM with g = 1 + c/L, the near-critical massive mode of lesson 2; and random bonds from the gauge quotient of the M = 1 sphere. Report L^2 delta_L/M against L for the worst cases found. A route whose constant grows with L on an explicit d = 2 chain is obstructed, and you should say which lesson it hits.
5. State exactly what is proved, for which class of chains and at which d, and the single remaining lemma if the general case is not closed.
