# C140 result

This run was done by a Claude Opus agent, not by Astra. 2026-11-07. Proofs in [PROOF.md](PROOF.md), log in [PROGRESS.md](PROGRESS.md). No referee pass yet.

## Bottom line

[Proved] (M) between the outermost roots, with margin, for 0 <= Delta <= 63/64, every L >= 2 and 1 <= n <= L/2: along the branch continued from Delta = 0 the ground-state roots satisfy sigma_L >= theta(gamma) r0 on [x_1, x_n], theta = 1 - 8(2-ln2)(pi-gamma)/pi^3 > 0, and the negative mass of sigma_L beyond the outermost roots is at most 1/L. The Yang-Yang Hessian is positive definite along the branch, so the identification hypothesis (I) of C139 holds at both parities on this range (PROOF Theorem 4.1).

[Proved] Consequence: C139 Theorem A is unconditional on 0 <= J_z/J <= 63/64, L(E_2L - 2E_L) <= 345600 M for every field, DM angle and L >= 2. With the sawtooth bound 3/(2L) that follows, C139 Section 8.2 gives (19935 + pi^2/4) M on the same range (degree-counting rescaling of C139's 2215, not redone line by line).

[Assumed, single remaining repulsive estimate] (L*): L r0(x_n) >= ((2-ln2)/pi^2) b = 0.1324 b along the branch. Everything else in the proof is uniform in Delta up to 1. Proved: L r0(x_n) >= (pi/(8(pi-gamma))) b, which suffices exactly for gamma > 0.176 and tends to 0.125 b as Delta -> 1. [Numerical] L r0(x_n)/b >= 0.442 in every case (L 4 to 1024, all sampled fillings, Delta to 0.9999). Given (L*), the bound (19935 + pi^2/4) M holds on all of [0,1], Heisenberg included.

[Not proved] (D) and the C138 continuum bounds. For Delta < 0 the hole kernel is negative with norm B-1, so the method does not transfer; (D) is an accuracy statement, not a sign statement. No common lemma found.

## Method

Dual identity sigma = r0 + J*(sigma dx) - (1/L) sum J(.-x_j) with J >= 0 unimodal. Pairing each gap (mass exactly 1/L) with its adjacent root leaves only the two roots flanking x, so sigma >= r0 - J(0)/L - (outer correction), with J(0) <= ln2/(pi gamma) and the outer correction <= 2(1-ln2)/(pi gamma L) from positive Neumann comparisons on {|y| > x_n}. The upper comparison sigma <= s_(x_n) there plus the exact tail mass 1/(2L) + c - nu/2 locates the outermost root. A continuity argument in gamma closes the loop.

## Numerics

sigma_L >= 0.952 r0 on all of R in every case (Delta 0.01 to 1, L 2 to 4096; minimum at L=2), ratio_scan_{A,B,C}.json. Lemma checks in check_lemma1.json (dual identity residual <= 1.5e-3, the J-grid error; all lemma margins positive). The naive continuum tail inequality fails far beyond x_n for Delta >= 0.5 (tail_cont.json), so strict positivity beyond the last root is not obtained this way, though only the 1/L negative-mass bound is needed downstream.

## Assembly (J_z/J in [-1,1])

- [0, 63/64]: proved, C <= 345600 (C139 Theorem A plus Theorem 4.1 here; C129 reductions, C136 seam gluing).
- (63/64, 1]: conditional on (L*), C <= 19938.
- [-1+delta, 0]: C136, 6e15/delta^(25/3). Near -1: conditional on (D) and C138's four continuum bounds, C <= 628. Delta = -1: zero defect (C129).
- Conditional on (L*) and (D): one constant C = 19938 on the whole class.

Files: ratio_scan.py, check_lemma1.py, location_scan.py, tail_cont.py, edge_scan.py, tail_profile.py, with .json outputs; solvers copied from C139 and C138. ED at L <= 12 not rerun (C139's check, max error 6.2e-14, covers the same solver).
