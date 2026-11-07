# C141 result: the fast-sector budget

This run was done by a Claude Opus agent (2026-11-07), not by Codex. Derivations are in PROOF.md. Scripts: xx_chi.py, fib_emergence.py, cauchy_sv.py (NumPy/SciPy, each under 10 s locally).

## Verdict

[derived, conditional] No unconditional physics-level derivation of C(d) < infinity was obtained, for d = 2 or general d. What was obtained is a reduction, plus a conditional derivation at d = 2.

1. Compactness reduction (PROOF 5.1). B_d (bonds modulo gauge, scale and on-site unitaries) is compact, of dimension d^4 - 2d^2. So sup A/M < infinity if and only if A is locally bounded at every bond, and every hidden tower is a perturbative problem around a single accumulation bond h*.
2. Regular bonds are harmless whatever their crossover length (weak-perturbation c-theorem, PROOF 5.2). C106 is a continuation of a regular limit bond.
3. Emergence lemma (PROOF 5.3, Schrieffer–Wolff plus Lieb–Robinson). Sectors emerging from a gapped degenerate manifold have v_em <= C|delta|M. Sectors from z > 1 bands have v ~ M k_F^{z-1}. Neither can produce a fast sector with a long crossover. [numerical] At the d = 2 Fibonacci point, v/delta tends to the projected-model velocity linearly in delta.
4. d = 2 (PROOF 6). There are three singular strata: the Fibonacci classical point, the ferromagnetic frustration-free points (Bravyi–Gosset classification), and the quadratic Lifshitz points. Each is single-species, with emergent c <= 1 and v -> 0. Given (IR), the weak c-theorem and completeness of this list, C(2)_asymptotic < infinity. With non-coexistence, c_tot <= 1 and A <= 4eM. Completeness of the strata at strong coupling is the weakest link. No qubit bond with c > 1 is known.
5. General d (PROOF 8). C(d) < infinity follows from (P1), a multi-velocity c-theorem for the central charge carried above each velocity, and (P2), a single-level emergence budget.

## Budget candidates

- (a) Ground energy [derived, negative]. A sector with crossover xi lowers the energy density by about c v/xi^2. The free bound A <= pi|e_hop| is really the Kohn f-sum with K = 1. Its interacting form, sum K v <= 4 pi M, misses neutral sectors and sectors with small K.
- (b) Energy-current f-sum [derived, negative]. The per-sector weights are S_eps = c v^2 q^3/(24 pi), m_1 = c v^3 q^4/(24 pi), chi_eps = c v q^2/(12 pi) (checked on XX to five digits, including the factor 1 - n^{-2}). The lattice side is O(q^2) and the CFT side O(q^4), so the bound weighs sectors by v^3/xi^2. A sum rule is crossover-free only for conserved densities of scaling dimension 1 (charges).
- (c) Entanglement rate [derived, negative]. It weighs by v (p = 1) and needs ell > every crossover.
- (d) Operator rank [derived, partial positive]. Rank bounds hold only at first order. At second order, Zolotarev decay of Cauchy port matrices gives a crossover-free budget A_em <= nu eta^2 n a [1 + (4/pi^2) log(4b/a)], with nu <= d^2 - 1 ports, n the cell length and a the smallest intermediate gap. The bound is verified numerically. A is bounded for gapped cells and grows at most like log n for critical cell interiors.

## Hidden towers and the sharpest obstruction

[derived] The equal-velocity hierarchy of C133 reduces to one of two things: a single bond with an infinite fast hierarchy (Casimir term ~ M log L), or a fast sector emerging late. Both require a z = 1, Lorentz-violating flow that increases the fast central charge. The sharpest remaining obstruction is (P1) together with (P2) beyond second order and quadratic payloads: the harmonic tower (v_k ~ M/k), produced as many sectors at one emergence level. At second order it is capped at log growth in the cell length. For the rigorous program, the d = 2 Fibonacci stratum is new to the rigorous list. It maps to free fermions on L - N sites.
