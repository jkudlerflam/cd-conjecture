# C201 result

2026-11-18. Manager Claude Fable 5.2, derivations and referee reports by Astra (gpt-6.2-astra via the shim, five calls). Proofs and tags are in [PROOF.md](PROOF.md). No human review. The general target, a finite K(d) with E_L ≥ L e₀ − K(d) M/L for all chains, was not reached.

Notation: g_d(L) = sup over bonds with M = 1 of (L e₀ − E_L). The trivial bound is g_d(L) ≤ 2, and Lemma P is g_d(L) ≤ K(d)/L.

## Theorems

A (all chains, no hypothesis on the chain). Fix d, n ≥ 1, 0 < δ ≤ 1, and set B₁ = 1, t_i = B_i + 1, q_i = d^(t_i), B_{i+1} = (⌈2(2q_i − 1)² δ⁻² ln q_i⌉ + 1) t_i, f = 3δ + 2√(2δ) + 4√(2 ln d/n). For every L ≥ 2B_{n+1}, every bond and every k ≥ 2,

    L e₀ − E_L ≤ f M,    E_{kL} − k E_L ≤ k f M.

So g_d(L) → 0 and sup_h (E_{2L} − 2E_L)/M → 0 as L → ∞ at fixed d, and every translation-invariant state on the L-ring has its two-site marginal within f/L = o(1/L) of the marginals of translation-invariant states on ℤ. Rigorous given the recovery theorem of Junge, Renner, Sutter, Wilde and Winter (Corollary 4.1 with Remark 2.4(2)).

Quantitative caveat. The bound beats the trivial bound 2M only when 4√(2 ln d/n) + 3δ + 2√(2δ) < 2. At d = 2 this needs n ≥ 6 with δ of order 5·10⁻⁴ or smaller, or n ≥ 8 with δ = 10⁻², and in either case L ≥ 2B₇ or more, an iterated-exponential scale (already B₂ is of order 10⁶ to 10⁹ for these δ, and at n = 1, δ = 1 the threshold is L ≥ 548 with f ≈ 10.5). So Theorem A is qualitative and has no content at any length a computer or a physical system reaches. The Astra referee call gives a construction at d = 4 (read, not verified independently) showing that the full-window decoupling step (2.2) of the proof cannot be carried out on polynomial scales. Other joins are not excluded by it.

Consequences. Rigorous: no chain has a ring energy a fixed amount below L e₀, and the positive part of the doubling defect vanishes uniformly. Rigorous for chains with a positive spatial transfer operator only: a third law in the crossed channel, lim [L E_{L+1} − (L+1) E_L] = 0 (false without positivity, as h = Z⊗Z shows). Physics-level only (Corollary A3): if the finite-size form E_L = L e₀ − (π/6L) Σ_{ξ_a ≤ L} c_a v_a + r_L holds with r_L uniformly small across a family, the velocity-weighted central charge visible at scale L is o(L) M, which excludes an extensive fast tower. That uniformity hypothesis does most of the work, and the statement is empty below the scale above. Theorem A does not bound A, since deficits M N/(L + N²) obey it with A = (6/π) M N.

B (exchange class, every d ≥ 2). Site space ℋ₀ ⊕ ℋ₁ of any dimensions, bond −t(W + W†) + V n⊗n + (μ/2)(n⊗1 + 1⊗n) with W any unitary exchanging the two species, entangling or not. Then E_L is at least the scalar (XXZ) ring energy at minimizing flux, e₀ is the scalar one, M ≥ t, and
- V = 0: L(L e₀ − E_L) ≤ (4 − 8/π) M for all L, sharp at d = 2, L = 2, and A ≤ 2M. Fully rigorous, no citation, uniform in the field, including a field tuned with L.
- V ≤ −2t: E_L = L e₀. Fully rigorous.
- 0 ≤ V/2t ≤ 1: L(L e₀ − E_L) ≤ 1682 M, given the project's XXZ theorem (REPORT §9.100, refereed in §9.104, stated there for Δ < 1, with the endpoint V = 2t by continuity at fixed L).
- −1 + δ ≤ V/2t ≤ 0: the same with 2·10¹⁵ δ^(−25/3), given REPORT §9.94.

The hypothesis is structural and does not assume the conclusion. It is restrictive: the only fast sector is the scalar charge mode, and the internal register enters ring spectra only as a boundary twist (proved for every W and L in PROOF §3.3, Step 3′). The class is not new. The reduction is the spin-charge factorization known from the U = ∞ Hubbard ring (Ogata and Shiba 1990) and the XXC models (Maassarani 1998), extended to entangling W. The new content is the bounds in units of M.

C (critical ferromagnetic Potts chains, d = 3, 4). For H = −Σ e_i in the Q-state Potts representation of the periodic Temperley–Lieb algebra, Q = 2, 3, 4: E_L equals the lowest level of the XXZ ring of 2L sites at Δ = √Q/2 in the sector Sᶻ = 0 with twist 2γ (√Q = 2 cos γ), e₀ = 2ε₀^XXZ, and the whole ring spectrum lies in an explicit finite list of real-twist XXZ sectors. These three statements (C1 to C3) are known: Alcaraz, Barber and Batchelor (1987–88) and the affine Temperley–Lieb literature (Martin–Saleur, Graham–Lehrer, Richard–Jacobsen, Morin-Duchesne–Saint-Aubin). PROOF §4.8 gives a self-contained rigorous proof of the ground-state identity. New here are only (C4) M = (√Q + 1)/2 for Q = 3, 4 and M = 1 for Q = 2, fully rigorous, and (C5)

    L(L e₀ − E_L) ≤ 841 = 841 M (Q = 2), 615.7 M (Q = 3), 560.7 M (Q = 4), for every L ≥ 2,

conditional on the same XXZ theorem (with Q = 4, Δ = 1, by continuity). At Q = 3 the infrared is a neutral, non-Gaussian sector (c = 4/5). Its bound is inherited from the XXZ host, of which the Potts ground state is one level. It is not a new method for neutral sectors and gives nothing for a neutral sector without such a host.

D, E, F, G. Decoupled stacks of the proved classes have constant at most K_max log₂ d. A band whose r-th hopping moment in lattice units is at most D M has Σ over its Fermi points of |v|^(1/(r−1)) ≤ 2π (D M/r)^(1/(r−1)), the velocity-tail estimate with exponent below one that REPORT §9.93 asked for, on that class. A superlattice whose inter-cell hopping has rank R has at most 2R fast crossings. Permutation-invariant ring states are within 30 d⁵/L² of Ω_∞.

## What is conditional

- A: the published recovery theorem (Corollary 4.1 and Remark 2.4(2) of arXiv:1509.07127), confirmed from the text of the paper by the C203 red team. Universality of the recovery map is not used.
- B at V ≠ 0 and the bound in C: the project's XXZ finite-size theorems, which are AI-produced and AI-refereed.
- The spectral inclusion in C: the cited intertwiner theorem, with theorem numbers as quoted by Astra.
- Every translation into a bound on A uses the infrared assumption.

## Targets that failed, and why

- General Lemma P (B1 in full, B2). Recovery-map joins pay the square root of an information quantity per join. For conformal ring states the best such join has error of order L^(−x/(2x+1)), x the lowest scaling dimension, and a linear-response version would give at most A ≲ M Σ c_a. The missing estimate is a join whose trace-norm error is O(1/L) (PROOF §7.3, §7.7).
- Fast tower (B3). None found. Excluded: an extensive tower (physics-level only, Corollary A3 of Theorem A, under a uniformity hypothesis), internal registers moved by unitary exchange (Theorem B), sawtooth bands with bounded hopping moment and superlattices of bounded hopping rank (Lemmas E, F), mean-field ring states (Proposition G). Open: neutral interacting factors, weakly charged modes with vanishing Luttinger weights, dressed strings. A tower need not show an order-M deficit at any accessible size (PROOF §7.4), which weakens earlier searches that looked for one.
- Reversed-sign Potts (Q = 2, 3). Spectral inclusion holds, but the bound does not transfer, because the Potts representation keeps only excited host levels (PROOF §7.2).
- Site insertion. E_{L+1} ≤ (L+1)E_L/L + KM/L² would imply Lemma P, and it is false on incompatible lengths (PROOF §7.5).
- B4. No invariant bounded by a function of d was found beyond g_d(L) itself.

## Referee

Theorems B, C and items D to G were refereed by an independent Claude Opus 5.5 agent with its own exact-diagonalization code, see [../work_C203/REFEREE_BC.md](../work_C203/REFEREE_BC.md). Verdict: B and C correct after one shared repair (the cited XXZ theorem stops at Δ < 1, so V = 2t and Q = 4 need a continuity step) and small fixes (the boundary-twist proof moved into PROOF §3.3, M = 1 at Q = 2, literature credit and a narrower novelty claim), D, E, F, G correct. All of these are applied in this version. Theorem A was refereed by a second Astra call and then by an independent Claude Opus 5.5 red team, see [../work_C203/REDTEAM_A.md](../work_C203/REDTEAM_A.md). Verdict: (A1) and (A2) correct given Corollary 4.1 with Remark 2.4(2) of the recovery paper, confirmed from its text, with no gap found in the scale selection, the decoupling, the join geometry, the gauge handling of M, the k-fold bookkeeping or the incommensurate lengths. The statement appears to be new and is assembled from standard tools (the literature found gives only O(1/n)). The red team's corrections are applied here: the positivity proviso on the third law, the physics-level status of the tower exclusion, the narrower scale obstruction, the positive part of the doubling defect, and the quantitative caveat.

## Checks

Exact diagonalization with my own scripts (PROOF §8): Potts inclusion to L = 5 and the ground-state identity to L = 7, the Potts Casimir ratio (0.46 M at Q = 3 against the proved 615.7 M), Haar-random exchange unitaries at d = 3, 4, the sharp V = 0 bound (saturated at d = 2, L = 2), the SU(3) swap chain as a control outside the class (2.67 M at L = 3, above the class constant 1.45 M), and the lemmas inside Theorem A (marginal correction, covariance identities, the chain-rule decoupling estimate). Theorem A itself cannot be tested at reachable sizes. Heisenberg, XX and XXZ enter as the scalar members of class B and the hosts of class C, and transverse Ising is the Q = 2 Potts chain.

## Provenance

Astra gpt-6.2-astra (trailer `model=gpt-6.2-astra` checked on all five outputs): first proof of A (`astra_A2_thirdlaw.out`), its referee report with the shorter proof and the doubling statement (`astra_A4_referee_T1.out`), C (`astra_A1_potts.out`) and its referee report with the second proof of the ground-state identity (`astra_A5_referee_potts.out`), the arbitrary-W form of B and Lemmas E, F (`astra_A3_adversary.out`). Claude Fable 5.2: choice of targets, all specs, the label-conserving precursor of B, Proposition D, PROOF §7.3, all numerics, and the line-by-line checks. No gpt-6-astra (previous generation) output was used. Files: PROOF.md, RESULT.md, `spec_A*.md`, `_common.md`, `astra_A*.out`, and the scripts `potts_xxz_inclusion.py`, `potts_table8_check.py`, `potts_casimir.py`, `xxc_check.py`, `exchange_class_check.py`, `b2_bound_check.py`, `t1_lemma_checks.py`, `t1_chainrule_check.py`.
