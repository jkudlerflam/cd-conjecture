# Referee report on the input-S runs C48, C58, C51, C52, C49

Referee: Claude Opus 5.5 agent, acting as an independent referee (mathematical physics, cluster expansions, Lieb-Robinson and quasi-adiabatic methods), 2026-10-28. This is an AI review. It is not a human review.

Scope. The fixed-d₀ = 4464 construction gives c = n + 1 for every n conditional on the Wilsonian input S (`contrib/expert_note/note.tex`, Section "S"). REPORT.md Sections 9.25 to 9.30 summarize Astra runs C47 to C63, which claim to prove several pieces of S. This report referees five of them, in order of importance: C48, C58, C51, C52, C49. Nothing outside this folder was modified.

## Summary table

| Claim | Status (this referee) | How it was established and checked |
|---|---|---|
| C48. Resummed massive-sector matching on the enlarged rooted cell model: all-degree retarded kernels with exponential space-time decay, uniform in cell number and in β ≥ 2T_n/Δ, correction to the retained second-order Hamiltonian linear in the C45 remainder | proved, for the enlarged model and under its stated hypotheses (bounded local fermionic interval representation, complete-gap heat trace, C45 block-diagonal form at the strengthened exponent). No fatal or serious gap found in the argument. Two fixable presentation gaps. | Read PROOF.md Sections 2 to 7 line by line, recomputed the key inequalities (8), (12) to (17), (21), (32), (33). Reran checks.py in a scratch copy (identical output). New multicell checks of the embedding, the pressure bound (33), the C46 resummed formula and the u-uniformity of the leakage-linear error. |
| C58. Transfer from the enlarged model to the physical fixed-length ring, pressure and free energy to activated accuracy C T e^(−δ/T), marker number summed or at fixed interior density | proved conditional on C48 and on imported collaborator lemmas. The new pieces (fixed-marker pinning, completion removal, canonical extraction, root bracketing) check out. Serious dependency gap: several imported lemmas have only AI review, and I re-derived only some of them. | Read PROOF.md fully; re-derived (4), (5), (8), (11), (13) to (17), (21), (24), (26), (30). Reran checks.py (identical output). New checks of the pinning lemma on spin-1 rings, the partial-trace lower bound (14), and the semigroup heat-trace constant used to feed C48. |
| C51. c_L = c_R = n + 1 for a complete compact Gaussian limit with unequal velocities and moving zero modes; thermal coefficient (π/12) Σ (1/u_R + 1/u_L) | proved, conditional on the completeness hypothesis G. Correct and standard; given G the central charge statement is close to definitional. | Read PROOF.md fully. Reran checks.py (identical output). New checks of the canonical diagonalization, the charged-sector Virasoro norm, and the T² coefficient for free fermions with unequal left and right Fermi speeds. |
| C52. No-go for exactly invariant product-cell codes with invertible nearest-cell logical hopping; leakage of the literal C44 marker move | proved. Correct. It is a negative result about one redesign class and does not bear on whether S holds. | Read PROOF.md fully. Reran checks.py (identical output). New Fock-space check of the three-cell identity with Jordan-Wigner signs and quartic terms, the head-leakage formula and monotone ratio, the Gram singular-value bound. |
| C49. Finite prethermal Schrieffer-Wolff dressing with 2^(−m) leakage, and the Schur-complement spectral theorem with error 8Mε_m²/δ | proved. Correct. At fixed coupling it certifies only a finite size window and does not give o(1/M). | Read PROOF.md fully, re-derived (4) to (13), (21) to (26), (35), (36). Reran checks.py (identical output). New random tests of theorem F on 3- and 4-cell chains and of the product-rotation example. |
| Input S as a whole | assumed. Not closed by these runs. | See the verdict. |

Tags used below. proved means I followed the argument and found it correct at the standard of a careful mathematical-physics paper. numerical means supported by a finite computation only. heuristic means plausible but not argued. assumed means an input taken from elsewhere.

## Reproduction

The five work directories were copied to a temporary scratch directory outside the repository, keeping the `codex/work_C*` layout, and each `checks.py` was rerun with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py`. All five pass, and every regenerated `checks.json` is identical to the committed one.

The independent checks are in `referee_checks.py` (numpy and scipy only, about one second). Run `OPENBLAS_NUM_THREADS=1 python3 -B referee_checks.py` in this folder. It writes `referee_checks.json`. None of it reuses code from the Astra check scripts. All assertions pass. These are small toy models that test the inequalities and identities the proofs rely on. They are not simulations of the 4464-state chain.

## C48. Resummed massive matching

### What is claimed

On an enlarged cell space E = C Ω ⊕ QK, the code is mapped to logical states tensor Ω and excited states to logical vacuum tensor the excited vector (eq. (3)). With B commuting with the local all-code projection P_X, the coupling K = J B J† − H_log ⊗ I annihilates the auxiliary vacuum from both sides. A coarse-time polymer expansion, keeping the whole excited cell Hamiltonian as the auxiliary semigroup and expanding only the finite logical blocks, then gives a convergent influence functional with kernels decaying as exp(α · spatial diameter + κ · temporal diameter), κ = ΔL/(16T_n), and ‖W_∞[K + R̃]‖ ≤ A₀ r with r = ‖R‖, uniform in cell number. Finite temperature adds a baseline B_* exp[−βΔ/(4T_n)]. Invalid embedding states cost at least Δ − w and change the pressure by at most β⁻¹ log[1 + C₀ e^(−β(Δ−w))].

### Line-by-line findings

- Embedding (3), (4) [proved]. The valid sector of the enlarged Hamiltonian is unitarily equivalent to the original one, since a_aux J = J A and H_log ⊗ I + K = J B J†. The vacuum-annihilation property of K holds because B_X reduces P_X, so J B J† maps |u⟩ ⊗ Ω_X to (P B P u) ⊗ Ω_X. The text says "(4) vanishes on every |u⟩ ⊗ Ω_X" right after stating that hypothesis, so this is correct, but it would read better with the dependence on [B_X, P_X] = 0 made explicit (cosmetic).
- Coefficient conversion (5) to (7) [proved]. 4^(mr) matrix units, at most 2^(mr) terms from the vacuum projector, degree at most 2mr, so (8h²)^(mr). With h = 4 this is 128^(mr), which is why a = L + m log 128 + 1. The fermionic sign structure is delegated to the "graded tensor product" and to a hypothesis that a local fermionic representation exists. C58 supplies that through the collaborator's SEPARATE_JW note.
- Block bound and activities (12) to (15) [proved]. The Schatten-Hölder step |Tr Π e^(−s_i A/2) O_i| ≤ Tr e^(−βA/2) Π‖O_i‖ is correct with exponents β/s_i. The passive cost e^(−T/2) ≤ e^(−L) follows from T ≥ 2L, and the trace factor Tr_E e^(−τ a) ≤ 1 + C₀ e^(−T) < e from T ≥ log C₀ + 2. In (13), e^x − 1 ≤ x e^x with x = 2τ v_X ≤ 1/8 gives η ≤ e^(1/8)/8. The animal bound 16^N, the weight bound exp(L/4) exp[(α + L/8)N], and q = 16 exp(−7L/8 + α + 6) all check. With L = 8(α + log 32 + 10), q is of order 10^(−40), so the Kotecký-Preiss-type criterion holds with a very large margin (C48's own script gives 9.3e-40 for the KP ratio). The rooted-tree and Penrose-type tree bounds are standard and correctly stated.
- Cancellation W_∞[K] = 0 and the leakage-linear bound (16) [proved]. At zero temperature every K-only vacuum word starts on Ω and vanishes. At finite temperature a K-only polymer must meet every time layer, which gives the B_* e^(−β/(4τ)) term. The Cauchy estimate along K + sR̃ gives (16). This cancellation is the load-bearing idea of C48. It is what makes the error linear in the C45 remainder rather than in the order-λ coupling, and it is correct.
- Moments (17), (18) [proved]. I re-derived 2^p(1 + Δ Σ|t_i − t_a|)² ≤ [2 + 16(Δ/κ)²] 4^p e^(κd) using p² ≤ 4 · 2^p and x² ≤ 2e^x. The frequency estimate is also correct.
- C46 resolution, Section 6 [proved]. Moving the internal energy u into the auxiliary semigroup removes the divergent u/Δ series. The operator identity (31) is correct. My check `C48_C46_resummed_formula` confirms that the four lowest levels of the one-cell problem agree with −λ²/Δ + λ²u/[Δ(Δ+u)] n₁n₂ up to λ⁴/Δ³, for u/Δ = 1, 3 and 20. C48 is careful to say that this does not make C46's raw action converge; the source definitions change.
- Low-energy and pressure comparison (32), (33) [proved]. Sign and constants check. My check `C48_enlarged_embedding_multicell` builds the original and enlarged Hamiltonians explicitly on a ring of three cells with random block-diagonal B, generic R and u ∈ {0, 1, 5}Δ. The spectra below E₀ + Δ − w coincide to 1e-15, and the pressure difference lies in [0, bound (33)] at β = 2, 5, 12. C48's own check is single-cell only.
- C42/C44 substitution, Section 5 [proved conditional on imports]. The heat trace (21) is correct given C44's legal wrong-length cost, and Δt₀ ≤ log(2d₀) follows from Δ ≤ γ. The parameter choices (26), (29) put λ last, so there is no circularity at fixed n. Constants are astronomically large (A_dec above e^120 and input exponent a above 108), which is irrelevant at fixed n but should be said whenever C48 is described as "explicit".

Gaps.

1. Fixable. The passage from the operator trace to a Grassmann influence functional with split insertion times ("continuous fields and kernels") is asserted, not proved. The bounds needed are uniform in the split, so I expect it to go through by a standard time-lattice argument, but it is a sketch.
2. Fixable. Several constant reserves are bundled rather than tracked ("a factor 64", "enlarging the numerical factor within B_* if needed" for periodization). Nothing depends on the precise value, only on finiteness at fixed n and independence from the wrapper ratios, and that independence holds.
3. Serious, acknowledged by C48 and handled in C58. The hypotheses (local fermionic interval representation preserving the reference, all-string heat trace, C45 block-diagonal form at exponent a) are not proved in C48.

My independent u-uniformity check (`C48_leakage_linear_uniform_in_u`, three cells, internal massive interaction up to u = 100Δ) finds that the low-spectrum deviation from H_log = P B P is zero to 1e-15 when R = 0. With R it is linear in r, with err/r between 0.42 and 0.44 for every u. That is consistent with C48's leakage-linear, u-independent claim, though it is a finite-size spectral proxy and not a test of the kernel norms.

Verdict on C48. A careful mathematical physicist should accept C48 as a correct large-gap polymer expansion for the enlarged rooted model. The result is the S(a) estimate for that model, at fixed n, with constants independent of the wrapper ratios. It is less than S(a) for the physical chain until the hypotheses are discharged (C58).

## C58. Transfer to the physical ring

### What is claimed

At fixed n and small fixed couplings, sup over |s| ≤ ε of |p_o(T,s) − p_m(T,s)| ≤ C T e^(−δ/T), where p_o is the original rooted pressure and p_m the pressure of the exact matched (enlarged) thermal action. Physical free energies with marker number summed, or at fixed interior marker density, then agree to the same activated accuracy, with the cyclic identification and periodic boundary conditions included. Under (a) and (b), f(T) = e_phys − (π/12) Σ (1/u_R + 1/u_L) T² + o(T²).

### Line-by-line findings

- Fixed-marker pinning (4) to (7) [proved]. The compression identity (5) is correct: Σ_i [p_i, [p_i, H]] = 2MH − 2Σ p_i H p_i on the M-marker sector. In fact ‖[p,[p,h]]‖ ≤ ‖h‖, so the constant 4 is generous. The Jensen trace inequality for the isometry V and the min-max comparison give the upper bound. My check `C58_fixed_marker_pinning` tests both inequalities on spin-1 rings with a random complex marker-conserving bond, L = 3 to 6, every M, three temperatures (162 cases, all pass). This is a different model from C58's own XXZ-type test.
- Completion removal, Section 3 [proved, conditional on imported C44 lower bounds]. The mismatch flags are invariant under the literal selector rewrites. C58's script checks 6400 rewrites, and I followed the case analysis by hand. The fragment bound A_bad ≥ U + |e|ℓ follows from K_sel = 4t_sel + 9U, since K_sel z − (2t_sel + 4U)(z+1) = U + (2t_sel + 5U)(z−1). I checked the same arithmetic in the collaborator's ILLEGAL_SECTORS.md (6). The heat bound (11) and the lower bound (14) are correct. (14) is the Peierls-Bogoliubov-type partial-trace inequality with a unital trace-preserving conditional expectation. It uses that the selected code has equal rank at lengths b and b+1, so the centered length averages to zero. My check `C58_partial_trace_lower_bound` tests (14) on 60 random instances. The combination (15) to (17) with δ_def ≥ U/4 is correct. The completion comparison is the one new idea in C58 that a reader should check against C44's rule table. It depends on the claim that the conditional pins are the only strong intercell terms and that a legal predecessor always enables them. I could not verify that against C44 here.
- Heat trace feeding C48, Sections 4 and 5 [proved]. At first sight C₀ = max(1, e^(t₀Δ)/(2d)) looks wrong, because the generic all-length bound decays only at rate κ = Δ/(b+2) < Δ/2. It is right by the semigroup factorization Tr_Q e^(−tA/2) ≤ ‖Q e^(−(t−2t₀)A/2)‖ Tr_Q e^(−t₀A) ≤ e^(−(t−2t₀)Δ/2)/(2d). My check `C58_semigroup_heat_constant` confirms the bound on synthetic spectra saturating A ≥ max(Δ, κ(ℓ+1)) with d^ℓ states at length ℓ. The tilt window (21) and the invalid-state bound (24) follow from TILTED_EMBEDDING.md Sections 3 and 4, which I re-derived, including A + s(q + N_len − 1/2) ≥ A/4 on invalid states.
- Root placement and monotonicity (25), (26) [proved]. Every cell and every enlarged state has length at least one, so p(s+t) ≤ p(s) − t. That gives |s_o − s_m| ≤ d_T directly.
- Canonical extraction, Section 8 [proved at the level of a sketch]. The blocking, Chernoff and multinomial argument is the standard proof of equivalence of ensembles for a one-dimensional system with finite cut strength. The saddle bracket using the strict curvature p_G,0'' = 1/[2πW j_len sin(πν)] and a uniform pressure error E_T is correct, given the order of choices (ε_stat, then λ, then T).
- Final conditional theorem (30) to (32) [proved, conditional on (a), (b)]. The expansion of p_o(T, s_o) = 0 around s_* is right and does not differentiate the remainder. The conversion u = ℓ v to physical speeds is right.

Gaps.

1. Serious, a dependency rather than an error. C58 rests on ILLEGAL_SECTORS, PINNED_MARKER, LENGTH_FUGACITY, TILTED_EMBEDDING, SEPARATE_JW and TILTED_MATCHING (collaborator notes, AI-audited only), on C44's local lower bounds and pin structure, and on C45's SW identity extended to the all-string completed space. I re-derived PINNED_MARKER Section 1, TILTED_EMBEDDING Sections 3 and 4, and the ILLEGAL_SECTORS fragment arithmetic. I did not re-derive SEPARATE_JW (the claim that the Jordan-Wigner reference extension keeps K_D vacuum-annihilating and the retained hoppings quadratic), LENGTH_FUGACITY Section 2, or TILTED_MATCHING. These carry the C48 hypotheses and need a human check.
2. Fixable, in the statement. REPORT.md 9.27 says that under (a) and (b) "C51 then gives c = n + 1". C51 needs the completeness hypothesis G for the Virasoro central charge. With (a) and (b) as C58 states them, one gets the thermal coefficient with n + 1 dressed modes and their speeds. That determines c in the operational thermal sense (each identified mode contributes c = 1 and the coefficient leaves no room for an extra sector with finite speed). It does not by itself construct the stress tensor. C58's own text is careful about this. The summary should name which c is meant.
3. Fixable. The canonical-ensemble lemma is a condensed sketch of a standard argument.
4. Cosmetic. The seam constant B = 2d⁴J is quoted without derivation. Only its finiteness matters.

Verdict on C58. The new transfer arguments are correct. A careful reader should accept that, given C48 and the imported finite and collaborator lemmas, the physical fixed-length free energy equals the matched pressure up to activated errors, so nothing in the ensemble, illegal-cell or cyclic-gluing steps can change a T² coefficient. This is a free-energy result only. It gives no o(1/L) spectrum and no tower multiplicities.

## C51. c from a complete Gaussian limit

Findings [proved].

- Canonical diagonalization (4), (5). R = V^(−1/2) O U^(1/2) gives RᵀVR = U and RᵀSR = S, and the signed speeds are the eigenvalues of SV. My check `C51_canonical_diagonalization` confirms this for 80 random positive V with N = 1 to 4, to 4e-14.
- Charged-sector Virasoro algebra (6) to (10). ‖L_(−m)|x⟩‖² = m x² + (m³ − m)/12 is correct, so c = 1 per chiral current in every charge sector. I recomputed the contraction count independently.
- Thermal coefficient (2), (16). ∫₀^∞ −log(1 − e^(−x)) dx = π²/6 gives π/(12βu) per chirality. My check `C51_asymmetric_free_fermion_T2` computes the grand potential of free fermions with dispersion −2cos k + a sin 2k, whose left and right Fermi speeds differ (for example 2.65 and 1.33). The measured coefficient of −T² agrees with (π/12) Σ 1/|v_F| to between 3e-7 and 1.5e-5 relative error at T = 0.005. This confirms the unequal-speed form of (2).
- Moving zero modes, theta bounds, scalar summation (Sections 4, 6, 7). Correct. The point that a raw Casimir fit need not converge at irrational density, and that the zero-mode sum contributes O(log L) to log Z, is right and useful for any numerical test of the construction.

Gaps. None in the mathematics. One observation for the reader. Given G, which assumes the infrared theory is the complete theory of N compact bosons with no further massless factor, c = N is close to definitional. The substantive content of C51 is the bookkeeping (unequal chiral speeds, moving offsets, the summation criterion P with p > 2). The spectator example and the degeneracy example (one boson of speed 1 versus two of speed 2) are correct and show why a thermal coefficient alone does not fix c without the velocities. Whether a multi-velocity Gaussian theory counts as "a collection of unitary CFTs" for the conjecture is the scope question in the expert note, and C51 does not, and cannot, settle it.

Verdict on C51. Accept, conditional on G (or on P plus a mode and velocity identification for the thermal version).

## C52. No-go for exact product-cell codes

Findings [proved].

- The commuting-bond theorem (Section 2) is C8's, correctly reproduced. The key step applies a traceless functional on cell j+1 to the leakage equation, so the right boundary lifts reduce P_j, and lifts with disjoint support commute after compression.
- The three-cell identity (5). ⟨1a|[K₁₂, K₂₃]|3b⟩ = (T₁T₂)_ab. My check `C52_three_cell_commutator` builds the Jordan-Wigner Fock space of three cells with n = 1, 2 flavours, random complex hoppings, random onsite terms, density-density quartics and a pair-hopping quartic, and finds agreement to 4e-16.
- Invertibility of C39's hopping (7). σ_min(−W(Γ₁ + iηΓ₂)) ≥ W λ_min(Γ₁) for positive Γ₁ and Hermitian Γ₂. Confirmed on 200 random instances (`C52_gram_singular_value`).
- Head leakage (8) and the monotone ratio (10). Confirmed for b from 4 to 200 (`C52_head_leakage`). The b = 4 value 0.1701573216 matches C52's report.
- Rank bound (12). m₁ ≤ ⌊d²/4⌋ = 4,981,824 at d = 4464, correct.

Gaps. None. Its scope is narrow, as C52 says. It excludes exact invariance with separated endpoints and finite-range exact dressing of range below (b − 2)/4. It does not exclude a dressing of range comparable to n, and it does not bear on whether S holds for the perturbative construction. Its use is to show that S cannot be avoided by a cheap redesign.

Verdict on C52. Accept as stated.

## C49. Prethermal dressing and the Schur-complement theorem

Findings [proved].

- Local inverse and weighted commutator bounds (4) to (6). I re-derived ‖[C,D]‖_(a−h) ≤ (4/h)‖C‖_a‖D‖_a and the BCH majorant.
- Iteration (10), (11). With h = 64eΛ/Δ one gets ζ_j = o_j/(8Λ) ≤ 1/8 and ‖T_j‖ ≤ o_j/7, so the off-diagonal part halves at each step while d_j + o_j ≤ Λ, and m = ⌊aΔ/(128eΛ)⌋ steps fit in half the locality weight. The coefficient k/(k+1)! is correct. This is a clean, conservative version of the Abanin, De Roeck, Ho and Huveneers prethermal scheme, specialized to inverting only code-to-complement transitions.
- Physical support of the generators (15). Correct and a nice observation. A P_X on one side plus length conservation confines each generator to physical length at most (b+1)|X|.
- Theorem F (16) to (26). The orthogonality P O_X O_Y P = 0 for disjoint supports gives ‖QOP‖² ≤ M ε'². The completely positive extension φ_j(y) = P_j y P_j + tr(ρy) Q_j gives (21). The Young-inequality bound (23) and the Feshbach count (25) are correct. My check `C49_theorem_F` tests (18) and (26) on random 3- and 4-cell chains (cell dimension 4, code dimension 2) with leakage up to 0.16. All levels in the window satisfy k_i − t_M ≤ h_i ≤ k_i. The bound is loose by one to two orders of magnitude in these examples.
- The product-rotation example (35), (36). Confirmed for M up to 8 (`C49_product_rotation`). The bare QHQ threshold falls toward the ground energy with M while the true gap stays √(Δ² + 4t²).

Gaps. Cosmetic only. In (9) the sum is bounded by Σ o_j/2 where the proof gives o_j/7. That is harmless.

Verdict on C49. Accept. Its spectral theorem certifies a size window of order λ⁻² exp[c/λ] at fixed n, not an o(1/M) statement at fixed coupling. C49 says this plainly and correctly identifies (31) as the missing state-dependent estimate.

## Overall verdict

A careful mathematical physicist should accept the following.

1. C48 proves the massive-elimination half of S, S(a), for the enlarged rooted cell model at fixed n, uniformly in volume and at all sufficiently low temperature, with constants independent of the wrapper ratios. The decisive point, that block-diagonality makes the influence functional linear in the C45 remainder, is correct.
2. C58 correctly transfers this to the physical fixed-length ring at the level of the free energy, with activated errors. The transfer is conditional on collaborator lemmas that have only had AI review, some of which I did not re-derive.
3. C51 correctly identifies c = n + 1 once a complete multichannel Gaussian limit is available, including unequal chiral speeds and moving zero modes.
4. C52 and C49 are correct, limited results. C52 rules out one exact redesign. C49 gives a finite prethermal window and isolates the spectral estimate that would be needed for o(1/M).

What remains open, so that S is still assumed.

- (a) The all-scale constructive RG for the matched retarded action of arbitrary degree, with counterterms and length-source control. This is Lemma N with the corrected source budget (C55, C56, C60, C61), and it is not closed.
- (b) The temperature-sensitive scalar estimate (P in C51, TP in C59 and C62), which delivers the T² coefficient. Not closed.
- For the Virasoro central charge rather than the thermal one, completeness of the Gaussian limit (C51's G), which neither (a) nor (b) as stated supplies.
- Human review of the imported finite construction (C39, C42, C44, C45) and of the collaborator lemmas that C58 imports.
- The scope question in the expert note (whether a multi-velocity Gaussian theory with moving zero modes is "a collection of unitary CFTs" in the sense of the conjecture).

No fatal gap was found in any of the five runs. The serious gaps are dependencies, not errors. The fixable gaps are presentation (the Grassmann time-split limit in C48, the condensed ensemble lemma in C58, and the wording of which c the conditional theorem delivers). In short, these runs move the massive elimination and the ensemble transfer from assumed to proved-conditional, and leave S resting on the two constructive-RG estimates for the logical fermion action. The counterexample remains conditional.

## Files

- `README.md`, this report.
- `referee_checks.py`, independent numerical checks (numpy and scipy).
- `referee_checks.json`, its output.
