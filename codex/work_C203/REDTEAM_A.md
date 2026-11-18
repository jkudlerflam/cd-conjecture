# C203 red team, report A: Theorem A of C201

2026-11-18. Referee: Claude Opus 5.5 (red-team subagent). No Astra call was made for this report. Every verdict below rests on my own line-by-line reading of `codex/work_C201/PROOF.md` §2, `astra_A4_referee_T1.out` §1 to §4, and the text of arXiv:1509.07127, which I extracted from the arXiv PDF and read. The numerical probe is my own script, `codex/work_C203/redteam_marginal_probe.py`. No human review.

Line references are to `codex/work_C201/PROOF.md` unless stated.

## 1. The citation

Verdict: correct. The inequality used is the one the paper proves, with one improvement to the attribution.

What the paper says. Corollary 4.1 of Junge, Renner, Sutter, Wilde and Winter (eq. (101) of the arXiv version) states, for any density operator ρ_ABC on finite-dimensional A ⊗ B ⊗ C, I(A:C|B)_ρ ≥ −2 log F(ρ_ABC, R_{B→BC}(ρ_AB)), with R = ∫ dt β₀(t) R^{t/2}_{ρ_BC, tr_C} the β₀-average of rotated Petz maps (Theorem 2.1, eqs. (18), (19), and Remark 2.2, eq. (21); the corollary's display drops the 1/2 in the rotation parameter, a typo that does not matter). The fidelity is F(ρ, σ) = ‖√ρ √σ‖₁, unsquared (eq. (8)). The same logarithm appears in the entropies and in −2 log F, so the inequality holds in nats as used. There is no full-support hypothesis: the corollary is Remark 2.2 applied with σ = 1_A ⊗ ρ_BC, and every ρ_ABC lies in S_σ automatically. The map is completely positive and trace non-increasing, trace preserving on the support of ρ_B, and Remark 2.2 says it may be taken to vanish off that support.

The exact-recovery property 𝓡(τ_B) = τ_BC, which PROOF §2.6 and §2.7 need for "the equality controls every bond inside the re-created interval" (line 98), is not in Corollary 4.1 itself. It is Remark 2.4, item 2 ("perfect reconstruction of σ from N(σ)"), applied to σ = 1_A ⊗ ρ_BC and N = id_A ⊗ tr_C, which gives 𝓡(1_A ⊗ ρ_B) = 1_A ⊗ ρ_BC and hence 𝓡(ρ_B) = ρ_BC. PROOF line 44 should cite Corollary 4.1 together with Remark 2.4(2). The completion off the support (prepare a fixed state on the complement of supp τ_B) is CPTP and does not change the action on any operator whose B-marginal is supported in supp τ_B, which is the case at both places it is used, since 𝒟(η_B) = τ_B exactly and τ_AB has marginal τ_B.

The conversion in (2.3) is right: Fuchs–van de Graaf gives ‖α − β‖₁ ≤ 2√(1 − F²), and 1 − F² ≤ −ln F² ≤ I in nats.

Universality of the map is never used. The proof applies the recovery theorem to the single state τ_ABC, so any rotated-Petz result with exact recovery of τ_BC would do. The dependence on JRSWW is therefore milder than RESULT suggests.

## 2. The construction and its error accounting

Verdict: correct. I found no gap.

Scale selection (§2.3, lines 54 to 58). SSA on two overlapping translates gives u_{j+1} ≤ u_j for j + 2 ≤ L. Weak monotonicity with two disjoint j-site intervals flanking one site gives u_j ≥ 0 for 2j + 1 ≤ L, which covers j ≤ B_{n+1} − 1 because L ≥ 2B_{n+1}. The drops telescope to at most u₀ = S₁ ≤ ln d. In (2.1), S(A|B) ≤ S(A_ℓ|B_ℓ) + S(A_r|B_r) uses subadditivity of conditional entropy and SSA, and each term is S_{b+1} − S_b = u_b by translation invariance (no reflection symmetry needed, since B_r is to the left of A_r and the pair is still a (b+1)-interval). S(A|BC) = S_m − S_{m−2} because removing the two exterior sites leaves a contiguous interval. Both bounds need m ≤ B′ ≤ L/2, which holds. Line 58 is correct.

Measurement (§2.4, lines 62 to 66). The 2-design identity and the dual X_y = (q + 1)P_y − 1 reconstruct every operator, ‖X_y‖₁ = 2q − 1, and the trace-norm amplification bound follows from the triangle inequality on Σ_y Z_y ⊗ X_y. Correct.

Decoupling (§2.5, lines 70 to 78). The chain rule and I(U:Y) ≤ S(U) for classical Y are correct. Pinsker in nats per record value, then Jensen, then reconstruction, gives (2.2). The averaged target Σ_v p_v ω^v_{UW_j} = ρ_{UW_j} holds because the record is taken on disjoint windows. With K = ⌈2(2q − 1)² δ⁻² ln q⌉ the right side is at most δ. This is the Brandão–Harrow local-measurement chain-rule step, with the IC-measurement cost made explicit.

Marginal correction (§2.7, line 92). I re-derived the lemma. T ≥ s gives supp s ⊆ supp T, so K T K† = s, K†K ≤ Π_T, and K r K† ≤ s because r ≤ T. The second branch has weight tr[(1 − K†K) r] ≤ tr[(1 − K†K) T] = z. The Hilbert–Schmidt estimate ‖(1 − K)√r‖₂² ≤ ‖√T − √s‖₂² ≤ tr T − tr s = z uses r ≤ T and tr(√T √s) ≥ tr s, and the factorization through √ξ gives 2√z. The bound 2√z + z is correct, including singular r and s and a vanishing residual.

The join bound (2.4), line 96: with z ≤ δ/2 the three terms are √(2δ) + δ/2, δ, and 2√(2 ln d/n), which sum to f/2. Correct.

Geometry and energy count (§2.8, lines 102 to 116). In ring order the interval is A_ℓ B_ℓ C B_r A_r, and the retained arc runs B_r A_r ⋯ A_ℓ B_ℓ, so concatenating arc i, a fresh C, and arc i + 1 reproduces the order A_ℓ B_ℓ C B_r A_r at every join. For k ≥ 2 the two windows at a join come from different copies, so their marginal is η_AB = Σ_v p_v ω^v_U ⊗ ω^v_{W_j}. (For k = 1 they would come from the same copy, and the construction would fail, which is why k ≥ 2 is required.) The join channels act on disjoint sites and are trace preserving, so each join's marginal is (id ⊗ 𝒥)(η_AB) and each internal arc bond keeps marginal ρ₂. The count N − 2b − 1 + (m − 3) + 2 = L bonds per unit is right, and so is the site count N + (j − 1)t = L. The case j = 1 (empty C, m = 2t) is covered: the two buffers become adjacent and the corrected B-marginal is the 2b-site interval state.

From (2.4) to energy. In the gauge where ‖h‖ = M, the gauge and constant shifts cancel on every ring state by telescoping, and they cancel in tr(hρ₂) up to the same constant because ρ₂ has equal left and right marginals. So the trial-state energy is at most kL e + k f M. Line 114 is correct.

From the ring inequality to (A1). E_{kL}/(kL) → e₀ as k → ∞, so e₀ ≤ tr(hρ₂) + fM/L for every translation-invariant ρ on the L-ring. The trace-norm form also follows without energies: taking k = ∞ copies and averaging over L translations gives a translation-invariant state on ℤ whose marginal is ((L − 2)ρ₂ + σ_ℓ + σ_r)/L with ‖σ_{ℓ,r} − ρ₂‖₁ ≤ f/2, at trace distance at most f/L from ρ₂. The minimax identity of `astra_A4_referee_T1.out` §4 is then not needed, though it is also correct (and uses M(h) ≤ ‖h‖).

One remark on (2.4). The connecting bonds are controlled only in the unconditional output state, not conditionally on the record. That suffices, because the energy is linear in the state. `astra_A4_referee_T1.out` makes the same point, and I agree.

## 3. The doubling and k-fold statement

Verdict: correct. No assumption about the infinite chain enters.

The trial state on the kL-ring is the one built in §2.8 from k conditionally independent copies of the retained arc of the translation-averaged L-ring ground state, with the shared record and the k join channels. It is a state on exactly kL sites (the second Astra proof makes the join length deterministic, which the first random-length proof did not, as `astra_A4_referee_T1.out` §4 notes). The reference energy is tr(hρ₂) = E_L/L of the L-ring itself. Nowhere is E_{kL} ≤ kL e₀ + c used, and no compatibility of L with the order of the infinite chain is needed.

Incommensurate lengths. For a classical bond, cutting a ring configuration and gluing k copies is exact, so E_{kL} ≤ k E_L, and the theorem adds nothing. For h = Z ⊗ Z on odd L, E_{2L} − 2E_L = −2L − 2(−L + 2) = −4 ≤ 2fM. For a period-3 bond and L ≢ 0 mod 3, E_L carries the incompatibility cost, and the inequality only bounds E_{kL} from above in terms of that already-raised E_L. Frustration makes (A2) easier, not harder. The statement is one-sided: E_{kL} − kE_L can be −Θ(M) (the Z ⊗ Z example), so "the doubling defect vanishes uniformly" (RESULT, "Consequences") is accurate only for its positive part. The formula in RESULT, sup_h (E_{2L} − 2E_L)/M → 0, is correct as written, since the supremum is nonnegative (a frustration-free bond gives zero).

## 4. The use of M

Verdict: correct. The constant is right.

The proof never needs the "equal one-site marginals" form of the dual-norm statement in `_common.md`. It fixes an optimal gauge representative h̃ with spectrum in [−M, M], so ‖h̃‖ = M, uses |tr h̃ (σ − ρ₂)| ≤ ‖h̃‖ ‖σ − ρ₂‖₁ ≤ M ‖σ − ρ₂‖₁ at the two connecting bonds, and then returns to h. That return is exact: Σ_j (X_j − X_{j+1}) vanishes on every ring, and tr[(X ⊗ 1 − 1 ⊗ X)ρ₂] = 0 because ρ₂ is the marginal of a translation-invariant state. The minimum defining M is attained (the spread is continuous and coercive in X modulo multiples of the identity), so the gauge choice exists. The bound in the trace-norm form uses M(h) ≤ ‖h‖, which holds with X = 0.

## 5. Counterexample search

Verdict: no counterexample, and none is possible at reachable sizes, because Theorem A is vacuous there.

Analytic checks of the listed families.
- GHZ ring states: the two-site marginal ½(|00⟩⟨00| + |11⟩⟨11|) is the marginal of the GHZ mixture on ℤ. Distance zero.
- Singlet coverings on even L: the translation average of the dimer covering is the marginal of the translation-averaged dimer state on ℤ. Distance zero.
- W and Dicke states: permutation invariant, so within 30d⁵/L² of Ω_∞ by Proposition G (PROOF §7.1), whose proof I checked. Explicitly for W, tiling ℤ with open W-blocks of length L gives a marginal at trace distance O(1/L²).
- Half-filled hopping ring: distance Θ(L⁻²), computed exactly in `astra_A4_referee_T1.out` §3 (x_L = 1/(L sin(π/L))). I rechecked the expansion.
- SU(N) swap chain at L = N (the antisymmetric singlet, the most "ring-adapted" state there is): E_N = −N, and with Sutherland's e₀(N) = 1 − (2/N)[ψ(1) − ψ(1/N)] one gets N e₀ − E_N = 2ζ(2)/N + O(N⁻²) → 0 (checked numerically: 0.307 at N = 10, 0.0327 at N = 100) with M = 1. This family has d = L, outside the theorem, and it still decays. At fixed d = 3 the RESULT checks give L(Le₀ − E_L) = 2.67, 2.32, 2.25 at L = 3, 6, 9 (I reproduced 2.67 at L = 3 from e₀(3) = −0.7032), the CFT 1/L law.

Numerical probe (`redteam_marginal_probe.py`, Claude Opus 5.5). For d = 2 and L = 4, 5, 6 I lower-bounded L·dist₁(ρ₂, Ω_∞) by the trace distance to the outer relaxation Ω^{(r)} (two-site marginals of r-site states with equal left and right (r − 1)-site marginals, r = L + 2), which contains Ω_∞. States: GHZ, W, Dicke(⌊L/2⌋), the dimer average, ring ground states of Heisenberg and XX, and four rounds of an alternating ascent in which the optimal witness of the distance SDP is fed back as the next bond. RESULTS_PLACEHOLDER

Why no test can bite. At δ = 1 and n = 1 the threshold is already L ≥ 2B₂ = 548 for d = 2, and there f = 3 + 2√2 + 4√(2 ln 2) ≈ 10.5, five times the trivial bound 2. The bound beats the trivial one only when 4√(2 ln d/n) + 3δ + 2√(2δ) < 2. At d = 2 this needs n ≥ 6 with δ ≲ 5·10⁻⁴ (then B₂ ≈ 10⁹), or n ≥ 8 with δ = 10⁻² (then B₂ ≈ 3·10⁶). Either way L ≥ 2B₇ or more, five or more iterated exponentiations above B₂. Theorem A has no content at any length a computer or a physical system reaches. It can only be refuted by finding an error in the proof.

## 6. Literature

Verdict: I found no prior statement of Theorem A. The ingredients are known.

- Fawzi, Fawzi and Scalet, "Entropy constraints for ground energy optimization", J. Math. Phys. 65, 032201 (2024), arXiv:2305.06855. They study translation-invariant marginal relaxations with entropy constraints and state convergence rates O(1/m) for both the relaxation and the finite-ring limit (their Appendix B builds a translation-invariant infinite state from a finite one with error decaying like 1/m). This is the trivial order, not o(1/m).
- Blakaj and Wolf, "On the set of reduced states of translation invariant, infinite quantum systems", Lett. Math. Phys. 114, 28 (2024), arXiv:2308.14585. They show the set of two-body reduced states of translation-invariant chains (Ω_∞ here) is not semialgebraic and approximate it from inside by MPS and from outside by finite-system marginals. No rate in system size.
- Wang and Navascués and coworkers (Kull, Schuch, Dive, Navascués, PRX 14, 021008 (2024), arXiv:2212.03014) relate ring energies to locally translation-invariant relaxations, at O(1/n).
- Zeng, Yang, Zhang and Wang, arXiv:2608.10365 (August 2026), on block-product periodic marginals. From the abstract, no ring-to-chain rate.
- Techniques: the universal recovery theorem (JRSWW, and the Sutter–Fawzi–Renner precursor arXiv:1504.07251), and the local-measurement chain-rule decoupling of Brandão and Harrow (arXiv:1210.6367; also their product-state approximation paper, arXiv:1310.0017). Gluing finite TI pieces with recovery maps is the Markov-entropy-decomposition circle of ideas (Poulin and Hastings) as used for Gibbs states by Fawzi, Fawzi and Scalet (arXiv:2209.14989).

So, as far as my search reached, the o(1/L) ring-to-chain inclusion with a d-only rate is new as a statement, and its proof is a careful but standard assembly of known tools. I would not call it known. I would also not call it deep, and a referee in quantum information would likely see it the same way. My search was web-based and not exhaustive. In particular I did not check the Carlen–Lebowitz–Lieb extension-problem papers or the Aharonov–Irani and Watson–Cubitt complexity papers in full text.

## 7. Consequences as written in RESULT.md

The theorem is correct, but RESULT overstates three of its consequences.

(i) The "third law". PROOF Corollary A1 (line 128) tags lim s̃_L = 0 as rigorous given theorem R under the positive-transfer reading, but the step from limsup(Le₀ − E_L) ≤ 0 to lim s̃_L = 0 is not written, and it is false without positivity. For h = Z ⊗ Z, s̃_L = L E_{L+1} − (L + 1) E_L equals 2L on even L and −2L − 2 on odd L. Under positivity the step is true, by an argument that should be added: (O2) gives concavity and E_{2L} ≥ 2E_L, hence g_L = Le₀ − E_L is nonincreasing and nonnegative, Theorem A forces g_L → 0, s̃_L is nonincreasing with s̃_L = g_L + L(g_L − g_{L+1}) ≥ g_L ≥ 0, and a positive limit of s̃_L would make Σ(g_L − g_{L+1}) diverge. In the crossed-channel picture the content is that the crossed-channel ground state has no residual degeneracy, which is a fair use of "third law". RESULT's phrase "no chain has a ring energy a fixed amount below L e₀ (a third law in the crossed channel)" is fine for the first half. The label should carry the positivity proviso. Verdict: correct after the stated repair.

(ii) "An extensive fast tower is excluded." In PROOF this is Corollary A3, tagged physics-level, with an explicit hypothesis that the remainder r_L is uniformly small across the family. RESULT's summary drops both the tag and the hypothesis. That hypothesis does most of the work, and without it the deficit family MN/(L + N²) of PROOF §2.10 shows that Theorem A says nothing about the coefficient. The rate is o(L) only in the sense of an inverse tower function, so at every length below 2B₇ (§5) the statement is consistent with an extensive tower. Verdict: overclaimed in RESULT as worded. The careful version in PROOF §2.10 is acceptable.

(iii) "Showing that this method cannot shorten the scales to polynomial ones." The d = 4 antisymmetric construction (`astra_A4_referee_T1.out` §2, which I read but did not verify independently) defeats full-window trace-norm decoupling (2.2) on spans polynomial in t. It does not exclude a join that decouples only the end sites, controls only the observables the two connecting bonds pull back to, or reconnects rather than decouples. Astra itself says the example gives no lower bound on the energy defect. "This method" should read "the full-window decoupling step (2.2)". Verdict: mild overclaim.

Smaller points. "The doubling defect vanishes uniformly" should say its positive part does (§3 above). PROOF line 44 should cite Remark 2.4(2) for exact recovery (§1 above), and can now drop "I did not open the paper": I did, and Corollary 4.1 is the right number.

## Overall verdict

Theorem A as stated, both (A1) and (A2), is correct given JRSWW Corollary 4.1 together with Remark 2.4(2), which I confirmed from the paper's text. I attacked the citation, the scale selection, the decoupling, the marginal correction, the join geometry, the bond count, the gauge handling of M, the k-fold length bookkeeping and the incommensurate-length cases, and found no gap. The constants and thresholds are as stated. The theorem is weaker than the trivial bound at every length below 2B₇ at d = 2, an iterated-exponential scale, so no numerical test can bear on it. It appears to be new, assembled from standard tools. Of RESULT's consequences, the "third law" needs the positive-transfer proviso and a short omitted argument, the exclusion of an extensive fast tower is physics-level and rests on a uniformity hypothesis that RESULT's summary omits, and the scale obstruction applies to the full-window decoupling step, not to the method in general.
