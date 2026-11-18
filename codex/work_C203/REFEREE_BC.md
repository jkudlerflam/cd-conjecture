# C203 referee report on C201 Theorems B and C (and D, E, F, G)

2026-11-18. Referee: Claude Opus 5.5 (independent agent, adversarial brief). Every numerical statement below comes from ED code I wrote for this report (not C201's scripts), kept in the session scratchpad as `ref_common.py`, `ref_exchange.py`, `ref_potts.py`, `ref_potts_incl.py`. No Astra call was made; every check here was done by me. No human review.

Files read: `work_C201/RESULT.md`, `work_C201/PROOF.md` (§§3 to 7), `astra_A3_adversary.out` §1, the parts of `astra_A1_potts.out` and `astra_A5_referee_potts.out` that state the twist list, REPORT §§9.68, 9.94, 9.100, 9.104, `work_C142/PROOF.md` Theorem I.9 and its corollary, `work_C136/PROOF.md` §1 and (6.13).

## Summary of verdicts

| Item | Verdict |
|---|---|
| B1 domination, e₀ equality, M ≥ t | correct |
| "internal register enters only as a boundary twist" (RESULT line 21, PROOF §3.8) | correct, but the proof is not in PROOF.md (it is in `astra_A3_adversary.out` §1); add it |
| B2 sharp V = 0 bound | correct, sharp, reproduced |
| B3 V ≤ −2t | correct |
| B4 0 ≤ V/2t ≤ 1 | correct after a stated repair: the cited X1 (C142 Theorem I.9) covers 0 ≤ Δ < 1 only, so Δ = 1 needs a continuity sentence |
| B5 −1 + δ ≤ V/2t ≤ 0 | correct given X2 (C136 (6.13)), which covers the flux-minimized sectors used |
| C1 inclusion | correct given theorem M, reproduced numerically |
| C2 ground-state identity, twist 2γ, Sᶻ = 0 | correct, reproduced to L = 9, 8, 7 for Q = 2, 3, 4 |
| C3 e₀ = 2ε₀ | correct, reproduced |
| C4 M = (√Q + 1)/2 | correct for Q = 3, 4 as PROOF states. RESULT line 23 overclaims it for Q = 2, where M = 1 |
| C5 L(L e₀ − E_L) ≤ 841 | correct after the same Δ = 1 repair (needed at Q = 4) |
| Novelty of C1 to C3 | not new: Alcaraz, Barber and Batchelor (1987–88), with the representation theory in Martin–Saleur, Graham–Lehrer and Morin-Duchesne–Saint-Aubin. C5 is new, conditional on the project's X1 |
| Novelty of B | the reduction is a mild generalization of the known spin-charge factorization (U = ∞ Hubbard, multispecies exclusion chains). The bounds in units of M are new |
| D, E, F, G | correct as stated, with minor remarks |

Overall: B correct after one repair (Δ = 1 in B4) and one relocation of a proof into PROOF.md. C correct after the same Δ = 1 repair and a one-word correction of RESULT line 23, with the novelty statement to be narrowed to C5.

## 1. Theorem B

### 1.1 The reduction (B1), every W and every L

Domination (PROOF §3.3 Step 2). Correct. The Cauchy–Schwarz step needs only that each hop is −t times a partial isometry between fibers, which holds for any unitary W. At L = 2 the two oriented bonds give two edges between the words 10 and 01, each bounded separately, and the comparison scalar problem has amplitude −2t, which is what the scalar ring at L = 2 has. Odd L is covered with no change, since no sublattice rotation is used anywhere. The statement that hopping −t is the flux minimizer follows from the same inequality with W replaced by a phase. Every sector m is covered because the inequality is over the whole Hilbert space.

Open chain (Step 3). Correct. Termination (inversion number) plus local confluence (two available 10 patterns cannot share a site) gives path independence by the diamond lemma, and every open edge is a first move of some sorting path from its higher-inversion end, so the gauge is consistent on every open edge, not only on the sorting tree. I checked this numerically by building U_s by breadth-first search from the sorted word and testing U_s = U_{s′} W_j on every open edge: deviation ≤ 2.4·10⁻¹⁶ in 36 cases, (d₀, d₁) ∈ {(1,2), (2,1), (2,2), (1,3), (3,1), (2,3)}, L = 2 to 5, Haar-random W, three (V, μ) pairs.

The boundary twist. RESULT line 21 and PROOF §3.8 assert that the internal register enters ring spectra only as a boundary twist, for every W. PROOF §3.3 does not prove this. It proves domination and the open-chain identity, which suffice for (B1) but not for the twist statement. The proof exists in `astra_A3_adversary.out` §1 ("the transformed seam unitary is the same matrix T_{L,m} for every allowed interior") and it is correct: seam source words have fixed occupations at sites L and 1, their interiors are connected by swaps on bonds (2,3), …, (L−2, L−1), which commute with the seam hop, so Θ = U_{s′} W_seam U_s⁻¹ is independent of s, and diagonalizing Θ splits sector m into scalar rings with seam phases equal to the eigenvalues of Θ. I verified numerically that (i) Θ is independent of the seam source word to 4.6·10⁻¹⁶ and (ii) the ring spectrum in each sector 1 ≤ m ≤ L − 1 equals, with multiplicities, the union over eigenphases of Θ of the twisted scalar spectra, to 1.6·10⁻¹⁴, in the same 36 cases, including V ≠ 0 and both signs of V. Repair: copy the five-sentence argument into PROOF §3.3 as Step 3′. A corollary worth stating: on the open chain every W, entangling or not, is unitarily equivalent to the flip, so the whole content of W on the ring is the single unitary Θ.

M ≥ t (Step 5). Correct against the definition in REPORT §9.68 (M is the gauge-optimized half-spread, the minimum over Hermitian X of half the spread of h + X⊗1 − 1⊗X, with the additive constant free). The witnesses ρ_± are positive (P_mixed ± (W + W†) is twice the projector on the ±1 eigenspace of W + W†), have trace one, equal left and right marginals, and tr h ρ_± = μ/2 ∓ t. By SDP (cvxpy, CLARABEL) I find M(h_W) = M(h_sc) in 16 cases, and in all of them M equals the half-spread at X = 0, which is what the four equal-marginal witnesses of the Astra proof give.

### 1.2 The sharp bound at V = 0 (B2)

Correct. The free-fermion formula E_L^sc(m) = μm − 2t sin(πm/L)/sin(π/L) has the right boundary phase (−1)^(m−1), the field drops out before any limit, and the monotonicity of (x/sin x − 1)/x² follows from Euler's product. My ED:

- Haar-random W, (d₀, d₁) ∈ {(1,1), (1,2), (2,1), (2,2), (1,3), (3,1)}, t ∈ {1, 0.63}, μ ∈ {−2.5, −0.9, 0, 0.4, 1.3, 2.6}, L from 2 up to 12 (d = 2), 8 (d = 3), 7 (d = 4): no case with E_L(W) < E_L(sc) − 10⁻⁸, no case above 2tL(csc(π/L) − L/π), no case with M < t. That is 516 (W, t, μ, L) cases in all. The smallest E_L(W) − E_L(sc) is −1.3·10⁻¹³ (rounding), and the largest value of L(L e₀ − E_L)/((4 − 8/π)M) is 0.99994, from a random W that nearly saturates the bound.
- Saturation. With the flip W and μ = 0 the ratio L(L e₀ − E_L)/((4 − 8/π)M) is 1.000000 at L = 2 for d = 2, 3 and 4, and on even L it equals 2tL(csc(π/L) − L/π)/((4 − 8/π)t) (0.775904 at L = 4, 0.744185 at L = 6), so the intermediate bound in (B2) is attained on every even ring. With a random phase at d = 2 the L = 2 value is −2t|cos θ| and the ratio is below one, as it must be. So "sharp at d = 2, L = 2" is right, and the constant is in fact attained at every d by the flip.

### 1.3 V ≤ −2t (B3)

Correct. The shifted bond h_W − ((V + μ)/2)(n⊗1 + 1⊗n) has eigenvalues 0 on the 00 and 11 blocks and −V/2 ± t on the mixed block. ED: E_L − L min(0, V + μ) between −5·10⁻¹⁴ and 0 in 36 cases ((1,2), (2,2), (1,3); V = −2, −2.4, −3.1; L = 2, 3, 5), including V = −2t exactly.

### 1.4 V ≠ 0 imports (B4, B5)

What is needed. Only the flux-minimized sector energies. By (B1), E_L(W) ≥ E_L(sc), and the scalar ring with hopping −t in sector m is the XXZ ring J Σ(SˣSˣ + SʸSʸ + Δ SᶻSᶻ), J = 2t, Δ = V/2t, at flux π per bond, which is stoquastic and hence the flux minimizer. That is the F_L(n) of C142 (stated at flux π per bond, "at both parities of L") and of C136 §1 (F_L(n) = min_Φ e_L(n, Φ)). So no fixed-twist input is needed in B, and the bookkeeping G_L(m) = F_L(m) + (V + μ)m − VL/4 is right (I rederived it from n⊗n = 1/4 + (Sᶻ⊗1 + 1⊗Sᶻ)/2 + Sᶻ⊗Sᶻ).

Gap at Δ = 1. PROOF §3.2 quotes X1 as "for 0 ≤ Δ ≤ 1". C142 Theorem I.9 is stated "for every 0 ≤ Delta < 1" (`work_C142/PROOF.md` line 141). Heisenberg is reached in C142 only for the doubling corollary I.10, "by finite-matrix continuity at fixed L". The same continuity repairs B4: at fixed L, E_L(W) is continuous in V, and e₀ is concave in V (an infimum of affine functions of the bond), hence continuous on the open interval containing V = 2t. So the inequality E_L(W) ≥ L e₀ − 2K t/L, proved for V < 2t, passes to V = 2t, and t ≤ M finishes. This is one sentence, but it has to be written, and X1 should be quoted with its actual range. Verdict for B4: correct after this repair.

B5 imports C136 (6.13), K_all ≤ 10¹⁵ δ^(−25/3), in the convention K_L = Σ[−SˣSˣ − SʸSʸ + Δ SᶻSᶻ] with −1 < Δ ≤ 0, which is the stoquastic form above. The range and the flux convention match. 2K M with J = 2t ≤ 2M gives the stated 2·10¹⁵. Correct given X2.

### 1.5 Novelty of B

The domination inequality is a diamagnetic (Perron–Frobenius) inequality, and the open-chain identification is the squeezed-space or spin-charge factorization that is classical for the U = ∞ Hubbard ring (Ogata and Shiba, Phys. Rev. B 41, 2326 (1990), where the ring spectrum is spinless fermions with a twist set by the eigenvalues of the cyclic spin permutation) and for multispecies exclusion and XXC chains (Maassarani, Phys. Lett. A 244, 160 (1998); Alcaraz and Bariev's multispecies models). The label-conserving case of B with V = 0, d₀ = 1, d₁ = 2 is in effect the U = ∞ Hubbard chain up to the fermion sign, which is itself a twist. The extension to entangling W is a modest generalization, since on the open chain every W is unitarily equivalent to the flip. The finite-size bounds in units of M (B2 with its sharp constant, B4, B5) are, as far as I can tell, new, and B2 is an elementary free-fermion computation once the reduction is in place. RESULT and PROOF make no explicit literature claim for B, so nothing is wrong, but a citation to Ogata–Shiba and Maassarani should be added so the class is not presented as new.

## 2. Theorem C

### 2.1 The correspondence and the twist

My ED (Potts ring −Σ(E_j + F_{j,j+1}) on Q^L against H_X(2L; a, φ) in the Sᶻ = 0 sector, built independently from the definition in PROOF §4.1):

| Q | γ | L range | max |E_P − min Spec(2L; a, 2γ, 0)| |
|---|---|---|---|
| 2 | π/4 | 2 to 9 | 8.5·10⁻¹⁴ |
| 3 | π/6 | 2 to 8 | 7.1·10⁻¹⁴ |
| 4 | 0 | 2 to 7 | 3.9·10⁻¹⁴ |

The twist is not interchangeable. At Q = 3, L = 8 the Sᶻ = 0 lowest levels at φ = 0, γ and π differ from E_P by −0.020, −0.015 and +0.159, and the lowest level over the sectors Sᶻ ≠ 0 at φ = 2γ lies 0.48 above E_P. So E_P is the lowest level of Sᶻ = 0 at φ = 2γ specifically, and not the lowest level of the twisted ring over all sectors (at φ = 2γ the Sᶻ = 0 level is the lowest, but the untwisted ring lies below E_P, as it must since the twist raises the energy). At Q = 4, γ = 0, the twist is trivial and E_P is the untwisted Heisenberg Sᶻ = 0 ground state.

Literature. In the convention of PROOF §4.1 (H_X = (1/2)Σ[σˣσˣ + σʸσʸ + Δ(σᶻσᶻ − 1)], seam e^{iφ}σ⁺_Nσ⁻_1), a seam twist φ gives non-contractible loops the weight 2cos(φ/2), so φ = 2γ is the value that gives them weight √Q, and the resulting c = 1 − 6γ²/(π(π − γ)) gives 4/5 at γ = π/6, which matches the Casimir ratio below. This is the correspondence of Alcaraz, Barber and Batchelor, Ann. Phys. 182, 280 (1988) and Phys. Rev. Lett. 58, 771 (1987), whose abstract states that the Potts ground-state sector coincides with an XXZ chain with a "defect seam", with Δ = −√Q/2 in their sign convention, which is +√Q/2 in C201's. The finite twist list of (C1) is the standard decomposition of the periodic Potts space into affine Temperley–Lieb standard modules (Martin and Saleur; Graham and Lehrer; Richard and Jacobsen, Nucl. Phys. B 750, 250 (2006); Jacobsen, Ribault and Saleur, SciPost Phys. 14, 092 (2023)). The intertwiner paper is correctly identified: Morin-Duchesne and Saint-Aubin, J. Phys. A 46, 285207 (2013), arXiv:1203.4996, whose abstract states an intertwiner from the link modules of the periodic TL algebra to the Sᶻ eigenspaces of the XXZ chain and the conditions under which it fails to be an isomorphism. I did not open the paper to confirm theorem numbers either.

(C1). Correct given theorem M. My own check (`ref_potts_incl.py`): every Potts eigenvalue, both signs, lies in the union of Spec(2L; s·a, φ, m) over m = 0 with φ ∈ {±2γ, π} and m = j ≥ 1 with φ = 2πk/j, to 2.5·10⁻¹⁴, for Q = 2, 3 at L = 2, 3, 4 and Q = 4 at L = 2, 3, with the sector m = 1 removed. This tests set inclusion, not multiplicities.

(C2). Correct. The second proof (§4.8) is the stronger one: Perron–Frobenius through the explicit map f_π, then the J(v) intertwiner, invertible for large v, then equality of characteristic polynomials as Laurent polynomials in v. I checked the four cases of the f_π action and the arc table in Step 3 by hand. The first proof (§4.3 to 4.5) depends on the planar-algebra composition rule, which I did not recheck, and on theorem M.

(C3). Correct, and reproduced: with ε₀ from the Bethe integral (ε₀ = −2/π at Δ = 0 and −1 at Δ = 1/2 as checks), E_P/L converges to 2ε₀ with Casimir values L(L e₀ − E_L) falling monotonically to 0.3716 (Q = 2, L = 9), 0.6318 (Q = 3, L = 8; CFT limit π(4/5)(3/2)/6 = 0.6283) and 0.8328 (Q = 4, L = 7). These agree with `potts_casimir.py`.

### 2.2 M

(C4) is right for Q = 3, 4: the SDP gives M = 1.366025 and 1.500000, and C₀ has spectrum [0, √Q + 1]. At Q = 2 the spectrum of C₀ is [√2 − 1, √2 + 1] (u^⊥ is one-dimensional, so the antisymmetric zero vector does not exist), and M = 1 by SDP, not (√2 + 1)/2 = 1.207. PROOF (C4) correctly restricts to Q = 3, 4, but RESULT line 23 lists "Q = 2, 3, 4 … M = (√Q + 1)/2, all fully rigorous". Repair: write "M = (√Q + 1)/2 for Q = 3, 4 and M = 1 for Q = 2". At Q = 2 the bound reads 841 M.

### 2.3 Does X1 apply at fixed twist and Sᶻ = 0?

Yes, since only a lower bound is used. min Spec(N; a, φ, m) ≥ 2F_N(n) − Na/2 for every real φ because F_N(n) is the minimum over flux (C142 identifies it with the stoquastic flux π per bond, which is the minimizer by the diamagnetic inequality). The normalization H_X = 2Σ(S·S)_Δ − NΔ/2 is right. With n = N/2 (Sᶻ = 0) the chain gives min Spec ≥ N ε₀(a) − 1682/N, and N = 2L gives 841/L. The passage ε₀(a) ≤ 2e_bulk(a, ν) − a/2 for every rational ν is correct. The arithmetic 841/1.366025 = 615.66 and 841/1.5 = 560.67 is right.

Gap at Q = 4. Here a = cos 0 = 1, which is outside C142 Theorem I.9 (0 ≤ Δ < 1), contrary to PROOF §4.6 ("X1 covers 0 ≤ Δ ≤ 1"). Repair as in B4: prove min Spec(2L; Δ, 0, 0) ≥ 2Lε₀(Δ) − 841/L for Δ < 1, then let Δ → 1 at fixed L, using continuity of the finite matrix and of ε₀ (concave in Δ). Q = 2, 3 are inside the theorem. After this, (C5) is correct given X1.

### 2.4 Novelty and framing of C

C1, C2 and C3 are not new results. The level correspondence and the ground-state identity with the seam twist are in Alcaraz–Barber–Batchelor 1987–88 (at the level of exact level matching, used for numerics), and the module-theoretic inclusion is in the affine TL literature above. What C201 adds is a self-contained rigorous proof of C2 (§4.8, Perron–Frobenius plus a generic-v intertwiner), the value of M, and the finite-size bound C5, which is new only through X1, itself AI-produced and AI-refereed. RESULT presents C without any of these references. It should cite ABB 1988 for C1 and C2 and say plainly that the new content is C4 and C5.

One framing overclaim. RESULT line 27 and PROOF §4.1 say the three-state Potts infrared is a neutral, non-Gaussian sector that "no stiffness or free-fermion argument covers" and that it "now has a crossover-free bound". The bound is true, but it is obtained by embedding the Potts ground state as one level of a U(1) Gaussian host and bounding the host. So it gives no technique for neutral sectors that lack such a host, which is the case the project needs (PROOF §7.6 lists neutral interacting factors as open). The sentence should say that the bound is inherited from the XXZ host.

## 3. Proposition D, Lemmas E and F, Proposition G

Proposition D. Correct. The compression argument for M(h) ≥ M_i is right (compression to (ℂ^{d₁} ⊗ χ)^{⊗2} turns any gauge representative of h into one of h₁ plus a constant, and compression does not increase the spread). The log₂ d form uses m ≤ log₂ d. Minor: Theorem C is a list of three chains, not a class, so "drawn from the classes of Theorems B and C" should read "from class B and the chains of Theorem C".

Lemma E. Correct. I rechecked the interpolation step (g has r + 1 zeros with multiplicity, so g^(r) = f^(r) − b r! vanishes somewhere), the cyclic extension when the number of zeros is below r, and the Hölder step with Σ_i (x_{i+j} − x_i) = 2πj, which gives 2π(‖f^(r)‖/r)^p via ((r − 1)!)^p. Two small points: the proof assumes finitely many zeros, which is harmless because the argument applies to every finite subset of zeros and the bound is uniform; and the hypothesis r ≥ 3 can be relaxed to r ≥ 2. The microscopic consequence is correctly tagged physics-level.

Lemma F. Correct. A term in the column expansion with more than R columns from zT vanishes since those columns are dependent, so z^R det is a polynomial of degree at most 2R, nonzero away from flat bands. The consequence is correctly restricted to non-interacting superlattices.

Proposition G. Correct. I rederived each step: the measurement (1 ± G_α)/(2d²) is a POVM because ‖G_α‖_op ≤ 1; the outcome matrix ±d²G_α is unbiased, also conditionally; E‖X̄_U − X̄_V‖₂² = 2(a − c)/R; the identity ρ₂ − Σ_v p_v γ(v)^{⊗2} = E[(X̄_U − γ)⊗(X̄_W − γ)]; the trace-to-Hilbert–Schmidt factor d; the periodized construction with one join per L − R sites; and ⌊L/3⌋ ≥ L/5 for L ≥ 3, giving 4·15/2 = 30. The aside on W states of nonzero momentum needs one more sentence (twirl over the global U(1) before applying the site-dependent rotation, so the rotated state stays translation invariant), but it is a remark.

## 4. Overall verdict

Theorem B as stated: correct after one repair. B1, B2 and B3 are rigorous and reproduced, B2's constant is sharp (attained at L = 2 for every d by the flip). B4 needs the Δ = 1 endpoint added by continuity, because the cited sector theorem stops at Δ < 1. The boundary-twist assertion is true for every W and every L, but its proof must be moved from the Astra output into PROOF.md. The class should be credited to the U = ∞ Hubbard and XXC literature.

Theorem C as stated: correct after the same Δ = 1 repair (needed at Q = 4) and the Q = 2 correction of M in RESULT. The correspondence, the twist 2γ, the sector Sᶻ = 0 and the normalization e₀ = 2ε₀ all check out, against ED to L = 7 to 9 and against the literature, but C1 to C3 are known from Alcaraz, Barber and Batchelor (1987–88) and the affine TL literature, so the novelty claim should be narrowed to C4 and the conditional bound C5.

Lemmas E, F and Propositions D, G: correct, with the minor remarks above.

## Sources

- [Alcaraz, Barber, Batchelor, Ann. Phys. 182, 280 (1988)](https://www.sciencedirect.com/science/article/abs/pii/0003491688900152)
- [Morin-Duchesne, Saint-Aubin, arXiv:1203.4996](https://arxiv.org/abs/1203.4996)
- [Jacobsen, Ribault, Saleur, arXiv:2208.14298](https://arxiv.org/abs/2208.14298)
- [Maassarani, The XXC models, solv-int/9712008](https://archive.org/details/arxiv-solv-int9712008)
- Ogata and Shiba, Phys. Rev. B 41, 2326 (1990), and Richard and Jacobsen, Nucl. Phys. B 750, 250 (2006), cited from memory, not opened.
