This report was written by a Claude Opus 5.5 agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C171: referee report on C159 and C163 (the T and K_I chain)

Tags: [proved], [derived], [numerical], [heuristic], [assumed], as in C159. Notation and section numbers are those of codex/work_C159/PROOF.md and codex/work_C163/PROOF.md. Derivations behind the new counts are in PROOF.md. The operator counts were redone independently by opcount.py and lemmaL_check.py (console output in opcount.out and lemmaL_check.out). They enumerate every relevant vertex, current and discrete-sector product, impose invariance under one-site translation T1 and under the antiunitary, and count invariant Hermitian directions. Runs take seconds, and no Typhon job was used.

## Summary of grades

| Claim | Grade |
|---|---|
| C159 Lemma 1 (normal form of K, sheets K_I, K_II) | correct |
| C159 2.1 (Ogata–Tasaki on T, none on K) | correct |
| C159 Lemma A (one field of a floating boson is commensurate) | correct for one boson with T1 in the identity component. Gap: T1 acting by charge conjugation, and multi-sector automorphisms, are not covered |
| C159 2.2 consequences (T liquid stable iff K > 1/2, K_I shift-0 branch K > 2, K_II K < 1/2) | correct |
| C159 Lemma B (no symmetric c < 1 critical point in T) | correct (Ising, tricritical Ising, Potts derived; D and E series assumed, as tagged) |
| C159 2.4 (SU(2)_2, SU(2)_3 excluded) | correct for translation acting as the centre. Chiral actions reduce to Ising × boson rows, so nothing is lost |
| C159 table 4.1, rows as written | correct (reproduced by opcount.py) |
| C159 table 4.1, as a complete list | wrong. Potts × half-floating boson (period-3 order) has codimension 1 at every K > 1/2 |
| C159 Proposition C, as a statement about the listed rows | correct |
| C159 Proposition C, as used ("the only codimension ≤ 1 window is attractive") | wrong, because of the Potts row and the θ-commensurate second boson |
| C159 dismissal of codimension ≥ 2 candidates, and the replacement H_T → A6_TK | wrong (non-sequitur; drops C154's pattern P2) |
| C163 Fact 1 | correct, but the tag should be [derived], not [proved] |
| C163 Fact 2 (zero-field attractive parent has only the Gaussian tower) | correct [literature] |
| C163 Fact 3 / §1.2 ("every soft wavevector of the U(1) parent lies in Zα + Zπ") | wrong at finite magnetization (2k_F is incommensurate). Harmless on T and K_I, where the magnetization along m vanishes |
| C163 Lemma L, branch Q_2 = n Q_1 | correct for n = 1, 2 at K_1 > 1. For n ≥ 3 it holds only where n²/(4K_1) + 1/(4K_2) < 2 |
| C163 Lemma L, branch Q_2 = n Q_1 + π | wrong. No spin-0 locking operator is relevant for n = 2 at any 1 < K_1 ≤ 2.5 tested |
| C163 Lemma U (uniform Θ-odd Ising × boson needs codimension 2) | correct count. "The tilt is the boson's own density channel" is derived only near h_F |
| C163 §1.4(b) energetics ("dimerization needs frustration") | wrong as stated, though not load-bearing |
| C163 Proposition P (one boson to all orders in the harmonics) | correct at the physics standard |
| C163 Proposition G (twisted curve) | correct, conditional as tagged. It covers only the large-K corner near γ |
| C163 numerics | consistent with REPORT §9.125–9.126. KI22, inside W, is not yet cleared |
| (A_W) as the single residual assumption on T and K_I | incomplete. The reduction is not complete (section 4) |

## 1. C159, sections 1 and 2

1a. Lemma 1 [proved, checked]. The antiunitary Θ R_π(z) acts on Pauli vectors as Q = diag(1, 1, −1), and invariance of the bond needs Q J Q = J and Q B = B. With Q A(D) Qᵀ = det(Q) A(Q D) this forces D ∥ z, no xz or yz symmetric part, and B in the plane. One rotation about z diagonalizes the xy block. Counting (J_x, J_y, J_z, D_z, B_x, B_y) modulo scale gives dimension 5. The sheet condition det M = J_x J_y + D_z² decides whether the planar polar factor is a rotation (spiral axis z, parent at zero field) or a reflection across u (O is the π rotation about u, parent in the field u·B). This is correct.

1b. Lemma A [derived, re-derived]. At generic radius the antiunitary automorphisms are T0 ∘ shift (θ reversed) or T0 ∘ C ∘ shift (φ reversed). Take A reversing φ, with A(φ) = −φ + c and A(θ) = θ + d. Then A V_{p,n} A⁻¹ ∝ V_{−p,n}. T1 multiplies V_{p,n} by e^{i(pb + 2na)}, and A T1 A⁻¹ multiplies it by e^{i(pb − 2na)}. Commutation therefore needs 4a = 0 mod 2π, which means the reversed field's shift is 0 or half a period, while b is free. This is correct. Two scope limits matter later.
- The proof assumes T1 acts in the identity component. If T1 = C ∘ shift, then T1² acts trivially, the boson is fully commensurate, and C159's table has no row for it (section 2d).
- For an infrared with two bosons the automorphism group is the larger group of the rank-2 Narain lattice. Θ can then exchange the two sectors (φ1 → −φ2, θ1 → θ2, and so on). In that case commutation forces equal φ shifts and opposite θ shifts, and both bosons are fully floating, not half-floating. A hand count gives cos and sin of 2(φ1 − φ2) and of (θ1 + θ2), whose dimensions 2K and 1/(2K) cannot both exceed 2, each contributing two real couplings. So this option has codimension ≥ 2 [derived] and falls under section 3c.

1c. Application to T [derived]. For the liquid connected to the zero-field parent, Θ reverses the spin density along m, hence φ. The step "shift 0 gives a symmetric unique gapped state at small K, contradicting LSM" is right: the anomaly is constant on the conformal manifold, and cos(2φ − c) with c the Θ shift pins φ at a Θ-fixed point. One refinement. The anomaly needs the φ shift π/2 and the Θ shift π of θ. An anomaly-free boson can still have V_c at momentum π, if Θ does not shift θ. C159 4.1 (iv) allows this, so the table is consistent.

What Lemma A does not decide is which field Θ reverses for an emergent second boson. If its density is Θ-even (a fluid of neutral composites), Θ reverses θ2. Its commensurate vertex is then e^{iθ2}, with x2 = 1/(4K2), and its density is free. C159 handles this in the x-language of the table, but the translation "x = K" in Proposition C, and C163's restriction to hosts with K > 1, lose it (section 3b).

1d. 2.2 consequences [derived, numerical check]. opcount.py, case S1, finds no invariant relevant operator of the T host for K > 1/2, and only cos 4φ for K < 1/2. The K_I and K_II statements follow in the same way. Correct.

1e. Lemma B [derived]. Every invertible automorphism of Ising fixes ε. The argument is that σ → ±σ, and an antiunitary map preserves the real coefficient C_σσε. The KW map ε → −ε is non-invertible and is not a symmetry of the local algebra. Hence ε is allowed, its disordered sign gives a unique symmetric gapped vacuum, and Ogata–Tasaki forbids that. The same structure holds for TCI (ε, dimension 1/5) and Potts (S_3-singlet ε, dimension 4/5). Gapped anomalous phases in 1+1 dimensions must break a symmetry, so "inside a broken background" is the only exit. I checked the free-locus corollary from the BdG band. Gaplessness at k = 0 or π needs t = 0, which closes both points. Interior zeros solve (4D² − δ²) sin²k = t² cos²k, a pocket with two Fermi points of one band, c = 1. At B = 0 every gapped point is in the ordered (Kitaev-topological) phase, which breaks Θ. Correct.

1f. 2.4 [derived]. The k mod 2 centre anomaly and the dimensions 1 (SU(2)_2) and 4/5 (SU(2)_3) of the spin-1 primary are right. Each spin-1 primary has nine translation- and Θ-even components when translation acts as the centre. If translation acts as an incommensurate chiral rotation, only Φ^{zz} survives, and SU(2)_2 becomes a point on the Ising × boson manifold (C154 4.2). Nothing new results. Correct, with this scope remark.

## 2. C159, sections 4 and 6: the codimension table and Proposition C

2a. Rows as written [numerical, opcount.out S1, S2, S2']. The single half-floating boson needs no tuning for K > 1/2. Ising (σ at π) × boson gives ε plus σ V_c at 1/8 + K, so k = 2 for K < 15/8 and k = 1 above, for both Θ parities of σ. Ising (uniform, Θ-odd) × boson gives ε at spin 0 plus σ ∂_xφ at dimension 9/8 and spin 1. The thresholds K = 15/8 at Δ = cos(11π/15) = −0.669, and the ordered-side pinning for K < 2, check out. The rows are correct.

2b. Missing row 1, Potts × half-floating boson [derived, numerical S3]. Let translation act on a three-state Potts sector by its Z3 generator, meaning period-3 order at momentum ±2π/3. Then σ, σ† (2/15) and Z, Z† (4/3) carry momentum ±2π/3. No vertex of a half-floating boson neutralizes that, because its vertex momenta are p(π + α) + nπ, which hit ∓2π/3 only on resonance surfaces of α. So the only invariant relevant operator is ε (4/5). Θ acts as σ → σ† (forced by commuting with the Z3 action) and fixes ε. opcount.py gives k = 1 at every K from 0.6 to 3.0.

So c = 9/5 is kinematically allowed on hypersurfaces of T and K at every K > 1/2, with no attraction. On the ordered side ⟨σ⟩ cannot pin the boson for the same momentum reason, so its neighbour is "liquid plus period-3 order". This is pattern P1 with Z3 in place of Z2, and in T it is less constrained than Ising, because the Z2 order is neutralized by V_c and Z3 order is not. The same omission sits in C154 table 4.2 for N, where Potts × floating boson also has k = 1. There C154's (A6) still covers it. Clock orders with N = 4 (Ashkin–Teller × boson, c = 2) are similar, except that the charge-2 polarization at π can pair with V_c, which brings back a K threshold.

2c. Missing row 2, a sector on which translation acts by charge conjugation [derived, numerical S4]. Such a boson always has at least one relevant invariant, cos 2θ2-type (1/K2) or cos 4φ2-type (4K2), so it alone needs one tuning. Its T1-odd operators pair with V_c1 of the host at momentum π. Minimizing over all admissible translation and Θ actions (opcount.py scans them), codimension 1 needs K1 > 2 − min(K2, 1/(4K2)) ≥ 3/2 (k = 1 first appears at K1 = 1.6 for K2 = 0.6). This stays in the attractive window, but the second sector sits at momentum 0 or π, which is inside Zθ + Zπ, so the concrete list of (A_W) misses it.

2d. Proposition C, the K-translation [derived, numerical S5]. In the two-boson row the requirement is x1 + x2 > 2 with both x > 1/2. If the second boson is θ-commensurate (1c), then x2 = 1/(4K2). Take K1 = 0.9 (a repulsive host) and K2 = 0.2. opcount.py finds no invariant relevant operator at all: k = 0, a codimension-0 c = 2 candidate whose host is repulsive. A6_TK covers it, since boson 2 has x2 = 1.25 > 1. The phrase "the second needs at least one attractive liquid (K > 1, Δ < 0)" does not, and neither does (A_W), which is stated for hosts with K > 1. In T, LSM in addition forbids boson 2 from being the only gapless sector, so the c = 2 region cannot border a boson-2-only liquid. Its c = 1 neighbours are the repulsive host or boson 2 inside a broken background.

C163's remark that a strongly repulsive composite (K2 < 1/4) "would be gapped" is true only for a commensurate density. Here the density is free, and the floating crystal stays gapless.

2e. Codimension ≥ 2 [derived]. Proposition C says candidates with k ≥ 2 "must lie near strata where extra symmetry removes operators". That does not follow. A k = 2 locus is the transverse intersection of two hypersurfaces, a dimension-3 set in the 5-dimensional stratum, and it is generic. C154 4.2 makes only the weaker statement "should be found at or near", and it is unsupported there too.

C159 §6 then replaces H_T by A6_TK by quoting C154 5.2 as "P1/P3-type sub-loci". C154 5.2 lists P1, P2 or P3. P2 (uniform order × boson, k = 2) has, in T, a Θ-odd σ coupled through σ ∂_xφ at any K. Its neighbours are c = 1 phases only: ε of either sign, with ⟨σ⟩ ∂_xφ acting as a field along m. So a P2 component in T does not need any attractive sector nearby, and A6_TK does not exclude it. The same holds for Ising (σ at π) × boson at K < 15/8 (k = 2) and SU(2)_1 × Ising (k ≥ 2). These are the codimension-2 multicritical points where an independent Ising mode meets the liquid. Grade: wrong. The replacement H_T → A6_TK is valid only for codimension ≤ 1 candidates.

## 3. C163

3a. Fact 1 [derived]. It rests on Lemma A, which is [derived], so the [proved] tag is too strong. Fact 2 [literature] is correct. The zero-field critical XXZ chain has the c = 1 Gaussian finite-size spectrum at every γ, and the strings at Δ < 0 belong to the same theory. The breather identification is a side remark and is not needed.

3b. §1.2 and Fact 3. The single-harmonic n-magnon bands above the polarized state, with minima at nα or nα + π, are right. The sentence "every soft wavevector of the U(1) parent lies in Λ" is false at finite magnetization, where the 2k_F density modes sit at incommensurate momenta set by the filling. On T and K_I the magnetization along m is zero (Lemma A), 2k_F = π lies in Λ, and nothing downstream changes. Restrict the sentence to zero field.

3c. Lemma L [derived, numerical lemmaL_check.out]. The dimension formula x_conv = n²/(4K1) + 1/(4K2) is right. For n = 2 and K1 > 1 the conversion vertex is relevant whenever K2 exceeds 1/(4(2 − 1/K1)), which is below 1/4, so the main case is correct.

The branch Q2 = nQ1 + π is not. The proposed lockers are O_conv cos 2φ1, which is always irrelevant since 1/K1 + K1 ≥ 2, and O_conv², of dimension n²/K1 + 1/K2, which for n = 2 is irrelevant whenever K1 ≤ 2. The product O_conv e^{2iφ2} is momentum-allowed, but it has conformal spin 1. With boson 2's V_c at π, the scan finds no spin-0 locking operator for n = 2 at any (K1, K2) tested between (1.2, 0.8) and (2.5, 2.0). Each of these is a codimension-0 two-boson point with K1 + K2 > 2. For Q2 = nQ1 with n = 3 the vertex is irrelevant at K1 = 1.2. For n = 4 it is irrelevant up to K1 = 1.8. The neutral case n = 0 (a composite at momentum 0 or π) is not treated.

Corrected statement: a charge-n composite at Q2 = nQ1 is locked if and only if n²/(4K1) + 1/(4K2) < 2. At Q2 = nQ1 + π it is locked only if n²/K1 + 1/K2 < 2. Otherwise two bosons whose wavevectors lie in Λ survive at codimension 0. The consequence in §1.5, "two bosons in W need two independent incommensurate wavevectors", is therefore wrong as stated.

3d. Lemma U. The count is correct (2a). The claim that the soft uniform Θ-odd mode is the boson's own density channel is derived only in the semiclassical regime near h_F (C160 M3, refereed in C161). Elsewhere it is heuristic. In any case, codimension 2 does not exclude the candidate (2e).

3e. §1.4(b) heuristic. "Dimerization needs frustration, which a nearest-neighbour bond does not have" is contradicted inside the programme. The cyclic stratum near XZ (in the closure of T) dimerizes spontaneously (C154 3.4). C159 3.2 gets a dimerized phase near the AF Heisenberg point from the sign of the marginal coupling. The statement is already absorbed into (A_W)(a), so nothing else changes. The continuity remark (a c = 3/2 segment at K > 15/8 cannot join a KT segment at K = 1/2) is correct.

3f. Proposition P [derived]. The resonance condition (2.1) and the tongue widths a^{k(Q)/(2 − Q²/(4K))} with k(Q) = ⌈|Q|/2⌉ are correct, as is relevance for Q² < 8K. Perturbation theory around a CFT with a lattice-scale gap cannot create a sector, so the all-orders statement holds at the physics standard. Two remarks. At finite K the staircase is finite, because only Q < √(8K) tongues open. The threshold "harmonic fraction about 0.3" has no derivation and is [heuristic].

3g. Proposition G [derived, conditional as tagged]. Identity (3.1) checks: R(jφ)(−R(φ))R(−(j+1)φ) = −1. In (3.2) the charge-q parts of X R(−φ) pick up e^{iqjφ}. The q = 0 part is XXZ + DM + field along m, and −R(π) = XX + YY − ZZ is the staggered image of h_F. The gap estimate at momentum qφ is of order ε^{1/2} J |sin qφ| in the easy-plane liquid, which is consistent with the spinon bandwidth near Δ = −1. At φ = π/2 the charge-2 term is resonant at first order with e^{2iθ} cos 2φ, which the text lists and which is irrelevant at K of order ε^{−1/2}. The proposition is correct. It covers only K of order ε^{−1/2}, the large-K corner, and not the intermediate-K window W.

3h. Numerics. The C159 job readings and the reading of dist_T_50 and hi_T_40 agree with REPORT §9.125 and with C166's Néel-cat resolution, which places both outside W. One first-stage pass inside W is not cleared. KI22 is on a K_I path with parent Δ = −0.75, K_S(16) = 3.2, c_two 1.76 and c_prof 1.31 at L = 64 and 128. It sits at λ = 0.436, where S_mid jumps from 0.50 to 0.12 between λ = 0.43 and 0.44 (c163 typhon_out/summary.txt), which suggests a level crossing. Its FES second stage is in C166's job 17602270, and codex/work_C166/typhon_out does not exist yet. Until it reports, "no c > 1 in W" is numerical with one open point.

## 4. Is the reduction to (A_W) complete?

No. Every route below reaches c > 1 on T or K_I without passing through (A_W) as stated, or through any derived exclusion in C159 or C163.

- R1. Potts (Z_N, N ≥ 3) order × half-floating boson, codimension 1 at every K > 1/2 (2b).
- R2. Every codimension ≥ 2 candidate at any K: P2-type uniform Ising × boson, Ising at π × boson for K < 15/8, SU(2)_1 × Ising, and Θ-exchanged boson pairs (1b, 2e).
- R3. Two bosons whose wavevectors lie in Λ but are not locked: the Q2 = 2Q1 + π branch at 1 < K1 ≲ 2.5, n ≥ 3 composites, and neutral composites (3c). (A_W)(b) covers only wavevectors outside Zθ + Zπ.
- R4. A repulsive host with a θ-commensurate second boson, x2 = 1/(4K2) > 2 − K1 (2d). (A_W) is stated for hosts with K > 1.
- R5. A second sector translated by charge conjugation, codimension 1 for K1 > 3/2 (2c). It is at momentum 0 or π, outside the concrete list of (A_W).

The general sentence of (A_W), "no lattice-scale excitation outside the boson's charge lattice reaches zero energy while the liquid survives", covers R1 and R5 for attractive hosts. It does not cover R3, where the composite carries charge inside the lattice, or R4 and R1 at K < 1, where the host is repulsive. R2 is outside it entirely.

Corrected residual assumption [assumed]:

    (A_W') In T and K_I, outside the controlled neighbourhoods (λ-neighbourhood of U(1) and XYZ, AF Heisenberg, cyclic SU(2)_1, h_F and γ, ferromagnetic onsets, free locus), no gapped sector becomes critical or forms an independent liquid while a half-floating liquid is present, at any K and at any codimension. This includes Ising modes at 0 or π independent of the boson's vertices, Z_N orders at 2πk/N with N ≥ 3, and second bosons whether their wavevector is outside Λ, inside Λ with an irrelevant conversion vertex, θ-commensurate, or translated by charge conjugation.

This is C154's (A6) restricted to T ∪ K_I. What survives of C163's narrowing is the restriction to harmonic amplitudes of order the exchange (Proposition P, Proposition G and C160 cover the rest). The restriction to attractive liquids does not survive. It holds only for the codimension-1 Ising-at-π row (K > 15/8), the charge-conjugation row (K1 > 3/2), and the composites that Lemma L does lock.

Two lemmas would restore part of the narrowing. (i) Exclude Z_N (N ≥ 3) order independent of the boson, which is the analogue for period 3 of C165's statement at momentum π. (ii) Show that a codimension-2 c > 1 point in T needs a codimension ≤ 1 c > 1 structure next to it. That is false for P2 kinematically, so it would have to come from dynamics.

Scope notes. K_II is not covered by (A_W) or by anything in C163. Its window needs x = 1/(4K) > 15/8 or x1 + x2 > 2, which XXZ parents in a field (K ≥ 1/4) do not reach, so it too lives at harmonics of order J. STATUS's "on the time-reversal strata it holds given one assumption" should therefore read "on T and K_I, given (A_W'), with K_II open". The Potts × floating boson omission in C154's table (k = 1 in N) should go to C174, since C154's statement that P1 and P3 are the only candidates on hypersurfaces or open sets is incomplete. The conditional theorem survives, because (A6) covers it.
