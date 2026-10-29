# Referee report on the collaborator lemmas imported by C58

Referee: Claude Opus 5.5 agent, acting as an independent referee (mathematical physics), 2026-10-29. This is an AI review. It is not a human review. Nothing outside this folder was modified.

## Claim

`contrib/referee_S/README.md` accepted C48 and C58 but left a dependency gap. C58 rests on lemmas written by a collaborator (GitHub: salpace) that referee_S did not re-derive, namely SEPARATE_JW, LENGTH_FUGACITY Section 2, TILTED_MATCHING, and the pin structure of C44. This report referees those four, line by line where they are load-bearing for C48 and C58, and asks whether the chain C48 → C58 is complete modulo the two open renormalization-group inputs, Lemma N (C58's estimate (a)) and lemma TP (C58's estimate (b)).

## Status

| Imported item | Where C58 uses it | Verdict | Worst gap |
|---|---|---|---|
| SEPARATE_JW (`contrib/thermal_ensemble_20261028/exact/SEPARATE_JW.md`) | Section 5, the fermionization paragraph. It supplies C48's hypothesis of a bounded local fermionic interval representation, the reference extension, and the bound ‖K_D‖ ≤ 2‖B‖ with two-sided auxiliary-vacuum annihilation | proved | cosmetic |
| LENGTH_FUGACITY Section 2 (`.../exact/LENGTH_FUGACITY.md`) | Section 7, eqs. (25), (26), the matched renewal paragraph, and Section 9 (the proof of (30) repeats Lemma 2.2) | proved, given hypothesis (4) | cosmetic |
| TILTED_MATCHING (`contrib/direct_attack_20261028/infrared/TILTED_MATCHING.md`) | Section 5, the derivatives paragraph, eq. (21), and Section 6, eq. (23) | proved, with one sketched limit | fixable |
| C44 pin structure (`codex/work_C44/PROOFS.md` Sections 1 and 3) | Section 3, eqs. (8), (9), the invariant flags and the set D, and Section 4, eq. (18) | proved from the C44 rule text | fixable (one residual dependency on the C8 rule table) |
| Chain C48 → C58 | | complete modulo (a) and (b), with the conditions listed in the verdict | none fatal, none serious within these four lemmas |

Tags. proved means I followed the argument and found it correct at the standard of a careful mathematical-physics paper. numerical means supported by a finite computation only. heuristic means plausible but not argued. assumed means an input taken from elsewhere.

## How it was established and checked

- I read C58 PROOF.md in full and marked every sentence that cites one of the four items. I then read SEPARATE_JW, LENGTH_FUGACITY, TILTED_MATCHING, FINITE_CELL_MATCHING, the relevant parts of INFINITE_CELL_PROGRESS and COMPLETION_AUDIT, and C44 PROOFS Sections 1 to 3, and re-derived the load-bearing equations by hand. Each finding below says which equations.
- Reproduction. I copied `contrib/thermal_ensemble_20261028/`, `contrib/direct_attack_20261028/`, `codex/work_C44/` and `codex/work_C58/` to a scratch directory outside the repository, keeping the repository-relative layout, and reran `exact/jw_checks.py`, `exact/checks.py`, `ensemble/checks.py`, `ensemble/illegal_checks.py`, `check_tilted_embedding.py`, C44's `bounds_checks.py` and C58's `checks.py` with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B`. All exit cleanly. `jw_checks.json`, `illegal_checks.json`, `bounds_results.json` and C58's `checks.json` are byte-identical to the committed files. `exact/checks.json`, `ensemble/checks.json` and `tilted_embedding_checks.json` differ only in floating-point roundoff (largest relative difference 2.4e-8, on log ratios of size 1e-6, absolute difference below 1e-13). Every assertion in those scripts passes.
- Independent checks. `referee_imported_checks.py` (numpy and scipy only, about ten minutes on one thread, almost all of it in the 8000-dimensional Jordan-Wigner check) writes `referee_imported_checks.json`. It imports no code from the collaborator or Astra scripts. All assertions pass. The checks are toy instances that test the identities and inequalities the proofs rely on. They are not simulations of the 4464-state chain.
- No Codex delegation. Each check is a short script whose specification would be as long as the code.

## SEPARATE_JW

### Statement (as used)

Per cell, take n payload CAR modes c_{j,a}, one hard-core length qubit b_j, and an auxiliary space E_j = ℂΩ ⊕ QK_j carrying integer payload and length charges with parities π_j, λ_j. With p_j, q_j the logical parities, P_j = p_jπ_j, Q_j = q_jλ_j, strings S_j = ∏_{k<j} P_k, T_j = ∏_{k<j} Q_k and Klein factor K = ∏_k P_k, set ψ_{j,a} = S_j c_{j,a}, χ_j = K T_j b_j, and Ξ_j(O_{r,s}) = (Kπ_j)^s S_j^r T_j^s p_j^r q_j^s O_{r,s} for auxiliary operators of parity type (r,s). The claims are: (i) ψ, χ are CAR and Ξ_j is an isometric unital *-homomorphism with the single-parity graded rules (6), (7). (ii) Every operator on an interval X that is neutral under both total charges is an element of the CAR/graded algebra on X with the same norm. (iii) The reference extension E_JW (write the code operator in code-only CAR variables, then substitute ψ, χ) is isometric, agrees with the code operator on the auxiliary vacuum, and K_D,X = J_X B_X J_X† − E_JW(P_X B_X P_X) annihilates the auxiliary vacuum on both sides with ‖K_D,X‖ ≤ 2‖B_X‖. (iv) The coefficient bound (13), Σ_I H^{deg I}‖O_{X,I}‖ ≤ (8H²)^{(n+1)r}‖O_X‖.

### Line-by-line findings

- (3) to (5), CAR [proved]. ψ and χ anticommute because c anticommutes with K (which contains p_j) and commutes with every other factor in χ. Within each family it is the standard string computation.
- Ξ_j [proved]. The prefactor of O_{r,s} is built from involutions that commute with B(E_j), since Kπ_j = p_j ∏_{k≠j} P_k no longer contains π_j. The prefactors multiply as characters of ℤ₂ × ℤ₂, so Ξ_j is a homomorphism, and it is *-preserving because each O_{r,s}† has the same type. On a joint eigenspace of the two commuting involutions U_j = S_j p_j and V_j = (Kπ_j)T_j q_j, Ξ_j acts as conjugation by a product of π_j and λ_j, so it is isometric. I re-derived the exponent ru + sv + rv + su in (7).
- Section 4, locality [proved]. Inverting (4), (5), a monomial in the original variables on X picks up S_a to the power (number of payload fields plus Σr), T_a and K to the power (number of length fields plus Σs). Neutrality makes both even. The local parities p_j, π_j at j ∈ X are themselves local in the new variables (π_j = Ξ_j(π_j) since it is even). The argument correctly needs neutrality under each charge separately, not only total parity.
- Section 1, termwise neutrality [proved]. Averaging a term over e^{iθQ_X} is the same as averaging over the global U(1), preserves support and does not increase the norm. It also preserves [B_X, P_X] = 0 because P_X commutes with both charges. The note asserts this without the last remark (cosmetic).
- (9) [proved]. ψ_j†ψ_{j+1} = c_j† p_jπ_j c_{j+1} and χ_j†χ_{j+1} = b_j† q_jλ_j b_{j+1} = b_j† λ_j b_{j+1}, using b†q = b†.
- (10), (11), (12) [proved]. In a fixed auxiliary parity sector the substitution multiplies ψ_j by a prefix sign of the α's and χ_j by a prefix sign of the β's times a global Klein sign. The global and exterior signs cancel by neutrality, and the interior signs are a diagonal gauge unitary, so each sector block is unitarily equivalent to O_X. On the auxiliary vacuum both J B J† and E_JW(PBP) act as PBP. J B J† kills states outside ran J, maps excited valid states to excited states by block diagonality, and E_JW is diagonal in the auxiliary basis, so the two-sided annihilation holds. The norm bound follows from ‖J B J†‖ ≤ ‖B‖ and (11).
- (13) [proved]. 4^m matrix units with blocks of norm at most ‖O‖, at most 2^m nonzero terms from the vacuum projector, and degree |I| + |J| + 2|S| ≤ 2m because S must avoid I ∪ J (otherwise c_k† n_k = 0).

Numerical checks (`separate_jw` in the JSON) [numerical]. Three cells, one payload mode, one length qubit, a five-state auxiliary space with charges (0,0), (0,0), (1,0), (0,1), (−1,1), total dimension 8000. CAR errors are zero. Graded relations (6), (7) and the homomorphism property hold to 9e-16, and Ξ is isometric to 2e-14. Random doubly neutral operators on X = {0,1} and X = {1,2}, and an explicit charge-exchange vertex (payload lowering times auxiliary raising), commute with every generator outside X to 2e-15. As a negative control, an operator that is even under combined parity but odd under each parity separately fails to commute with ψ outside its support (commutator norm 2), which confirms that both charges are needed, as Section 7 of the note says. For random block-diagonal neutral B on a physical cell with four code and four excited states, E_JW was built from an independent Majorana-basis expansion. Its coefficients live only on modes in X, every auxiliary-configuration block has the spectrum of O_X to 1e-14, ‖E_JW(O_X)‖ = ‖O_X‖, K_D annihilates the auxiliary vacuum on both sides to 9e-16, and ‖K_D‖/(2‖B‖) is 0.56 and 0.60. Identity (9) holds exactly. The naive extension O ⊗ I_aux differs from the CAR quadratic hopping by an operator of norm 4, so the note's warning about the naive extension is correct. The coefficient bound (13) on two cells holds with ratio at most 0.013 of the bound for H = 1, 2, 3.

Gaps.

1. Cosmetic. The note should say that U(1) averaging preserves local block diagonality.
2. Dependency, not a gap. SEPARATE_JW shows that the retained code Hamiltonian, written in the separate code-only CAR variables, keeps its form after the reference extension. That this form is n payload bands plus one length band (that is, quadratic) is C44/C45's statement about the physical transport rules, refereed separately in `contrib/referee_C44_C45/`.

Verdict on SEPARATE_JW. Correct. It discharges C48's local fermionic representation hypothesis and the vacuum-annihilation and norm hypotheses for K_D, as C58 Section 5 states.

## LENGTH_FUGACITY Section 2

### Statement (as used)

For rooted cells with length-graded finite-dimensional K_ℓ, ℓ ≥ ℓ_min ≥ 1, on-cell A commuting with L_cell, bounded length-conserving interactions with finite incidence and cut strength, a trace-class tilt interval I, and the physical comparison (4) between Z_phys,L and Σ_M Z_{M,L} with constants (or polynomial factors), the rooted pressure p_cell(T,s) exists, is convex and strictly decreasing with slope at most −ℓ_min (8), and if it changes sign in I its unique zero equals the physical pressure (6). Section 2.1 gives an actual renewal limit lim (T/L) log S_L via superadditivity (9b). Lemma 2.2 gives s_T = s_0 + [A(s_0)/ℓ_0]T² + o(T²) from a uniform expansion and one derivative of p_0.

### Line-by-line findings

- (7) [proved]. Deleting the crossing terms is a bounded perturbation of norm at most C_cut, and the min-max principle moves each eigenvalue by at most that much, which is valid for semibounded operators with trace-class heat kernel. Finiteness of Z_M(s) follows from Tr e^{−β(H_0 + B)} ≤ e^{β‖B‖}(Tr e^{−β(A+sL)})^M. The Fekete step with remainder r < k is standard.
- Convexity and (8) [proved]. Z_{M,L} ≥ 0 because H_M commutes with L_M, so log Z_M(s) is a log-sum-exp in s. L_M ≥ ℓ_min M gives (8).
- (9) and the abscissa argument [proved]. Tonelli gives the generating identity. Strict monotonicity and existence of the limit identify the convergence abscissa of the right side with s_*. Hypothesis (4) transfers the abscissa to the physical series, whose abscissa is p_phys by (5). Polynomial factors in L do not change an abscissa, as the note says.
- Section 2.1 [proved]. The subspaces E_M (prefix of M cells with total length L_1) are mutually orthogonal because cell lengths are positive. Local length conservation lets the compressed Hamiltonian differ from the decoupled one by the cut terms only. The pinching inequality Tr e^{−βH} ≥ Σ_M Tr e^{−βE_M H E_M} (scalar convexity in an eigenbasis of each compressed block) gives (9b). Superadditivity plus eventual positivity and an exponential upper bound give the limit.
- Lemma 2.2 [proved]. |s_T − s_0| ≤ |p(T,s_0)|/ℓ_min = O(T²) from (8), then one derivative of p_0 at s_0 and continuity of A. No derivative of the remainder is used. C58 (30) uses this verbatim, with f = e − s.

How C58 uses it. In Section 7 the original bond's conditional pins are large but bounded intercell terms, so they enter B with finite cut strength (smallness is not needed), and the unpinned on-cell operator still has the trace-class window from C58 (8), (9). The matched renewal paragraph applies Section 2.1 to the enlarged space. Its hypotheses (positive integer lengths including invalid states, finite fixed-length dimension, local length conservation) hold there because enlarged states have original length ℓ or ℓ + 1 with a finite logical factor. I agree with C58's check of these.

Numerical checks [numerical]. (a) An interacting classical renewal model with cell lengths 1 to 40, length-dependent nearest-cell couplings up to 1.5 and β up to 6 (`length_fugacity_classical`). The fixed-L open sum over all cell numbers, computed to L = 4000, grows as exp(β s_* L) where s_* is the transfer-matrix pressure zero (residual growth rate below 1e-15). Rescaling with s_* + 0.01 instead leaves the expected residual rate −0.0100, so the test has resolution. (8) holds with positive slack. (b) A quantum model with cell lengths 1, 2, 3, internal multiplicities 1, 2, 2, on-cell blocks, and a random nearest-cell interaction that transfers length between cells (`length_fugacity_quantum`). The cut inequality (7) holds for all (M,N) with M + N ≤ 6 at three tilts and three temperatures, and (9b) holds for all L_1 + L_2 ≤ 11. (c) Lemma 2.2 on a model with a non-smooth p_0 and an oscillating remainder, where (s_T − prediction)/T² decreases with T.

Gaps.

1. Cosmetic. Hypothesis (4) is written with open rooted traces, while PINNED_MARKER and C58 (6) compare with periodic rooted traces and a factor L. The seam is a bounded constant factor and the L is polynomial, both allowed by the note's own remark, but one sentence would make the match explicit.
2. Dependency. Hypothesis (4) itself is PINNED_MARKER (7), whose seam Section 3 neither referee_S nor I re-derived. C58's fixed-M statement (6) uses the same seam constant B.

Verdict on LENGTH_FUGACITY Section 2. Correct. It is a clean pressure-zero theorem for interacting all-length rooted cells, and it is used correctly in C58.

## TILTED_MATCHING

### Statement (as used)

If the auxiliary operator a and the auxiliary length charge q strongly commute and |q| ≤ C_a a (C58 uses |q_aux| ≤ bQ + A_λ/κ with A_λ ≥ ΔQ), then on the strip |Re s| ≤ s_* = 1/(2C_a) one has Re a_s ≥ a/2, the semigroup e^{−t a_s} is a contraction with ‖e^{−t a_s}Q‖ ≤ e^{−tg/2} and trace norm Tr e^{−t(a + Re s q)} ≤ Tr e^{−ta/2}, the massive matching map is Banach-holomorphic in s, and Cauchy's formula bounds ‖∂_s^m E_{s,β}‖ ≤ m!(2/s_*)^m C_match r on the inner strip, with the analogous bound for the thermal periodization error.

### Line-by-line findings

- C58's input [proved]. From (18), A_λ ≥ κ(ℓ + 1) on Q_ℓ, so q = ℓ − b satisfies |q| ≤ bQ + A_λ/κ. With |Re s| ≤ s_a = min(κ/4, Δ/(4b)) and A_λ ≥ Δ on Q, |Re s||q| ≤ Δ/4 + A_λ/4 ≤ A_λ/2. Since min(1/(4x), 1/(4y)) ≤ 1/(2(x+y)), C58's s_a is inside TILTED_MATCHING's strip with C_a = b/Δ + 1/κ.
- (3), (4) [proved]. a_s is normal, |e^{−t a_s}| = e^{−t Re a_s} in the joint spectral representation, and the unitary factor e^{−it Im s q} commutes.
- Section 2, holomorphy [proved, with a sketched limit]. Domain equality, closedness and normality are correct. The derivative (−tq)^k e^{−t a_s} is bounded on smaller strips since t|q| ≤ C_a t a and Re a_s ≥ c a with c > 1/2. The passage from finite spectral cutoff to the full map uses strong convergence of compressed vertices, Schatten-class semigroups and dominated convergence. That is the right argument, but it is written as a sketch.
- Transfer to C48 [proved, by inspection of C48 and INFINITE_CELL_PROGRESS]. TILTED_MATCHING is phrased for the abstract matching theorem of INFINITE_CELL_PROGRESS. C48's estimates use only that free segments are contractions (INFINITE_CELL_PROGRESS Section 3), the passive cost ‖Q e^{−τa}‖, and the Schatten-Hölder trace bound with exponents β/s_i. For normal a_s each of these is controlled by Re a_s, so all of C48's constants hold uniformly on the strip with gap Δ/2 and heat trace Tr e^{−tA_λ/2}. C58 states this in one sentence. It is correct.
- (10), (11) [proved]. Circles of radius ρ < s_*/2 about real centers |s| ≤ s_*/2 lie in the open strip, and the supremum on the closed strip is C_match r.
- Section 4 [proved]. Valid embedded length is b + N_len + q_aux, the scalar and the logical number term stay in the logical Hamiltonian, and a length-preserving SW unitary leaves D and R independent of s. With SEPARATE_JW, N_len is unchanged by the Jordan-Wigner map, so the tilt adds the quadratic term s χ†χ.

What it is load-bearing for. TILTED_MATCHING is not used in the activated comparison (1) to (3) of C58. Those compare two exact traces, p_o and p_m, through (17) and (24), with real tilts only. TILTED_MATCHING (together with SEPARATE_JW) is load-bearing for the interface with (a) and (b). It is what makes the matched action S_{β,s}, to which (a) and (b) are applied, a C48 output with all-degree kernel bounds and controlled source derivatives uniformly in s.

Numerical checks (`tilted_matching`) [numerical]. For spectra saturating C58 (18) with b = 4, 5, 9, 30, the bound on |q| holds, min Re(a_s)/a ≥ 0.755 at |Re s| = s_a, and accretivity is lost only at a tilt about four times s_a. The complex-tilt Schatten-Hölder bound |Tr ∏ e^{−t_i a_s} O_i| ≤ Tr e^{−β Re a_s} ∏‖O_i‖ ≤ Tr e^{−βa/2} ∏‖O_i‖ holds in 200 random trials with up to four noncommuting insertions and |Im s| up to 5. For a Feshbach-type influence map W(s) = PKQ(a_s + QKQ)^{−1}QKP the supremum over the strip (|Im s| up to 50) stays below the uniform resolvent bound ‖K‖²/(Δ/2 − ‖QKQ‖), and the exact derivatives up to order five obey the Cauchy bounds with ratio at most 0.003.

Gaps.

1. Fixable. The cutoff removal in Section 2 should be written as a lemma with the dominated-convergence majorant displayed.
2. Fixable. C58 should say in one line which C48 estimates are re-used under complex tilt (contraction, passive cost, Schatten-Hölder) and that each depends only on Re a_s. I checked that this is so.

Verdict on TILTED_MATCHING. Correct, with one sketched limit. It does what C58 says.

## C44 pin structure

### Statement (as used)

C58 Section 3 needs: (i) the selector arrows #:(1,0), L:(0,0), X:(0,1), R:(1,1) and legality of a bond as equality of outgoing and incoming arrows; (ii) pins on the cut projector X# + R#, track B on the left and track A on the right; (iii) the only unperturbed intercell terms that depend on a neighbor's state are these conditional pins (wrapper bonds are cut at markers, the K_sel penalty is diagonal in selector letters, ports and marker motion are weak); (iv) the per-cell mismatch counts z_j are invariant under head hopping LX ↔ XR and marker motion R# ↔ #L; (v) the fragment bound giving A_bad ≥ U + |e|ℓ, and the legal-cell bound (9) with γ = t_sel η_b − U; (vi) the all-string bound (18) A_λ ≥ κ(ℓ + 1)Q_ℓ with κ = Δ/(b + 2).

### Line-by-line findings

- (i), (ii) [proved from C44 Section 1 text]. A cell is legal exactly when its selector word is # L^a X R^c, so there is no legal length-one cell. A legal bond into # requires an outgoing arrow 1, so the predecessor's last letter is X or R and both pins are enabled. Hence cells outside D = S ∪ (S+1) carry both pins independently of their neighbors.
- (iv) [proved]. The outer arrows of LX and XR are both (0,1), and those of R# and #L are both (1,0), and all four internal bonds are legal. For R# → #L, the old aR mismatch becomes a#, with the same incoming arrow, and the old #c becomes Lc, with the same outgoing arrow. I re-derived the case analysis.
- (v) [proved, given C44's fragment bound]. K_sel z − (2t_sel + 4U)(z + 1) = U + (2t_sel + 5U)(z − 1) ≥ U for z ≥ 1. C58 adds μ_sel on top. Whether μ_sel is already inside C44's "fragment costs at worst −2t_sel − 4U" is not stated, but the bound ≥ U holds without it, so nothing depends on that. Subtracting eℓ with e < 0 gives (8). COMPLETION_AUDIT (3) states the same inequality, with μ_sel, for the completed bond. Dropping positive pins and helix terms gives a lower operator bound, so (9) is C44 (1), which does not use the pins.
- (vi) [proved]. With Δ ≤ min(complete selected gaps, γ, U, |e|) and κ = Δ/(b+2): at ℓ ∈ {b, b+1}, κ(ℓ + 1) ≤ Δ. On legal wrong lengths, γ(1 + |ℓ − b|) ≥ κ(ℓ + 1) because γ ≥ Δ. On illegal cells, U + |e|ℓ ≥ Δ(1 + ℓ) ≥ κ(ℓ + 1). This avoids COMPLETION_AUDIT's γ_all = min(γ, U/b, |e|), which it does not need, because the illegal bound is U + |e|ℓ and not just γ_all(1 + |ℓ − b|). The heat trace (18) is then Σ_ℓ d^ℓ e^{−tκ(ℓ+1)} with dim K_ℓ ≤ 4464^ℓ, equal to 1/(2d) at t_0 = log(2d)/κ.

Numerical checks (`c44_pin_structure`) [numerical]. All 320,000 ring words of length 2 to 9 over {#, L, X, R} containing a marker, and all 701,536 literal selector rewrites on them, with ordered marker labels tracked through marker motion. Every per-cell mismatch count is preserved. In every word, every good cell ends in X or R and every cell after a good predecessor sees X or R before its marker (zero pin-enable failures), and every good cell has the shape # L^a X R^c. This is an independent implementation, separate from C58's 6400-rewrite test. The fragment arithmetic and the κ bound hold on 20,000 random parameter sets (slack at roundoff, since both are tight at their defining points).

Gaps.

1. Fixable, a residual dependency. C58's claim that the conditional pins are the only strong neighbor-dependent intercell terms is read off C44 Section 1's description of the C8 selector_wrapper bond. I did not reopen the C8 rule table. Everything C44 states is consistent with it, and referee_C44_C45 rebuilt the selector bounds.
2. Cosmetic. The μ_sel term in C58 (8) should be dropped or justified. It is not needed.
3. Cosmetic. C58 (18) is stated with "the imported legal and illegal estimates imply". The three-case verification above should be written out, since it is the reason C58 can use Δ/(b + 2) rather than γ_all.

Verdict on the C44 pin structure. Correct as C58 uses it.

## Verdict

None of the four imported lemmas has a fatal or serious gap. SEPARATE_JW and LENGTH_FUGACITY Section 2 are correct as stated. TILTED_MATCHING is correct with one limit argument written as a sketch. The C44 pin structure supports C58 Section 3 with one residual dependency on the C8 rule table. Every supplied script reruns with identical or roundoff-identical output, and the independent checks pass.

On the question asked. The chain C48 → C58 is complete modulo the two open inputs, (a) the all-scale constructive RG with the length-source component (Lemma N) and (b) the uniform subtracted scalar estimate (lemma TP), provided one also accepts:

- (i) the imported finite construction C39, C42, C44, C45 and C49's physical-support lift, refereed separately (`contrib/referee_C44_C45/`, `contrib/referee_C39/`, and referee_S for C49), with no fatal gap reported;
- (ii) collaborator pieces outside these four that C58 also imports: TILTED_EMBEDDING Sections 3 and 4 (re-derived by referee_S), PINNED_MARKER Section 1 (re-derived by referee_S), PINNED_MARKER Section 3 (the seam constant, not re-derived by anyone), and ILLEGAL_SECTORS with COMPLETION_AUDIT (3) (the illegal-cell bound A ≥ μ_sel + U + |e|ℓ for the completed bond, whose arithmetic agrees with C44 (2) but whose open-wrapper input I did not re-derive);
- (iii) the fixable presentation items above.

The activated free-energy comparison (1) to (3) of C58 depends on LENGTH_FUGACITY and the pin structure, and not on SEPARATE_JW or TILTED_MATCHING. Those two are what make the matched action a legitimate C48 output, uniformly in the length tilt, so that (a) and (b) can be posed on it. The scope of the conclusion is unchanged from referee_S. C58 is a free-energy result with activated errors. It gives no o(1/L) spectrum, and the Virasoro central charge still needs C51's completeness hypothesis G on top of (a) and (b).

Forward directions, in order of value: (i) a human check of PINNED_MARKER Section 3 and of the ILLEGAL_SECTORS open-wrapper bound, the two imported pieces nobody has re-derived; (ii) writing TILTED_MATCHING's cutoff removal as a short lemma; (iii) one sentence in C58 listing which C48 estimates survive complex tilt.

## Files

- `README.md`, this report.
- `referee_imported_checks.py`, independent numerical checks (numpy and scipy).
- `referee_imported_checks.json`, its output.
