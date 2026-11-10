This report was written by a Claude Opus 5.5 agent (2026-11-09), not by Codex. Astra was not used (out of credits).

# C179: referee report on C173 PROOF sections 3 to 5, with the hidden-symmetry check

Tags: [proved], [derived], [numerical], [heuristic], [assumed], as in C173. Section numbers in the left column are C173's. Derivations are in PROOF.md of this folder (cited as P1 to P5), scripts and outputs are in this folder. Everything was re-derived independently. The numerics use my own codes (sw_exact.py, sw_ring.py, ringk.py, ringv.py, twobody_pi2.py), validated against exact limits and against C173's codes where they overlap. All local runs, no Typhon job.

## Summary of grades

| C173 claim | Grade |
|---|---|
| 3.1 t_2 = −Δ²/(2V), sign-free NNN hop | correct |
| 3.1 shift 3Δ²/V, V_2 = 2t_2, V_3 = (2/3)t_2, nothing at distance ≥ 4 | correct (at second order; V_4 appears at relative ε²) |
| 3.1 ΔN = 2 only at wΔ/V and Δ³/V² | correct, with a refinement: the wΔ/V part is ∝ J_x + J_y and vanishes to all orders on the stratum |
| 3.2 the closed form agrees with the full chain, error ∝ V⁻³ | correct (confirmed by an exact extraction) |
| 3.2 "ΔN = 2 terms enter energies only at sixth order" | gap (wording): they enter at relative (Δ/V)², the same order as the N-conserving corrections |
| 3.3 legs are ferromagnetic Heisenberg at w = 0, C153's large-b box, N_A − N_B exact | correct |
| 3.3 band, valleys, offsets, curvatures | correct |
| 3.3 the symmetry U on J_x + J_y = 0 | correct, and it is the only one within inversion × z-rotations |
| 3.3 every inter-valley conversion vertex carries momentum | correct vertex by vertex, gap as used: two products are resonant on the whole stratum |
| 4(a) same-valley pairs are identical fermions | correct |
| 4(b) inter-valley pairs are impenetrable (no zero-energy resonance) | wrong on the stratum: there is an inter-valley bound state at K = π |
| 4(c) and the wedge |δ| ≲ n³, |η| ≲ n^{1/2} | correct off the stratum for δ ≫ E_b, wrong near it (the window and the competing phase are set by the binding) |
| 5.1 no two-body bound state, "a gas, not a pair fluid" | correct at φ = 1.3 where it was checked, wrong at the candidate φ = π/2 |
| 5.2 DMRG table | not re-run (C176's task) |
| 5.3 valley-balanced rings imply easy-plane or antiferromagnetic | gap: inter-valley pairing also gives a balanced ground state and the observed odd-N doublet |
| 5.3 the gauge map to the J1–J2 chain | correct |
| 5.4 open inside the parity class | gap: needs c = 2 off the stratum, supported by one L = 48 point |
| 5.4 parity-odd terms are ΔN = 1 of order Δ²/V and off resonance | off-resonance correct, size wrong as stated (they are first order in the perturbation, PXP-type) but harmless |
| 5.4 open in B_2 | gap (fixable): O_4 and the valley backscattering are resonant on the stratum, the openness survives only through a commensurate-incommensurate argument off the stratum |
| Candidate bond numbers | correct |
| Hidden symmetry (C163 §4.2) | none found: not free, not integrable, not decoupled, and D breaks C153's U(1) at leading order |

## 1. The Schrieffer-Wolff model (C173 3.1, 3.2)

1a. [proved, numerical] Conventions. I re-derived the Jordan-Wigner form (0.1) from the Pauli bond: w = J_x + J_y + 2iD, Δ = J_x − J_y, V = 4J_z, μ = 4J_z + 2h_z, particle = down spin. The candidate bond maps to Δ = 1, V = 40, |w| = |t_2| = 0.0125, φ = π/2, μ_eff = −2.2. Correct.

1b. [derived] Second order. t_2 = −Δ²/(2V) from the single path (create (j + 1, j + 2) at cost 2V, remove (j, j + 1)), with product Δ² c†_{j+2}c_j and no string sign. Diagonal: −Δ²/(V m_b) per empty bond, giving the bulk shift +3Δ²/V, V_2 = −Δ²/V and V_3 = −Δ²/(3V), V_4 = 0. All as in C173 (P1.1).

1c. [numerical, sw_exact.out, sw_ring.out] I extracted the all-orders effective Hamiltonian on P from the full chain by the des Cloizeaux construction (the canonical Schrieffer-Wolff Hamiltonian), on an open chain of L = 11 and a ring of L = 12. At V/Δ = 40: t_2, the shift, V_2 and V_3 agree to 0.04%, 0.05%, 0.13% and 0.4%, and the deviations fall as (Δ/V)². V_4 = 2.7e-5 at V/Δ = 40, relative ε² × 3.4. This is stronger than C173's spectral check, since it tests each coupling separately. Correct.

1d. [derived, numerical] ΔN = 2. The second-order piece is −(2Δ Re w/V) a†_{b+2}a†_b (spin language), proportional to J_x + J_y. The third-order piece is +(2Δ³/(3V²)) a†_{b+3}a†_b, of relative size (4/3)Δ/V. Measured: the separation-3 element agrees to 0.07% at V/Δ = 40, and the separation-2 element agrees to 0.1% at Re w ≠ 0 and is below 1e-14 at Re w = 0. So C173's statement is correct. Refinement [proved, P2]: on J_x + J_y = 0 the symmetry U makes every even-separation pair creation U-odd, so (1.1) vanishes to all orders there. This also explains the C153 limit: at w = 0 only A–B (odd-separation) pairs can be created, which is C153's conversion vertex.

1e. [derived] "Odd-order ΔN = 2 terms enter energies only at sixth order." As an energy, (Δ³/V²)² divided by an intra-P denominator of order Δ²/V is Δ⁴/V³, relative ε², the same as the fourth-order N-conserving corrections, so the V⁻³ error law C173 measured is unaffected. The sentence should say that, and it should add that for the infrared the amplitude (relative ε) is what matters, not the energy shift. Gap (wording).

1f. [derived] Omitted terms, listed with their order. N-conserving corrections relative ε². ΔN = 2 relative ε (only odd separations on the stratum). ΔN = 4 directly in H_eff only at relative ε⁴ (no path (2, 0) or (3, 1) of pair creations and annihilations lands in P), measured below 3e-10 at V/Δ = 20. Parity-odd perturbations of size λ act at first order inside P as ΔN = ±1 PXP-type terms of size λ (P1.5).

1g. [derived] Exact consistency check I added. At w = 0 the empty state stops being the ground state at μ_eff = −2. In the full chain that must be C153's factorization (Kurmann–Thomas–Müller) field h_f = √(4J_z² − Δ²), and expanding, 4J_z − 2h_f − 3Δ²/V = −Δ²/V = −2|t_2|. It is.

## 2. The structure of (3.2) (C173 3.3)

2a. [proved] Each leg alone is −2 S·S, ferromagnetic Heisenberg. At w = 0 and J_x + J_y = 0 the full chain conserves Q = Σ (−1)^j Z_j, i.e. N_A − N_B, exactly. Correct.

2b. [proved] The band ε(k) = 2|w| cos(k − φ) − 2 cos 2k, its minima near −s and π + s with s ≈ |w| sin φ/4, offset and curvature difference 4|w| cos φ: correct. At φ = π/2 exactly, the minima are at sin k = −|w|/4, so s = arcsin(|w|/4), ε_min = −2 − |w|²/4, ε'' = 8 − |w|²/2, and two minima exist iff |w| < 4.

2c. [proved, P2] U = I ∘ Π_j R_z(jπ + π/2) is a symmetry iff J_x + J_y = 0, for any D, J_z, h_z. Within inversion times z-rotations the only other possibility is plain inversion at D = 0. Correct.

2d. [derived, P5.1] Vertex list. Each listed conversion vertex does carry nonzero momentum (π + 2s, 4s, π, −2s, 2s). But the list is not closed under products, and on the stratum two composites are resonant at every density: the two-pair vertex O_4 = (ψ_a†ψ_b†)², with momentum 2(k_a + k_b) = 2π ≡ 0, and the valley backscattering O_s = cos(2√2 φ_s), with momentum 2π(ρ_a − ρ_b) = 0 for the balanced (U-symmetric) state. A general vertex has P = n_a k_a + n_b k_b + 2π(m_a ρ_a + m_b ρ_b), and on the stratum P ≡ 0 iff n_a = n_b is even and m_a = −m_b. So the infrared symmetry on the stratum is U(1)_s × Z_4 (charge), not U(1)_s × Z_2. Gap (fixable): PROOF 5.4 must treat these two vertices, and section 4 below does.

## 3. The dilute cusp (C173 4)

3a. [derived] 4(a) is correct.

3b. [numerical, twobody_pi2.out, twobody_check.out] 4(b) assumes no zero-energy resonance in the inter-valley channel. On the stratum there is a bound state. At φ = π/2 the lowest two-body level at K = π = k_a + k_b lies below 2ε_min by E_b = 0.0017, 0.0195, 0.1226 for |w| = 0.5, 1, 2 (L = 400 and 800 agree to five digits). C173's own twobody.py at φ = π/2 gives the same −0.01953. There is no bound state at the intra-valley momenta. The bound state is about 20 sites across at |w| = 1. So the premise of 4(b) is wrong at the candidate point. Off the stratum the valley offset 4|w| cos φ lifts the inter-valley pair above 2ε_a once it exceeds E_b, which is why C173's check at φ = 1.3 (offset 1.07) saw nothing. Wrong on the stratum, correct for δ ≫ E_b.

3c. [derived] Minor. In a two-valley band the relative motion at K = π has two low-energy points ±p_0, so the threshold analysis is a two-channel S-matrix in the valley labels, not the "even channel" of a single-minimum problem. The generic threshold law is then S → −1 in that 2 × 2 space, which is the impenetrable statement C173 wants. The language should be fixed, the conclusion off the stratum is unaffected.

3d. [derived, heuristic where marked] Consequences. With E_b > 0 the dilute limit on the stratum is a gas of inter-valley molecules at momentum π, with a pseudospin gap E_b/2 as n → 0, c = 1 [derived, assuming molecule-molecule repulsion]. At n ≈ 0.17 the Fermi scale exceeds E_b, a BCS-like regime where a 1D pseudospin gap is expected to be much smaller than E_b and the pair coherence length much longer than 20 sites [heuristic, by analogy with the attractive Gaudin–Yang gas]. A coherence length beyond 96 sites would read c ≈ 2 at C173's sizes. So the "wedge" picture of PROOF 4 (impenetrable valleys, a pseudospin XXZ chain with couplings of order n³, c = 2 iff easy-plane) does not describe the stratum. The relevant scales are the binding and the pseudospin gap, and the competitor C173 did not list is inter-valley pairing (c = 1).

3e. [numerical, manybody.out, krho_control.out] Ring ED of (3.2) shows that the valleys are far from impenetrable. The total-density Luttinger parameter is K_ρ ≈ 3.2 to 3.6 at n = 0.13 to 0.2 on the stratum (3.6 at L = 36, n = 1/6), against 1.3 to 1.8 for the one-valley control at φ = 1.3, with the estimator calibrated on hard rods to 1% (0.686 and 0.690 against the exact 0.694). Impenetrable valleys would give about the one-valley value. Two weakly coupled valleys give about twice it, which is what is seen. The parity gap L·Δ_pb at n = 1/6 reads 2.76 (L = 24) and 1.91 (L = 36), so no pseudospin gap is resolved at L ≤ 36, which is what a coherence length of several tens of sites or more predicts and what a gapless pseudospin predicts too (registered as R5, indicative only).

## 4. Numerics and the status of the candidate (C173 5)

4a. 5.1 [numerical]. The two-body statement is correct where checked (φ = 1.3) and wrong at the candidate (3b). Its use, "the dilute gas is a gas, not a pair fluid", does not hold at φ = π/2.

4b. 5.3 [derived]. "The pseudospin is not ferromagnetic, so it is easy-plane or antiferromagnetic" omits pairing. A paired (Luther–Emery-type) valley sector also gives a valley-balanced ground state at even N and an odd-N ground doublet with imbalance ±1, which is what valley_ed.out shows. The ring data do not discriminate. The gauge map c_j → i^j c_j to real NN hop |w| and NNN hop +1 is correct [proved].

4c. 5.4, parity-odd directions [derived]. The ΔN = 1 terms off the parity class are first order in the perturbation λ (PXP-type), not "of order Δ²/V and smaller". Their vertices are off resonance as C173 says (on the stratum P ≡ 0 needs n_a + n_b even). So this part of the openness argument is correct once the size statement is replaced by "for λ ≪ |t_2|".

4d. 5.4, the momentum-2π two-pair vertex O_4 [derived, numerical input]. Status at φ = π/2:
- Lowest order. It is absent from (3.2), which conserves N. In the qubit chain the dominant source is second order in the separation-3 ΔN = 2 vertex A_2 = (4/3)(Δ/V)|t_2|, each insertion carrying momentum π or ∓2s and fusing to P = 0. So y_4 ≈ c_4 (16/9)(Δ/V)²|t_2| = c_4 (16/9) Δ⁴/V³, with c_4 of order one, the same order as the first corrections to (3.2). The direct ΔN = 4 term in H_eff is only relative (Δ/V)⁴. The larger separation-2 contribution, ∝ (Δ Re w/V)², is removed by U.
- Dimension. On the stratum U decouples charge and pseudospin, and x(O_4) = 2/K_c = 4/K_ρ, in either the two-component phase or the paired phase (there it is the squared pair field, 1/K_p with K_p = K_ρ/4).
- Relevance. With K_ρ ≈ 3.2 to 3.6 from the rings, x_4 ≈ 1.2, so O_4 is relevant. It can lock a mode, as C153's resonant conversion vertex locks the density at w = 0. It pins θ_c, removing the charge mode (two-component phase → c = 1, paired phase → c = 0). The locking scale is Δ_4 ~ |t_2|(1.8 c_4 ε²)^{1/(2 − x_4)}, i.e. a correlation length of a few thousand sites at V/Δ = 40 and of order 100 at V/Δ = 8 [heuristic]. So it is invisible in effective-model DMRG and in the V/Δ = 40 qubit chain at L ≤ 256, and may be visible in C177's small-V/Δ qubit chains on the stratum.

4e. 5.4, openness [derived, conditional]. Off the stratum both O_4 and O_s carry momentum (2η and 2π(ρ_a − ρ_b)) and are cut off by a commensurate-incommensurate transition once their momentum exceeds their gap. So the openness argument survives, in the form "a c = 2 phase of (3.2) at a nonzero valley offset is open in B_2", but the specific candidate bond, which sits on the stratum, is the one point where both dangerous vertices are resonant. C173's evidence for c = 2 off the stratum is one L = 48 entropy at φ = 1.54 and the unconverged L ≤ 96 runs at φ = π/2 and 1.56. Gap (fixable): the claim should be restated as conditional on c = 2 at nonzero offset, and the candidate bond should be moved off J_x + J_y = 0 (by more than the C-IC slab, of order Δ_4/v in the momentum mismatch).

4f. [heuristic, P5.4] An alternative c = 2 mechanism the binding suggests. With an inter-valley bound state, the valley offset acts like a Zeeman field on an attractive two-component gas. In the 1D attractive Gaudin–Yang gas, a field above the spin gap gives a partially polarized two-component (FFLO-like) liquid with c = 2 on an open set. If the same holds here, the c = 2 region would sit at nonzero offset, between the pseudospin gap and the offset where the second valley empties, not on the stratum. That would explain why the reading persists to φ = 1.54 (offset 0.12) and is gone by 1.5 (offset 0.28). It is untested.

## 5. Hidden-symmetry check (C163 §4.2)

5a. [proved] Site unitaries. The in-plane exchange block has singular values Δ/2 ± D, unequal for D ≠ 0, so the bond is in neither the uniform U(1) class (rotation block) nor C153's staggered-charge class (reflection block), and a z-twist cannot remove D without making the pairing wind. The symmetric part of J has eigenvalues (0.5, −0.5, 10), none zero, so it is not on the XY-in-a-field free-fermion locus, and the Jordan-Wigner form keeps the quartic V n n with V = 40.

5b. [derived] D = 0 neighbour. Particle-hole on the odd sublattice (π-rotation about x on odd sites) maps the D = 0 bond to C153's model at (u, b) = (20, 39.95), g = b³(r − u) = 80.6, deep in C153's polarized gapped phase (liquid needs g below about 1). In effective terms μ_eff = −2.2 lies below the w = 0 band bottom −2.

5c. [proved] Does D break the structure at the order that matters? Yes. Under the same map D becomes a bond-alternating ΔQ = ±2 term, and inside P it is the first-order hop with |w| = 2D/|t_2| = 1, equal to the leg hop t_2 itself. Q is broken at the leading order of the effective model.

5d. [numerical] Level statistics. (3.2) at φ = π/2, |w| = 1, rings L = 30, N = 5, four resolved momentum sectors of dimension 2125: <r> = 0.535 (GOE 0.531, Poisson 0.386). Controls: φ = 1.3 gives 0.524, |w| = 2 gives 0.521, and w = 0 with N_A − N_B conserved but unresolved gives 0.406, so the test detects a hidden conservation law when one exists (levelstat_L30N5.out). The full qubit ring in its lowest band (L = 16) gives 0.509 ± 0.007 at V/Δ = 6 and 0.493 ± 0.008 at V/Δ = 12, against 0.407 ± 0.008 for the D = 0 control with C153's hidden Q (full_lowband_L16.out).

5e. Conclusion [derived, numerical]. The candidate is not equivalent to a free, integrable or decoupled model. Its only extra structures are U on J_x + J_y = 0 and the approximate U(1)_N of the cusp (broken at relative order Δ/V). Neither makes it solvable, and both enter the infrared through the two resonant vertices of 2d. The candidate passes the hidden-symmetry check, but it fails something more mundane: on its own stratum the valleys bind (3b), and the two-pair vertex appears relevant in the qubit chain (4d).

## 6. Overall verdict

C173's Schrieffer-Wolff model is right in every sign and coefficient, and the lattice symmetry is right. The model has no hidden integrable, free or decoupled structure. What fails is the dilute-limit analysis at the candidate itself. On J_x + J_y = 0 the two valleys form an inter-valley bound state (E_b = 0.0195|t_2| at |w| = 1), which C173's two-body check at φ = 1.3 could not see. That makes inter-valley pairing (c = 1) the leading competitor on the stratum, ahead of the easy-plane versus antiferromagnetic question C173 posed, and it makes the total-density stiffness large (K_ρ ≈ 3.3 at n ≈ 1/6 on rings), so that in the qubit chain the momentum-2π two-pair vertex, generated at order Δ⁴/V³, appears relevant on the stratum and removes a further mode there. Neither effect is visible at the sizes C173 reached, so its c ≈ 2 readings are not contradicted, only reinterpreted as a possibly long crossover.

Concrete repairs:
- (a) PROOF 3.1/3.2: add the U selection rule, the explicit separation-3 vertex (1.2), and the correct order counting for energies versus amplitudes.
- (b) PROOF 3.3 and 5.4: add O_4 and O_s to the vertex list, with x(O_4) = 4/K_ρ and x(O_s) = 2K_s, and state that both are resonant only on the stratum.
- (c) PROOF 4 and 5.1: restrict the impenetrability argument to valley offsets well above E_b, and add inter-valley pairing to the list of competitors.
- (d) PROOF 5.4: restate openness as conditional on c = 2 at nonzero valley offset, and move the candidate bond off J_x + J_y = 0.
- (e) Tests to add for C176 to C178: odd-N parity gap and pair correlations at momentum π on long open chains at φ = π/2, the same at φ = 1.54 to 1.52 to test the partially polarized (FFLO-like) alternative, and, for the qubit chain at small V/Δ, the difference between J_x + J_y = 0 and a small nonzero J_x + J_y (O_4 on and off resonance).
