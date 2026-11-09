This work was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C180 result: field theory and literature for the Fibonacci-cusp c = 2 candidate

Derivations are in PROOF.md and the log in PROGRESS.md. One script, optable.py (output optable.out, seconds of CPU), checks the band facts, enumerates the vertex operators and evaluates the scaling dimensions. No Typhon job.

## Verdict

The field theory gives no kinematic or symmetry obstruction to c = 2. At the Σ point (J_x + J_y = 0) the N-conserving two-valley liquid has exactly one nonlinear operator at generic parameters, the inter-valley backscattering cos 2√2φ_−. In the parity class there is one more, a ΔN = 4 vertex that is irrelevant in the repulsive regime. Off Σ there is none. So c = 2 is stable on an open set of B₂ exactly when the valley pseudospin is easy-plane, and in the dilute limit that is a statement about three two-body excluded lengths.

The literature does not settle the sign, but it narrows it. The closest published model, hard-core bosons on the frustrated triangular ladder at t = |t′| (Mishra, Greschner and Santos, arXiv:1411.7270), is the candidate's effective model without the two generated attractions V₂ and V₃. It is two-component (c = 2) at low density, both in its dilute analysis (on the c = 2 side for nearest-neighbour repulsion above about 3.6t) and in DMRG at large repulsion (c = 2 for ρ < 0.4). What is left is whether the generated attractions −2n_jn_{j+2} − (2/3)n_jn_{j+3} flip the pseudospin to ferromagnetic. A Born-level sign indicator, which is not reliable in the hard-core regime, says they push in that direction. C175's exact two-body computation decides it.

No c = 2 phase of a translation-invariant nearest-neighbour spin-1/2 chain is claimed in the literature I could check, and none is excluded. The nearest-neighbour studies closest to this corner (the Rydberg-constrained Kitaev chain, the interacting Kitaev chain with repulsion) report c = 1, but none enters the cusp. If confirmed, the candidate violates Latorre and Sierra's c ≤ d − 1 at d = 2 (c_T = 2), within their own definition, which explicitly allows several CFTs with different velocities.

## Key claims

1. [proved] Symmetries. Θ = I·K is an antiunitary symmetry of every parity-class bond and preserves each valley. Σ = I × Π_odd R_z(π) × Π R_z(π/2) maps (J_x, J_y, D) → (−J_y, −J_x, D), so it is exact iff J_x + J_y = 0, and it exchanges the valleys (k → π − k). At φ = π/2, c_j = i^j d_j makes (3.2) real: a J1–J2 chain (J1/J2 = |w|) with the blockade, J2^z = −J2 and a third-neighbour attraction. Valleys k_{a,b} = π/2 ∓ k₀ with cos k₀ = −|w|/4, ε_min = −2 − |w|²/4, 1/m = 8 − |w|²/2 (checked numerically, optable.out).
2. [derived] Bosonized form (PROOF 1). Fields φ_±, θ_± from the two valleys. Σ allows (∂φ_±)², (∂θ_±)² and two drift terms ∂φ_+∂θ_−, ∂θ_+∂φ_− (from the band asymmetry ε‴ ≠ 0 about each valley). It forbids ∂φ_+∂φ_− and ∂θ_+∂θ_−. Dimensions without drift: Δ = [(p_a+p_b)²K_+ + (p_a−p_b)²K_− + (n_a+n_b)²/K_+ + (n_a−n_b)²/K_−]/8. With drift, Δ = vᵀGv/4 from the Gaussian covariance G (checked against the closed form, optable.out).
3. [derived, enumeration checked numerically] Operator table (PROOF 2). The momentum is P = n_a k_a + n_b k_b + π(p_aρ_a + p_bρ_b). On an open set at Σ only n_a = n_b ∈ 2ℤ, p_a = −p_b ∈ 2ℤ survive. Off Σ nothing survives.

   | operator | momentum | dimension | allowed | relevant if |
   | :-- | :-- | :-- | :-- | :-- |
   | cos 2√2φ_− (inter-valley backscattering, g_1⊥) | 2π(ρ_a − ρ_b), zero when balanced | 2K_− | Σ point, all classes | K_− < 1 |
   | e^{2√2iθ_+} (two ab pairs, ΔN = 4) | 2(k_a + k_b) = 2π | 2/K_+ | Σ point, parity class only | K_+ > 1 |
   | ∂φ_− (valley offset δ) | 0 | 1 | off Σ | always (a pseudospin field) |
   | ∂φ_+∂φ_−, ∂θ_+∂θ_− (mass imbalance) | 0 | 2 | off Σ | marginal |
   | resonant only: valley transfer, pair transfer, ΔN = 1, 2, commensurate locking | zero on codimension-one sets | see PROOF 2.2 | on hypersurfaces | tongues of c ≤ 1 |

4. [derived] Stability conditions for c = 2 at Σ (PROOF 3):
   - (S1) The quadratic form is positive. The binding condition is the true inverse valley susceptibility, u_−/K_− − λ²/(u_+K_+) > 0. Its vanishing is the vector-chiral instability, with ⟨∂θ_+⟩ ≠ 0 and valley imbalance appearing together.
   - (S2) K_− > 1, or K_− = 1 on the marginally irrelevant side.
   - (S3) K_+ < 1 in the parity class.
   - (S4) Away from the resonance tongues.

   Off Σ, c = 2 persists while |δ| < h_sat = J_⊥(1 + Δ_σ) (order n³) and |η| ≲ n^{1/2}. This reproduces C173's wedge from the field-theory side. Off the parity class, ΔN = 1 vertices are resonant only, so the region is open in B₂ if (S1) to (S4) hold with margin. What protects the phase is incommensurability (emergent U(1)_a × U(1)_b from translation), not a symmetry.
5. [derived] Dilute limit (PROOF 4). The impenetrable two-valley gas maps to a pseudospin XXZ chain on the particle index with J_⊥ = (2π²n³/3m)|A_o|, J_z = (2π²n³/3m)(a_s − A_d). Here a_s is the same-valley excluded length, and [[A_d, A_o], [A_o*, A_d]] is the inter-valley excluded-length matrix on the ordered channels |ab⟩, |ba⟩ (equal diagonals by Θ). Then K_− = 2K_σ(Δ_σ), Δ_σ = (a_s − A_d)/|A_o|, and K_+ → 1/2.
   - c = 2 iff a_− < a_s < a_+ (a_± = A_d ± |A_o|).
   - Vector chiral if a_s < a_−. This is KHGV's g₁₁ = g₁₂ boundary.
   - Valley crystal if a_s > a_+. This side is not contained in KHGV's two-coupling criterion.

   Check: the δ-gas gives the standard antiferromagnetic Yang–Gaudin exchange 4π²n³/(3m²g). The drift is of order n² and subleading as n → 0, but at n = 0.17 its correction to the valley susceptibility is of the same order as the susceptibility itself [heuristic estimate].
6. [derived] Competing c = 1 phases (PROOF 5):
   - vector chiral (valley ferromagnet, Σ broken, twofold);
   - valley-crystal liquid (φ_− locked, unique ground state, a Luther–Emery analogue whose dominant correlation is the valley density wave for K_+ < 1 and the ab pair for K_+ > 1);
   - bound-state liquids (ab pairs give the previous phase, aa pairs give a two-valley pair gas that repeats the question);
   - charge-locked commensurate phases;
   - ΔN = 4 θ_+ locking (needs K_+ > 1).
7. [derived] A resonance line inside the candidate's window (PROOF 2.2). At |w| = 1 the 4k_F-dressed pair transfer (aa ↔ bb with p = ±2) is momentum-zero at n = 1.011/(2π) ≈ 0.161, close to the candidate density n ≈ 0.17. Its dimension is 2/K_− + 2K_+, relevant only if K_− > 1/(1 − K_+) ≈ 1.5. If the pseudospin is strongly easy-plane, a c = 1 tongue sits near n ≈ 0.16.
8. [literature check] C173's second-order coefficients agree with Verresen–Vishwanath–Pollmann's independent 1/U expansion (arXiv:1903.09179, appendix C) at t = Δ: next-nearest hop −Δ²/(2V), n_jn_{j+3} term −Δ²/(3V), and the n_jn_{j+2} term after removing the −2t²/U piece that vanishes at the cusp.

## What the literature says (all arXiv numbers below were opened and checked; "abstract" means only the abstract page)

- Two-valley magnon gases near saturation. Kolezhuk, Heidrich-Meisner, Greschner, Vekua, arXiv:1112.5720 (read): two-component Bose gas theory with couplings from the exact two-magnon problem, TLL2 vs vector chiral at g₁₂ = g₁₁, isotropic antiferromagnetic spin-1/2 boundary predicted at J1/J2 = 1.682 (VC below, TLL2 above), with an unexplained quantitative disagreement with spin-1/2 DMRG on the antiferromagnetic side. Hikihara, Momoi, Furusaki, Kawamura, arXiv:1004.0622 (read): TLL2 (c = 1 + 1) and vector chiral phases in the antiferromagnetic J1–J2 chain, with each vector chiral region adjacent to a TLL2 region, and the backscattering g_1⊥ whose irrelevance defines TLL2. Hikihara, Kecke, Momoi, Furusaki, arXiv:0807.0858 (abstract): ferromagnetic J1, vector chiral and multipolar bound-magnon phases. Okunishi, Tonegawa, cond-mat/0211458 (abstract): cusps and 1/3 plateau.
- Constrained (blockaded) frustrated hopping, the closest calibration. Mishra, Greschner, Santos, arXiv:1411.7270 (read): hard-core bosons with t, t′ = −t, nearest-neighbour V and three-body W. Dilute theory puts V = ∞ on the c = 2 side (critical V ≈ 3.6t). The blockaded squeezed model (their model 8) shows SF and 2SF, with "no region with finite chirality". DMRG at V = W = 8 gives c = 2 for ρ < 0.4, with the caveat that open-chain c fits at low filling drift upward. Related, abstract only: Greschner–Santos–Vekua arXiv:1202.5386, Mishra–Pai–Mukerjee arXiv:1309.4412.
- Nearest-neighbour B₂ chains near the corner. Verresen, Vishwanath, Pollmann, arXiv:1903.09179 (read): Rydberg-constrained Kitaev chain, emergent U(1), c = 1, K = (1 − ρ)², at t = Δ, real hopping, no cusp. Chepiga, Mila, arXiv:2206.11754 (read in part): interacting Kitaev chain with repulsion, floating phases and an eight-vertex line, c ≈ 1, real hopping only, so the Σ degeneracy (which needs DM) is absent. C156's list (XYZ in fields, with DM, Kitaev–Gamma family) is all c ≤ 1 and is not repeated.
- Integrability. Corcoran, de Leeuw, Pozsgay, arXiv:2405.15848 (read in part): reflection-symmetric integrable Rydberg-constrained chains at range 3 (RSOS) and range 4 (constrained XXZ, double golden chain), with a partial range-5 list. The candidate's real-gauge model (0.2) is range 5 in their counting and matches none of the listed models. No exact solution is known, though their range-5 list is partial.
- Other constrained models, abstract only: Lesanovsky–Katsura arXiv:1204.0903 (golden chain physics, c < 1), Chen–Iadecola arXiv:2103.09773 (PXP in a staggered field, two-site cell), Sohal–Bernien–Verresen arXiv:2610.40771 (blockaded ladders, c = 2 only with three legs), Keselman–Berg arXiv:1502.02037 (spin-gapped c = 1 gapless SPT, spinful), Ruhman–Altman arXiv:1705.00635 (paired spinless fermions).
- Latorre and Sierra, arXiv:2403.17242 (read). c_T = Σ_a c_a ≤ d − 1 for one-site-translation-invariant nearest-neighbour chains whose low-energy spectrum is "a collection of unitary CFTs". Different velocities are explicitly allowed. At d = 2 they state that the Heisenberg chain saturates the bound and that no spin model gives more.

## Known or excluded?

Neither. Every published c = 2 two-component phase in spin-1/2 chains has range ≥ 2 or is a ladder. The nearest-neighbour studies near this corner stay in the strong-rung regime or have real hopping. No theorem reaches a neighbourhood of the Fibonacci point. If confirmed, the candidate is new, and it lies inside Latorre–Sierra's hypotheses with c_T = 2 > d − 1 = 1. The project's replacement bound Σ c_a v_a ≤ π M log₂ d is not threatened (A is of order n|t₂| against M of order J_z).

## Predictions registered before any C175 to C179 output exists (those folders were absent when this was written)

These are what this field theory says the other agents' data should show if the c = 2 reading is right.
- (R1) Ring spectrum at φ = π/2, |w| = 1, n ≈ 0.17. There are two velocities, u_+ ≈ πn/m ≈ 4.0|t₂| (within about 30 percent, interaction corrections of order n·a), and u_− well below u_+.
- (R2) The lowest state with N_a − N_b changed by 2 at fixed N lies at π(u_−/K_−)_eff/L, scaling as 1/L. The subscript means the drift-corrected susceptibility of claim 4. A saturating gap means the valley-crystal phase. Valley-imbalanced levels falling below the balanced ground state mean vector chiral.
- (R3) K_+ between 0.3 and 0.5 at n ≈ 0.17.
- (R4) n(k) has two equal-weight valley peaks at k = π/2 ∓ k₀.
- (R5) If c = 2 holds with K_− > 1.5, a c = 1 dip should appear in a density scan near n ≈ 0.161 at |w| = 1 (claim 7). Its absence means either K_− < 1.5 or a small bare amplitude.
- (R6) The exact two-body lengths should satisfy a_− < a_s < a_+ (c = 2). As a calibration, the same computation with V₂ = V₃ = 0 must give a_− < a_s (two-component, as MGS and KHGV find).

## What is left open

1. The sign of the pseudospin anisotropy with the generated attractions included. This is the deciding input, C175's two-body problem. The literature covers J2^z/J2 ∈ [0, 1] (KHGV) and J2^z = 0 with the blockade (MGS), not J2^z/J2 = −1 plus a third-neighbour attraction.
2. Bound states and near-threshold resonances at φ = π/2 (C173 saw a two-body level 0.002 to 0.008 above threshold at φ = 1.3). A large |a| would break the dilute expansion at n = 0.17.
3. The drift correction to the valley susceptibility at n = 0.17, which is not small by the estimate here, and the possible pair-transfer tongue near n ≈ 0.161.
4. The constants in the Σ action, which fix the sign of cos 2√2φ_− and matter only at the marginal point Δ_σ = 1.
5. Whether the pseudospin-gapped competitor has a gapless-SPT refinement.
