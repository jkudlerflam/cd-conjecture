This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C179 derivations: referee and hidden-symmetry check of C173 sections 3 to 5

Tags [proved], [derived], [numerical], [heuristic], [assumed] open the paragraph they apply to. Every numerical statement names its script and output in this folder. Equation numbers (0.1), (3.1), (3.2) refer to codex/work_C173/PROOF.md. Units |t_2| = Δ²/(2V) = 1 inside the effective model, and ε = Δ/V.

## 0. Conventions, checked

[proved] With the bond h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z), write σ^± = (X ± iY)/2. Then J_x XX + J_y YY = (J_x + J_y)(σ⁺_j σ⁻_{j+1} + h.c.) + (J_x − J_y)(σ⁺σ⁺ + σ⁻σ⁻), and X_jY_{j+1} − Y_jX_{j+1} = 2i(σ⁺_jσ⁻_{j+1} − σ⁻_jσ⁺_{j+1}). Take the particle to be a down spin, n = (1 − Z)/2, a† = σ⁻, and c_j = (−1)^{Σ_{l<j} n_l} a_j. Then c†_{j+1}c_j = a†_{j+1}a_j = σ⁺_jσ⁻_{j+1}, c_jc_{j+1} = −a_ja_{j+1}, J_z ZZ = J_z(4n_jn_{j+1} − 2n_j − 2n_{j+1} + 1) and h_zZ_j = h_z(1 − 2n_j). This reproduces C173's (0.1) with w = J_x + J_y + 2iD, Δ = J_x − J_y, V = 4J_z and μ = 4J_z + 2h_z (bulk sites, each in two bonds). In the spin (hard-core boson) language that I use for all numerics, the pairing is +Δ(a_ja_{j+1} + h.c.), with no signs.

[proved] The candidate bond J_x = 0.5, J_y = −0.5, D = 0.00625, J_z = 10, h_z = −19.97625 has Δ = 1, V = 40, w = 0.0125 i, t_2 = −0.0125, so |w|/|t_2| = 1 and φ = π/2. The chemical potential of (3.2) is μ_eff = μ − 3Δ²/V = 0.0475 − 0.075 = −0.0275 = −2.2|t_2|. All of C173's numbers for the candidate are consistent.

## 1. The second-order model and everything it omits

### 1.1 The N-conserving terms

[derived] Unperturbed H_0 = V Σ n_jn_{j+1}, with the blockaded space P at energy 0 and a configuration with m adjacent occupied pairs at energy mV. The perturbation is the hop w, the pairing Δ and −μn. The hop and −μn act inside P at first order. The pairing always leaves P (it creates an adjacent pair, or needs one to annihilate). At second order:

- NNN hop. The only path from a particle at j to j + 2 creates the pair (j + 1, j + 2) (three in a row, cost 2V) and then removes (j, j + 1). In fermions the product of matrix elements is (−Δ c_jc_{j+1})(−Δ c†_{j+2}c†_{j+1}) = Δ² c_j c_{j+1}c†_{j+2}c†_{j+1} = Δ² c†_{j+2}c_j on states with j + 1 empty. So t_2 = −Δ²/(2V), as in (3.1). The middle site is empty by the blockade, so there is no Jordan-Wigner sign, and the amplitude does not depend on particles at distance ≥ 2 from the path (the intermediate cost stays 2V). A particle at j + 3 makes the final state leave P, so the hop is blockaded there. No other N-conserving off-diagonal second-order path exists.
- Diagonal. Every empty bond (b, b + 1) contributes −Δ²/(V m_b), with m_b the number of adjacent pairs after creating the pair at b. A particle blocks two bonds (+2Δ²/V) and raises its two next bonds from V to 2V (+Δ²/(2V) each), a bulk shift of +3Δ²/V. Distance 2 shares a blocked-and-raised bond on each side, an overcount of Δ²/(2V) twice, so V_2 = −Δ²/V = 2t_2. Distance 3 has one bond at 3V instead of the additive estimate, −Δ²/V + 2·Δ²/(2V) = 0, so V_3 = −Δ²/(3V) = (2/3)t_2. Distance 4 shares nothing, V_4 = 0 at this order.

All signs and magnitudes in C173 (3.1) and (3.2) are correct.

[derived] Consistency check against C153. At w = 0 the band of (3.2) is −2 cos 2k with bottom −2, so the empty state stops being the ground state at μ_eff = −2. In the full chain at D = 0 and J_x + J_y = 0 this should be C153's factorization line u = r, whose field is h_f = √(4J_z² − Δ²) (the Kurmann–Thomas–Müller field of an XYZ chain with J_x = −J_y). Expanding, μ = 4J_z − 2h_f = 2Δ²/V + O(Δ⁴/V³), so μ_eff = μ − 3Δ²/V = −Δ²/V = −2|t_2|. This tests t_2 and the shift 3Δ²/V together against an exact result.

### 1.2 The ΔN = 2 terms

[derived] At second order the only path from P to P with N + 2 is pair creation followed by a hop that separates the pair. From a local vacuum: create (b, b + 1) at cost V, hop b + 1 → b + 2 with amplitude w; or create (b + 1, b + 2), hop b + 1 → b with amplitude w*. In spin language the sum is

    −(Δ/V)(w + w*) a†_{b+2}a†_b = −(2Δ Re w/V) a†_{b+2}a†_b,                  (1.1)

i.e. +(2Δ Re w/V) c†_{b+2}c†_b in C173's fermions. It is proportional to J_x + J_y and vanishes on the candidate's stratum.

[derived] At third order, with three pairings, the only paths into P from a local vacuum create (b, b + 1) and (b + 2, b + 3) in either order (costs V, then 3V) and then remove (b + 1, b + 2). Both orders contribute with the same sign (disjoint pair creators commute), so

    + (2Δ³/(3V²)) a†_{b+3}a†_b  =  −(2Δ³/(3V²)) c†_{b+3}c†_b,               (1.2)

of relative size (4/3)ε compared with |t_2|. Paths using w or μ at third order are O(Δ⁵/V³) or smaller in the cusp. A ΔN = 2 path cannot start next to an existing particle in a way that avoids these two shapes, so the separation-2 and separation-3 vertices, modified near other particles, are the complete list at this order.

[derived] C173's sentence "ΔN = 2 inside P appears only at order wΔ/V and Δ³/V²" is correct. Two refinements. (i) The wΔ/V piece is (1.1), proportional to J_x + J_y, and on J_x + J_y = 0 it vanishes to all orders (section 2). (ii) C173 says the odd-order ΔN = 2 terms enter energies "only at sixth order". As an energy shift, (1.2) enters at second order through denominators of order |t_2| inside P, giving (Δ³/V²)²/(Δ²/V) = Δ⁴/V³, the same order as the fourth-order N-conserving corrections. That is consistent with the V⁻³ error law C173 measured, so nothing numerical changes, but the right counting is "relative (Δ/V)² in energies, relative Δ/V as an amplitude". For the infrared it is the amplitude that matters (section 5).

### 1.3 Higher orders and the ΔN = 4 term

[derived] N-conserving corrections to (3.2) are O(Δ⁴/V³), i.e. relative ε², including a nonzero V_4 and corrections to t_2, V_2, V_3. A direct ΔN = 4 term inside P needs net two pair creations. With pairings only, (k creations, l annihilations) = (2, 0) and (3, 1) cannot land in P (every arrangement leaves an adjacent pair), so the first connected ΔN = 4 term is (4, 2), O(Δ⁶/V⁵) = |t_2| O(ε⁴). Paths with two w's and two pairings are also O(Δ⁶/V⁵) in the cusp.

### 1.4 Numerical check by exact extraction

[numerical, sw_exact.py → sw_exact.out, open chain L = 11; sw_ring.py → sw_ring.out, ring L = 12; sw_dn4.py → sw_dn4.out] I diagonalize the full qubit chain and build the all-orders effective Hamiltonian on P by the des Cloizeaux (direct-rotation) formula H_eff = U E U†, U the unitary polar part of the overlap ⟨P|low eigenvectors⟩. This is the canonical Schrieffer-Wolff effective Hamiltonian (Bravyi, DiVincenzo and Loss, Ann. Phys. 326, 2793 (2011)). It is independent of C173's blockade.py, and it uses the Pauli bond directly, so it also tests section 0. Ratios measured/predicted at w = i|t_2|, μ_eff = −2.2:

| V/Δ | t_2 | 3Δ²/V shift | V_2 | V_3 | (1.2), separation 3 | (1.1), separation 2 |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 20 | 0.9985 | 0.9981 | 0.9948 | 0.9842 | 0.9974 | 0 (3e-16) |
| 40 | 0.9996 | 0.9995 | 0.9987 | 0.9960 | 0.9993 | 0 (9e-15) |
| 80 | 0.9999 | 0.9999 | 0.9997 | 0.9990 | 0.9998 | 0 (1e-14) |
| 160 | 1.0000 | 1.0000 | 0.9999 | 0.9998 | 1.0000 | 0 (3e-15) |

The deviations fall by 4 per doubling of V, i.e. relative ε². V_4 is 2.7e-5 at V/Δ = 40, which is ε² × 3.4 in units of |t_2|. At w = (0.6 + 0.8i)|t_2| the separation-2 element is −1.492e-3, −3.745e-4, −9.372e-5, −2.344e-5 at V/Δ = 20, 40, 80, 160 against −2ΔRe w/V = −1.5e-3, −3.75e-4, −9.375e-5, −2.344e-5. On the ring the separation-2 element from the vacuum is below 3e-15 at every b when Re w = 0, the separation-3 element is translation-invariant to 5e-15 and equals 0.9974 and 0.9993 of (1.2) at V/Δ = 20 and 40, and the largest ΔN = 4 element is 2.4e-10 and 1.9e-12 at V/Δ = 20 and 40 (sw_ring.out). That falls like V⁻⁷, faster than the Δ⁶/V⁵ bound of 1.3, so the leading connected ΔN = 4 terms apparently cancel. Either way the direct term is negligible. The pre-registered criterion R1 is passed.

[numerical, sw_dn4.out] On the open chain, ΔN = 2 elements whose position sum changes by an even number are not exactly zero (5.0e-6, 3.1e-7, 2.0e-8 at V/Δ = 20, 40, 80, ∝ V⁻⁴). These are edge terms. The open chain is inversion-symmetric about its centre, so U relates the element at b to the mirrored element with a sign (it forbids only the uniform bulk combination), and on the ring the bulk element is zero. This is a reminder for C176 and C177: an open chain on the stratum carries U-odd edge terms of relative size ε³.

### 1.5 Off the parity class

[derived] A parity-odd perturbation of size λ (a transverse field B_x, B_y, or J_xz-type couplings) acts inside P at first order. For example B_x X_j restricted to P is B_x P(a_j + a†_j)P with the blockade projectors, the PXP term, and X_jZ_{j+1} restricted to P acts the same way with a sign set by n_{j+1}, which is 0 whenever site j is occupied in P. So the leading parity-odd terms are ΔN = ±1 terms of size λ, not of size Δ²/V. C173 PROOF 5.4 writes "ΔN = 1 terms of order Δ²/V (and smaller)", which is right only if λ is itself taken of order Δ²/V. The correct statement is that the controlled model survives for λ ≪ |t_2|, so the open set has a width of order Δ²/V in the parity-odd directions. That is still open. Second-order mixed terms λΔ/V carry ΔN = ±1 and ±3.

## 2. The lattice symmetry on J_x + J_y = 0

[proved] Rotations by α_j about z send the hop amplitude to w e^{i(α_j − α_{j+1})} and the pairing to Δ e^{i(α_j + α_{j+1})} (for one sign convention of the phases, the other gives the same conditions). Site inversion j → −j sends w → w* and leaves Δ, J_z and h_z unchanged. For U = I ∘ Π_j R_z(α_j) to be a symmetry of the translation-invariant chain, α_j − α_{j+1} = −β must be independent of j and α_j + α_{j+1} ≡ 0 mod 2π for every j. Subtracting consecutive conditions gives 2β ≡ 0, so β ∈ {0, π}. For β = 0, α_j = α_0 ∈ {0, π} and the hop condition is w* = w, i.e. D = 0 (plain inversion, the D = 0 case). For β = π, α_j = jπ ± π/2 and the hop condition is −w* = w, i.e. Re w = J_x + J_y = 0, with D, J_z, h_z arbitrary. So within inversion times z-rotations the candidate's only extra symmetry is U = I ∘ Π_j R_z(jπ + π/2), present exactly on J_x + J_y = 0. In momentum space U maps k → π − k (inversion k → −k, the staggered π-rotation k → k + π) and multiplies each particle by a phase ±i. This confirms C173 PROOF 3.3.

[proved] Selection rule. Under U the uniform ΔN = 2 operator Σ_b a†_{b+r}a†_b maps to itself (inversion maps b → −b − r and hard-core boson creators commute) times e^{−i(α_b + α_{b+r})} = e^{−i((2b + r)π + π)} = (−1)^{r+1}. So even separations are U-odd and are forbidden in a U-symmetric translation-invariant chain to all orders, and odd separations are allowed. This is why (1.1) vanishes identically at Re w = 0 while (1.2) survives, and it is what sw_ring.out shows at machine precision.

## 3. Hidden-symmetry check (C163 §4.2)

The question is whether the candidate bond is equivalent, by an exact site unitary, a duality or a known map, to a free, integrable or decoupled model, and whether D breaks such a structure at the order that matters.

3.1 [proved] Site unitaries. The field and the dominant coupling J_z fix the z axis (a rotation tilting it produces an O(J_z) transverse field term that is absent), so the relevant uniform site unitaries are rotations about z, possibly with a π-rotation about an in-plane axis. Under these the in-plane block M = [[J_x, D], [−D, J_y]] transforms as M → R^T M R', and its singular values are invariant for uniform maps. They are Δ/2 + D and Δ/2 − D (M^TM = [[1/4 + D², D], [D, 1/4 + D²]]), unequal for D ≠ 0. The uniform U(1) class (XXZ + z-DM + z-field) needs M ∝ rotation, and C153's staggered-charge class needs M ∝ reflection, both with equal singular values. So the candidate is in neither class [proved]. A uniform twist α per site sends w → w e^{iα} and Δ → Δ e^{i(2j+1)α}, and keeping the pairing uniform requires α ∈ {0, π}, so the DM term cannot be twisted away either [proved].

3.2 [proved] Free fermions. In the Jordan-Wigner frame the bond has the quartic term V n_jn_{j+1} with V = 40, the largest coupling. The symmetric part of J has eigenvalues (0.5, −0.5, 10), none zero, so no rotation brings the bond into the XY-in-a-field free-fermion form (one exchange eigenvalue zero along the field). The candidate is not on the free-fermion locus.

3.3 [derived] The hidden U(1) at D = 0. A π-rotation about x on odd sites (particle-hole on the odd sublattice) maps XX − YY → XX + YY, ZZ → −ZZ and the uniform field to a staggered one, so at D = 0 the bond is C153's staggered-charge model with Q = Σ_j (−1)^j Z_j conserved. In C153's units (H = Σ XX − YY + uZZ + bZ) the candidate at D = 0 sits at u = 20, b = 39.9525, r = √(1 + b²/4) = 20.00126, so g = b³(r − u) = 80.6. C153's liquid strip is g < g_c ≈ 0.6 to 1, so the D = 0 neighbour is deep in the polarized, gapped phase. In effective-model terms the w = 0 band bottom is −2, above μ_eff = −2.2, so the band is empty. The gas exists only because of D.

3.4 [proved] D breaks Q at the order that matters. The same odd-sublattice rotation maps D(XY − YX) to the bond-alternating term ∓D(XY + YX), which changes Q by ±2. Inside P, D is the first-order hop w = 2iD, and |w|/|t_2| = 2D/|t_2| = 1. So Q is broken at the leading order of the effective model, by an amount equal to the leg hop t_2 itself, not perturbatively. The w = 0 structure of C173 PROOF 3.3 (two ferromagnetic legs, N_A − N_B conserved, C153's resonant conversion vertex) is an O(1) distance away.

3.5 [numerical, levelstat.py → levelstat_L30N5.out] Hidden integrability of (3.2). On rings of L = 30 with N = 5 (n = 1/6), hard-core bosons, resolved momentum sectors that U does not map into themselves (dimension 2125 each), the mean adjacent-gap ratio is <r> = 0.539, 0.531, 0.529, 0.542 in four sectors, mean 0.535, against GOE 0.531 and Poisson 0.386. Each K sector carries the antiunitary inversion × conjugation, so GOE is the chaotic expectation. Controls: φ = 1.3 gives 0.524, |w| = 2 gives 0.521, and w = 0, where N_A − N_B is conserved but not resolved, gives 0.406. So the test does detect a hidden conservation law when there is one, and (3.2) at the candidate point has none. The code ringk.py is validated against C173's real-space ring_ed.py to 5e-14 (validate_ringk.out). Criterion R2 is passed for the effective model.

3.6 [numerical, full_levelstat.py, full_lowband.py] The full qubit ring. Over the whole spectrum the candidate family reads <r> = 0.489, 0.485, 0.462 at V/Δ = 4, 8, 40 (L = 14), against 0.411, 0.405, 0.399 for the D = 0 control (full_levelstat_L14.out). The whole spectrum is the wrong place to look, because at high energy D is a tiny perturbation of C153's Q-conserving model, and intermediate values are expected. In the lowest band the D term is O(1). The low band of the L = 16 ring (the Lucas(16) = 2207 blockaded states, separated from the next band by more than 55|t_2|) gives <r> = 0.509 ± 0.007 at V/Δ = 6 (mean Var(N) per eigenstate 0.93) and 0.493 ± 0.008 at V/Δ = 12 (Var(N) = 0.50), against 0.407 ± 0.008 for the D = 0 control with its hidden Q (full_lowband_L16.out). The candidate's values sit close to GOE, slightly below it, as expected from sectors of only about 70 states with an approximate N that is mixed at relative order Δ/V. The D = 0 control shows that a hidden U(1) would be seen. R2 is passed for the effective model and, at this resolution, for the full ring.

3.7 Duality and decoupling. [derived] The legs decouple (in the sense of separately conserved N_A, N_B) only at w = 0, which is 3.3. Bond-site and Kramers–Wannier type dualities of this chain are nonlocal maps that preserve level statistics, so a duality onto an integrable or free model would show up in 3.5. I found no such map, and the level statistics exclude one for (3.2).

Conclusion of the check. The candidate is not equivalent to a free, integrable or decoupled model, and the D term breaks C153's hidden U(1) at the leading order of the effective model. The only extra structures are the exact symmetry U on J_x + J_y = 0 (section 2) and the approximate U(1)_N of the cusp, broken at relative order ε by (1.2). Neither makes the model solvable. Both matter for the infrared, through the two resonant vertices of section 5.

## 4. The two-body problem on the stratum: the valleys bind

[numerical, twobody_pi2.py → twobody_pi2.out; twobody_check.py → twobody_check.out] C173 PROOF 5.1 states that no two-body bound state lies below 2ε_min, from twobody.py at φ = 1.3. At φ = π/2 this fails. On rings of L = 400 and 800 (identical to five digits), the lowest two-body level at total momentum K = π lies below 2ε_min by

    E_b = 0.0017, 0.0195, 0.1226  for |w| = 0.5, 1, 2.

K = π = k_a + k_b is the momentum of one particle in each valley. At the intra-valley momenta K = 2k_a, 2k_b there is no bound state (the lowest level is 1e-4 above the continuum at L = 800, a finite-size offset). C173's own twobody.py, run at φ = π/2, gives the same −0.01953 (twobody_check.out). The bound state's size is about 1/√(m E_b) ≈ 20 sites at |w| = 1 (m = 1/ε''(k_a) = 1/7.5).

[derived] Why φ = 1.3 missed it. At φ = 1.3 the valley offset is 4|w| cos φ ≈ 1.07, much larger than E_b, so the inter-valley pair sits above 2ε_a. The bound state is a property of the degenerate valleys.

[derived] Consequences for C173 PROOF 4. The argument there assumes that the even inter-valley channel has no zero-energy resonance, so that inter-valley pairs are impenetrable at low density. On the stratum the inter-valley channel is past resonance, with a shallow bound state. So:
- (b) and (c) of PROOF 4 do not apply at φ = π/2. In the dilute limit, μ just above E_2/2 = ε_min − E_b/2, particles enter as inter-valley molecules at momentum π. Single particles and any valley imbalance cost E_b/2. That is a one-component liquid of molecules, c = 1, with a pseudospin (valley) gap [derived from the two-body data, with the assumption that the bosonic molecules repel each other, which is the generic case and is what hard-core molecules do].
- The dichotomy "easy-plane pseudospin (c = 2) or a polarized or ordered valley state (c = 1)" in PROOF 4 and 5.3 is incomplete. The third option, inter-valley pairing with a pseudospin gap (a Luther–Emery-type liquid of momentum-π pairs), is the one the two-body problem points to. It is c = 1, and in 1D it gives a valley-balanced ground state and an odd-N ground doublet with a single unpaired particle, which is what C173's valley_ed.out shows (N = 3: imbalance ±1; N = 4: imbalance 0). So the ring data do not discriminate between easy-plane and paired.
- At n ≈ 0.17 the Fermi scale (μ − ε_min ≈ 0.05 to 0.15) exceeds E_b, so the gas is in a BCS-like regime, where the 1D pseudospin gap is expected to be much smaller than E_b and the pair coherence length much longer than the two-body size of 20 sites [heuristic, by analogy with the attractive Gaudin–Yang gas]. A coherence length beyond L = 96 would read c ≈ 2 in C173's DMRG, so C173's numerics do not exclude this.

[numerical, manybody.py → manybody.out; krho_control.py → krho_control.out] The many-body rings agree with a strongly attractive inter-valley channel. The total-density Luttinger parameter (boson convention, Tonks = 1, two decoupled Tonks gases = 2) from K_ρ = √(vK / (v/K)), vK = πL E''(Φ), v/K = (L/π)[E(N + 2) + E(N − 2) − 2E(N)]/4, is

| L, N (n) | φ = π/2 | φ = 1.3 (one valley) |
| :-- | :-- | :-- |
| 24, 4 (0.167) | 3.16 | 1.46 |
| 30, 4 (0.133) | 3.35 | 1.34 |
| 30, 6 (0.200) | 3.42 | 1.81 |
| 36, 6 (0.167) | 3.6 ± 0.15 (a) | 1.57 |

(a) At L = 36, N = 6 the φ = π/2 ground state is an exactly degenerate ±K doublet, so the lowest level has a cusp at Φ = 0 and manybody.py's flux step returns a negative curvature. stiff36.py takes the flux curvature of the doublet mean (each branch is E_0 ± αΦ + E''Φ²/2, flux conserves K so the branches do not mix) at Φ = ±0.002 and ±0.004, which gives vK ≈ 2.6 with about 7% rounding uncertainty from the printed levels (stiff36.out). With v/K = 0.2035 from manybody.out this is K_ρ ≈ 3.6.

The estimator reads 0.686 and 0.690 on hard rods at L = 24 and 36, against the exact (1 − n)² = 0.694. K_ρ ≈ 3.2 to 3.6 at φ = π/2 is above 2, the value for two decoupled Tonks valleys. For impenetrable valleys (PROOF 4's premise) the total density would behave like one species, with K_ρ near the one-valley value. The one-valley control at the same n reads 1.3 to 1.8, so the intra-valley channel is itself attractive (V_2, V_3 < 0), and K_ρ ≈ 3.3 at φ = π/2 is roughly twice the one-valley value, which is what two weakly coupled valleys give. Impenetrable valleys would instead give about the one-valley value. So the rings show that the valleys are far from the impenetrable regime of PROOF 4. The sign of the inter-valley interaction itself comes from the two-body bound state, not from K_ρ. Parity gaps L·Δ_pb, with Δ_pb = E(N + 1) + E(N − 1) − 2E(N): 2.76 (L = 24, n = 1/6), 2.53 (L = 30, n = 0.133), 2.89 (L = 30, n = 0.2), 1.91 (L = 36, n = 1/6). The φ = 1.3 control reads 5.43, 5.30, 3.53, 3.98. So at n = 1/6 the stratum's L·Δ_pb falls from 2.76 to 1.91 between L = 24 and 36, while the one-valley control stays at 4 to 5.4. No pseudospin gap is resolved at L ≤ 36. Registered expectation R5: this test is indicative only at these sizes.

## 5. Vertices, the momentum-2π two-pair vertex, and openness

### 5.1 Which vertices are resonant on the stratum

[derived] Write c_j ≈ e^{ik_a j}ψ_a(j) + e^{ik_b j}ψ_b(j) and bosonize each valley, ψ_ν ∝ e^{iθ_ν} Σ_m e^{2im(φ_ν + πρ_ν x)} (with half-integer m allowed alongside odd charge for fermionic valleys, which changes nothing below). A vertex e^{i(n_aθ_a + n_bθ_b)} e^{2i(m_aφ_a + m_bφ_b)} has charge ΔN = n_a + n_b and lattice momentum

    P = n_a k_a + n_b k_b + 2π(m_a ρ_a + m_b ρ_b)   (mod 2π).

A uniform perturbation contributes at first order only through P ≡ 0 vertices. On the stratum U forces k_a + k_b = π and, in a U-symmetric state, ρ_a = ρ_b = n/2. With k_a = −s,

    P = (n_b − n_a)s + n_b π + π n (m_a + m_b).

For generic s and n (s = arcsin(|w|/4) at the band level, renormalized but U-preserving), P ≡ 0 iff n_a = n_b, n_b even and m_a + m_b = 0. So on the stratum:
- N_a − N_b is conserved by every non-oscillating vertex (an emergent U(1)_s), and the total charge is conserved mod 4 (U(1)_c broken to Z_4), not mod 2.
- The non-derivative resonant vertices are e^{2√2 i p θ_c} e^{2√2 i m φ_s}, θ_c = (θ_a + θ_b)/√2, φ_s = (φ_a − φ_b)/√2. The leading ones are O_s = cos(2√2 φ_s) (p = 0, m = 1), the valley backscattering ψ†_{aR}ψ_{aL}ψ†_{bL}ψ_{bR}, which is already in (3.2) and is what decides the pseudospin, and O_4 = cos(2√2 θ_c + α) (p = 1, m = 0), the two-pair vertex (ψ_a†ψ_b†)², which is absent from (3.2) (it has ΔN = 4) and present in the qubit chain.

C173 PROOF 3.3 lists the conversion vertices (single transfer, pair transfer, inter- and intra-valley pairing) and states that each carries momentum. That is correct for each of them, but the list is not closed under products. The product of two inter-valley pairings has momentum 2π ≡ 0 on the whole stratum, and the valley backscattering has zero momentum whenever the valleys are balanced. Both are resonant on the stratum itself, which is where the candidate sits. This is a gap in PROOF 5.4 (fixable, below).

Off the stratum, k_a + k_b = π + η and ρ_a ≠ ρ_b, so O_4 carries 2η and O_s carries 2π(ρ_a − ρ_b). Both oscillate, and the only P ≡ 0 vertices left are N-conserving and valley-conserving (pattern P3), up to resonance hypersurfaces of high-order vertices.

### 5.2 Order and dimension of O_4

[derived] Order. In the effective Hamiltonian on P, O_4 has a direct ΔN = 4 term only at O(Δ⁶/V⁵) = |t_2| O(ε⁴) (section 1.3; the largest such element of the exact H_eff is 2.4e-10 at V/Δ = 20 and 1.9e-12 at V/Δ = 40, sw_ring.out). The dominant source is second order in the ΔN = 2 vertex (1.2), A_2 = (4/3)ε|t_2| per bond. Each A_2 insertion carries lattice momentum π (inter-valley) or ∓2s (intra-valley), so it averages out at first order, and two insertions fuse into the P = 0 vertex O_4 with a coefficient set at the lattice scale of (3.2),

    y_4 ≈ c_4 A_2²/|t_2| = c_4 (16/9) ε² |t_2| = c_4 (16/9) Δ⁴/V³,

with c_4 a nonuniversal number of order one. So O_4 enters the infrared at the same order as the first corrections to (3.2), Δ⁴/V³. The separation-2 vertex (1.1) would have given a larger contribution, ∝ (Δ Re w/V)², but it vanishes on the stratum by U. At V/Δ = 40 the bare y_4 is of order 1e-3 |t_2|.

[derived] Dimension. On the stratum U exchanges the valleys, so the charge and pseudospin sectors decouple at quadratic order (c is U-even, s is U-odd). Then

    x(O_4) = (2√2)²/(4K_c) = 2/K_c = 4/K_ρ,

where K_c is the charge parameter in the two-component normalization and K_ρ = 2K_c is the total-density Luttinger parameter in the boson convention (Tonks = 1, two decoupled Tonks gases = 2). The same number follows in the paired phase: there O_4 is the square of the pair field, e^{2iθ_p}, with dimension 1/K_p and K_p = K_ρ/4 for pairs of charge 2. So O_4 is relevant iff K_ρ > 2, in either phase. For comparison x(O_s) = 2K_s.

[numerical, manybody.out] At n ≈ 0.13 to 0.2 the ring estimate is K_ρ ≈ 3.2 to 3.6 (table in section 4), so x_4 ≈ 1.1 to 1.25 and O_4 is relevant. Two decoupled Tonks valleys would make it marginal. Impenetrable valleys (PROOF 4's premise) would give K_ρ near the one-valley value 1.3 to 1.8 and x_4 ≈ 2.2 to 3, irrelevant. The measured stiffness puts it on the relevant side, because at this density the valleys are far from impenetrable and each is attractive.

### 5.3 Can O_4 lock a mode

[derived, given K_ρ > 2] Yes. A relevant cos(2√2θ_c) pins θ_c, as C153's conversion vertex pins the sum of the leg phases at w = 0. The difference is the charge it carries (ΔN = 4 here, ΔN = 2 in C153) and its origin (second order in ε here, a direct third-order term in C153). If the pseudospin is gapless, pinning θ_c leaves the pseudospin mode, c = 1. If the pseudospin is paired (section 4), pinning the pair phase gaps the last mode, c = 0. Either way O_4 removes one gapless mode on the stratum of the qubit chain, a mode that (3.2), which conserves N, keeps. The gap scale is

    Δ_4 ~ |t_2| (y_4/|t_2|)^{1/(2 − x_4)} ≈ |t_2| (1.8 c_4 ε²)^{1.25}   (x_4 ≈ 1.2),

i.e. ξ_4 ~ v/Δ_4 of a few thousand sites at V/Δ = 40, but of order 100 sites at V/Δ = 8 and a few hundred at V/Δ = 12 to 16 [heuristic, c_4 ~ 1 and the ring K_ρ]. So effective-model DMRG cannot see O_4 at all, the V/Δ = 40 qubit chain cannot see it at L ≤ 256, and the small-V/Δ qubit chains of C177 may.

### 5.4 Openness

[derived, conditional on a c = 2 phase of (3.2) somewhere near the stratum] O_4 and O_s are resonant only on the stratum (and the codimension-1 hypersurface k_a + k_b = π that continues it off the parity class). A small departure η, of either sign, gives them momentum 2η and 2π(ρ_a − ρ_b), and a relevant vertex with small momentum q is cut off once vq exceeds its gap (commensurate-incommensurate transition). So O_4 can at most remove a slab of width of order Δ_4/v around the stratum, and the openness argument of PROOF 5.4 survives O_4. It does not survive in the form written, because it needs the c = 2 phase to exist off the stratum, at a valley offset δ = 4|w| cos φ ≠ 0. C173's DMRG reading at φ = 1.54 (δ = 0.12) is the evidence for that, and it is a single L = 48 run.

[heuristic] What the binding implies off the stratum. With an inter-valley bound state, a valley offset δ acts like a Zeeman field on an attractive two-component gas. In the 1D attractive Gaudin–Yang gas a field above the spin gap produces a partially polarized two-component liquid (FFLO-like) with c = 2, on an open set of (μ, h). If the same happens here, the c = 2 region would sit off the stratum, at |δ| between the pseudospin gap and the saturation offset where the second valley empties, with O_4 and O_s oscillating there. That is a different mechanism from the easy-plane pseudospin of PROOF 4, and it would make the stratum the worst place to look (paired, and in the qubit chain also O_4-locked). The φ = 1.54 reading is consistent with it. I have not tested it, and the competing outcome (first-order jump from the paired state to one valley, as C173 PROOF 4 argues for the impenetrable case) is not excluded.
