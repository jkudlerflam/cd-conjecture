# C186 derivations: the general second-order reduction at the Fibonacci cusp, and what it allows

This run was done by a Claude Opus 5.5 agent (C186, 2026-11-09), not by Codex (Astra out of credits).

Standard: theoretical physics. Tags [proved], [derived], [numerical], [heuristic], [assumed], [literature] open the paragraph they apply to. Every numerical statement names its output file (all in this folder).

## 0. Conventions and the general bond

Spin language throughout: particle = spin down, n = (1 − Z)/2, b† = σ⁻ = (X − iY)/2, b = σ⁺. No Jordan–Wigner signs are used, which matters below because the bond will not conserve parity.

[proved, algebra] Up to telescoping single-site terms and a constant, every nearest-neighbour qubit bond can be written in boson form as

    h = V n₁n₂ − (μ/2)(n₁ + n₂)
        + [ w b₂†b₁ + p b₁†b₂† + b₁†(a + β n₂) + b₂†(a′ + β′ n₁) + h.c. ].                      (0.1)

Counting: V, μ (2 real), w, p (4), and of a, a′ only λ = a + a′ survives the telescoping (2), β, β′ (4), which is 12 = 9 + 3, the dimension of (J, B) in C154's form. Summed over bonds, the creation operator at site j carries the amplitude

    c_j = λ + β′ n_{j−1} + β n_{j+1},                                                          (0.2)

so f_L ≡ λ + β′ is the amplitude when only the left neighbour is occupied, f_R ≡ λ + β when only the right one is, and f_L + f_R − λ when both are. A nearest-neighbour bond cannot condition on both neighbours independently, which is the structural fact behind section 3.

Dictionary to the Pauli convention h = J_x XX + J_y YY + D(XY − YX) + K(XY + YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) + (h_x/2)(X⊗1 + 1⊗X) + (h_y/2)(Y⊗1 + 1⊗Y) + Γ_xz XZ + Γ_zx ZX + Γ_yz YZ + Γ_zy ZY [proved, algebra, checked numerically through sw_general.py's direct boson build]:

    V = 4J_z,  μ = 4J_z + 2h_z,  w = J_x + J_y + 2iD,  p = J_x − J_y + 2iK,
    c_j = (h_x + ih_y) + (Γ_xz + iΓ_yz)(1 − 2n_{j+1}) + (Γ_zx + iΓ_zy)(1 − 2n_{j−1}),
    so λ = (h_x + ih_y) + (Γ_xz + Γ_zx) + i(Γ_yz + Γ_zy),
       d ≡ f_L − f_R = 2[(Γ_xz − Γ_zx) + i(Γ_yz − Γ_zy)],
       s ≡ f_L + f_R = 2λ − 2(Γ_xz + Γ_zx) − 2i(Γ_yz + Γ_zy).                                 (0.3)

So the inversion-odd facilitated amplitude d is a transverse Dzyaloshinskii–Moriya vector (along y for real d, along x for imaginary d), and s is the symmetric off-diagonal exchange (Γ-type) together with the transverse field.

The parity class (bonds commuting with Π_j Z_j) is exactly λ = β = β′ = 0.

## 1. The cusp regime

[assumed, definition] The cusp regime is V ≫ |p|, |β|, |β′| (order ε), with every first-order process inside the blockaded space P (no two adjacent particles) tuned to order ε²/V: the hop w, the unconditional creation λ, and μ − (blockade value). This is the regime of C173 extended to the full bond. With H₀ = V Σ n_j n_{j+1} and H₁ the rest, the effective Hamiltonian to order ε²/V is

    H_eff = P H₁ P − P H₁ Q H₀⁻¹ Q H₁ P,                                                       (1.1)

with P H₀ = 0. Cross terms involving w or λ are of order ε³/V² or smaller and are not part of the second-order model (they are seen in the all-orders check of section 4 as deviations falling like 1/V).

## 2. The second-order model

[derived; enumeration of virtual processes; checked to machine precision in section 4] The only operators that leave P are the pairing p (always, cost V(1 + n_{j−1} + n_{j+2}) for a pair on (j, j+1)) and the facilitated creation (creation next to an occupied site, cost V times the number of occupied neighbours). Write pp = |p|². With λ = 0 at this order:

(a) Pairing alone (C173): NNN hop −pp/(2V) (path: create the pair (j+1, j+2) next to a particle at j, cost 2V, annihilate (j, j+1)); diagonal −pp/(V(1 + n_{k−1} + n_{k+2})) per empty bond (k, k+1). This gives the on-site shift +3pp/V, V₂ = −pp/V, V₃ = −pp/(3V), V₄ = 0 (C173 PROOF 3.1, C179 P1.1).

(b) Facilitated terms twice. An empty site j with occupied neighbours contributes the diagonal −|β′n_{j−1} + βn_{j+1}|²/(V(n_{j−1} + n_{j+1})). Per isolated particle this is −(|f_L|² + |f_R|²)/V. For two particles at distance 2 the shared empty site has both neighbours occupied, amplitude f_L + f_R and cost 2V, so

    δV₂ = (|f_L|² + |f_R|²)/V − |f_L + f_R|²/(2V) = |f_L − f_R|²/(2V) = |d|²/(2V) ≥ 0.          (2.1)

No other separation is affected, so δV₃ = δV₄ = 0. The same pair of operators also moves a particle by one site (create at j next to j + 1, then remove j + 1), a NN hop −f_R f_L*/V for j+1 → j, which is P-projected and uniform, so it only renormalizes w.

(c) Pairing and facilitated term together, ΔN = ±1. Create the pair (j, j+1) (cost V(1 + n_{j+2})), remove j + 1 (amplitude conj(f_L + f_R n_{j+2}) with λ = 0); or create (j−1, j) and remove j − 1. The sum is

    P b_j† [ −p s*/V + (p d*/(2V)) (n_{j+2} − n_{j−2}) ] P.                                   (2.2)

(d) Nothing else. Two facilitated creations cannot return to P with ΔN = 2 (the first leaves an adjacent pair the second cannot remove), pairing twice gives no ΔN = ±2 inside P at this order (C173), and a facilitated annihilation from P is first order (it stays in P) and equals conj(λ) there.

Collecting, in units |t₂| = pp/(2V) and with g ≡ −V₃/|t₂|, v ≡ V₂/|t₂| + 2,

    NNN hop −1,  V₃ = −2/3,  V₂ = −2 + v,  v = |d|²/|p|²,
    NN hop w_eff = w − f_L f_R*/V,  on-site −μ + 3pp/V − (|f_L|² + |f_R|²)/V,
    ΔN = +1:  P b_j† [ λ_eff + r (n_{j+2} − n_{j−2}) ] P + h.c.,
              λ_eff = λ − p s*/V,  r = p d*/(2V),  |r| = √v |t₂|.                              (2.3)

## 3. The obstruction at second order

[derived from (2.3)] Three rigidities follow.

1. g = 2/3 for every nearest-neighbour bond in the cusp regime. V₃ comes only from the pairing (a facilitated creation never sees a particle at distance 3, since the intermediate state with particles at j, j+1, j+3 costs V, the same as with j + 3 empty).
2. v = |d|²/|p|² ≥ 0, and v > 0 requires d ≠ 0, a transverse DM component (0.3). In the parity class v = 0 identically. So inside the parity class the collaborator's unbound side v > g u²/[2(2 − g)] is out of reach at second order, for every hop modulus u and phase.
3. Whenever v > 0 there is a conditioned ΔN = 1 vertex r b_j†(n_{j+2} − n_{j−2}) of modulus exactly √v |t₂|. It cannot be cancelled. The only first-order ΔN = 1 operator a nearest-neighbour bond puts inside P is λ P b_j† P (a creation inside P needs both neighbours empty, so its conditioning on n_{j±1} is invisible), and at second order the only parity-odd source is (2.2). So λ can cancel the uniform part p s*/V and nothing can cancel r. The phase of r is a gauge choice (a global z-rotation, once λ_eff = 0), so the second-order cusp family with v > 0 is the one-parameter extension

       H_v = H_C173(u, φ, μ; V₂ = −2 + v, V₃ = −2/3) + √v Σ_j P[ b_j†(n_{j+2} − n_{j−2}) + h.c. ]P.   (3.1)

At the collaborator's boundary, v_c = u²/4 (g = 2/3), the vertex has r = u/2: 0.5|t₂| at u = 1 and |t₂| at u = 2. Number non-conservation at that order is not a small perturbation of the two-valley gas. Whether it rebinds the pair, clusters the gas or leaves an unbound window is a dynamical question, section 6.

Symmetry [proved]. The vertex breaks C173's Σ (it maps r b_j†(n_{j+2} − n_{j−2}) to a staggered term), so the valley degeneracy is no longer protected. An antiunitary Θ′ = inversion × complex conjugation × R_z(π) survives (with r real), and it preserves each valley. A substitute for Σ that exchanges the valleys does not exist within inversions × z-rotations, so degenerate dressed valleys need φ retuned, a codimension-one condition. That is allowed by the spec ("a fine-tuned point counts").

## 4. Checks against the full chain (sw_general.py)

[numerical, sw_partA.out] Part A builds the full spin chain from the 4 × 4 bond (0.1) and computes (1.1) exactly as matrices, at w = λ = 0, random complex p, β, β′, μ, V ∈ [3, 30], on open chains and rings, L = 8 to 12 (18 cases). The closed form of section 2, built independently from the local rules (edges included), agrees entry by entry to at most 3.6e-15. This checks every second-order coefficient including the absence of ΔN = ±2 and of V₃, V₄ corrections.

[numerical, sw_partB.out (open L = 11) and sw_partB_ring.out (ring L = 12)] Part B computes the all-orders des Cloizeaux effective Hamiltonian on P for a generic bond (p = e^{0.3i}, β = 0.35 + 0.2i, β′ = −0.25 + 0.4i, w = (0.4 + 0.9i)|t₂|, λ = (0.3 − 0.5i)|t₂|) at V = 20, 40, 80, 160 and reads the bulk coefficients. Every deviation from (2.3), in units of |t₂|, falls like 1/V (for example the conditioned ΔN = 1 element: 1.3e-2, 6.8e-3, 3.5e-3, 1.8e-3; the NN hop: 2.6e-2, 1.3e-2, 6.6e-3, 3.3e-3). The ΔN = 2 elements are nonzero at relative order ε/V (6.6e-2 → 8.3e-3 for separation 3), third order as expected, and V₄ ∝ V⁻² relative. Ring and open chain agree. R1 passed.

## 5. The parity class beyond second order

[numerical, sw_parity.out; my run, open L = 11, u = 0.5, φ = π/2; C179 sw_exact.out agrees at u = 1] All-orders des Cloizeaux coefficients in units of the measured NNN hop:

| V/Δ | v = V₂/\|t₂\| + 2 | g = −V₃/\|t₂\| | V₄/\|t₂\| | y₃ (ΔN = 2, sep. 3) | v_c = g u²/[2(2 − g)] at u = 0.5 |
| --: | --: | --: | --: | --: | --: |
| 8 | 0.0536 | 0.6102 | 0.0513 | 0.1657 | 0.0549 |
| 10 | 0.0347 | 0.6301 | 0.0333 | 0.1328 | 0.0575 |
| 14 | 0.0179 | 0.6478 | 0.0172 | 0.0950 | 0.0598 |
| 20 | 0.0088 | 0.6574 | 0.0085 | 0.0666 | 0.0606 |
| 40 | 0.0022 | 0.6643 | 0.0021 | 0.0333 | 0.0612 |

So v ≈ 3.5(Δ/V)², 2/3 − g ≈ 3.7(Δ/V)², V₄ ≈ 3.3(Δ/V)² (repulsive), and y₃ = (4/3)(Δ/V) to 0.5 percent. The fourth order moves in the unbinding direction, and at u = 0.5 the number-conserving part reaches the collaborator's boundary near V/Δ ≈ 8 (more easily with the repulsive V₄, which the boundary formula does not include). In general unbinding needs u ≲ 3.7 Δ/V.

[derived] At the same order the parity class carries the third-order ΔN = 2 vertex at separation 3, of amplitude y = (4/3)(Δ/V)|t₂|, at total momentum zero. Creating a pair at K = 0 out of the vacuum costs min_q[ε(q) + ε(−q)] − 2ε_min = u²/2 (at φ = π/2, ε(q) + ε(−q) = −4 cos 2q). Its second-order effect is of relative size y²/(u²/2) = (32/9)(Δ/V)²/u². At the unbinding boundary u² ≈ 14 (Δ/V)² this is about 0.25 |t₂|, independent of Δ/V, while the boundary itself sits at v_c ≈ 3.5(Δ/V)² → 0. At the tabulated V/Δ = 8, u = 0.5 point: y₃²/(u²/2) = 0.22 against v_c = 0.055. So in the parity class the pair can only unbind where the momentum-zero pair-creation vertex is near resonance and dominates every scale of the two-body problem. As u → 0 this vertex becomes C153's resonant conversion vertex, which locks the total density (c = 1). The number-conserving effective model is not a controlled description there [derived, scaling; no numerics at the required V/Δ ≲ 10 and u ≲ 0.3 were run].

### 5.1 The chemical-potential dependence of the fourth-order coefficients

[derived; numerical check sw_mu.out] Every second-order process goes through a virtual pair whose energy is kV − 2μ_bare (k adjacencies, two extra particles), with μ_bare = (μ_eff + 6)|t₂| + O(ε⁴/V³). So the bond energies of 2(a) are e_k = −(p²/kV)[1 + a/k + …] with a = (μ_eff + 6)(Δ/V)², and the NNN hop carries the factor 1 + a/2. Then E{m} − E₀ = −4e₁ + 2e₂ and E{m, m+2} − E₀ = −6e₁ + 2e₂ give V₂ = 2e₁ − 2e₂ = −(p²/V)(1 + 3a/2), and V₂/|t₂| + 2 = −2a + O(a²). The other fourth-order terms do not depend on μ_eff at this order. So v(μ_eff) = v(−6) − 2(μ_eff + 6)(Δ/V)². Measured at u = 1: v(Δ/V)² = 17.8, 21.1, 22.3 at μ_eff = −12 and 2.9, 3.0, 3.1 at μ_eff = −2.25 for V/Δ = 8, 16, 40, i.e. slopes −1.53, −1.85, −1.98 (approaching −2). The μ_eff = −12 values reproduce C177's two-body fit (twob_fit.out model B: V₂/|t₂| = −1.723, −1.917, −1.986, against −1.722, −1.917, −1.986 here). So that fit and the des Cloizeaux coefficients agree once both are taken at the same μ. The physically relevant μ for the dilute cusp gas is the onset, μ_eff ≈ −2.2, where v ≈ 3(Δ/V)².

## 6. Dynamics of the mixed model (3.1)

### 6.1 Calibration

[numerical, calib_r0.out] mixed_fewbody.py (infinite line, fixed total momentum, relative coordinates, hard wall at extent R) at r = 0 reproduces the collaborator's theorem: E_b(v = 0) = 0.019527 (u = 1) and 0.12255 (u = 2), equal to their exact values; the bound state disappears between 0.9 v_c and 1.1 v_c at u = 0.5, 1, 2, with the unbound side showing the −1/R² approach to threshold.

### 6.2 Dressing

[numerical, dressing.out] The single particle at a bare valley minimum, μ at the bare band bottom, with particle number truncated at N_max = 1 to 5:
- u = 1, v = 0.25 (r = 0.5): valley a is lowered by 0.044 (93 percent N = 1), valley b by 0.11 (about 50 percent N = 1). The vertex splits the valleys by about 0.066 |t₂| at φ = π/2, three times the r = 0 binding energy.
- u = 1, v = 0.5: valley b is lowered by 0.21 with 4 to 9 percent N = 1 weight, the rest in N = 2 to 4.
- u = 2, v = 1 and 2: at the bare band bottom the lowest state at the valley momentum has 3 to 5 particles and keeps falling with N_max. This only says that the bare band bottom lies far above the dressed onset (the ring onsets of 6.3 sit 0.47 and 0.74 below it), so particles proliferate there. It is not by itself evidence of bound clusters, and I do not use it as such.
- u = 0.5, v = 0.0625 (r = 0.25): converged, lowering 0.0065 (a) and 0.015 (b).

[derived, interpretation] Why valley b is dressed more. The vertex turns a particle at momentum q into two particles of total momentum q with relative amplitude ∝ 2r sin q at separation 2. A particle at k_b ≈ π + s lands next to the unlike-pair channel at total momentum π, where the bound ab pair (or, unbound, the near-resonant pair continuum) lives. A particle at k_a ≈ −s lands next to the like-pair channel at 2k_a, which is the weak odd channel. So the vertex hybridizes valley b with the ab pair band, and the valley degeneracy of Σ is lost at order r².

The N_max truncation is not size-consistent for the pair question, so the binding test moved to exact ring ED (6.3).

### 6.3 Exact ring ED and the density onset

[numerical, check_ring.out, check_trunc.out] ring_mixed.py (numba, Bloch states over translation representatives, all particle numbers or a cap N_max) reproduces a direct real-space construction of (3.1) with λ ≠ 0 to 5e-14 at L = 12, 13, and the capped basis equals the full one restricted to N ≤ N_max exactly.

[numerical, onset.out, L = 20 only] Onset μ_on where the lowest non-vacuum level crosses zero. Controls: r = 0, u = 1, v = 0, 0.1, 0.2 give a K = π object with ⟨N⟩ = 2 first, v = 0.3 a single particle at K = 3.456 (the grid point next to k_b). With r = √v at u = 1 the first object is a dressed single particle at K = 3.456 (⟨N⟩ = 1.16 to 1.20) for v = 0.1 to 0.6, with μ_on falling from −2.292 to −2.396 (bare band bottom −2.25), and at v = 1 a dressed particle with ⟨N⟩ = 1.45 followed within 0.02 by a three-particle object (⟨N⟩ = 2.9). At u = 2 the onset lies 0.31 to 0.74 below the bare band bottom, with ⟨N⟩ = 1.27 to 1.37 at onset and a ⟨N⟩ = 2.7 object right above it at v = 1. At u = 0.5 the first object has ⟨N⟩ = 1.65 to 2.0 and a ⟨N⟩ = 3.1 to 3.7 state lies within 0.02 above. L = 20 does not resolve binding energies below about 0.02, so these runs only show the scale of the dressing (the onset moves by 0.04 to 0.15 at u = 1, comparable to or larger than v_c-scale binding energies) and that the first object at u = 1 sits in valley b, the more strongly dressed one, at φ = π/2.

### 6.4 Dressed band at u = 2, v = 1 (the naive boundary, r = 1)

[numerical, band_u2v1_N3.out, ring L = 48, N_max = 3, φ = 1.6043 and μ = −3.4033 tuned on the ring] The dressed one-particle branch (⟨N⟩ = 1.0 to 1.6) is flat to within 0.06 from q = 3.80 to 5.63, with shallow minima at q = 4.06 (−0.053) and 5.50 (−0.016), against a bare band that rises to +1.4 in between. The two-valley structure of the bare band is gone, replaced by a heavy hybrid branch between the bare valleys, and the K = π two-particle states sit 0.19 above. At this coupling the second-order cusp model is no longer a two-valley gas at all. This N_max = 3 band is not converged in N_max (the hybrid has ⟨N⟩ up to 1.6), so only the qualitative statement is used.

### 6.5 Pair binding at u = 1.0353 (valleys on the L = 48 grid)

Method [numerical, pairbind2.py]: ring L = 48, particle number capped at N_max, φ tuned on the ring until the two dressed valley minima are degenerate, μ tuned so they sit at zero (the dilute onset of quasiparticles), then B = T(K*) − E_low(K*) at K* = q_a + q_b = π. Calibration at r = 0 (pb2_r0.out): B_ring = 0.0278, 0.0148, 0.0051, −0.0021, −0.0075 at v = 0, 0.1, 0.2, 0.3, 0.4, against the exact infinite-line E_b = 0.0218, 0.0088, 0.0015, 0, 0 (boundary 0.268). So B_ring tracks E_b with an upward offset of about 0.005 for bound states, and crosses zero within about 0.03 of the exact boundary.

[numerical, pb2_raw.out (tuning at N_max = 3, 4), pairconv_v*.out (N_max = 3 to 6 at fixed tuning), true_threshold.out] Two corrections to the naive reading turned out to matter.

(i) Particle number. A pair needs about twice the extra particles of one quasiparticle, so at equal N_max the pair is under-dressed. At fixed (φ, μ) the quasiparticle levels converge by N_max = 5, and the K = π level converges geometrically by N_max = 6 (at v = 0.2: 0.0285, 0.0041, −0.0005, −0.0010; at v = 0.3: 0.0553, 0.0155, 0.0060, 0.0046; at v = 0.1 and 0.05 the N_max = 5 and 6 values agree to 6e-5 and 1e-6).

(ii) The grid. The vertex moves the dressed valley minima off the ring momenta (k_a = −0.273, −0.285, −0.308, −0.331, −0.354 at v = 0.05 to 0.4, against the grid point −π/12 = −0.262), so tuning the grid minimum to zero leaves the true infinite-volume threshold T_true below zero. T_true is computed from the infinite-line dressed minima (true_threshold.py, N_max = 5, box 26), whose grid values reproduce the ring E1 to about 1e-4.

| v | r | E_pair(N_max = 6) | T_true | E_pair − T_true | r = 0 on the same ring: E_pair − T | exact r = 0 E_b |
| --: | --: | --: | --: | --: | --: | --: |
| 0.05 | 0.224 | −0.0053 | −0.0007 | −0.0046 | −0.0208 | about 0.015 (interpolated) |
| 0.1 | 0.316 | +0.0027 | −0.0039 | +0.0065 | −0.0148 | 0.0088 |
| 0.2 | 0.447 | −0.0010 | −0.0156 | +0.0145 | −0.0051 | 0.0015 |
| 0.3 | 0.548 | +0.0046 | −0.0108 | +0.0154 | +0.0021 | unbound |
| 0.4 | 0.632 | +0.0117 | −0.0033 | +0.0150 | +0.0075 | unbound |

Reading [numerical, with the r = 0 calibration as the yardstick]. On the same ring a bound pair sits below the true threshold by about E_b + 0.005 (r = 0 rows), and an unbound one above it. With the vertex, the K = π level lies below T_true only at v = 0.05 (by 0.0046, against 0.0208 at r = 0), and above it by 0.0065 to 0.015 for v ≥ 0.1, where the number-conserving model is still bound by 0.0088 and 0.0015. Binding energies below about 0.002 (pair sizes beyond the ring) cannot be resolved at L = 48. So within the ring resolution the vertex moves the binding boundary at u = 1.035 from v_c = 0.268 down into 0.05 < v* ≲ 0.1. The third addendum's registered reading (unbound if B < −0.002 with the grid threshold) gives "unbound" at v = 0.3 and 0.4 (B = −0.0049, −0.0122) and "marginal" at v = 0.2 (+0.0011); the T_true reading above is the corrected one, and it agrees at v = 0.3, 0.4.

[heuristic, mechanism] The vertex lowers every quasiparticle by its dressing (0.0007 to 0.016 here, valley b more). Two quasiparticles close together lose part of it: a virtual particle between them at distance 2 from both has amplitude r(n_{j+2} − n_{j−2}) = 0, and the blockade removes others. That costs energy at short distance, a repulsion in every channel, larger than the binding scale already at r ≈ 0.3.

[numerical, band_u1035v01_N3.out, band_u1035v03_N5.out] At the tuned (φ, μ), the lowest level over all 48 momenta is the dressed quasiparticle at the two valleys (⟨N⟩ = 1.03 to 1.36), with nothing below it: no bound pair or cluster precedes the single-quasiparticle onset at v = 0.1 (N_max = 3) and v = 0.3 (N_max = 5). The two dressed valleys are degenerate after retuning φ (1.5762, 1.5793, 1.5831, 1.5881, 1.5895 at v = 0.05 to 0.4), a codimension-one condition.

[numerical, channels_v005_N4.out, channels_v01_N3.out, channels_r0.out] Two-quasiparticle channels at v = 0.1 (N_max = 3): ab +0.0116, aa +0.0238, bb +0.0467 above the grid threshold, so the lowest two-quasiparticle state is the balanced ab one. Three-quasiparticle channels are not converged at these N_max and are not used.

### 6.6 The number-conserving baseline past the boundary

[numerical, pseudospin_r0.out, C175's scat.py imported read-only; trimer_r0.out] For the model with r set to zero and v > v_c, C175's dilute criterion gives the two-component window at every point tested (u = 1, 1.035, 2, v/v_c = 1.2 to 3): the unlike length a₋ has passed through −∞ (−88 to −2.7) and a_s lies inside (a₋, a₊) with Δ_σ = 0.60 to 0.95 (easy-plane, approaching the isotropic point at the boundary). No three-body bound state at u = 1, v = 0.3 and 0.4 (the lowest N = 3 level falls as 1/R², +0.0187 → +0.0070 and +0.0245 → +0.0095 from R = 80 to 140). So the number-conserving extension past the boundary is a dilute two-component liquid by C175's and C180's criteria. The vertex is present in every real bond, though, and the pseudospin lengths of the mixed model were not computed.

## 7. What the vertex does in the infrared

[derived, C180's operator table extended] With Σ broken and φ retuned, the dressed valleys sit at k_a + k_b ≠ π (−3.119 at v = 0.3, a mismatch of 0.022). So C179's two-pair vertex (momentum 2(k_a + k_b)) carries momentum 0.045 and is off resonance over lengths beyond about 20 sites, far shorter than C179's locking length. The vertex itself and every ΔN-odd product carry n_a k_a + n_b k_b with n_a + n_b odd, which is nonzero for incommensurate valleys, so they are irrelevant away from resonance tongues at special densities. The backscattering cos 2√2φ₋ stays the deciding operator, as at Σ. What the vertex changes is the microscopic input (dressing, valley splitting, the binding boundary), not the list of resonant operators.
