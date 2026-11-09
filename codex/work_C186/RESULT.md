This run was done by a Claude Opus 5.5 agent (C186, 2026-11-09), not by Codex (Astra out of credits).

# C186 result: can a nearest-neighbour bond unbind the Fibonacci-cusp pair?

Derivations are in PROOF.md, the log in PROGRESS.md. Inputs read first: the collaborator's (GitHub: salpace) contrib/cusp_binding_certificate_20261109/, contrib/cusp_scattering_20261109/ and contrib/cusp_validation_20261109/ (read-only). The validation README cites a file FACILITATED_REPAIR.md (a transverse-DM repair and its density-assisted creation term), which is not present in that folder, so it could not be read or compared.

## Pre-registered predictions and pass criteria (written before any code in this folder was run)

Conventions (PROOF 0). Spin language, particle = spin down, b† = σ⁻. Near the Fibonacci point a general nearest-neighbour bond is, up to telescoping terms,

    h = V n₁n₂ − (μ/2)(n₁ + n₂) + [w b₂†b₁ + p b₁†b₂† + b₁†(a + β n₂) + b₂†(a′ + β′ n₁) + h.c.],

with V large, p (pairing) of order ε, the facilitated amplitudes β, β′ of order ε, and the first-order dynamics inside the blockaded space P (the hop w, the unconditional creation λ = a + a′ and μ − blockade value) of order ε²/V. Write f_L = λ + β′ and f_R = λ + β (creation amplitude at j when only the left, or only the right, neighbour is occupied), d = f_L − f_R, s = f_L + f_R.

R1 (second-order Schrieffer–Wolff, my derivation). On P, to order ε²/V,
- NNN hop ⟨j+2|H|j⟩ = −|p|²/(2V), unchanged by any other term;
- V₃ = −|p|²/(3V), unchanged; V₄ = 0;
- V₂ = −|p|²/V + |d|²/(2V);
- NN hop (coefficient of b_{j+1}†b_j) w − f_L f_R*/V;
- on-site −μ + 3|p|²/V − (|f_L|² + |f_R|²)/V;
- ΔN = +1: P b_j† [ λ − p s*/V + (p d*/(2V)) (n_{j+2} − n_{j−2}) ] P, and its conjugate;
- no ΔN = ±2 and no other term at this order.

Pass: the exact second-order matrix −P H₁ Q H₀⁻¹ Q H₁ P + P H₁ P of the full open chain (L = 8 to 12, random complex couplings, also rings) agrees with the closed form entry by entry to 1e-10. The des Cloizeaux (all-orders) effective Hamiltonian of the full chain reproduces each coefficient with relative deviation falling at least as fast as 1/V.

R2 (obstruction at second order, derived from R1). In units |t₂| = |p|²/(2V) every cusp bond has g = 2/3 exactly, v ≡ V₂ + 2 = |d|²/|p|² ≥ 0, and an uncancellable conditioned ΔN = 1 vertex of modulus |r| = √v (the unconditional part is cancelled by tuning λ). In the parity class d = 0, so v = 0. So no bond reaches the collaborator's unbound side v > g u²/[2(2 − g)] without a number-changing vertex of strength √v |t₂|.

R3 (two-body test of the mixed model, numerical). The second-order cusp model with v > 0 and r = √v (number not conserved, vacuum exact) is solved on the infinite line at fixed total momentum K = π in relative coordinates, with particle number truncated at N_max, and the dressed single-particle minimum set to zero by μ.
- Calibration: at r = 0 the binding energy vanishes at v_c = u²/4 (u = 1: 0.25, u = 2: 1.0) within the box resolution. If not, the code is wrong and nothing from it is used.
- Criterion: the pair is bound if the lowest K = π level lies below the dressed two-particle threshold by more than 1e-4, stable between the two largest boxes, and changing by less than 30 percent from N_max = 3 to 4. It is unbound if the gap closes within box resolution at the two largest boxes.
- Expectation [heuristic, no sign registered]. With N_max = 2 the vertex only lowers the threshold (the K = π pair cannot mix with one particle, the amplitude is ∝ sin K), which helps unbinding. At N_max ≥ 3 the pair has its own three-body virtual processes, with a repulsive part (the created third particle interferes destructively between the two, and the self-energy is blocked) and an attractive part (V₂, V₃ with the partner in the intermediate state). I do not know which wins.

R4 (parity class at moderate V/Δ). Prediction: the fourth-order corrections push v up and g down (C179's ratios suggest so), so the N-conserving part unbinds only at u ≲ κ Δ/V with κ of order one. There the momentum-zero ΔN = 2 separation-3 vertex (amplitude (4/3)(Δ/V)|t₂|) is detuned by only about u²/2 (in |t₂|), so the number-conserving description is not controlled exactly where unbinding happens.

R5 (only if R3 finds a robust unbound window). Pseudospin criterion with C175's lengths (a₋ < a_s < a₊), then the full qubit chain. Registered pass for a c > 1 candidate is the spec's criterion (open-chain fits above 1.2 at truncation below 1e-7 that do not fall with L, the parity gap against L, controls). Not reached means no candidate is claimed.


### Addendum to R3, registered after dressing.out and before any ring run

The truncated few-body comparison of R3 is not size-consistent (a single particle uses the particle-number budget for its own dressing, a pair cannot dress both constituents), so B from mixed_fewbody.py at N_max = 2, 3 is not used. dressing.out (dressed one-particle states, N_max up to 5) showed strong dressing, unequal between the valleys, and at u = 2, v = 1 and 2 states dominated by 3 to 5 particles far below the bare band. Replacement test, exact in particle number: momentum-resolved ring ED of the mixed model (ring_mixed.py, validated against a real-space build in check_ring.out), and the density onset μ_on(L) where the lowest non-vacuum level crosses zero (onset.py). The object that appears first is identified by its momentum and ⟨N⟩.
- Controls at r = 0. u = 1, v = 0: first object at K ≈ π with ⟨N⟩ = 2 (the known pair). u = 1, v = 0.35 > v_c: first object is one particle (⟨N⟩ = 1) at a valley momentum.
- Prediction for r = √v: at u = 2 and v ≥ 1, the first object has ⟨N⟩ ≥ 3 and μ_on lies well below the bare band bottom (clustering). At u = 1 and v ≥ v_c = 0.25 I expect the first object to be a multi-particle state (⟨N⟩ ≥ 2) as well, with no unbound window. Pass criterion for an unbound window: some v at fixed u where the first object has ⟨N⟩ < 1.6 at K near a valley at the two largest L, and the dressed valleys can be made degenerate by retuning φ.

### Second addendum, registered after the L = 20 onset runs and the r = 0 ring calibration, before any locked ring-binding output

The L = 20 onsets (onset.out) cannot resolve binding energies below about 0.02 (at r = 0 they show the pair first at v = 0.1 and 0.2 only because the one-particle grid misses the valley). The binding test is therefore pairbind2.py: rings with a particle-number cap N_max (ring_mixed.Ring(L, N_max), validated in check_trunc.out), φ tuned on the ring so the dressed valleys are degenerate, μ tuned so they sit at zero, and B = T(K*) − E_low(K*) at K* = q_a + q_b. Calibration at r = 0 (pb2_r0.out, u = 1.0353 where L = 48 contains the valley momenta): B = 0.0278, 0.0148, 0.0051, −0.0021, −0.0075 at v = 0, 0.1, 0.2, 0.3, 0.4, against the exact 0.0218, 0.0088, 0.0015, unbound, unbound (boundary 0.268). So L = 48 overbinds by about 0.005 and places the boundary within about ±0.05 in v.
- Registered reading for r = √v at u = 1.0353, L = 48, N_max = 3 and 4: the vertex prevents unbinding if B > +0.005 (beyond the calibration offset) at every v ≤ 0.5 with N_max = 3 and 4 agreeing in sign. An unbound window is indicated if B < −0.002 at some v at both N_max, with the dressed one-particle band keeping two separated minima.
- Prediction: B stays positive for all v up to 0.5 (rebinding by the vertex). I hold this weakly; the dressing lowers the threshold, which works the other way.

### Third addendum, registered after v = 0.05 to 0.2 at N_max ≤ 6 and before the v = 0.3, 0.4 convergence runs

N_max = 3 and 4 agree in sign at v = 0.1 and 0.2 (B < 0), which formally meets the second addendum's "unbound window indicated". But B moves by +0.008 and +0.022 from N_max = 3 to 4, so that criterion is not trustworthy. pair_conv.py holds (φ, μ) from the N_max = 4 tuning and raises N_max to 6 at K = π only. At v = 0.2 the pair level converges geometrically (0.0285, 0.0041, −0.0005, −0.0010 for N_max = 3 to 6) while the quasiparticle levels converged by N_max = 5. The decisive number is therefore the N_max-extrapolated B_ring. Registered reading: unbound if it is below −0.002 (the r = 0, v = 0.3 value), bound if above +0.002, marginal in between. Prediction for v = 0.3 and 0.4: unbound.

### Typhon job, registered before submission

(A) pairbind2.py on L = 48 rings at N_max = 4 and 5 for u = 0.522, 1.035, 1.531, 2 (bare valley momenta on the grid), v/v_c = 0.25, 0.5, 1, 1.5, 2 with r = √v, and r = 0 calibrations. Reading as in the second addendum, with N_max = 4 against 5 as the convergence test (sign of B must agree).
(B) bm_driver.py: open-chain block DMRG of the mixed model at u = 1.035, L = 96 and 192 (3-site blocks, χ ≤ 128), μ = μ_onset + 0.03 and + 0.06, for r = √v at v = 0.1, 0.2, 0.3, and the r = 0 controls v = 0 (bound pair) and v = 0.4 (unbound, N conserved). Two-size c = 6[S̄(192) − S̄(96)]/ln 2 with S̄ the mean of the three middle block cuts. Registered reading: a locked point reads "two-component at these sizes" if c ≥ 1.6, truncation below 1e-6, and within 0.3 of the v = 0.4 control at the same μ offset; "one mode" if c ≤ 1.3. This is a first-stage screen only; it is not the spec's candidate criterion, which also needs FES and the full chain.

## Results

## Verdict

A nearest-neighbour bond can move the Fibonacci cusp across the collaborator's binding boundary, but in one way only, and the way it does so changes the question.

- [derived, checked to 3.6e-15 against the full chain] At second order every cusp bond has V₃ = −(2/3)|t₂| exactly (g = 2/3 is rigid), and v = V₂/|t₂| + 2 = |d|²/|p|² ≥ 0, where d is the inversion-odd facilitated amplitude, a transverse Dzyaloshinskii–Moriya term Γ(XZ − ZX). In the parity class v = 0 identically, so no parity-class bond unbinds the pair at second order. Whenever v > 0 the reduction also contains an uncancellable number-changing vertex r Σ_j P b_j†(n_{j+2} − n_{j−2}) P + h.c. with |r| = √v |t₂| exactly. So the cusp family with v > 0 is a one-parameter extension of C173's model in which particle number is not conserved.
- [numerical, exact ring ED of that model, N-converged, calibrated] The vertex weakens the binding instead of restoring it. At u = 1.035 the K = π pair, bound by about 0.015 at v = 0.05 without the vertex, is bound by much less with it, and for v ≥ 0.1 no pair bound by more than about 0.002 survives (the number-conserving model is still bound at v = 0.1 and 0.2, boundary 0.268). The dressed valleys stay two and can be made degenerate by retuning φ, and no pair or cluster lies below the single-quasiparticle onset (checked at v = 0.1 and 0.3). The registered prediction (that the vertex rebinds) failed.
- [numerical, C175's criterion] The number-conserving baseline past the boundary is a dilute two-component liquid by C175's criterion (Δ_σ = 0.60 to 0.95, easy-plane), with no trimer at u = 1, v = 0.3, 0.4.
- Not established: that the resulting gas is c = 2. The pseudospin lengths of the model with the vertex, the many-body central charge, and the full qubit chain were not computed in this session. Two Typhon jobs are pending: backfill 17604628 (the decisive v = 0.3 block DMRG at L = 96, 192 with C173's c = 1 and c = 2 controls) and 17604527 (the rest of the block DMRG and the larger-N_max ring runs at four values of u).
- [derived, scaling] In the parity class beyond second order, v grows as 3.5 (Δ/V)² and g drops, so the number-conserving part does reach the boundary (near V/Δ ≈ 8 at u = 0.5). But that only happens where the momentum-zero ΔN = 2 vertex of amplitude (4/3)(Δ/V) is detuned by only u²/2. Its second-order size, about 0.25 |t₂|, dwarfs the boundary scale v_c ≈ 3.5(Δ/V)², so the parity-class route is not controlled. As u → 0 it turns into C153's density-locking vertex.
- At u = 2 the naive boundary needs r = 1. There the vertex flattens the dressed band between the bare valleys and the two-valley structure is lost (N_max = 3, not converged), so the useful window is u ≲ 1.5.

So the obstruction asked about in the spec does not exist at the two-body level. The exact obstruction at second order is a different one, the forced √v vertex. It turns out to help, and it is what any further test has to include.

## Best candidate bond

The transverse-DM cusp bond (pauli_bond.py, candidate_bond.out), in h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) + Γ_xz XZ + Γ_zx ZX, at Δ = J_x − J_y = 1, V/Δ = 40, u = 1.035, v = 0.3 (r = 0.55), φ = 1.5881 (ring-tuned valley degeneracy), μ_eff = onset + 0.03 = −2.3244:

    J_x = 0.49895031, J_y = −0.50104969, D = 0.0064695034, J_z = 10,
    h_z = −19.978903, Γ_xz = 0.13693064, Γ_zx = −0.13693064.

Energy scale |t₂| = 0.0125. The same at V/Δ = 16: J_x = 0.49737578, J_y = −0.50262422, D = 0.016173759, J_z = 4, h_z = −7.9472566, Γ_xz = −Γ_zx = 0.13693064. Γ_xz = −Γ_zx = Δ√v/4 in general. The φ and h_z values are second-order tunings. A full-chain run must retune them (h_z sets the density in a window of order 1e-3, and higher orders shift the valley degeneracy). This bond is a candidate for the premise "unbound two-valley gas", not a demonstrated c = 2 point.

## Key claims

1. [proved, algebra; numerical check sw_partA.out] The second-order effective Hamiltonian of the most general nearest-neighbour bond near the Fibonacci point is PROOF eq. (2.3). Which terms do what: the pairing (XX − YY, XY + YX) gives the NNN hop, V₂ = −2|t₂| and V₃ = −(2/3)|t₂|. The symmetric off-diagonal exchange XZ + ZX, YZ + ZY and the transverse fields (amplitude s, λ) only renormalize the NN hop, the on-site energy and a uniform ΔN = 1 term, which λ cancels. The antisymmetric (transverse DM) part d adds δV₂ = |d|²/(2V) and nothing at distance 3. The ΔN = 1 vertex has r = p d*/(2V). There is no ΔN = 2 term at this order. The exact second-order matrix of the full chain agrees with the closed form to 3.6e-15 on 18 open chains and rings, L = 8 to 12.
2. [numerical, sw_partB.out, sw_partB_ring.out] All-orders des Cloizeaux coefficients of a generic bond at V = 20 to 160 deviate from (2.3) by amounts falling as 1/V. The ΔN = 2 elements appear at relative order ε/V (third order).
3. [derived] Rigidity: g = 2/3, v = |d|²/|p|², |r| = √v|t₂|, v = 0 in the parity class (PROOF 3).
4. [proved] The vertex breaks C173's Σ. An antiunitary Θ′ = I·K·R_z(π) survives and preserves each valley, so valley degeneracy is a codimension-one tuning of φ (PROOF 3).
5. [numerical, sw_parity.out] Parity class beyond second order: v ≈ 3.5(Δ/V)², 2/3 − g ≈ 3.7(Δ/V)², V₄ ≈ +3.3(Δ/V)², y₃ = (4/3)(Δ/V) (PROOF 5 table).
6. [derived, scaling] Parity class: unbinding only at u ≲ 3.7 Δ/V, where the ΔN = 2 vertex is near resonance with second-order size 0.25|t₂| ≫ v_c (PROOF 5).
7. [numerical, calib_r0.out] The few-body code reproduces the collaborator's E_b exactly and the boundary v_c = u²/4 at u = 0.5, 1, 2.
8. [numerical, dressing.out, onset.out, band_u2v1_N3.out] The vertex dresses valley b more than valley a (it hybridizes b with the ab-pair channel at total momentum near π), splitting the valleys at φ = π/2. At u = 2, v = 1 the dressed band loses its two valleys.
9. [numerical, pb2_r0.out, pairconv_v*.out, true_threshold.out, PROOF 6.5 table] The pair unbinding with the vertex at u = 1.035, L = 48 rings, N-converged: E_pair − T_true = −0.0046, +0.0065, +0.0145, +0.0154, +0.0150 at v = 0.05, 0.1, 0.2, 0.3, 0.4, against −0.0208, −0.0148, −0.0051, +0.0021, +0.0075 without it. Boundary 0.05 < v* ≲ 0.1, against 0.268, within the ring resolution of about 0.002 in E_b.
10. [numerical, band_u1035v01_N3.out, band_u1035v03_N5.out] No state below the quasiparticle onset at any momentum at v = 0.1, 0.3. Two degenerate dressed valleys after retuning φ.
11. [numerical, pseudospin_r0.out, trimer_r0.out] The number-conserving baseline past the boundary is two-component by C175's criterion, with no trimer (u = 1, v = 0.3, 0.4).
12. [derived] With Σ broken the dressed valleys have k_a + k_b ≠ π (mismatch 0.022 at v = 0.3), so C179's two-pair vertex carries momentum 0.045 and is off resonance. All ΔN-odd vertices carry momentum for incommensurate valleys (PROOF 7).
13. [heuristic] The mechanism of the weakening is dressing blockade: two nearby quasiparticles lose part of their dressing energy (PROOF 6.5).

## Report against the registered predictions

- R1: passed (3.6e-15; deviations ∝ 1/V).
- R2: derived as registered.
- R3: calibration passed. The N_max-truncated infinite-line criterion was abandoned as registered in the first addendum (not size-consistent). The ring onset (first addendum): at u = 1 the first object was a dressed single quasiparticle for v = 0.1 to 0.6 at L = 20, not the multi-particle state I predicted. The prediction failed, though L = 20 is coarse. The second addendum's grid-threshold B met "unbound indicated" at N_max = 3, 4 for v = 0.1, 0.2, but it moved strongly with N_max. The third addendum's N-extrapolated reading gives unbound at v = 0.3 (−0.0049) and 0.4 (−0.0122) and marginal at v = 0.2 (+0.0011). That reading still used the grid threshold. Comparing with the true off-grid threshold (added after seeing the non-monotonic v = 0.1 value, so not registered) gives the table of claim 9, and it agrees with the registered reading at v = 0.3, 0.4. My registered prediction (the vertex rebinds) failed.
- R4: confirmed in the form predicted (claim 6).
- R5: reached only partly (C175's criterion for the baseline). The many-body and full-chain tests are in the pending Typhon job or left open.

## C177's two-body fit against these coefficients (coordinator's question)

C177 fitted V₂, V₃ of a number-conserving two-body model to full-chain two-magnon levels at K = π (twob_fit.out) and found v ≈ 20(Δ/V)² (v = 0.267 at V/Δ = 8, pair unbound there). My des Cloizeaux coefficients gave v ≈ 3.5(Δ/V)². The factor of about 6 is not the ΔN = 2 mixing. It comes from the chemical potential at which each number was taken.

- [numerical, sw_mu.out; open L = 11, u = 1, φ = π/2] Rerunning the des Cloizeaux extraction at C177's μ_eff = −12 reproduces their fit. At V/Δ = 8, 16, 40 I get V₂/|t₂| = −1.722, −1.917, −1.986 (C177 model B: −1.723, −1.917, −1.986) and V₃/|t₂| = −0.519, −0.621, −0.659 (C177: −0.496, −0.616, −0.658). At the dilute onset μ_eff = −2.25 the same code gives V₂/|t₂| = −1.954, −1.988, −1.998. So v(Δ/V)² = 17.8, 21.1, 22.3 at μ_eff = −12 against 2.9, 3.0, 3.1 at μ_eff = −2.25. The resonance function D is +0.22 (unbound side) at V/Δ = 8, μ_eff = −12, and −0.49 (bound side) at V/Δ = 8, μ_eff = −2.25.
- [derived, matches the numbers] Mechanism. The virtual pair that generates every second-order term costs kV − 2μ_bare, with μ_bare = (μ_eff + 6)|t₂| at leading order, and k = 1, 2, 3 for the vacuum bond, the NNN hop and the separation-3 bond. Expanding 1/(kV − 2μ_bare) gives each bond energy −(p²/kV)[1 + (μ_eff + 6)(Δ/V)²/k], and the counting of PROOF 2(a) then gives V₂/|t₂| + 2 → −2(μ_eff + 6)(Δ/V)² from this term. The predicted slope is dv/dμ_eff = −2(Δ/V)². Measured: −1.98(Δ/V)² at V/Δ = 40 and −1.84 at 16. This is an energy (Brillouin–Wigner) dependence of the number-conserving coefficients, present at fourth order whenever μ_bare is of order Δ²/V.
- [derived] The ΔN = 2 mixing cannot account for it. At μ_eff = −12 the four-particle states that the ΔN = 2 vertex couples to sit about 20|t₂| above the two-particle ones (C177 chose that μ for this reason), so its second-order effect is about y₃²/20 ≈ 1e-3|t₂| at V/Δ = 8, against a v shift of 0.23. The des Cloizeaux numbers are the number-conserving block itself and do not contain it either.
- [derived] Which one governs the physical binding. The full chain's h_z enters the two-excitation energies, so the binding is a property of the chain at its actual h_z. For the dilute cusp gas μ_eff sits at the onset, about −2.2 (C177's own many-body target). There the number-conserving coefficients are the μ_eff = −2.25 ones: v ≈ 3(Δ/V)², and D(V/Δ = 8) = −0.49, on the bound side. C177's unbinding at V/Δ = 8 is correct for μ_eff = −12, a chain whose band lies 10|t₂| above zero density, and does not transfer to the gas. At the onset the momentum-zero ΔN = 2 vertex is no longer far off resonance (detuning about u²/2 = 0.5|t₂|, y₃ = 0.17 at V/Δ = 8), and its effect on the binding was not computed. If it acts like the ΔN = 1 vertex of the transverse-DM bond (dressing blockade), it would push toward unbinding [heuristic]. A full-chain two-magnon calculation at the onset h_z, with the four-particle sector included, is the test that decides the parity-class route at V/Δ ≈ 8. It is added to "What is left open".
- My conclusions for the transverse-DM bond are unaffected. They live at second order, where μ only shifts energies uniformly. The fourth-order μ dependence is a relative (Δ/V)²|μ_eff + 6| correction there.

## What is left open

1. The central charge of the dilute unbound gas with the vertex. Collect Typhon jobs 17604628 and 17604527 and read (B) as registered. Then go beyond: converged open chains to L ≥ 256 at the candidate, FES, and the parity-gap analogue (the model has no particle-number symmetry, so the test has to be the valley-flip or entanglement scaling).
2. The pseudospin lengths a_s, a_± of the model with the vertex, which needs a size-consistent infinite-line method (N_max ≥ 6 per pair). Only the r = 0 baseline is known to be easy-plane.
3. The exact boundary v*(u) with the vertex. L = 48 resolves E_b ≳ 0.002 only, and the dressed valleys sit off the ring momenta. The Typhon (A) runs at u = 0.522, 1.531, 2 should be re-read with true_threshold.py rather than the grid B.
4. Clusters with three or more quasiparticles at N_max ≥ 7.
5. The full qubit chain with the transverse-DM bond at V/Δ = 12 to 16, with φ and h_z retuned at that V/Δ. Third-order terms (ΔN = 2 at separations 2 and 3, from p·F·F, measured at relative order ε/V in sw_partB.out) enter there.
6. The parity-class route at V/Δ ≈ 8 with its near-resonant ΔN = 2 vertex was not run. The decisive test is the full-chain two-magnon spectrum at the onset h_z (μ_eff ≈ −2.2) with the four-particle sector kept, not at μ_eff = −12 (see the C177 section).
7. The collaborator's FACILITATED_REPAIR.md (cited in contrib/cusp_validation_20261109/README.md) is not in the folder. From its description it treats the same transverse-DM repair and its density-assisted creation term, and should be compared with PROOF 2 to 3 once available.

## Typhon jobs

Two jobs, both module anaconda3/2025.06 (numba works there), partition typhon. Reading rules as registered above (Typhon (B) block-DMRG reading; the ring (A) outputs should be re-read against the true off-grid threshold with true_threshold.py).

Backfill job 17604628 (split off at the coordinator's request, 2026-11-09): 16 cores (8 slots of 2 threads), 12 h, directory /scratch/lustre/jkudlerflam/c18x/c186b/, tasks in typhon_bf/tasks.txt. It holds the decisive runs: the mixed model with the vertex at u = 1.035, v = 0.3, μ = onset + 0.03 and + 0.06, L = 96 and 192 (bm_driver.py), and C173's estimator controls at the same L, run with C173's block_dmrg.py unmodified (copied as block_dmrg_C173.py) and C173's own T_c1 and T_c2 configurations (hard rods, c = 1; two decoupled XX legs, c = 2; C173 read 0.90 and 1.88 for them at 48 → 96). Status PENDING at submission. Expected outputs: out/bm_*_Lb{32,64}.log, out/ctrl_c{1,2}_Lb{32,64}.log, logs/, and in the job directory bm_u1.035_v0.3_locked_dmu0.03.jsonl, bm_u1.035_v0.3_locked_dmu0.06.jsonl, bd_ctrl_c1_32.jsonl, bd_ctrl_c1_64.jsonl, bd_ctrl_c2_32.jsonl, bd_ctrl_c2_64.jsonl. Collect, from this folder, with

    mkdir -p typhon_out_bf && scp -r typhon:/scratch/lustre/jkudlerflam/c18x/c186b/out typhon:/scratch/lustre/jkudlerflam/c18x/c186b/logs typhon:/scratch/lustre/jkudlerflam/c18x/c186b/bm_u1.035_v0.3_locked_dmu0.03.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186b/bm_u1.035_v0.3_locked_dmu0.06.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186b/bd_ctrl_c1_32.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186b/bd_ctrl_c1_64.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186b/bd_ctrl_c2_32.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186b/bd_ctrl_c2_64.jsonl typhon_out_bf/

Main job 17604527: one 48-core node (24 slots of 2 threads), 24 h, directory /scratch/lustre/jkudlerflam/c18x/c186/. Trimmed after the split: the four v = 0.3 block-DMRG tasks were removed from its task list (/scratch/lustre/jkudlerflam/c18x/c186/typhon/tasks.txt, now 60 tasks; the local copy typhon/tasks.txt matches), so nothing is duplicated. It keeps (A), the ring binding at N_max = 4 and 5 for u = 0.522, 1.035, 1.531, 2, and the remaining block-DMRG runs (v = 0.1, 0.2 with the vertex, and the r = 0 points v = 0 bound and v = 0.4 unbound). Status PENDING (Priority), estimated start 2026-11-09 16:21 at the last check. Expected outputs: out/pb_u{0.522,1.035,1.531,2.0}_f{0.25,0.5,1.0,1.5,2.0}_N{4,5}.jsonl, out/pb_u*_zero.jsonl, out/bm_*.log, logs/, and in the job directory bm_u1.035_v0.0_zero_dmu0.03.jsonl, bm_u1.035_v0.0_zero_dmu0.06.jsonl, bm_u1.035_v0.4_zero_dmu0.03.jsonl, bm_u1.035_v0.4_zero_dmu0.06.jsonl, bm_u1.035_v0.1_locked_dmu0.03.jsonl, bm_u1.035_v0.1_locked_dmu0.06.jsonl, bm_u1.035_v0.2_locked_dmu0.03.jsonl, bm_u1.035_v0.2_locked_dmu0.06.jsonl. Collect, from this folder, with

    mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c18x/c186/out typhon:/scratch/lustre/jkudlerflam/c18x/c186/logs typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.0_zero_dmu0.03.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.0_zero_dmu0.06.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.4_zero_dmu0.03.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.4_zero_dmu0.06.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.1_locked_dmu0.03.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.1_locked_dmu0.06.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.2_locked_dmu0.03.jsonl typhon:/scratch/lustre/jkudlerflam/c18x/c186/bm_u1.035_v0.2_locked_dmu0.06.jsonl typhon_out/

Files of this run: sw_general.py, sw_parity.py, pauli_bond.py, mixed_fewbody.py, calib_r0.py, dressing.py, ring_mixed.py, check_ring.py, check_trunc.py, onset.py, pairbind.py (superseded), pairbind2.py, pairbind5.py (unused), pair_conv.py, true_threshold.py, band_diag.py, channels.py, trimer_r0.py, pseudospin_r0.py, blockmix.py, check_blockmix.py, bm_driver.py, dmrg_d.py (C173's, unmodified), and their .out files.
