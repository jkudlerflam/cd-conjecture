# C175 derivations: the valley-pseudospin sign of the Fibonacci-cusp gas

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

Standard: theoretical physics. Tags [proved], [derived], [numerical], [heuristic], [assumed], [literature] open the paragraph or displayed claim they apply to. Every numerical claim names the file it came from (all in this folder).

## 0. Setting and conventions

Model (3.2) of C173 in units |t₂| = 1: hard-core particles with no two adjacent, hop amplitudes a_d for a move by +d, a_{+1} = w, a_{−1} = w̄, a_{±2} = −1 (the move by 2 is over an empty site, so particles never pass), interactions V₂ = −2 at distance 2 and V₃ = −2/3 at distance 3. Band ε(k) = Σ_d a_d e^{−ikd} = 2|w| cos(k − φ) − 2 cos 2k. At φ = π/2 the minima are at sin k = −|w|/4, k_a + k_b = π (mod 2π), ε_min = −2 − |w|²/4, 1/m ≡ ε″ = 8 − |w|²/2, and ε‴(k_a) = −ε‴(k_b) ≠ 0 [proved, algebra].

Two particles x₁ < x₂, r = x₂ − x₁ ≥ 2, ψ(x₁, x₂) = e^{iKx₁} f(r):

    (h f)(r) = Σ_d a_d [ e^{−iKd} f(r + d) + f(r − d) ] + V(r) f(r),   f(r) = 0 for r < 2.      (0.1)

A plane wave f = e^{ipr} has right momentum p, left momentum K − p and energy ε(K − p) + ε(p). Code: scat.py. The N-body generalization in gap coordinates is nbody_box.py, the ring ED with momentum resolution is ring_k.py.

## 1. Impenetrable two-valley gas and the pseudospin chain

[derived; inputs: Girardeau, the excluded-volume (squeezed-space) form of the strong-coupling 1D mixture (C162 PROOF 1.2), generalized here to two channels] At low density, if no two-body bound state exists, every pair is impenetrable at threshold and the order of particles is fixed. Give the particle at position x_j the fast phase e^{ik_{σ_j}x_j}, σ_j ∈ {a, b}, times a slow envelope. The envelope is a Tonks gas of total density n in the length L − Σ_bonds a_bond, so E = π²N³/(6m(L − Σ a)²) and each bond contributes c·a_bond, c = π²n³/(3m). For a bond between valleys (σ, σ′) the excluded length is an operator on the ordered channels. For a like pair it is the scalar a_s (a_aa = a_bb by Σ). For an unlike pair it is the 2 × 2 matrix A on (|ab⟩, |ba⟩), defined from the zero-energy solution in the outer region

    f(r) = e^{ik_b r}(α_ab + β_ab r) + e^{ik_a r}(α_ba + β_ba r),   α = −A β.      (1.1)

The amplitudes multiply exactly the fast phases e^{ik_a x₁ + ik_b x₂} and e^{ik_b x₁ + ik_a x₂}, so A acts directly on the pseudospin of the ordered pair. Writing the bond as an XXZ term,

    H_ps = Σ_k [ J_⊥ (S^x_k S^x_{k+1} + S^y_k S^y_{k+1}) + J_z S^z_k S^z_{k+1} ],
    J_z = 2c (a_s − A_d),   J_⊥ = 2c |A_o|,   Δ_σ = (a_s − A_d)/|A_o|,      (1.2)

with A_d = A₁₁ = A₂₂ (by the antiunitary Θ that maps ab ↔ ba) and A_o = A₁₂. The phase of A_o and the sign of J_⊥ are removed by a rotation about z on the squeezed (open) chain. Then c = 2 iff −1 < Δ_σ < 1, i.e. A_d − |A_o| < a_s < A_d + |A_o|. Below that range the gas is valley ferromagnetic (vector chiral, one valley), above it valley-Néel (valley crystal). This coincides with C180's criterion (C180 RESULT claim 5) [derived independently].

Kinematic anti-Hermitian part [derived, checked numerically]. The matrix A obtained from (1.1) is not Hermitian. Its anti-Hermitian part is diag(−iγ, +iγ) with γ = ε‴/(6ε″) (0.0630, 0.1291, 0.2887 at |w| = 0.5, 1, 2, equal to ε‴/(6ε″) at the valley), the same with and without V₂, V₃ (scat.py run, recorded in PROGRESS.md). It comes from the cubic term of the relative dispersion 2ε(p), which enters the conserved discrete Wronskian as a term D|β|² (so A − A† = −D/C is fixed by kinematics). It is valley-odd, cancels against the valley-odd cubic dispersion term of the envelope, and drops out of (1.2). The physical matrix is A_H = (A + A†)/2.

Checks [numerical]. (i) Pure nearest-neighbour hard rods (no NNN hop, no V's): a = 1 exactly (excluded length of rods of length 2). (ii) The V₂ = V₃ = 0 limit at |w| = 1, which is the V → ∞ limit of Mishra, Greschner and Santos (arXiv:1411.7270, per C180), gives a_s = 0.629, A_d = 0.500, |A_o| = 0.516, window (−0.016, 1.016), Δ_σ = 0.25, two-component, as the coordinator and C180 expect (calib.out). The same holds at |w| = 0.5 and 2 (Δ_σ = 0.125, 0.50).

## 2. Exact threshold data at φ = π/2

### 2.1 The unlike pair in closed form

[proved, algebra; analytic_A.py, analytic_res.py] At K = π the gauge f(r) = i^r g(r) makes (0.1) real,

    2g(r+2) + 2u g(r+1) + 2u g(r−1) + 2g(r−2) + V(r) g(r) = E g(r),   u = |w|.      (2.1)

At E = 2ε_min = −4 − u²/2 the characteristic equation is (2x + u)²/2 = 0 in x = z + 1/z, so the free modes are e^{±iκr} and r e^{±iκr} with cos κ = −u/4, a double root at each valley. The boundary rows (r = 2, 3) are −2u g(1) − 2g(0) + V₂ g(2) = 0 and −2g(1) + V₃ g(3) = 0 with g continued by the free form. Solving them for the 4 mode coefficients gives, for general V₂, V₃,

    det(α block) ∝ √(16 − u²) · D,   D = 2V₂V₃ + 4V₂ + V₃u² + 4V₃ + 8,      (2.2)

and the eigenvalues of A are

    a_∓ = [10V₂V₃ + 12V₂ + 5V₃u² + 12V₃ + 8 ∓ √(...)] / (2D).      (2.3)

At the Schrieffer-Wolff values V₂ = −2, V₃ = −2/3,

    A_d = 5/2 + 8/u²,   |A_o| = √(u⁴ + 160u² + 256)/(2u²),   a_∓ = [5u² + 16 ∓ √(u⁴ + 160u² + 256)]/(2u²),      (2.4)

matching the numerical solution of scat.py to all printed digits (scan_w.out: at u = 1, A_d = 10.5, |A_o| = 10.211, a_− = 0.2889, a_+ = 20.711). Small u: a_+ ≈ 16/u², a_− ≈ 3u²/8.

### 2.2 An exact zero-energy resonance at V₂ = 2t₂, detuned by V₃

[proved, algebra] With V₃ = 0, D = 4(V₂ + 2), so D = 0 exactly at the Schrieffer-Wolff value V₂ = 2t₂ = −2, for every |w|. The unlike pair at K = k_a + k_b then has a zero-energy resonance (A has a pole, calib.out shows A_d and |A_o| of order 1e13 to 1e16 there). With V₃ switched on, D = V₃u² + (V₂ + 2)(2V₃ + 4), which at V₃ = −2/3 equals −2u²/3 ≠ 0.

[heuristic, interpretation] At w = 0 each leg of (3.2) is the ferromagnetic Heisenberg chain 2t₂ S·S exactly (C173 PROOF 3.3), whose two-magnon state at zero leg momentum is an SU(2) descendant, a threshold state with infinite scattering length. The unlike pair at K = π is built from such same-leg pairs. The calculation shows that the resonance survives the inter-leg hop w at any |w| as long as V₃ = 0, and that the attractive V₃ = (2/3)t₂ moves the pair to the bound side.

[derived] Robustness. The detuning −2u²/3 is of order one in units of t₂. Higher Schrieffer-Wolff orders change V₂, V₃ and the hops by relative amounts of order (Δ/V)² (C173 PROOF 3.2: error ∝ V⁻³), 6e-4 at V/Δ = 40, and the resonance moves only when D changes sign, i.e. when V₂ + 2 reaches u²/4 (0.25 at u = 1). So the binding is a robust feature of the cusp, not an artifact of the leading order.

### 2.3 The bound pair

[numerical, scan_w.out, phi_scan.out, ring2.out] At φ = π/2 the lowest two-body level at K = π lies below the absolute two-particle threshold 2ε_min for every |w| in (0, 4):

| \|w\| | a_s | A_d | \|A_o\| | a_− | a_+ | Δ_σ (formal) | E_b = 2ε_min − E₂ | ε″/a_+² |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 0.25 | −9.376 | 130.50 | 130.48 | 0.023 | 260.98 | −1.072 | 1.18e-4 | 1.17e-4 |
| 0.5 | −5.282 | 34.50 | 34.41 | 0.087 | 68.91 | −1.156 | 1.71e-3 | 1.66e-3 |
| 1.0 | −3.173 | 10.50 | 10.21 | 0.289 | 20.71 | −1.339 | 1.95e-2 | 1.75e-2 |
| 1.5 | −2.501 | 6.056 | 5.542 | 0.514 | 11.60 | −1.544 | 6.36e-2 | 5.11e-2 |
| 2.0 | −2.295 | 4.500 | 3.786 | 0.714 | 8.286 | −1.795 | 0.1225 | 0.0874 |
| 3.0 | −3.567 | 3.389 | 2.410 | 0.979 | 5.799 | −2.887 | 0.195 | 0.104 |
| 3.5 | −15.42 | 3.153 | 2.181 | 0.972 | 5.335 | −8.51 | 0.164 | 0.066 |

- The bound state is shallow and universal at small |w|: E_b ≈ ε″/a_+² ≈ |w|⁴/32, i.e. it is the bound state that the large positive a_+ of (2.4) announces. Its size is a_+, about 21 sites at |w| = 1 and 8 at |w| = 2.
- It sits at total momentum K = k_a + k_b = π, the pair of one particle in each valley. Over all K the minimum of the two-body energy is at K = π (bound_scan in scat.py, table in PROGRESS.md).
- Independent check [numerical, ring2.py, ring2.out]: brute-force two-particle ED on rings (full two-particle basis, separate code) gives E − 2ε_min = −0.01958, −0.01953 at L = 120, 240 (|w| = 1), −0.12255 at |w| = 2, and −0.00175 at |w| = 0.5, L = 240, in agreement with the relative-coordinate value.
- [derived from (2.2)] The bound state exists on the whole interval 0 < |w| < 4, not only on the grid. A 1D bound state can leave the spectrum only through threshold, where the scattering length has a pole (D = 0 in (2.2)). At the Schrieffer-Wolff values D = −2u²/3, which vanishes only at u = 0 (the w = 0 resonance), and the valley structure exists up to u = 4. So the bound state found on the grid 0.25 to 3.9 persists for all 0 < |w| < 4.
- The like pair (K = 2k_a) has a negative scattering length a_s (an attractive p-wave-like channel, no shallow bound state) for |w| < 3.6, and the all-K scan at |w| = 0.5, 1, 2 finds no like-pair bound state. a_s diverges at |w| ≈ 3.63 (a_s = −73.4 at 3.6, +26.2 at 3.7), above which like pairs bind too.

### 2.4 Which interaction does what (coordinator's calibration)

[numerical, calib.out] Switching V₂ and V₃ on separately:

| \|w\| | V₂, V₃ | a_s | window (a_−, a_+) | Δ_σ | Tonks verdict | lowest two-body level − 2ε_min |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | 0, 0 | 0.629 | (−0.016, 1.016) | 0.25 | c = 2 | +3e-5 (box) |
| 1 | −2, 0 | −1.437 | resonance (a_+ = ∞) | 1.00 | boundary | +1e-5 (box) |
| 1 | 0, −2/3 | 0.023 | (−1.358, 0.644) | 0.38 | c = 2 | +3e-5 (box) |
| 1 | −2, −2/3 | −3.173 | (0.289, 20.71) | −1.34 | vector chiral | −0.0195 (bound) |
| 2 | 0, 0 | 0.789 | (−0.077, 1.077) | 0.50 | c = 2 | +3e-5 |
| 2 | −2, −2/3 | −2.295 | (0.714, 8.286) | −1.80 | vector chiral | −0.1225 (bound) |

So V₂ alone puts the unlike channel at its exact resonance (2.2), V₃ alone leaves the gas two-component, and the two together bind the unlike pair and, among scattering states, put the same-valley excluded length far below the window (vector-chiral side). C180's Born-level indicator had the right direction for the combination, but the exact answer is dominated by the bound state, which no Born estimate sees.

## 3. Consequence in the dilute limit

[derived] Because a two-body bound state lies below 2ε_min at K = π, the premise of section 1 (no bound state, every pair impenetrable) fails at φ = π/2, and the threshold law argument of C173 PROOF 4(b) is incomplete there: inter-valley scattering states are impenetrable, but the ground state is not made of scattering states. As n → 0 at fixed μ − ε_min, the particles bind into ab pairs (molecules of size a_+ with momentum π), and a pair costs ε_min − E_b/2 per particle against ε_min for a single-valley particle. The resulting state is a liquid of valley-singlet pairs: the relative valley mode φ_− is locked, valley-imbalance excitations cost a pair-breaking gap of order E_b, and c = 1. In C180's language this is the bound-state (ab pair) liquid, the Luther-Emery analogue of the valley crystal (C180 RESULT claim 6).

[numerical, formal] Even if pairing were ignored, the scattering-state excluded lengths put the gas on the vector-chiral side, Δ_σ < −1, for every |w| < 3.6 (table in 2.3). So neither instability of the dilute gas is easy-plane. The answer to "is |J_z| < J_⊥?" is no, with J_z < 0 (ferromagnetic) among scattering states, and the bound state pre-empts both by pairing.

## 4. Larger clusters

[numerical, nbody3.out, nbody3b.out, nbody45.out, nbody4_w1_R110.out] Exact few-body ground states on the infinite line (gap coordinates with hard walls at gap R, converged in R where stated):

| \|w\| | N | K | E − Nε_min | nearest threshold | bound? |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 2 | 2 | π | −0.12255 | 0 | yes, E_b = 0.1225 |
| 2 | 3 | π + k_a | −0.14082 (R = 100, 150, 220 identical) | dimer + atom −0.12255 | yes, by 0.0183 |
| 2 | 4 | 0 | −0.25408 (R = 60) | two dimers −0.2451 | yes, by ≥ 0.0090 |
| 2 | 5 | k_a | ≤ −0.28482 (R = 36, still falling) | dimer + trimer −0.2634, tetramer + atom −0.2541 | yes, by ≥ 0.021 |
| 1 | 3 | π + k_a | −0.020123 (R = 400, converged to 1e-5) | dimer + atom −0.019527 | yes, by 6.0e-4 |
| 1 | 4 | 0 | −0.03689 (R = 110, still falling: −0.0322, −0.0351, −0.0369 at R = 60, 80, 110) | two dimers −0.03905 | not resolved |

So atoms bind to dimers (trimers) at |w| = 1 and 2, and dimers bind to each other at |w| = 2 (tetramer, pentamer). [heuristic] At |w| = 2 the dilute dimer gas is unstable to clustering, so the onset from the vacuum is probably a weak first-order jump into a self-bound pair liquid. The binding per particle grows only slowly (0.061, 0.047, 0.064, ≥ 0.057 for N = 2 to 5), with a strong even-odd effect (adding the particle that completes a pair gains about 0.11, an odd one 0.02 to 0.03), the signature of a liquid of pairs. A ring ED of the same clusters (ring_clusters.out, ring_k.py) confirms the binding qualitatively, with finite-ring oscillations of a few 1e-3 at L ≤ 62 (cluster sizes of order 20 sites).

## 5. Off the symmetric line, and the c = 2 window at small n

[numerical, phi_scan.out] Off φ = π/2 the unlike pair stays bound relative to ε_a + ε_b (by 0.0195 at |w| = 1 for every φ in [1.54, π/2]), but it lies below 2ε_a only while the valley offset δ = ε_b − ε_a is below E_b: for |φ − π/2| ≲ 0.005 at |w| = 1 (δ = 0.011 at φ = 1.568, 0.022 at 1.565) and ≲ 0.017 at |w| = 2. No like-pair bound state appears there (lowest level 8e-4 to 1e-3 above 2ε_a).

[derived at leading order in n, with the cluster data of section 4 and C164's mass-imbalance result as inputs] Dilute phase diagram near φ = π/2:
- |δ| < E_b: pair (and cluster) liquid, valley pseudospin gapped, c = 1.
- |δ| > E_b: single-valley gas, c = 1 (C173 PROOF 4).
- The intermediate regime, a mixture of dimers and excess lower-valley atoms (the analogue of the partially polarized phase of the attractive Yang–Gaudin gas, which is two-component there), is not a two-component gas here. Atoms bind to dimers (trimers bound), and an impenetrable mixture of composites of unequal mass is immiscible as n → 0 (C164 finding 3, interface cost of order η²n² against exchange of order n³). Yang–Gaudin evades this only because it is integrable (reflectionless atom-dimer scattering), which (3.2) is not (C180's integrability check).

So the predicted c = 2 window in (φ, μ) at small n is empty. This replaces C173's wedge |δ| ≲ J_⊥(n) ~ n³, which assumed a gapless easy-plane pseudospin at δ = 0.

## 6. The density of the candidate, n ≈ 0.17

The dilute expansion does not reach n = 0.17 at |w| = 1: n a_+ = 3.5 (pairs overlap strongly) and n|a_s| = 0.54. What can be said:

(a) [numerical, ring_eos_w1.jsonl, ring_eos_w1b.jsonl; momentum-resolved ring ED of (3.2), L = 12 to 36, N = 1 to 8] The vector-chiral (valley-ferromagnetic) state is excluded at n ≥ 0.13. For every even N the ground state is at K = Nπ/2 exactly (the balanced sector, a condensate of pairs at momentum π), and the lowest states in the fully polarized sectors lie 0.46 to 0.66 higher (N = 4 to 8 at L = 24, 30). That is far more than the ring twist cost of a polarized state (order 0.01), so the scattering-state indicator Δ_σ < −1 does not control this density. The grand-canonical density at μ = −2.2 is N = 4 at L = 24 (n = 0.167), reproducing C173's n ≈ 0.17.

(b) [numerical, same files] The single-valley-flip energy (the sector N_a − N_b = 2, identified by momentum) at n = 1/6 is 0.210, 0.157, 0.117 at L = 12, 24, 36, so L·E_flip = 2.53, 3.77, 4.21. It grows with L, which is what a pseudospin gap does, but the sizes (N = 2, 4, 6) are far too small to separate a gap from 1/L² corrections; two-point fits give gaps from 0.04 to 0.10 depending on which sizes are used. The odd-even energy E(N_odd) − [E(N−1) + E(N+1)]/2 is 0.013 to 0.027 (plus about 0.018 from the curvature of E(N)) for L = 18 to 36, of the order of E_b. Both are consistent with pairing and do not exclude c = 2. [numerical, not decisive]

(c) [heuristic] If pairing is the leading instability at n = 0.17, the regime is the BCS side of a 1D pair crossover (crossover density n_* ≈ 1/a_+ ≈ 0.05 at |w| = 1). Using the attractive Yang–Gaudin gap as a template with |γ| = 2/(n a_+) = 0.57, Δ_s ≈ (16/π) E_F √(|γ|/π) exp(−π²/(2|γ|)) ≈ 4e-4 E_F ≈ 1e-4 |t₂|, with pseudospin correlation length ξ_s ≈ v_F/Δ_s of order 1e4 sites. On that estimate C173's L ≤ 96 reading, c ≈ 2, would be a crossover, with c = 1 at L ≫ ξ_s. The template has SU(2) and (3.2) does not. The like-valley channel is strongly attractive (a_s = −3.2), which pushes the pseudospin toward the ferromagnetic side and away from the pairing gap (C180 PROOF 3, larger K_−). So a gapless two-component window between the pair liquid and the vector-chiral state at intermediate density cannot be excluded by anything derived here.

(d) [numerical, smatrix.out; illustrative only] The energy-dependent eigenphases of the two-body S-matrix at relative momentum q give "excluded lengths" a(q) = −θ(q)/(2q). At |w| = 1, q = 0.25 to 0.55 (the Fermi momenta of n = 0.08 to 0.17), the bound channel falls from a_+ ≈ 6.6 to 3.9, a_− rises from 0.38 to 0.65 and a_s stays between −2.5 and −1.8. Inserted naively into (1.2), these give Δ_σ between about −1.8 and −3. The ring ED in (a) contradicts the vector-chiral prediction this would imply, which shows that the finite-q Tonks extrapolation is not reliable at this density.

## 7. Pre-registered checks, reported against

- P1 (criterion). The criterion is met in its second part (no bound state) nowhere at φ = π/2. Verdict as registered: the dilute gas is not two-component.
- P2(i) (Hermiticity of A to 1e-8). Failed as stated. The anti-Hermitian part is diag(∓iγ), γ = ε‴/(6ε″), independent of the potential to all printed digits. Its origin is kinematic (section 1). The Hermitian part is used.
- P2(ii) (hard rods a = 1). Passed (1.000 to 1e-12).
- P2(iii) (ring two-particle levels). Passed for the bound state: ring2.py, an independent code, reproduces E_b to 5e-5 at L = 120 and 2e-9 at L = 240 (|w| = 1).
- P3 (ring ED N = 3 to 6 against the Tonks XXZ prediction). Not testable as registered. The registered effective theory assumed no two-body bound state, and section 2.3 shows there is one, so its predicted spectra do not apply. The ED instead confirms the pair structure (ground states at K = Nπ/2 for even N, polarized sectors far above, odd-even energies of order E_b, N = 2 and 3 cluster energies reproduced on large rings). No other ring comparison was attempted.

## 8. What would settle n ≈ 0.17

Predictions of the pairing picture, registered for C176 and C178 (and for this run's Typhon job, RESULT.md):
1. At |w| = 2, where the pairs are compact (a_+ = 8.3, E_b = 0.12), the φ = π/2 effective model should read c ≈ 1 at L = 96 to 192 and at densities n ≈ 0.1 to 0.2, with the two-size estimate falling with L. The same model with V₃ = 0 (exact resonance, no binding) and with V₂ = V₃ = 0 (two-component by section 2.4 and MGS) are the controls.
2. At |w| = 1, close to threshold (n ≲ 0.05, below n_*), the reading should be c = 1 already at L ≈ 200 (pair size 21, pair gap of order E_b = 0.02).
3. At n = 0.17, |w| = 1: open-chain pair correlations ⟨(c_j c_{j+r})† c_0 c_r⟩ at momentum π should decay more slowly than the single-particle Green's function, which should eventually decay exponentially, and the parity gap E(N+1) + E(N−1) − 2E(N) should stay finite as L grows. If instead the valley pseudospin is gapless, the parity gap scales as 1/L and the two-valley single-particle correlations are power laws.
