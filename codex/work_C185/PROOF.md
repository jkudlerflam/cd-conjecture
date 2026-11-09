# C185 PROOF: degenerate points where second-order dynamics can win

Conventions. A bond of B₂ is h = Σ_ab J_ab σ^a ⊗ σ^b + Σ_a B_a (σ^a ⊗ 1 + 1 ⊗ σ^a)/2, so the field per site is B. For the Fibonacci point I use the bond convention of the spec, h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z), particle = spin down, n = (1 − Z)/2, Δ = J_x − J_y, V = 4J_z, |t₂| = Δ²/(2V), as in C173 and C179.

## 1. What is being classified

The mechanism of the Fibonacci cusp has three ingredients.
- (a) A dominant bond h₀ ∈ B₂ whose chain Hamiltonian H₀ has a ground space G_L of dimension growing exponentially in L, separated from the rest of the spectrum by a gap of order ‖h₀‖ (or at least by a scale much larger than the dynamics inside G_L).
- (b) Perturbations δh whose first-order projection P δH P onto G_L can be tuned to vanish or to be small, on a set of finite codimension.
- (c) The remaining perturbations generate, at second order, H₂ = −P δH Q (H₀ − E₀)⁻¹ Q δH P, a nontrivial many-body model on G_L.

If dim G_L grows only polynomially, (c) produces a model on a space with zero entropy density, and low-energy extensive physics must come from excitations of H₀ itself (section 2.5). So the classification target is the set of bonds with exponentially degenerate ground spaces. I call these extensive points.

## 2. Classification of extensive points in B₂

### 2.1 Classical bonds [proved]

A classical bond is diagonal in a product basis that is the same on every site. After the common rotation that diagonalizes it, h = J ZZ + b (Z ⊗ 1 + 1 ⊗ Z)/2 plus a constant (a translation-invariant chain cannot carry a staggered field, and Z ⊗ 1 − 1 ⊗ Z telescopes). Its four bond energies are E(++) = J + b, E(+−) = E(−+) = −J, E(−−) = J − b. The open-chain ground states are the minimal-weight paths of the 2 × 2 transfer graph whose edges are the minimal bond configurations, provided the minimum of each bond can be attained simultaneously (it can, since the minimal edge set always contains a closed walk). The ground space is exponentially degenerate if and only if the minimal edge set contains two distinct cycles through a common vertex. Enumerating:
- J > 0, |b| < 2J: minimal set {+−, −+}, Néel, 2 states.
- J > 0, |b| = 2J (say b = 2J): minimal set {+−, −+, −−}, all strings with no two adjacent +, Fibonacci many (F_{L+2}).
- J > 0, |b| > 2J: one of {++}, {−−}, unique.
- J < 0: minimal set {++, −−} at b = 0 (2 states), otherwise one of them (unique).
- J = 0, b = 0 is the zero bond, excluded by normalization.

So the only classical extensive point is the Fibonacci orbit J > 0, |b| = 2J. Checked by exact counting (classify.out, R1): F_{L+2} = 8, 55, 377, 2584 at L = 4, 8, 12, 16 at b = ±2, and 1 or 2 elsewhere.

"Classical after a local rotation." A uniform rotation gives the same point of B₂. A sublattice rotation of a classical bond gives a translation-invariant chain only if it also maps the field to a uniform field, which for b ≠ 0 forces the rotation to fix the quantization axis up to a sign flip, and a sign flip maps the field to a staggered one. A period-two chain with a twisted translation is unitarily a B₂ chain (C169), and here that chain is again J ZZ + b Z. So no new classical point arises this way.

### 2.2 Frustration-free bonds [proved for products, numerical otherwise, consistent with Bravyi–Gosset]

If h ≥ 0 and the open chain has zero-energy states, the ground space depends only on the kernel S ⊂ ℂ² ⊗ ℂ² of h, G_{L+1} = (G_L ⊗ ℂ²) ∩ (ℂ^{2^{L−1}} ⊗ S). Classification by rank r = 4 − dim S (classify.py, exact rank computations to L = 14, classify.out):
- r = 1, forbidden state ψ = |aa⟩. In the basis {a, a^⊥} this is the classical projector onto |11⟩, so dim G_L = F_{L+2} [proved, 2.1]. Numerically 3, 5, 8, …, 987 for L = 2 to 14, for a = |1⟩ and for random a.
- r = 1, every other ψ (entangled, or |ab⟩ with b not parallel to a): dim G_L = L + 1 exactly for every L tested, for 8 structured and 40 random ψ. This is Bravyi's quantum 2-SAT count on a chain (the product-state ground states are parametrized by one domain-wall position) [numerical, literature].
- r = 2: dim G_L ≤ 2 for every L, and 0 for random kernels (frustrated).
- r = 3: dim G_L ≤ 1.

So the frustration-free extensive points are the Fibonacci orbit only. This reproduces C143's statement, now with an explicit check of every rank.

### 2.3 Degeneracy from local integrals of motion [C143, proved there]

A strictly local integral of motion forces rank J ≤ 1, and at range two such operators exist only at classical bonds and at the XZ point. The XZ point h = X ⊗ Z has the commuting Majorana algebra b_j = Z_j X_{j+1} and an exactly degenerate manifold of dimension 2^{L/2−1} (even rings), on top of a critical Majorana chain. The free locus has a flat band only at XZ.

### 2.4 Unbiased numerical completeness check (R3)

degen_opt.py minimizes, over all 12 couplings (J, B) at unit norm, the width per site of the lowest K_L states on rings of L = 8 and 10, for twelve target counts (K₈, K₁₀) between (4, 6) and (60, 150), including (8, 16) (the XZ count), (47, 123) (the Fibonacci count, Lucas numbers) and (9, 11) (the L + 1 count, a calibration on the polynomial points). Calibration (degen_test.py): W = 0 to 1e-14 at Fibonacci for every K ≤ Lucas, at XZ for (8, 16), and at the FM Heisenberg point for (9, 11). Results: section 6 (Typhon job).

### 2.5 Polynomially degenerate points and why they differ [derived]

The ferromagnetic multiplet (dim L + 1), the rest of the Bravyi–Gosset rank-one set (L + 1), the product projectors |ab⟩ (L + 1) and the Kurmann–Thomas–Müller points (2) are not extensive.
- The degenerate space has zero entropy density. For the ferromagnetic multiplet it is the spin-L/2 representation, and the first-order projection of any two-site perturbation is a quadratic polynomial in the total spin. That is a single large spin (a Lipkin–Meshkov–Glick model), with a unique ground state for generic perturbations and a level spacing of order δ. It carries no extensive low-energy sector.
- The extensive low-energy states at these points are the excitations of H₀ itself, magnons with dispersion of order J k² that is fixed by H₀ at zeroth order. A second-order term of size δ²/J competes with this kinetic energy only at k ≲ δ/J, so it renormalizes the interactions of a dilute single-band gas (C143 stratum 2, one magnon band with a single minimum for every frustration-free bond) and cannot create a second valley.
- To make the magnon band of H₀ flat one must remove its in-plane exchange, which turns h₀ into a classical bond, and then 2.1 applies.

So the second-order mechanism needs exponential degeneracy, and in B₂ that means the Fibonacci orbit or the XZ orbit.

## 3. The Fibonacci point: first-order dynamics and the complete tuned locus

### 3.1 First-order projection [proved, fib_first.out]

Projecting every bond Pauli sum onto the blockaded space P (ring L = 10, residuals below 1e-14):
- XX and YY → T (real NN hop), XY → −T_i, YX → +T_i (imaginary hop), so the hop is w = (J_x + J_y) + 2iD in spin-basis matrix elements;
- X ⊗ 1, 1 ⊗ X, X ⊗ Z, Z ⊗ X → Σ_j PX_jP; the same with Y → Σ_j PY_jP (PXP creation, amplitude Ω);
- Z ⊗ 1, 1 ⊗ Z → −2N, ZZ → −4N (chemical potential).
The reason for X ⊗ Z → PXP is that X_j can act inside P only when both neighbours are empty, where Z = +1.

The 12 couplings (J, B) map onto 5 first-order parameters (w ∈ ℂ, Ω ∈ ℂ, μ). The kernel is 7-dimensional:
- the overall Ising scale along the Fibonacci line (ZZ together with its blockade field);
- the two rotations about x and y (exact symmetries of the spectrum, first-order invisible because P[G, H₀]P = 0);
- the pairing Δ = (J_x − J_y) + i(J_xy + J_yx) ∈ ℂ (two real directions);
- the transverse Dzyaloshinskii–Moriya vector, λ_y (X ⊗ Z − Z ⊗ X) and λ_x (Y ⊗ Z − Z ⊗ Y), written ℓ ∈ ℂ (two real directions).

So [proved] modulo rotations and scale, the first-order-invisible normal space at the Fibonacci point is four-dimensional, spanned by the pairing Δ and the transverse DM vector ℓ. The tuned locus is codimension 5 (w = Ω = 0 and μ near the band), and on it the second-order model is built from Δ and ℓ only. Visible terms tuned to size Δ²/V contribute at second order only at Δ⁴/V³.

### 3.2 Parity-class rigidity [derived]

The parity class (bonds commuting with Π_j Z_j) contains Δ but not ℓ, since X ⊗ Z − Z ⊗ X is odd under Π Z. So on the parity-class tuned locus the only second-order source is Δ, whose phase is removed by the U(1)_z rotation that fixes H₀. The second-order model is therefore unique, C173's eq. (3.2) with (V₂, V₃) = (2t₂, (2/3)t₂), and the only free parameters are w and μ. By C175 and the collaborator's exact boundary (v < g u²/[2(2 − g)] with v = 0, g = 2/3) the opposite-valley pair is bound at φ = π/2 for every 0 < u < 4. So within the parity class the cusp cannot be unbound at second order by any choice of bond. Unbinding needs either ℓ (leaving the parity class) or terms of order Δ⁴/V³ (third and higher order), which are small.

### 3.3 The transverse DM term at second order [derived; checked in 3.4]

Take δh = λ (X ⊗ Z − Z ⊗ X) with λ real (λ_x is the same after a z-rotation). Acting on site j it flips σ^x_j with amplitude A_j = λ(Z_{j+1} − Z_{j−1}), collecting the two bonds that contain j. Write z = 1 − 2n. In general a sum of a field and XZ-type terms gives A_j = h + λ_L z_{j−1} + λ_R z_{j+1}; here h = 0, λ_R = λ, λ_L = −λ. Define its values on configurations of the two neighbours: Ω = A(0,0) = h + λ_L + λ_R, g_L = A(1,0) (left neighbour occupied), g_R = A(0,1), g₂ = A(1,1).
- First order inside P: only A(0,0) acts (a particle can be created only between two empty sites), and Ω = 0 here. So the term is first-order invisible, as 3.1 says.
- Off P, a flip next to one particle creates a blockade violation at cost V, next to two particles at cost 2V. Second order, with the standard sign −⟨f|δH|m⟩⟨m|δH|i⟩/(E_m − E₀):
  - one particle at i loses −(g_L² + g_R²)/V from the two flips at i ± 1, an on-site shift;
  - two particles at distance 2 share the site between them, with amplitude g₂ and cost 2V, so V₂′ = −g₂²/(2V) + (g_L² + g_R²)/V;
  - a NN hop j + 1 → j goes through the virtual pair (j, j + 1): create j (amplitude g_R), then remove j + 1, whose left neighbour is now occupied (amplitude g_L), so hop′ = −g_L g_R/V;
  - distance 3 and NNN hop: no process at this order (the flipped site is adjacent to at most one of two particles at distance 3, exactly as for isolated particles, and a single flip moves a particle by one site).
- For the general linear amplitude with h fixed by Ω = 0 (h = −λ_L − λ_R), g_L = −2λ_L, g_R = −2λ_R and g₂ = −2(λ_L + λ_R), so V₂′ = [4λ_L² + 4λ_R² − 2(λ_L + λ_R)²]/V = 2(λ_L − λ_R)²/V ≥ 0. The symmetric case λ_L = λ_R (an infinitesimal rotation about y) gives V₂′ = 0, as it must, since a rotated Fibonacci point stays exactly degenerate.
- For the transverse DM, g_L = 2λ, g_R = −2λ, g₂ = 0:

  on-site −8λ²/V, V₂′ = +8λ²/V (repulsive), hop′ = +4λ²/V (real), V₃′ = 0, NNN hop unchanged.

In units |t₂| = Δ²/(2V) the repulsion is v = 16λ²/Δ², the hop shift is 8λ²/Δ² (real, so it moves φ off π/2 and must be cancelled by J_x + J_y = −4λ²/V), and the on-site shift is absorbed into h_z.

The mechanism, in words. The pairing-generated V₂ = 2t₂ is attractive because two particles at distance 2 share their virtual-pair channels sub-additively. The DM term creates a particle next to an existing one with an amplitude that is odd under left-right exchange, so between two particles at distance 2 the two contributions cancel and the shared site loses its second-order energy. Losing a negative energy is a repulsion.

- ΔN = ±1 cross terms. The pair Δ σ⁺σ⁺ followed by a DM removal gives a single-particle creation. For an isolated new particle at s the two routes, via the pair (s, s + 1) and via (s − 1, s), carry amplitudes ∝ g_L and ∝ g_R = −g_L and cancel. With another particle at s ± 2 one route is blocked (g₂ = 0 and the higher cost), so a conditioned creation survives, of the form (Δλ/V) Σ_s σ⁻_s (n_{s−2} − n_{s+2}) up to a numerical factor. It is odd in the conditioning distance, so it vanishes on long-wavelength density, and it carries no momentum to make a single valley resonant. Its size relative to |t₂| is of order λ/Δ, which is not small at the λ needed below.

### 3.4 Check against the full chain [numerical]

See sw_dm.py and section 6.

## 4. The XZ point [C143, C144, C157, summarized, no new numerics]

- Degenerate space: the b-Majorana algebra (b_j = Z_j X_{j+1}), 2^{L/2−1} states, on top of the critical A chain (a_j = X_j Z_{j+1}). Not gapped from the rest, so the elimination is controlled only through the kernel of the A chain (C144).
- First order: dJ_xz E₀^A + (dJ_zx + ⟨a⟩ dJ_yy) Σ_j b_j, ⟨a⟩ → −2/π (C144, checked to 1e-13). Tuned locus: the hyperplane dJ_zx = (2/π) dJ_yy, codimension one.
- Second order on the tuned locus: H_B = −2β² Σ_{r≥2} K(r) b_j b_{j+r}, purely quartic, with K positive semidefinite (Lehmann) and peaked at q = π. The induced interaction is attractive in every channel. NN Majorana hopping terms b_j b_{j+1} cancel because the A-chain kernel is real symmetric (C144). Time-reversal-odd sources that could make K complex are the P-odd directions YX and ZY, which are strongly relevant at the A Ising point (dimension 1/8) and leave no perturbative B model (C157).
- Result: no c > 1. The c = 3/2 phase of the interacting Majorana chain (RZFA) needs a repulsive staggered channel, which second order cannot give (C157). This is the XZ analogue of the Fibonacci binding: in both cases the second-order channel that would separate the two species is attractive.

So the structural question at both extensive points is the sign of a second-order induced interaction. At XZ it is fixed by positivity of the kernel (no free parameter). At the Fibonacci point the parity class fixes it (3.2), and the transverse DM term is the one first-order-invisible direction that changes it.

### 3.4 Check against the full chain [numerical, sw_dm.out]

Exact des Cloizeaux extraction on the open L = 11 chain, built on C179's sw_exact.py.
- Δ = 0, λ = 1, w = 0: hop 4λ²/V, on-site −8λ²/V and V₂ = +8λ²/V agree to 0.5 percent at V = 40 and 0.01 percent at V = 320 (deviations ∝ 1/V²). V₃ and the NNN hop are O(λ⁴/V³).
- Δ = 1, λ = 0.25: the pairing terms of C173 add (hop2 ratio 0.986 → 0.998, V₂ = 8λ²/V − Δ²/V to 1.6 → 0.3 percent, V₃ = −Δ²/(3V) to 3 → 0.5 percent, from V = 40 to 320). The uniform single-particle creation ⟨c|H|0⟩ vanishes (7e-8 at V = 40, 2e-11 at V = 320, falling as V⁻³ or faster), and the conditioned creation ⟨c, c+2|H|c+2⟩ = 2Δλ/V to 1e-4 at V = 320. Creation conditioned on c + 3 or c + 4 is O(V⁻²) or smaller.

So the second-order model of the transverse-DM cusp, in units |t₂| = Δ²/(2V), is

  H = (3.2 with V₂ = −2 + v) + κ Σ_s P[σ⁻_s (n_{s+2} − n_{s−2}) + h.c.]P,  v = 16λ²/Δ², κ = 4λ/Δ,     (3.5)

after the real hop 8λ²/Δ² is cancelled by J_x + J_y = −4λ²/V and the on-site shift by h_z → h_z − 4λ²/V. Note v = κ² [derived]. The repulsion cannot be switched on without the ΔN = 1 term at the same order.

## 5. The transverse-DM cusp

### 5.1 Two-body and three-body problems of the U(1) part [numerical, with the proved boundary]

The U(1) part of (3.5) is the symmetric model of the collaborator's theorem with g = 2/3, so the opposite-valley pair is bound exactly when v < u²/4 [proved, contrib]. Above the boundary (pseudo_v.py, C175's scat.py unmodified):
- the dilute pseudospin parameter Δ_σ = (a_s − A_d)/|A_o| lies in (−1, 1) for every v > u²/4 tested, at u = 0.5, 1, 1.5, 2 (v up to 8). It approaches +1 just above the boundary, where A diverges, and settles at 0.27 to 0.6 for v of order one. By the C175/C180 criterion the dilute gas of the U(1) part is two-component (easy-plane pseudospin, c = 2).
- below the boundary Δ_σ < −1 (bound pair, then vector chiral among scattering states), reproducing C175 at v = 0.
- trimers (trimer_v.py, C175's nbody_box.py unmodified): at u = 1 and v = 0.5, 1, 2, E₃ − 3ε_min is positive and falls with the box (0.017 → 0.008 from R = 90 to 140), so no three-body bound state. At v = 0, E₃ − 3ε_min = −0.016, the known trimer.
So the U(1) part of the transverse-DM cusp has the repulsive inter-valley channel the spec asks for, with no dimers or trimers, for |λ|/Δ > u/8.

### 5.2 The ΔN = 1 term [derived, numerical]

κ = √v > u/2 on the unbound side, so the ΔN = 1 term is of the order of the hops wherever the pair is unbound. Its properties:
- it vanishes on an isolated particle and on long-wavelength density (form factor 2i sin 2q), so it only acts through short-range two-particle configurations;
- it breaks U(1) completely (no parity is left, since the transverse DM is odd under Π Z);
- it is odd under inversion and becomes staggered under the valley exchange S∘I (c_j → (−1)^j c_{−j}), so it breaks the valley symmetry and shifts the valley-balance point away from φ = π/2. The antiunitary K∘I∘U(π) (conjugation, inversion, c → −c) is still a symmetry and fixes each valley;
- all its matrix elements carry momentum, so no single-particle vertex is resonant at incommensurate valleys, but at κ of order one it dresses the particles strongly: a particle at the valley bottom mixes with two-particle states at the same total momentum that are only ε(0) − ε_min = 0.25|t₂| higher at u = 1.

### 5.3 Ring ED of the full model (3.5) [numerical, ring_kappa.py]

Blockaded rings, no number conservation, u = 1, v = 1, κ = 1, μ = −2.1. n(k) is the hard-core-boson momentum distribution, W_a and W_b its weights on cos k > 0 and cos k < 0.
- κ = 0 at the same μ: N = 2 at L = 20 and 24 (n = 0.10, 0.083), W_a/W_b = 1.000 at φ = π/2 and fully polarized by |φ − π/2| = 0.05 (L = 20). The U(1) part is dilute here.
- κ = 1: the density rises to 0.18 to 0.25, the κ term pulls particles in. At φ = π/2, W_a/W_b = 0.59 (L = 20) and 0.48 (L = 24). The balance point moves to φ* ≈ 1.645 (L = 20) and 1.635 to 1.645 (L = 24, W_a/W_b = 0.85 at 1.631 and 3.8 at 1.661). Away from φ* the minority valley keeps a sizeable weight (W_b/W_a ≈ 0.2 to 0.26 at φ = 1.66 to 1.72), unlike the κ = 0 case, which is the dressing of 5.2.
- The valley offset induced by κ = 1 is therefore of order δφ ≈ 0.07, i.e. a band-offset of order 2u·δφ ≈ 0.14 |t₂|, much larger than the dilute pseudospin scales. It is tunable by one more parameter (the phase of the first-order hop), so it costs one extra codimension, not a closure.
