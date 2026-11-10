This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C191 derivations: the low-energy theory of the two-valley gas from the exact two-body data

Standard: theoretical physics (RG, perturbation theory in the phase shifts, calibrated numerics). Tags [proved], [derived], [numerical], [heuristic], [assumed], [literature] open the paragraph they apply to. Every number names the script and output in this folder.

## 0. Setting

Model (3.2) of C173 in units |t₂| = 1: particles on a chain, no two adjacent, hop amplitude a_d for a move by +d with a_{+1} = w, a_{−1} = w̄, a_{±2} = −1 (over an empty site), V(2) = V₂ = −2, V(3) = V₃ = −2/3. We take |w| = 1, φ = π/2, so w = i. On open chains particles never pass, so hard-core bosons b_j and Jordan–Wigner fermions c_j = b_j e^{iπN(<j)} give the same spectrum. Total densities n = 1/8, 5/32, 3/16, with ρ = n/2 per valley in the valley-balanced state.

The plan is the standard one for a 1D multicomponent gas at intermediate coupling. The two-body S-matrix at the Fermi-surface kinematics fixes the forward (Landau) parameters and the backscattering amplitude at the scale E_F. The Gaussian theory built from them gives the Luttinger parameters, and the KT flow of g₁ is run from there. What makes this case special is (i) the choice of frame (section 3), (ii) the band's cubic term, which splits the two kinds of inter-valley R–L pairs, and (iii) the strength of the like-valley attraction, which is treated beyond first order in the Bethe ansatz (section 5).

## 1. Kinematics

[proved, algebra] ε(k) = 2 sin k − 2 cos 2k. The minima are at sin k = −1/4, k_a = −π + s and k_b = −s (in the code's labels, s = arcsin(1/4) = 0.2527), with ε_min = −2.25, ε″ = 7.5 and k_a + k_b ≡ π. The symmetry U of C179 (inversion times a staggered rotation) acts on momenta as k → π − k. It exchanges the valleys and maps a right mover of one valley to a left mover of the other. The third derivative ε‴ is valley-odd, so inside one valley v_R ≠ |v_L|.

[derived] Fermi points. By Luttinger's theorem each valley's Fermi points k_v + p_L < k_v < k_v + p_R satisfy p_R − p_L = 2πρ = πn and ε(k_v + p_R) = ε(k_v + p_L). By U the other valley has the mirror pair. Solved by bisection in scat191.py (fermi_points), output fs_data.py:

| n | p_L | p_R | v_R | \|v_L\| | ε_F − ε_min |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1/8 | −0.1913 | 0.2014 | 1.356 | 1.506 | 0.142 |
| 5/32 | −0.2375 | 0.2534 | 1.643 | 1.878 | 0.221 |
| 3/16 | −0.2829 | 0.3062 | 1.901 | 2.241 | 0.314 |

So the chiral velocities differ by 10 to 16 percent. Under U, aR ↔ bL (both speed v_R) and aL ↔ bR (both speed |v_L|).

## 2. The two-body S-matrix at the Fermi-surface kinematics

### 2.1 Set-up

[derived] Two particles x₁ < x₂, r = x₂ − x₁ ≥ 2, ψ = e^{iKx₁} f(r):

    (h f)(r) = Σ_d a_d [ e^{−iKd} f(r + d) + f(r − d) ] + V(r) f(r),   f(r) = 0 for r < 2.      (2.1)

A mode f = e^{ipr} has right-particle momentum p and left-particle momentum K − p. Its relative velocity is v_rel = ε′(p) − ε′(K − p), positive for an outgoing (separating) pair. At energy E the free modes are the roots of a quartic in z = e^{ip}. Open channels have |z| = 1, decaying ones |z| < 1. The true equation differs from the free one only at r = 2, 3, which gives two boundary rows. The flux-normalized S-matrix is S_{out,in} = A_out √|v_out| / (A_in √|v_in|). scat191.py implements this independently of C175's code.

Processes at the Fermi surface (labels: ordered channel ab = left particle in valley a):
- Like pair, (aR, aL): K = 2k_a + p_R + p_L, E = 2ε_F. One open channel. S = −e^{2iδ_s}, with the fermion reference S = −1 (f ∝ sin), so δ_s > 0 is attraction.
- Unlike pair, K = k_a + k_b, E = 2ε_F. Two incoming states: in_ab = (aR on the left, bL on the right) and in_ba = (bR on the left, aL on the right). Reflection r = S(out_ab ← in_ab) is the backscattering aR + bL → aL + bR, momentum transfer 2k_F inside each valley. Transmission t = S(out_ba ← in_ab) leaves both momenta unchanged, so it is forward scattering. r′, t′ are the same for in_ba.

### 2.2 Checks of the code

[numerical, check1.py] (i) Unitarity to 1e-13. (ii) The eigenphases of −S at E = 2ε(k_b + q) reproduce C175's energy-dependent eigen-lengths a_±(q) (smatrix.out) to all printed digits, e.g. 20.491 and 0.289 at q = 0.01, and 10.373 and 0.314 at q = 0.13. (iii) Hard rods (nearest-neighbour hop only, no V): δ = −q to 1e-15 at ρ = 0.05, 0.1, 0.2, the exact excluded-length phase with a = 1.

### 2.3 Results

[numerical, fs_data.py, lowE.out]

| n | δ_s | \|r\| = \|r′\| | arg t | arg t′ | \|t\| |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1/16 | +0.291 | 0.390 | +0.203 | +0.480 | 0.921 |
| 1/8 | +0.496 | 0.103 | −0.191 | +0.137 | 0.995 |
| 5/32 | +0.564 | 0.022 | −0.335 | +0.031 | 1.000 |
| 3/16 | +0.616 | 0.039 | −0.470 | −0.058 | 0.999 |
| 1/4 | +0.694 | 0.126 | −0.736 | −0.211 | 0.992 |

|r| passes through zero near n ≈ 0.167 (0.0011 at n = 0.1672, scan.out).

[derived] Why the reflection is small. In C175's eigen-channel language S_± = −e^{−2iq a_±(q)}, with r = (S₊ + S₋)/2 up to a phase, so |r| = |cos q(a₊ − a₋)|. The even channel carries the shallow bound state (a₊ = 20.7 at threshold). Its phase π/2 − q a₊(q) unwinds through zero near q ≈ 0.19, because a₊(q) falls quickly with q (8.05 at q = 0.19), the effect of the hard core and the finite range. The odd channel stays near its hard-core value (a₋ ≈ 0.3 to 0.4). At q = k_F ≈ πn/2 the two eigenphases therefore coincide modulo π to within 0.02 to 0.1 rad, and the backscattering is nearly absent. A zero-range model with the same a₊ would give δ_e = arctan(1/(q a₊)) ≈ 0.2 at k_F, attractive, which is the Yang–Gaudin template of C175. The difference is the whole story for g₁.

## 3. The frame

### 3.1 The free references

[derived] For two particles with valley labels, the free (non-interacting) S-matrix in the ordered basis is pure transmission. Its sign is fixed by statistics, because a free pair state is e^{i(p_a x₁ + p_b x₂)} ± e^{i(p_b x₁ + p_a x₂)}. Free fermions give t = −1, free bosons t = +1, and impenetrable particles give r = −1, t = 0.

[numerical] The computed transmission is t ≈ e^{iτ} with |τ| ≤ 0.47 at the candidate densities (table 2.3). So relative to free bosons the inter-valley pair is weakly interacting, while relative to free fermions it carries a phase near π. The like pair, by contrast, is close to free fermions (δ_s ≈ 0.5 to 0.6, against a phase near π/2 relative to bosons). The weakly coupled description therefore uses valley fields that are fermionic within a valley and mutually bosonic between valleys. This is the frame in which two decoupled XX legs are free.

[heuristic] Origin. The blockade puts the effective contact at odd separation (r = 1 is the first forbidden distance), and the two valleys differ by momentum ≈ π. A hard wall at r = 1 multiplies the transmitted wave by e^{i(k_b − k_a)} ≈ −1 relative to a wall at r = 0, which converts the fermionic sign into the bosonic one. At w → 0 this is the statement that particles on different legs of the zigzag ladder are mutually bosonic.

### 3.2 Verification on a two-particle ring

[numerical, ring2_frame.py] Two hard-core bosons on a ring of L = 600 sites at K = π, exact diagonalization of the relative problem in the bosonic sector f(L − r) = e^{−iKr} f(r). Since |r| = 0.04 at the n = 3/16 kinematics, the levels near E = 2ε_F should follow the reflectionless Bethe–Yang condition (k_a + p) L = 2πJ + Φ with Φ = −τ in the boson frame, or Φ = −τ + π in the fermion frame. Mean deviation of the predicted levels from the exact ones: 2.7e-6 for the boson frame with Φ = −τ, against 5.8e-3, 8.6e-3 and 9.0e-3 for the three alternatives (opposite sign, fermion frame with either sign). The level spacing is 2e-2. So the frame and the sign convention of the phase used in section 4 are fixed by an exact computation.

### 3.3 Operator content

[derived] In this frame the hard-core boson is b_j ≈ Σ_ν e^{ik_ν j} Ψ_ν(j), with Ψ_ν a hard-core boson of valley ν (a Jordan–Wigner fermion with a string over its own valley only). Bosonize with ρ_ν = ρ − ∂φ_ν/π and Ψ_ν ∝ e^{iθ_ν} Σ_m e^{2im(φ_ν − πρx)}. The fermion is c_j = b_j e^{iπN(<j)}, with e^{iπN(<j)} ∝ Σ_{p odd} e^{ip(πnj − φ_a − φ_b)}. The leading operators are:
- the inter-valley pair b_j b_{j+2} at momentum k_a + k_b = π, ∝ e^{i(θ_a + θ_b)} = e^{i√2θ_c}, with dimension 1/(2K_c) and correlation exponent 1/K_c (the separation-2 and -3 pair functions carry no string, since the site in between is empty);
- the boson one-body function, e^{iθ_a}, with exponent (1/4)(1/K_c + 1/K_s) in the symmetric decoupled case;
- the fermion one-body function, e^{iθ_a ± i(φ_a − φ_b)} or e^{iθ_a ± i(φ_a + φ_b)}, with exponent (1/4)(1/K_c + 1/K_s) + min(K_s, K_c) in the symmetric decoupled case;
- the valley-imbalance density ρ_a − ρ_b = −(√2/π)∂φ_s, with ⟨(N_a − N_b)²⟩_ℓ = (2K_s/π²) ln ℓ in a block of length ℓ.

Here φ_{c,s} = (φ_a ± φ_b)/√2 and θ_{c,s} likewise. Free two-leg XX gives exponents 1, 1/2 and 3/2 for pair, boson and fermion, which gauss191.py reproduces. With the cubic asymmetry the charge and spin sectors mix slightly at quadratic order, and the exponents are computed from the full Gaussian theory (section 4.4).

## 4. The low-energy theory

### 4.1 Resonant operators on the symmetric line

[derived, following C179 P5.1 with the frame of section 3] (3.2) conserves N. A valley-conserving local operator is a function of ∂φ_ν, ∂θ_ν and the density harmonics e^{2i(m_aφ_a + m_bφ_b)}, which carry lattice momentum 2π(m_aρ_a + m_bρ_b) = πn(m_a + m_b) in the balanced state. It is resonant iff m_a + m_b = 0 (commensurability would need m_a + m_b to be a multiple of 16, 64 and 32 at n = 1/8, 5/32, 3/16, i.e. cos(2√2 m φ_c) with m ≥ 8 and dimension 2m²K_c > 100). This leaves cos 2m(φ_a − φ_b) = cos √8 m φ_s. The m = 1 term is the inter-valley backscattering g₁, with dimension 2K_s. The m ≥ 2 terms have dimension 8K_s and higher and are irrelevant. Valley-changing operators carry k_b − k_a ≈ π + 2s (single conversion) or 4s (pair hopping aa → bb), which density harmonics cannot cancel at generic n. The closest case, pair hopping at n = 5/32, misses by 0.03 rad and has dimension about 2/K_s + 2K_c ≈ 4, so it is irrelevant anyway. So the only candidate for a pseudospin gap in (3.2) is g₁, and pseudospin conservation (U(1)_s) is exact at low energy.

### 4.2 Landau parameters from the phases

[derived] Bethe–Yang form. To first order in the phases, the momentum of particle j on a ring is fixed by p_j L = 2πI_j + Σ_l Φ_jl. For a pair whose left particle j approaches the right particle l (v_j > v_l), Φ_jl = −τ_jl and Φ_lj = +τ_jl, with τ the transmission phase relative to the free value. Hard rods give τ = −2qa and Φ_jl = a(p_j − p_l), the exact excluded-volume form, and the two-particle ring of section 3.2 confirms the sign for the inter-valley pair. The energy shift of the pair is then ΔE_jl = v_jΔp_j + v_lΔp_l = −(v_j − v_l)τ_jl/L. Summing over pairs in the Fermi seas and differentiating twice in the chiral occupations N_P gives the zero-mode energy

    E = (π/L) Σ_P |v_P| N_P² + (1/L) Σ_{P<P′} f_{PP′} N_P N_{P′},   f_{PP′} = −(v_P − v_{P′}) τ_{PP′},      (4.1)

with the pair (P, P′) taken in its approaching orientation. Explicitly, in the basis (aR, aL, bR, bL):
- f_∥ = f(aR, aL) = f(bR, bL) = −(v_R + |v_L|) · 2δ_s (like pair),
- f₊₊ = f(aR, bL) = −2v_R arg t, and f₋₋ = f(aL, bR) = −2|v_L| arg t′ (inter-valley R–L, the analogue of g₂⊥),
- f₄ = f(aR, bR) = f(aL, bL) (inter-valley R–R, the analogue of g₄⊥), discussed in 4.3,
- f(aR, aR) and the like are absent (zero relative velocity, and no exchange term for identical fermions at the same point).

U maps aR ↔ bL and aL ↔ bR, so these are all the independent entries. The band asymmetry makes f₊₊ ≠ f₋₋ (arg t ≠ arg t′, an eigenvector phase χ = (arg t − arg t′)/2 ≈ −0.2 that C175 also found at threshold).

[derived] Backscattering. In 1D, |T_fi| = |S_fi| √(v_rel,i v_rel,f) for flux-normalized states (in the symmetric Born limit this is r = −iŨ(2q)/v_rel). So

    g₁ = 2√(v_R |v_L|) |r|,   y₁ = g₁/(π v̄),   v̄ = (v_R + |v_L|)/2,      (4.2)

and only |y₁| matters, since without SU(2) the sign of g₁ can be changed by φ_s → φ_s + π/√8.

### 4.3 The R–R inter-valley parameter

[derived, with an unavoidable ambiguity] For an R–R pair the relative velocity is small (|v_L| − v_R ≈ 0.15 to 0.34) and the relative momentum is near zero. At zero relative momentum the 1D threshold law makes the vacuum pair impenetrable (|r| → 1, here |r| = 0.98 at the actual R–R kinematics, rr_check.py), so the vacuum phase is not a perturbative quantity there. In weak-coupling g-ology the right object is the principal-value (K-matrix) amplitude, f₄ = −(ε″₁ + ε″₂)/a₊(K) for a contact-like even channel, which reproduces Born g exactly for the Gaudin–Yang gas (section 7.2). Here a₊ at the R–R total momentum K = π + πn is 12.0, 9.7, 8.0 at n = 1/8, 5/32, 3/16, giving f₄ = −1.15, −1.34, −1.53, while at K = π it gives −0.72. But the vacuum K-matrix at threshold is dominated by the shallow bound-state pole, whose relative momenta (q ≲ 1/a₊ ≈ 0.05 to 0.1) lie inside the Fermi seas and are Pauli-blocked in the medium. The in-medium value is therefore closer to the background amplitude at q ~ k_F, which is weakly repulsive. I carry four values, 0, the K-matrix at K = π, the K-matrix at the R–R kinematics, and +0.4, and quote the spread. It moves K_s by about ±0.1 at n = 3/16.

### 4.4 Scaling dimensions

[derived, gauss191.py] With ρ_P = s_P ∂ϕ_P/2π (s = +1 for R, −1 for L) the field Hamiltonian is (1/4π)∂ϕᵀ W_f ∂ϕ with W_f = S W S, where W_PP = |v_P| and W_PP′ = f_PP′/(2π) from (4.1). Normal modes are the left eigenvectors of S W_f, normalized to l S lᵀ = ±1 (rows of R). A vertex e^{iα·ϕ} has dimension x = (1/2)|R^{−T}α|². In the per-valley language e^{i(aθ + bφ)} has chiral charges Q_R = (a − b)/2, Q_L = (a + b)/2. Operational definitions: K_s ≡ x(cos √8φ_s)/2 and K_c ≡ 1/(2x(e^{i√2θ_c})). Tests: free theory gives x = 2, 1/2, 1/4, 3/4 for g₁, pair, boson, fermion. One valley with a pure g₂ gives K = √((1 − y/2)/(1 + y/2)) exactly. Pure drift (v_R ≠ v_L, no interaction) leaves x(g₁) = 2. (A sign error in an earlier version, W instead of S W S, was caught by the one-valley test.)

[derived, symmetric closed form] With v_R = |v_L| = v and f₊₊ = f₋₋ = f_⊥, (4.1) gives

    u_c/K_c = v + (f_∥ + f_⊥ + f₄)/2π,   u_cK_c = v − (f_∥ + f_⊥ − f₄)/2π,
    u_s/K_s = v + (f_∥ − f_⊥ − f₄)/2π,   u_sK_s = v − (f_∥ − f_⊥ + f₄)/2π.      (4.3)

At first order K_s ≈ 1 + (2/π)(δ_s − τ̄) + f₄/(2πv) with τ̄ the mean inter-valley transmission phase, so the like-valley attraction (δ_s > 0) and the inter-valley forward repulsion at k_F (τ̄ < 0) both raise K_s, and the threshold attraction (f₄ < 0) lowers it. The open-chain parity gap is L·Δ_P = (π/2)(u_c/K_c + u_s/K_s), in which the inter-valley parameters cancel.

## 5. The Luttinger parameters

### 5.1 Three treatments

- Landau-δ: (4.1) with the phases as they are.
- Landau-tan: phases replaced by their tangents (the K-matrix form). The difference is third order in the phases and measures the truncation error (largest for δ_s ≈ 0.6).
- Hybrid: the intra-valley block replaced by one valley solved in the asymptotic Bethe ansatz (5.2), the inter-valley parameters as in Landau-δ.

### 5.2 One valley in the Bethe ansatz

[derived, aba_valley.py] For N particles of one valley on a ring of L sites, p_j L = 2πI_j + Σ_l Φ(p_j, p_l), with Φ from the exact like-pair phase τ(p₁, p₂) computed at the actual (K, E) of each pair (two-momentum dependence kept, interpolated on a 121 × 121 grid), and E = Σ ε(k_a + p_j). Then u/K = (L/π)∂²E/∂N² from N ± 2 and uK = (L/4π)∂²E/∂s² from a shift s of all I_j. This treats the intra-valley phase to all orders. It is exact for hard rods and for integrable kernels.

[numerical, aba_valley.out] Checks: hard rods give K = 0.80994 and 0.63997 at ρ = 0.1, 0.2, against (1 − ρ)² = 0.81, 0.64. The cusp valley at ρ = n/2 gives K_∥ = 1.362, 1.463, 1.492 and u_∥ = 1.031, 1.153, 1.304 at n = 1/8, 5/32, 3/16 (L = 512, stable against L = 256 to 0.03). The first-order values are K = (1 + δ_s/π)² = 1.43 (dressed charge) and 1.51 (Landau) at n = 3/16, so first order is adequate for K_∥ but overestimates u_∥/K_∥ by 40 percent. That is why the hybrid reproduces L·Δ_P and the pure Landau versions do not (7.4).

[numerical, aba_phi13.out, against C179 manybody.out] At φ = 1.3, where only one valley is occupied, the same code gives K = 1.32 to 1.37 at ρ = 2/15 and 1.34 at ρ = 1/6, against C179's ring values 1.34 (L = 30) and 1.46 to 1.57 (L = 24 to 36). The ρ = 0.2 point is not stable in L (1.47 and 1.31) and is not used. So the one-valley treatment is at the right level and, if anything, underestimates the like-valley attraction.

### 5.3 Results

[numerical, predict.out] Ranges over the twelve variants (three treatments, four values of f₄), and the hybrid with f₄ from the K-matrix at K = π as the central value:

| n | K_s central (band) | K_c central (band) | pair exp. | boson exp. | fermion exp. | u_s/K_s | L·Δ_P hybrid | dp/dδ |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 1/8 | 1.34 (1.32 to 1.48) | 1.39 (1.32 to 1.50) | 0.72 | 0.37 | 1.64 | 0.86 | 2.38 | 2.95 (2.2 to 3.7) |
| 5/32 | 1.53 (1.47 to 1.67) | 1.39 (1.32 to 1.53) | 0.72 | 0.35 | 1.78 | 0.83 | 2.48 | 2.47 (1.6 to 3.2) |
| 3/16 | 1.67 (1.58 to 1.89) | 1.32 (1.27 to 1.54) | 0.76 | 0.34 | 1.77 | 0.83 | 2.75 | 2.05 (1.3 to 2.6) |

dp/dδ = 1/(π n u_s/K_s) is the slope of the valley imbalance p = (N_a − N_b)/N against the offset δ = ε_b − ε_a (from E = (πu_s/4K_s)Lm² − (δ/2)Lm, m = (N_a − N_b)/L). The registered bands in RESULT.md are these, slightly widened.

What sets the numbers [derived]. One valley alone has K_∥ ≈ 1.5 at n = 3/16, from the like-valley attraction (δ_s = 0.62, from V₂ = 2t₂ and V₃, the same couplings that make each leg a ferromagnetic Heisenberg chain at w = 0). If the valleys were decoupled, K_c = K_s = K_∥. The inter-valley R–L forward scattering at k_F is repulsive in the boson frame (τ < 0) and splits them, pushing K_s up and K_c down. The R–R parameter pushes the other way. Pseudospin is far from SU(2). SU(2) would need 2(K_s − 1) = ±|y₁|, whereas here 2(K_s − 1) ≈ 0.7 to 1.3 and |y₁| ≈ 0.01 to 0.07.

## 6. The KT flow

[literature, Giamarchi's normalization] With H_s = H₀(u_s, K_s) + (2g₁/(2πα)²)∫cos √8φ_s and y = g₁/(πu_s),

    dK_s/dl = −K_s² y²/2,   dy/dl = (2 − 2K_s) y.      (6.1)

Near K_s = 1, with y_∥ = 2(K_s − 1), this is dy_∥/dl = −y², dy/dl = −y_∥y, whose invariant is y_∥² − y². The flow ends on the fixed line (gapless) iff y_∥ ≥ |y|, i.e. 2(K_s − 1) ≥ |y₁|, and otherwise runs to strong coupling. For y > |y_∥| the runaway scale is l* = [arctan(y_∥/A) + π/2]/A, with A = √(y² − y_∥²).

[derived] Verdict at the candidate densities (predict.out, ktflow.out). The margin 2(K_s − 1) − |y₁| is +0.57 to +0.89 (n = 1/8), +0.92 to +1.33 (5/32) and +1.13 to +1.75 (3/16) over all twelve variants. All three start deep on the gapless side. Integrating the full flow (6.1) from l = 0 at a₀ = 1/k_F, K_s moves by at most 3e-3 (n = 1/8), 1e-4 (5/32) and 4e-4 (3/16), and its change between L = 128 and 256 is at most 3e-5. g₁ is irrelevant from the start, with dimension 2K_s ≈ 2.6 to 3.8, so its effect decays as L^{−2(2K_s − 2)}. If y₁ is normalized by the dressed spin velocity (up to 1.9 times smaller than v̄) instead of the bare v̄, y₁ grows to at most 0.13, 0.03, 0.05 and the margins stay at +0.5 or more. So there is no pseudospin gap and no logarithmic drift of K_s at n = 1/8, 5/32, 3/16.

[derived] The separatrix in density (scan.out). The margin changes sign between n = 0.066 and 0.074 over the variants (and up to about 0.085 with the dressed normalization of y₁). Below it, y₁ grows quickly (0.25 at n = 1/16, 0.31 at 0.05, because |r| → 1 at threshold) while K_s − 1 shrinks (δ_s → 0 and the forward phase turns attractive), and the flow runs to strong coupling. ktflow.out gives ξ_s ≈ a₀e^{l*} ≈ 750 sites at n = 0.05, 2e4 to 6e4 at n = 1/16, and 2e7 to 1e15 at n = 0.0703, the KT essential singularity. This is the paired (Luther–Emery) liquid of C175 at low density, which crosses over to the c = 2 liquid near n ≈ 0.07 to 0.085. At much lower density (E_F ≲ E_b, n ≲ 0.03) the g-ology inputs fail and the dilute pair picture of C175 takes over.

## 7. Checks

### 7.1 Free two-leg XX

[derived, gauss191.py] Two decoupled XX legs are, in the frame of section 3, free valley fermions with mutual bosonic statistics: all phases vanish, W is diagonal, K_c = K_s = 1, y₁ = 0 (N_A and N_B are separately conserved). The flow sits on the fixed line at K_s = 1, c = 2. The pair, boson and fermion exponents come out as 1, 1/2 and 3/2. The fermion value is the known string result for two independent free-fermion legs with a Jordan–Wigner string through both (1 from its own leg, 1/2 from ⟨e^{iπN}⟩ on the other).

### 7.2 Attractive Gaudin–Yang gas

[numerical, gy_check.py] Here the particles are spin-1/2 fermions, so the free reference is t = −1 (fermion frame). The inputs are the contact phases δ_e(q) = arctan(1/(qa)), a = −2/(mg) > 0, δ_o = 0, δ_s = 0, and f₄ = g. The construction gives, for γ = mg/n = −0.05 to −0.4, 2(K_s − 1) = −|y₁| to first order (−0.01005 against y₁ = 0.01013 at γ = −0.05), i.e. the SU(2) separatrix on the attractive side, which flows to the gap for any attraction. It also gives 2(K_c − 1) = +|y₁|, and l* = 1/|y₁| = 98.71, 49.37, 24.72, 12.44 against the exact Bethe-ansatz exponent π²/(2|γ|) = 98.70, 49.35, 24.67, 12.34 in Δ_s ∝ exp(−π²/(2|γ|)) (Krivnov–Ovchinnikov, Fuchs–Recati–Zwerger, recalled from the literature). So the construction reproduces the known gap for any attraction.

### 7.3 Hard rods

[numerical] The S-matrix code gives the exact hard-rod phase −q (2.2), and the Bethe-ansatz valley code gives K = (1 − ρ)² to 1e-4 (5.2). The hard-rod chain is c = 1 with one valley, so no two-valley quantities arise.

### 7.4 Consistency with numbers known before this run (not blind)

[numerical, REPORT §9.158, C176 at L = 256] Pair exponent 0.78 at n = 3/16, against 0.73 to 0.79 (hybrid). Fermion envelope exponent near 2.0, against 1.72 to 1.82 (hybrid), with the fermion-frame alternative not able to reach 2 at these K's. L·Δ_P = 2.37 (n = 3/16), 2.24 (5/32), 2.18 (1/8), against 2.75, 2.48, 2.38 (hybrid) and 3.0 to 3.95 (Landau). The hybrid is 10 to 16 percent high, which again suggests that the true intra-valley compressibility is a little higher than the Bethe-ansatz value, i.e. K_∥ is if anything larger.

### 7.5 Tension with C178's small rings

[numerical, REPORT §9.151] C178's ratio R = D_Néel/D_ps, which in this convention is K_s², falls from 1.68 to 1.57 at n = 1/6 and from 1.67 to 1.35 at n = 1/5 over N = 4 to 8. That reads K_s ≈ 1.16 to 1.30, below the prediction here and falling with N. With three or four particles per valley on rings of at most 40 sites the corrections to the Luttinger form are of relative order 1/N_valley, so these sizes cannot test the prediction. The fermion exponent at L = 256 (7.4) favours the larger K_s. This tension is not resolved here.

### 7.6 The open boundary

[numerical, wall.py] A single particle reflecting from an open end of (3.2) at the Fermi energy (half-line, ψ(0) = ψ(−1) = 0 for the free extension, flux-normalized) is mostly reflected inside its valley: |R_intra| = 0.979 and 0.954, |R_inter| = 0.202 and 0.301 at n = 1/8 and 3/16 (unitary to 1e-15). At w → 0 this is the statement that each leg reflects on its own wall and the bonding and antibonding combinations mix only through the half-site offset of the two walls. So at short scales the wall pins φ_a and φ_b, i.e. φ_c and φ_s.

[derived, heuristic for the length] The inter-valley part is a boundary operator e^{±i(θ_a − θ_b)}(0) = e^{±i√2θ_s}(0). At a φ_s-Dirichlet boundary it has dimension 1/K_s (twice the bulk 1/(2K_s)), relevant for K_s > 1 and exactly marginal for free fermions, where any unitary wall matrix is allowed. With K_s ≈ 1.3 to 1.7 the far-infrared boundary condition can therefore flow to θ_s pinned (valley-swapping wall). A crude estimate ℓ_b ≈ k_F^{−1}|R_inter|^{−1/(1 − 1/K_s)} gives about 70 sites at n = 3/16 and 3000 at n = 1/8, with the amplitude normalization uncertain by a factor of order one in the exponent base.

[numerical, stiff.py] The open-chain parity gap distinguishes the two. φ_s pinned: L·Δ_P = (π/2)(u_c/K_c + u_s/K_s) = 2.38, 2.48, 2.75 at n = 1/8, 5/32, 3/16. θ_s pinned: (π/2)(u_c/K_c + u_sK_s) = 3.42, 4.13, 4.94 (with u_sK_s = 1.53, 1.88, 2.22). C176 measured 2.18, 2.24, 2.37, flat from L = 64 to 256, so at these sizes the wall behaves as φ_s pinned. One consequence for C189: a susceptibility estimator that combines L·Δ_P with χ_s on the assumption that θ_s is pinned returns (u_s/K_s)(K_s/u_s) = 1 if φ_s is in fact pinned, independently of K_s.

### 7.7 Frame of the operator relations used in C189's registration

[derived] C189 registered η_D = K_s + 1/K_c for the pair and η_G = (K_c + 1/K_c + K_s + 1/K_s)/4 for the fermion, the relations for mutually fermionic valley fields (correct for its free-fermion ff2v control). In the frame of section 3 the relations are η_D = 1/K_c and η_G ≈ (1/K_c + 1/K_s)/4 + min(K_c, K_s), as for two XX legs. With the predicted K's, the C189 relations would return K_s^D ≈ 0 and K_s + 1/K_s ≈ 4 at n = 3/16, both artifacts of the frame. C189's two-leg XX control (pair exponent 1, fermion exponent 3/2) discriminates the two sets of relations.

## 8. Other resonant operators, and the qubit chain

[derived] As shown in 4.1, nothing besides g₁ is resonant in (3.2) at these densities. In the qubit chain C179's two-pair vertex O₄ = (Ψ_aΨ_b)² ∝ e^{2i√2θ_c} is resonant on the symmetric line, with dimension 2/K_c ≈ 1.45 to 1.55 in this frame (C179 had 1.1 to 1.25 from K_ρ ≈ 3.2 to 3.6), so it is relevant. With C179's bare amplitude y₄ ≈ c₄(16/9)(Δ/V)² ≈ 1.1e-3 c₄ at V/Δ = 40, the locking scale is Δ₄ ~ y₄^{1/(2 − x₄)} ≈ 1e-6 |t₂| [heuristic, c₄ ~ 1], a locking length of order 1e6 sites. At V/Δ = 8 it would be of order 1e3 sites. This does not enter the effective-model question.

## 9. Status

| statement | tag |
| :-- | :-- |
| kinematics, Fermi points, U action | proved / derived |
| two-body r, t, t′, δ_s at k_F | numerical, exact, checked against C175 and hard rods |
| frame (mutually bosonic valley fields) | numerical, exact two-body ring, 2.7e-6 |
| only g₁ resonant in (3.2) | derived |
| Landau parameters from phases, g₁ from \|r\| | derived (first order in the inter-valley phases) |
| one-valley K_∥ in the Bethe ansatz | derived (asymptotic Bethe ansatz, exact for hard rods) |
| K_c, K_s values | derived, physics standard, band ±0.15 to 0.2 from the R–R parameter and truncation |
| gapless verdict at n = 1/8, 5/32, 3/16 | derived, robust (needs only K_s > 1.03) |
| n_KT ≈ 0.07 to 0.085, ξ_s below it | derived, semi-quantitative |
| O₄ locking length in the qubit chain | heuristic |
