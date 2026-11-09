This work was done by a Claude Opus 5.5 agent (C170, 2026-11-08), not by Codex. Astra was not used (out of credits).

# C170 derivations: re-derivation of C165 PROOF §1.2, §2 and §4

Conventions are C165's. H = Σ_j [σ_j^T J σ_{j+1} + B·σ_j] with Pauli matrices. Parent: in-plane block −R(a), J_zz = Dz, B = (0, 0, hz), pitch q = −a. Frames: σ_j = R_z(−aj) τ_j, then μ_j = τ_j rotated by π about z on odd sites (quantum frame), and classically ν_j = (−1)^j μ_j. All scripts named here are in this folder and were written independently of C165's code (only C165's JSON outputs and its family definition were read).

## R1. The exact frames (C165 §1.2)

[proved, re-derived] σ_j^T R(a) σ_{j+1} with σ_j = R(−aj)τ_j gives τ_j^T R(aj) R(a) R(−a(j+1)) τ_{j+1} = τ_j·τ_{j+1} (in-plane), so the parent becomes −(τ^xτ^x + τ^yτ^y) + Dz τ^zτ^z + hz τ^z, and the π rotation on odd sites gives the antiferromagnetic XXZ chain μ^xμ^x + μ^yμ^y + Dz μ^zμ^z + hz μ^z. Both maps are products of on-site unitaries.

[numerical, frame.py] Open chain L = 8, a = 0.77, Dz = 1.3, hz = 0.4: the parent and the XXZ chain have equal spectra to 5.5e-14.

[proved] Translation in the μ frame. Writing σ^±_j = e^{±iqj}(−1)^j μ^±_j and σ^z_j = μ^z_j, a σ-frame term of U(1)_z charge Q acquires the modulation e^{iQ(q+π)j} in the μ frame. Hence the full bond (parent plus any harmonics) commutes with

    T_σ = T_μ exp(−i(a − π) S^z_tot)      (R1.1)

and the parent alone commutes with T_μ and S^z separately. [numerical, frame.py] On a ring L = 8, a = 2π·3/8, with random harmonics: [H_μ, T_μ e^{−i(a−π)S^z}] = 4.4e-15, while [H_μ, T_μ] = 0.44, [H_μ, S^z] = 0.57 and the other sign of the phase gives 0.81. For the parent all of these vanish to 4e-15.

C165's ν-frame bookkeeping (charge-2 at 2q, charge-1 bond terms at q + π, transverse fields at q, hz at π) is correct, but the ν frame is a classical device (flipping all three components of a spin is antiunitary). In the quantum μ frame charge-1 bond terms and transverse fields carry the same modulation q + π and hz is uniform. This is why products such as hx·Γ or hx·D_x generate a uniform μ^z field at second order (R4).

## R2. Proposition 1, rigid cone linearity (C165 §2.1)

[proved, re-derived] For n_j = sinθ(cos qj u + sin qj v) + (−1)^j cosθ e, the bond products have momenta 0, ±2q (in-plane × in-plane), ±q + π (mixed), and 0 (staggered × staggered, which carries (−1)^j(−1)^{j+1} = −1). Fields carry ±q and π. For q ∉ {0, π} only momentum 0 survives the average, which gives sin²θ e_s + cos²θ e_N with e_s and e_N as in C165 (2.1). I checked the momentum-0 coefficients of the in-plane products: ⟨cos qj cos q(j+1)⟩ = ⟨sin qj sin q(j+1)⟩ = cos q/2 and ⟨cos qj sin q(j+1)⟩ = −⟨sin qj cos q(j+1)⟩ = sin q/2. The statement holds for every q ∉ {0, π}, including q = π/2.

## R3. Proposition 2, the second-order quartic (C165 §2.2), corrected

### R3.1 The relaxation formula

[derived, re-derived] At Dz = 1, hz = 0 the ν-frame parent is the ferromagnet −Σ ν_j·ν_{j+1}, degenerate on the sphere of uniform states L. With ν_j = L√(1 − |δ_j|²) + δ_j, δ_j ⊥ L, the parent costs (1/2)Σ|δ_{j+1} − δ_j|², i.e. (1/2)(1 − cos k)|d|² per site for δ_j = Re(d e^{−ikj}). A harmonic exerts F_j = Re(B_k e^{ikj}) with k ∉ {0, π}, whose energy is Re(B_k·d)/2 per site. Minimizing over d ⊥ L gives

    δE = −|P_L B_k|² / (8(1 − cos k)),   |P_L B_k|² = |B_k|² − |L·B_k|²,      (R3.1)

as in C165 (2.2). For a single charge |B_k|² is a U(1)-invariant quadratic form of L, so it is linear in u = (L^z)². The quartic comes only from |L·B_k|².

### R3.2 A general identity for L·B_k

[proved] For bilinear bond terms with ν-frame matrices M_j, F_j = M_{j−1}^T ν_{j−1} + M_j ν_{j+1}, so at the uniform state

    L·F_j = e_{j−1}(L) + e_j(L),      (R3.2)

where e_j(L) is the harmonic's energy on bond (j, j+1) evaluated on the rigid σ-frame cone n_j = Rot_j L. With e_j = Re(ε_k e^{ikj}) this gives L·B_k = ε_k(1 + e^{−ik}) and a quartic contribution |ε_k|² cot²(k/2)/4. For a field, L·F_j is the field energy itself.

The point C165 missed: the ν-frame matrix of a σ-frame antisymmetric (DM) term is Rot_j^T A Rot_{j+1}, which is not antisymmetric, because Rot_j ≠ Rot_{j+1}. So "L^T A L = 0" does not apply. Physically, the DM energy (n_j × n_{j+1})·d on the staggered cone has a component sinθ cosθ at momentum q + π, so DM drives a relaxation that depends on the cone angle.

### R3.3 The coefficients

[derived] Evaluating e_j on the cone (s = sinθ, c = cosθ, n_j = (s cos qj, s sin qj, (−1)^j c)):

- charge 2 (J_xx − J_yy = 2g, J_xy + J_yx = 2g_xy): ε ∝ (L^+)², |ε|² ∝ (1 − u)², quartic +c_2 u² with

      c_2 = (g² + g_xy²) cot²q / 4,      (R3.3)

  as C165 has for g alone (g and g_xy add in quadrature).
- charge-1 bond terms (J_xz, J_zx, J_yz, J_zy): e_j = sc(−1)^j Re[e^{iqj} Z] with Z = −J_xz + J_zx e^{iq} + iJ_yz − iJ_zy e^{iq}. So |ε|² = u(1 − u)|Z|², k = q + π, and the u² coefficient is −c_1 with c_1 = |Z|² tan²(q/2)/4. In terms of Γ_xz = (J_xz + J_zx)/2, D_y = (J_zx − J_xz)/2, Γ_yz = (J_yz + J_zy)/2, D_x = (J_yz − J_zy)/2,

      c_1 = sin²(q/2) |tan(q/2) Γ_c − i D_c|²,   Γ_c = Γ_xz − iΓ_yz,   D_c = D_y + iD_x.      (R3.4)

  For pure Γ this is C165's G²(1 − cos q)²/(2(1 + cos q)). For pure DM it is D² sin²(q/2), concave, which C165 says is zero. Γ_xz interferes with D_x (and Γ_yz with D_y), so the charge-1 coefficients do not add term by term.
- transverse fields: ε ∝ L^+, |ε|² ∝ 1 − u, linear only. Uniform hz: ε ∝ L^z, linear only (−hz²(1 − u)/8).

So E(u) = E_s + (e_N − e_s + linear shifts) u + (c_2 − c_1) u² + O(harmonics³), with c_1 from (R3.4). The qualitative rule survives and becomes sharper: every charge-2 term is convex, every charge-1 bond term (symmetric or DM) is concave, fields are linear. Coexistence at second order requires c_2 > c_1 with c_1 including the DM vector.

Worked example (pure DM_x): δE = −(D_x²(1 − cos q)/2)(1 − u)², so DM_x alone gives a first-order switch at Dz* = 1 + D_x²(1 − cos q)/2.

### R3.4 Checks

[numerical, lr2.py] Exact linear response on a ring L = 100, q = −2π·13/100 (Hessian of the parent on the sphere, pseudo-inverse on the force). First-order energies vanish (≤ 3e-15), the Hessian has two zero modes and then 1 − cos(2π/100), and dE2 is independent of the in-plane angle of L to 1e-17. Fits of dE2(u) (amplitude 0.1):

| harmonic | u² coefficient | prediction |
|---|---|---|
| g | +2.2046e-3 | c_2 = 2.2046e-3 |
| g_xy | +2.2046e-3 | c_2 = 2.2046e-3 |
| Γ_xz (G) | −2.954e-4 | −c_1 = −2.954e-4 |
| D_x | −1.5773e-3 | −D² sin²(q/2) = −1.5773e-3 (C165: 0) |
| D_y | −1.5773e-3 | same |
| G + D_x | −5.075e-4 | (R3.4): −5.07e-4 (C165: −2.95e-4) |
| G + D_y | −1.8726e-3 | (R3.4): −1.8726e-3 |
| hx, hz | 0 (1e-16) | 0 |

[numerical, check_c1.py] Eight random harmonic vectors (all nine charged entries, amplitude 0.1) at a = 0.817 and 1.948: the fitted u² coefficient equals c_2 − c_1 from (R3.3) and (R3.4) to all printed digits. C165's formula gets the sign wrong in two of the eight (trials 0 and 1 at a = 0.817, where it predicts a coexistence window that does not exist).

[numerical, clmin.py, cl_g02*.out] Full nonlinear minimization (independent L-BFGS, 32 starts per point, ring L = 100):
- g = 0.2: continuous coexistence for 1.000 < Dz < 1.019, N² linear with slope 52 (c_2^eff = 0.0097, against 0.0088 and 0.0113 at pitches 13/100 and 12/100).
- g = 0.2, D_x = 0.12: window shrinks to 1.004 < Dz < 1.019, slope 67 (ratio of effective coefficients 0.77, predicted (c_2 − c_1)/c_2 = 0.74 to 0.83).
- g = 0.2, D_x = 0.25: no window. N² jumps from 0 to 1 between Dz = 1.0200 and 1.0225, with the spiral at pitch 13/100, where (R3.3) and (R3.4) give c_2 = 0.0088 < c_1 = 0.0099. C165's formula predicts the full window here.

## R4. The quantum argument (C165 §4.1, §4.2)

### R4.1 Selection rule, exact form

[derived] At (Δ_eff, h_eff) = (1, 0) the parent is the spin-1/2 Heisenberg chain in the μ frame, with SU(2)_1 infrared. The lattice symmetry (R1.1) acts on the infrared as T_μ (the chiral Z_2, −1 on the (1/2, 1/2) primary, +1 on currents) times exp(−iα S^z_IR), with α = a − π + O(harmonics²) (the pitch renormalizes, as C165's relaxed classical pitch shows). An infrared operator of momentum P ∈ {0, π} and charge Q can appear in the effective Hamiltonian only if P − Qα ∈ 2πZ. Relevant and marginal SU(2)_1 operators, and their resonance conditions:

| operator | dim | P | Q | allowed when |
|---|---|---|---|---|
| tr g (dimerization), N^z | 1/2 | π | 0 | never |
| N^± | 1/2 | π | ±1 | α ≡ π, i.e. a ≡ 0 (q ≡ 0) |
| J^± + J̄^± (transverse field) | 1 | 0 | ±1 | a ≡ π |
| J^z + J̄^z (uniform field) | 1 | 0 | 0 | always |
| J^z − J̄^z (twist, DM_z) | 1 | 0 | 0 | always |
| J^zJ̄^z, J^+J̄^− + h.c. | 2 | 0 | 0 | always |
| J^±J̄^± | 2 | 0 | ±2 | a ∈ {0, π} |
| charge-2 operators at P = π | ≥ 5/2 | π | ±2 | irrelevant |

The bound used is h ≥ m² for a chiral state of J^z_0 = m (boson at the self-dual radius), so a charge-Q operator has dimension at least m_L² + m_R² ≥ Q²/2, and a diagonal spin-j operator needs j_L + j_R ≥ j, hence dimension at least j_L² + j_R² ≥ j²/2 (equality at j_L = j_R = j/2). For j = 4 that is 8, as C165 states.

Two consequences. (i) The eigenvalue of T_σ on a Q = 0, P = π operator is −1 whatever α is, so N^z and dimerization cannot be generated at any order, and not only perturbatively: the argument needs only that the infrared stays in the SU(2)_1 basin with the image of T_σ continuously connected to the parent's. (ii) A relevant Q = 0, P = 0 operator is generated, the uniform field J^z + J̄^z. It splits the triplet N^a and drives the flow to the magnetized liquid. It is a tuning parameter, not an obstruction, but it means the SU(2)_1 point of a bond with harmonics sits at h_eff = hz + O(harmonics²) = 0, not at hz = 0.

[numerical, heff.py] Ring L = 12, a = 2π/3, Dz = 0.9, hz = 0: the ground-state magnetization m = ⟨Σσ^z⟩/L is 0 (to 1e-5) for the parent, for hx = 0.15 alone, G = 0.15 alone and g, hx together, and 1.61e-3 for hx = G = 0.15, 4.1e-4 at half the amplitudes (ratio 3.9, quadratic), 2.85e-3 for hx = D_x = 0.15. So h_eff ∝ hx·Γ, hx·D is generated at second order, as expected from the shared μ-frame modulation q + π (R1).

So the answer to the spec's question is: no harmonic generates a relevant operator at momentum π, and the only relevant one at momentum 0 is the uniform field, which moves the bond along the switch set rather than splitting the multiplet into an independent order. Off a ∈ {0, π} (more precisely, outside resonance tongues around them of width of order harmonic^{2/3} for the dimension-1/2 operators and of order harmonic for the dimension-1 ones) the infrared at the switch is the XXZ dichotomy.

### R4.2 The spin-flop line is not SU(2)_1

[derived] Classically the switch at hz ≠ 0 is the spin-flop field hz_sf = 2√(Δ² − 1) (Pauli units). Quantum mechanically the corresponding transition is the edge of the Néel plateau, hz_c1(Δ) = 2 × (π sinh γ/γ) Σ_n sech((2n + 1)π²/(2γ)), Δ = cosh γ (des Cloizeaux–Gaudin spin gap, J = 1, times 2 for Pauli units). [numerical, plateau.py]

| Δ | 1.1 | 1.25 | 1.3 | 1.4 | 1.6 | 2.0 |
|---|---|---|---|---|---|---|
| hz_c1 | 0.0004 | 0.022 | 0.0405 | 0.096 | 0.269 | 0.780 |
| classical hz_sf | 0.92 | 1.50 | 1.66 | 1.96 | 2.50 | 3.46 |

The plateau edge is a Pokrovsky–Talapov (z = 2) transition between gapped Néel and a magnetized Luttinger liquid, not an SU(2)_1 point. Away from Δ − 1 → 0 the "one SU(2)_1 primary" statement does not apply, and the dimension-8 count has no bearing there. The conclusion still holds by a different argument:

1. Off resonance, the effective μ-frame Hamiltonian is U(1)_z and T_μ invariant (R4.1 with momenta P ∈ {0, π} + 2n k_F, whose resonances are C154's hypersurfaces). So near the plateau edge the problem is a spin-1/2 XXZ chain in a field with weak additional U(1)-symmetric couplings of range a few sites and size O(harmonics²/J).
2. In a single boson, Néel_z order (cos 2φ locked) and in-plane quasi-order (θ gapless) exclude each other. A coexisting (biconical) phase therefore needs a second low-energy sector, which at the plateau edge means an S^z = 1 bound state of two spinons (a magnon) below the two-spinon threshold, condensing before the spinons do.
3. The parent has no such bound state (massive spin-1/2 XXZ, Bethe ansatz), and its spinons at the band bottom scatter with S(0) = −1, with no zero-energy resonance (in the Ising limit they are hard-core domain walls). A weak attraction then cannot bind them, since in 1D binding by a weak potential needs a zero-energy resonance of the unperturbed pair. So weak harmonics create no second sector at the plateau edge.

Tag: derived, with the standard inputs in 3. This replaces C165's sentence "perturbation theory in the harmonics does not create [a second sector]" by its mechanism.

### R4.3 Where the classical window sits relative to the quantum transition

[derived] At hz = 0 the classical switch and the quantum KT point coincide at Dz = 1 for the parent, so with harmonics both move by O(harmonics²). C165's "the quantum shift of the Néel onset is O(1)" is wrong at hz = 0. At hz ≠ 0 there is an O(1) separation, but it comes from the field, not the harmonics. At hz = 0.05 the classical spin-flop is at Δ = √(1 + hz²/4) = 1.0003 and the parent's quantum Néel onset (hz_c1(Δ) = 0.05) at Δ = 1.321 [numerical, plateau.py]. On C165's line Q2 the classical window (Dz ≈ 1.0) is therefore deep in the parent's magnetized liquid, as C165 says, but for this reason.

## R5. The plateau-edge numerics (C165 §4.3, last bullet)

[numerical, pfull.py, plateau2.py] C165 reports gapped Néel (N = 0.73, ξ = 4) up to hz = 2.2 at Dz = 2.0 with the harmonics of Q2. The parent's plateau ends at hz_c1 = 0.78. Energies per site:

| hz | C165 VUMPS (Néel) | parent ring ED L = 20, best sector | full bond, open-chain ED [E(14) − E(10)]/4 |
|---|---|---|---|
| 0 | | −2.4727 (Néel) | −2.4571 |
| 1.0 | −2.4760 | −2.4727 (Néel, finite-size gap) | −2.5108 |
| 2.2 | −2.4781 | −2.5707 (m = −0.1) | −2.5748 |
| 2.6 | −2.4787 | −2.6271 (m = −0.2) | −2.5762 |

The open-chain estimator reads 0.02 high at hz = 0 (against the VUMPS Néel energy), so its error is about 0.02. At hz = 2.2 and 2.6 both independent estimates lie 0.1 below the VUMPS energy, which is variational. The VUMPS states at hz ≥ 1.4 (and probably at 1.0) are metastable Néel states, not ground states. The two-site, χ ≤ 40 ansatz evidently cannot nucleate a finite spinon density from the Néel seed.

Consequence for the Typhon screen (job 17601089, family A2): its fields are hz = 0.3 to 2.0 at Dz = 1.6 (parent edge 0.27) and hz = 0.6 to 2.7 at Dz = 2.0 (parent edge 0.78). Only the lowest one or two fields of each list are near the edge, and any A2 point reading "gapped Néel" above the edge must be compared in energy with a magnetized seed before it is read.

On line Q2 (hz = 0.05) the parent's quantum onset is Δ = 1.32, so C165's points at Dz = 1.25 and 1.3, read as gapped Néel at χ ≤ 40, lie on the parent's liquid side of its PT edge. Harmonics can shift the edge by O(0.1) and the χ ≤ 40 reading is not decisive either way. Neither changes c ≤ 1.
