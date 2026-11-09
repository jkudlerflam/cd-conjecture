This run was done by Claude Opus 5.5 (claude-opus-5-5), not by Codex. Astra was not used (out of credits).

# C169 derivations: constructive search for a c > 1 nearest-neighbour qubit bond

Notation. B_2 is the space of bonds h = σ_j^T J σ_{j+1} + B·σ_j (C154). A "family" of a Pauli Hamiltonian is a Pauli string P together with all its translates T^j P. Claims are tagged proved, derived, numerical, heuristic, assumed or literature. Scripts: frustration.py, claw.py, interface.py, freecut.py, with outputs frustration.out, claw.out, interface.out, freecut.out.

## 0. What an exact map can and cannot do

[proved, elementary] Let α be any locality-preserving map under which the spectrum and the low-energy theory are invariant (a quantum cellular automaton, a product of on-site unitaries, or a Kramers–Wannier-type duality on a symmetry sector). Then c(α(H)) = c(H). So the image α^{-1}(B_2) is a relabelling of B_2 with identical spectra, and "does α^{-1}(B_2) contain c > 1" is literally the question "does B_2 contain c > 1". Exact maps cannot produce a counterexample that was not already there. What they can do is recognize one: a chain whose c > 1 is known in some other representation might be secretly in B_2. The constructive task therefore reduces to the following, for each class of maps, either find a known c > 1 model inside α^{-1}(B_2) or exhibit an invariant of the class that separates every known c > 1 model from B_2.

Known c > 1 qubit models used as targets (C156 list, extended).
- Cluster-Ising, H = −Σ X_{j−1}Z_jX_{j+1} + λ Σ Y_jY_{j+1}, c = 3/2 at λ = 1. [derived here] By Jordan–Wigner both terms are range-3 Majorana hoppings, b_{j−1}a_{j+1} and a_jb_{j+1}, i.e. γ_n γ_{n+3} for every n. At λ = 1 the chain is three decoupled uniform Majorana chains (n mod 3), c = 3 × 1/2.
- Rahmani–Zhu–Franz–Affleck spin form, Z, XX, ZZ, X_jX_{j+2}, c = 3/2.
- O'Brien–Fendley, XX, Z, XXZ, ZXX (tricritical point c = 7/10, included as the standard three-site benchmark).
- Titvinidze–Japaridze, XX + YY + J*(XZX + YZY), free, c = 2 for J* above a threshold (§5.3).
- J1–J2 XXZ in a field (c = 2 two-component phases), two XX legs in zigzag form (c = 2).
- d = 4 ladders, Ising × XX (c = 3/2) and two XX legs (c = 2).
- Integrable staggered six-vertex and Bariev-type chains (c = 2), all with range-2 or three-site terms [literature, C156].

## 1. Products of single-site unitaries, possibly site-dependent

[proved] Write H = Σ_S H_S, where H_S is the part of H that acts nontrivially (traceless) on every site of S and trivially elsewhere. This decomposition is unique. A product U = ⊗_j u_j of arbitrary single-site unitaries (any period, or none) maps H_S to an operator with exactly the support S, because u_j maps traceless operators to traceless operators. So the set of supports {S : H_S ≠ 0} is invariant, and U H U† is nearest-neighbour iff H is. Every target in §0 has a nonzero H_S with S of diameter 2 or more, so none is mapped into B_2. The obstruction is range, in the strict sense of Pauli support.

[proved, generalizing C156's screw lemma] A nearest-neighbour chain with bond data (J_j, B_j) is equivalent under on-site rotations to a one-site chain iff it has a twisted translation, i.e. there are W_j ∈ SO(3) with W_j^T J_{j+1} W_{j+1} = J_j and W_j^T B_{j+1} = B_j. (If R_j makes the chain uniform, W_j = R_{j+1} R_j^{-1} works, and conversely R_j = W_{j−1} ⋯ W_0 builds the gauge.) The bond singular values, |B_j|, and the relative orientation of J_{j−1}^T J_{j−1}, J_j J_j^T and B_j at each site are the local invariants that must be j-independent. So site-dependent rotations only identify which period-p nearest-neighbour chains are B_2 in disguise. The question "does a period-p nearest-neighbour qubit chain reach c > 1" is relevant here only as a test of whether one-site translation matters. Sections 5 to 7 show that in every controlled setting it does not.

## 2. Clifford dualities: Kramers–Wannier, Jordan–Wigner, cluster entanglers, any Clifford QCA

A Clifford map sends each Pauli string to a single Pauli string (up to sign) and preserves commutation. A translation-invariant Pauli Hamiltonian has a unique family decomposition, and its image is nearest-neighbour iff every family is sent to a family of support at most two adjacent sites.

### 2.1 Two invariants

For families F, G define the anticommutation offset set D(F, G) = {k ∈ Z : P_F anticommutes with T^k P_G}.

- [proved] If α commutes with translation (every translation-invariant Clifford QCA, and KW on the Z_2-even algebra with its half-shift relabelling), then α(P_F) = T^{s_F} Q_F with Q_F a B_2 string, and D(F, G) = D(Q_F, Q_G) + s_G − s_F. So the diameter of D(F, G) is invariant.
- [proved] For any Clifford map, covariant or not, the frustration graph (vertices the terms, edges the anticommuting pairs) is preserved up to isomorphism, so the number of infinite connected components is invariant.

### 2.2 The values on B_2 and on the targets

[proved by support, checked by frustration.py] For two B_2 families, two-site strings anticommute only when their supports overlap, so D(F, G) lies inside three consecutive integers, and every pair has diameter at most 2. Over all 4095 subsets of the 12 Pauli-aligned B_2 families, the frustration graph has at most 2 infinite components (attained by XX + YY, the two Majorana chains of the XY chain; this is also §7).

[proved, frustration.out] Target values.

| model | largest diameter of D | pair | infinite components |
|---|---|---|---|
| cluster-Ising | 3 | XZX vs YY, D = {−1, 2} | 3 |
| RZFA spin form | 3 | ZZ vs XIX, D = {−2, −1, 0, 1} | 1 |
| O'Brien–Fendley | 4 | XXZ vs XXZ | 1 |
| Titvinidze–Japaridze | 4 | XZX vs YZY, D = {−2, 2} | 1 |
| J1–J2 XXZ | 4 | XIX vs YIY | 1 |
| two XX legs, zigzag | 4 | XIX vs YIY | 4 |
| Ising × XX (two-site cell) | 2 in cells | | 3 |

So no translation-covariant Clifford map sends any one-site target into B_2 (diameter), and no Clifford map of any kind sends cluster-Ising, two XX legs or Ising × XX into B_2 (three or more components). RZFA, O'Brien–Fendley, Titvinidze–Japaridze and J1–J2 are excluded for translation-covariant maps. A non-covariant Clifford image of a one-site target would not be translation invariant in general, and the question is then the one of §3.

### 2.3 Which B_2 bonds survive KW or the cluster entangler

[proved] Kramers–Wannier with respect to ∏Z sends Z_j → μ^x μ^x, X_jX_{j+1} → μ^z, X_jY_{j+1} → μ^y μ^x, Y_jX_{j+1} → μ^x μ^y (all nearest-neighbour), but Y_jY_{j+1} → μ^x μ^z μ^x (three-site) and Z_jZ_{j+1} → μ^x_{j−1/2} μ^x_{j+3/2} (next-nearest). After a common rotation, the bonds whose KW image is again in B_2 are therefore those with J_zz = 0, J_yy = 0 and no Z_2-odd entries, which are quadratic after Jordan–Wigner. The KW-stable part of B_2 lies inside the free locus, where c ≤ 1 is proved (C143, and §5 for any period). Conversely, the KW preimages of B_2 contain XZX and next-nearest XX terms, a class whose c is the c of B_2's parity class (C152).

[proved] The cluster entangler ∏CZ_{j,j+1} sends X_j → Z_{j−1}X_jZ_{j+1} and fixes Z_j. On B_2 families the nearest-neighbour-preserving subspace is span{Z, ZZ, XZ, ZX, YZ, ZY}, on which it acts as XZ ↔ ZX and YZ ↔ ZY (with a shift), i.e. as the site reflection J → J^T. It produces no new B_2 points. Jordan–Wigner composed with an inverse Jordan–Wigner in a different axis is again Clifford and is covered by §2.1.

## 3. Blocking

[proved] If a chain on (C^{2^b})^{⊗L} is a blocked B_2 chain up to arbitrary block unitaries U(2^b)^{⊗L}, the part of H acting nontrivially on both sides of a block boundary is σ_last^T J σ_first conjugated by the block unitaries. Its left operator-Schmidt space V_L (a canonical subspace, independent of bases) is therefore unitarily a subspace of span{σ^a ⊗ 1}. Two consequences are invariant under block unitaries:
(i) dim V_L ≤ 3;
(ii) an orthonormal basis of V_L obeys the Clifford relations v_a v_b + v_b v_a = 2δ_ab, and for dim V_L = 3 the triple product v_1 v_2 v_3 is ±i times the identity.

[proved, interface.out] Every target fails at b = 2 and b = 3. Titvinidze–Japaridze, two XX legs (zigzag and ladder), J1–J2 and O'Brien–Fendley have dim V_L ≥ 4. RZFA and the Ising × XX ladder have dim V_L = 3 but break the Clifford relations (for example X_1 + X_2 has spectrum {2, 0, 0, −2}). Cluster-Ising passes the Clifford relations and fails only the triple product, which equals iX_1 rather than i. Its cut algebra is Cl_3 ≅ M_2 ⊕ M_2, not a single qubit factor. Two B_2 controls pass. So blocking never brings a known c > 1 model into B_2, and the precise invariant is "the cut is crossed by one qubit–qubit coupling".

## 4. General (non-Clifford) quantum cellular automata

[proved, elementary] For a general finite-depth circuit α there is no invariant beyond the spectrum itself, and the question is equivalent to the original one (§0).

[derived, with literature input] Integrable targets are excluded anyway. A QCA sends local conserved charges of range n to local conserved charges of range at most n + 2R. If α mapped a free or Bethe-solvable c > 1 model (Titvinidze–Japaridze, cluster-Ising, staggered six-vertex, Bariev) onto a B_2 bond h′, then h′ would carry nontrivial local charges of every range between 3 and L/2. For symmetric J (no DM part) with arbitrary fields, Yamaguchi, Chiba and Shiraishi (arXiv:2411.02162 and its companion 2411.02163) prove that only the known integrable chains have such charges, namely the XYZ, XXZ-with-longitudinal-field and free-fermion families. All of these have c ≤ 1 (C156 and de Leeuw–Pribytok–Ryan). Scope: bonds with a DM component are not covered by their theorem. For those this step is assumed.

## 5. Free fermions at any period: the obstruction is range, not translation

### 5.1 Lemma (cut rank) [proved]

Take any nearest-neighbour qubit chain of period p (one-site translation not assumed) that is quadratic after Jordan–Wigner, i.e. built from X_jX_{j+1}, Y_jY_{j+1}, X_jY_{j+1}, Y_jX_{j+1} and Z_j. With Majoranas a_j, b_j per site, the terms are a_jb_j, b_ja_{j+1}, a_jb_{j+1}, b_jb_{j+1}, a_ja_{j+1}. The single-particle matrix is block tridiagonal in 2 × 2 site blocks: on-site blocks D_j and inter-site blocks T_j (the in-plane part of J_j). The Bloch matrix h(z), z = e^{iK}, of size 2p has z only in the corner block T_p z and its partner −T_p^T z^{−1}. In the Leibniz expansion every term takes at most two entries from the rows of the last site and at most two from the rows of the first site, so z² det h(z) is a polynomial of degree at most 4. If det h ≡ 0 the chain has a flat zero band (a singular point, not conformal). Otherwise h(K) has at most 4 zero crossings on the circle, counted with multiplicity. Each linear crossing of an analytic band branch is one chiral Majorana, right and left movers are equal in number (spectral flow on the circle), and so

    c = (number of zero crossings)/4 ≤ 1.

If one bond in the cell has rank-1 in-plane block, the same expansion gives degree at most 2 and c ≤ 1/2. Checks: XX chain 4 crossings (c = 1), critical Ising 2 (c = 1/2). Numerical check (freecut.py): over 800 random chains with p ≤ 4 the Laurent span of det h(z) never exceeds 4, the rank-1 case has span 2, and 3000 random chains never show more than 4 unit-circle roots.

Consequences.
- [proved] c ≤ 1 on the whole free locus of B_2, including bonds that are free only after a site-dependent rotation with any period (the Kitaev chain form of XZ is period 2). This extends C143 §2.4.
- [proved] At range r the corner block is 2r × 2r and the same argument gives c ≤ r, attained by Titvinidze–Japaridze at r = 2 (c = 2). Cluster-Ising (c = 3/2) also needs r = 2.
- So for free fermions the bound is set by the Majorana cut rank (two Majoranas per site and nearest-neighbour coupling), and holds with or without one-site translation.

### 5.2 Free fermions in disguise [proved, claw.py]

Fendley's hidden free fermions and Chapman–Flammia's generalized Jordan–Wigner both require a claw-free frustration graph. Over all subsets of the 12 Pauli-aligned B_2 families, the claw-free non-commuting ones fall into 20 classes up to axis permutation and reflection (claw.out). Every one is either single-site or a J supported in a 2 × 2 block (after a uniform rotation) with the field perpendicular to that plane, i.e. ordinary Jordan–Wigner-free and covered by §5.1. So the Pauli-aligned part of B_2 contains no free-fermion-in-disguise point beyond the known free locus. Scope: bonds not aligned with any Pauli frame are not enumerated.

### 5.3 How far the free c = 2 phase is from nearest neighbour [derived]

For the Titvinidze–Japaridze band ε(k) = −cos k − J* cos 2k − μ, four Fermi points need ε′ to have four zeros, i.e. −sin k (1 + 4J* cos k) = 0 with |4J*| > 1. The c = 2 region requires a second-harmonic fraction above 1/4 of the first. Along any path from it to the nearest-neighbour locus, c drops to 1 at a Lifshitz point at finite distance.

## 6. Linear spin waves at any period: at most one soft species [proved algebra, derived physics]

Take a classical product state, periodic with any period p (spirals become periodic with p = 1 in the twisted frame). Linear spin-wave theory gives a quadratic boson problem coupled only between neighbouring sites, with 2 × 2 blocks in (a_j, a_j†). The dynamical matrix M(z) is again block tridiagonal with z only in the corner, so z² det M(z) has degree at most 4.
- Without anomalous terms (number-conserving case, e.g. any saturation onset) the problem is p × p and z det(M(z) − ω) has degree at most 2. A band minimum is a double root, so each band has a single minimum. Two degenerate magnon minima, the J1–J2 mechanism of the c = 2 two-component gas (Sudan et al., Kecke et al.), are impossible at every period. This extends C160's one-site statement.
- With anomalous terms, a stable state has M(K) ≥ 0, and at an instability each soft momentum is a root of even multiplicity of det M. So there are at most two soft momenta, a single K_0, a pair ±q bound into one complex boson, or {0, π}. Two independent incommensurate soft pairs would need four. This extends C154 Lemma 3 to every periodic reference state.
- At range r the degrees are 2r (conserving) and 4r, which is how J1–J2 gets ±q at saturation.

So the magnon route to c > 1 is blocked by range at every period, not by translation.

## 7. Exact decoupling into commuting sectors [proved]

Suppose H = H_A + H_B with H_A = Σ_j T^j a, H_B = Σ_j T^j b, a and b nearest-neighbour two-site operators (fields allowed), [T^i a, T^k b] = 0 for all i, k, and neither H_A nor H_B a sum of single-site terms or of mutually commuting terms.

Proof. Write the parts of a acting on its right site as the span S_a^R of {a_R, R_μ} (single-site part and right Schmidt factors), and similarly S_a^L, S_b^L, S_b^R. Expanding [a, T b] on three sites and using linear independence of the outer factors gives [S_a^R, S_b^L] = 0, and likewise [S_a^L, S_b^R] = 0. On a qubit, two nonzero spaces of traceless operators commute elementwise only if both lie in span{n·σ} for one n. Neither space can vanish, since otherwise the family is a sum of single-site fields. So a = α(m·σ)(n·σ) + β m·σ ⊗ 1 + γ 1 ⊗ n·σ and b = α′(n·σ)(m·σ) + β′ n·σ ⊗ 1 + γ′ 1 ⊗ m·σ. H_A must not be commuting, which needs m ∦ n and α ≠ 0. The same-bond condition [a, b] = 0 then forces m ⊥ n ([mσ⊗nσ, nσ⊗mσ] ∝ (m·n)|m × n| (kσ ⊗ 1 − 1 ⊗ kσ)), and its remaining six terms, kσ ⊗ nσ, mσ ⊗ kσ, kσ ⊗ mσ, kσ ⊗ 1, nσ ⊗ kσ and 1 ⊗ kσ with k = m × n, are linearly independent, so β = γ = β′ = γ′ = 0. ∎

Result. The only such bonds are α XZ + α′ ZX up to rotation: two decoupled uniform Majorana chains, each critical, c = 1/2 + 1/2 = 1 (α′ = 0 is the XZ point with a frozen second sector). A third commuting dynamical family would have to pair with a as b does, so it is proportional to b. The same argument with site-dependent axes and amplitudes covers every period p and again gives two Majorana chains, c ≤ 1. This proves C157's statement that XZ is the unique bond with an exactly decoupled second sector, in the strong form "any period, term-wise commuting". Scope: decompositions into commuting local Hamiltonians whose terms are longer than nearest-neighbour are not covered.

## 8. Deformations of c = 2 and c = 3/2 two-site-cell models

Is the obstruction range, translation by one site, or dynamics?

- Two XX legs (c = 2). In qubit form the legs are next-nearest-neighbour (zigzag, one-site cell) or range 3 (ladder ordering). Free version (Titvinidze–Japaridze): range, quantified in §5.3, and §5.1 shows period does not help. Interacting version (J1–J2 in a field): the two-minimum magnon band needs range (§6) at every period. Obstruction: range.
- Ising × XX (c = 3/2). As a free model it is three Majorana chains, i.e. three infinite frustration components, and B_2 allows two (§2.2, §7). The only place B_2 holds a second decoupled Majorana sector is XZ (§7). There the B chain's needed term is b_jb_{j+2} = Z_jX_{j+1}Z_{j+2}X_{j+3}, a four-site Pauli string [proved], so no nearest-neighbour bond supplies it at first order. The only direct B term available is Σ b_j = Σ Z_jX_{j+1}, which turns XZ into XZ + ZX (c = 1). The range-2 and range-3 B couplings arise only by virtual A processes, and C157 shows by Lehmann positivity that their second-order kernel is attractive, the wrong sign for RZFA's c = 3/2 phase. Obstruction: range at first order, then dynamics (sign) at second order. Higher orders at finite β are not controlled, but C144's DMRG to β = 0.5 saw no c above 1.09.
- Classical coexistence (C165): exists classically and fails at spin 1/2 because the staggered and spiral orders are one SU(2)_1 multiplet. Obstruction: dynamics (quantum, spin 1/2).
- Translation by one site appears in none of these obstructions. In the free, spin-wave and exact-decoupling settings the bounds hold at every period (§5 to §7), and kinematically one-site translation helps c > 1, since its incommensurate infrared image is what makes floating bosons codimension 0 (C154).

## 9. Targets added by the coordinator: irreducible c > 1 points and their dynamical obstructions

### 9.1 Emergent Z_m with m ≥ 5: a Fateev–Zamolodchikov point at a melting tongue (C174 R5)

Setup. B_2 has no on-site Z_m with m ≥ 3 beyond U(1). [proved] Under a rotation by 2π/m about an axis, the spherical components of J and B carry charges q ∈ {0, ±1, ±2}, and a component survives only if m divides q. For m ≥ 3 only q = 0 survives, so the bond is U(1) symmetric. Any Z_5 is therefore lattice translation acting on a period-5 locked phase (a Z_5 clock tongue of a floating liquid, C154 §2.3, C163 Proposition P).

(a) Chirality cannot be removed by symmetry [derived]. The chiral (Huse–Fisher) parameter is the asymmetry between a kink (the ordering label increases by one from left to right) and an antikink. Every lattice symmetry maps kinks to kinks. Translation and on-site unitaries commute with translation, so they shift labels uniformly. Time reversal commutes with translation, so it also shifts labels. Spatial reflection reverses both the spatial order and the labels (it inverts T), so a kink becomes a kink at momentum −k. In field-theory terms reflection acts as CP and time reversal acts in the σ basis as complex conjugation K, and the chiral clock Hamiltonian is CP- and K-invariant for every chiral angle. So chirality is generic in every stratum (N, T, K, parity, and the reflection-symmetric bonds) and vanishes only by tuning, as in Rydberg chains, where the Potts point is an isolated point on the Z_3 boundary.

(b) Counting [derived, dimensions heuristic for the chiral operator]. At the Z_5 parafermion the Z_5-neutral scalar relevant operators are ε_1 (4/7) and ε_2 (12/7), and the chiral operator (dimension about 11/7, spin 1) is relevant. Z_5-charged operators carry momentum 2πq/5 and are forbidden by translation. So the FZ locus has codimension 3 in any stratum: a 5-dimensional locus in N, 2-dimensional in T or K. This is not a kinematic obstruction.

(c) Gaussian obstruction [derived]. Near the tongue the theory is one compact boson. The Z_m locked phase pins V_m = O^m, where O is the order parameter, and the adjacent translation-invariant gapped phase pins a zero-momentum vertex W mutually nonlocal with O (index ℓ ≥ 1). For vertices e^{i(nθ + 2wφ)} with Δ = n²/(4K) + w²K, write Δ = a² + b² with a = n/(2√K), b = w√K. Lagrange's identity gives

    Δ(V_m) Δ(W) ≥ (a_1 b_2 − a_2 b_1)² = (mℓ)²/4 ≥ m²/4.

A direct continuous transition inside the Gaussian theory needs both to be relevant, Δ < 2 each, so m²/4 < 4, i.e. m ≤ 3 (m = 4 marginal, the Ashkin–Teller case). For m ≥ 5 there is a window where both are irrelevant, the intermediate c = 1 floating phase bounded by two KT lines on the commensurate line and by Pokrovsky–Talapov lines off it, as in José–Kadanoff–Kirkpatrick–Nelson. An FZ point needs this window to pinch shut non-perturbatively, which requires order-one vortex fugacity and order-one clock coupling together.

(d) Spin-1/2 nearest-neighbour specifics [derived, with (c)].
- The Z_5 clock vertex on the θ side is e^{5iθ}. Bond harmonics carry U(1)_m charge at most 2, so it appears only at third order (2 + 2 + 1 or higher) in the harmonic amplitude. It is relevant at weak coupling only for K > 25/8, i.e. strongly attractive liquids near h_F, where C160's semiclassics apply (ħ_eff small) and a quantum order–disorder transition of a locked spiral is excluded. There the locked spiral melts by soliton proliferation (Pokrovsky–Talapov, c = 1) or first order (C160, C161).
- On the φ side the period-m density wave e^{2imφ} has Δ = m²K, and every nearest-neighbour U(1) liquid has K > 1/4 (XXZ in a longitudinal field, Bethe ansatz [literature]). So Δ > 25/4 at m = 5, which is never relevant. To all orders in weak harmonics there are no φ-side Z_5 tongues at all.
- The disordering vortex's momentum must also lie in (2π/5)Z, an extra commensurability of the second wavevector.
So in every controlled neighbourhood (weak harmonics to all orders, C163 Proposition P; the semiclassical neighbourhood of h_F, C160) m ≥ 5 tongues melt through c = 1. An FZ point could only sit at order-one amplitude, out of perturbative and semiclassical control, which is the same region as (A6)/(A_W).

(e) Literature [literature]. In Rydberg chains, the qubit chains whose Z_p order also comes from translation, period-p crystals with p ≥ 5 always melt in two steps through a floating phase (Samajdar, Chepiga and Mila, PRR 3, 023049 (2021); Chepiga, PRB 108, 184425 (2023)). Sine-Gordon fits along the period-5 commensurate line give two BKT transitions (arXiv:2606.16128). The only direct z = 1 transition out of a p ≥ 5 crystal I found is the Z_6 transition in a two-leg Rydberg ladder (arXiv:2401.08087), a two-species geometry.

(f) Consistency with C154's closure lemma [derived]. An FZ locus would be a codimension-3 critical component on the frontier of the c = 1 floating component, with c = 8/7 ≥ 1. The closure lemma ("c at a frontier point ≥ c of the component") does not exclude it. It is excluded only if its own closure reaches the controlled set, which is an instance of (H) for C174.

Verdict on 9.1. No construction. Precise obstruction: dynamics in the form of (c). Inside every controlled region the Lagrange bound forces an intermediate c = 1 phase, and the clock coupling the FZ point needs (third order in harmonics, at K > 25/8) is out of reach there. Kinematics do not obstruct, and chirality can only be tuned away, never symmetry-forbidden.

### 9.2 Three-state Potts at momentum 2π/3 × a half-floating boson, c = 9/5 (C171 R1)

(a) Chirality [derived, from 9.1(a)]. With Z_3 acting through translation, no symmetry of any stratum maps kinks to antikinks. In T and K_I the antiunitary Θ maps momentum k → −k, so it sends σ to σ†. That is the action of K in the σ basis, under which the chiral clock chain is invariant at every chiral angle. So Θ does not forbid the chiral term. The Potts point needs ε and the chirality tuned, codimension 2 rather than 1, giving a 3-dimensional locus in T or K_I. In N a spin-1 coupling ε ∂_xφ to the boson (dimension 9/5, relevant) is also allowed. In T it is forbidden, because ∂_xφ, the magnetization along the spiral axis, is Θ-odd. So the codimension in N is at least 3.

(b) Where would the Potts sector come from? [derived] With one qubit per site, a period-3 order is one of the following.
- The host's own spiral angle locked by e^{3iθ}. This is the host boson itself, and locking it gaps the host (a Z_3 tongue), so there is no coexistence.
- A φ-side period-3 density wave (an up-up-down plateau). Its umklapp e^{6iφ} has Δ = 9K > 9/4 > 2 for every nearest-neighbour U(1) liquid (K > 1/4), so it does not exist to all orders in weak harmonics. This is the known absence of a 1/3 plateau in nearest-neighbour XXZ.
- Classical nearest-neighbour Ising-type ordering has period at most 2 (minimum mean cycle in a two-state transfer graph) [proved]. Continuous-spin classical ground states are uniform or single-period spirals (C160 M1), so period 3 classically is again the locked spiral.
So the Potts sector would have to be a second species, an independent gapped mode at ±2π/3 outside Zθ + Zπ that softens while the host liquid survives. That is (A_W)(b) in T and K_I and (A6) in N.

Verdict on 9.2. No construction. The obstruction is the one-species dynamics already isolated by C163 and C154. The kinematic count of R1 should read codimension 2 (T, K_I) and at least 3 (N), not 1.

### 9.3 A codimension-0 c = 2 with a repulsive host, K_1 = 0.9, K_2 = 0.2 (C171 R4)

(a) Kinematics [derived]. Agreed. With independent incommensurate momenta no vertex coupling the two bosons survives, and the only question is whether a second gapless boson exists.

(b) What K_2 = 0.2 needs [derived]. Every single-species nearest-neighbour U(1) liquid has K > 1/4 [literature], so K_2 = 0.2 cannot be a magnon fluid of the qubits themselves. It must be a fluid of composites with an extended hard core. For an n-quantum composite, C163 Lemma L's conversion vertex e^{i(θ_2 − nθ_1)} has dimension n²/(4K_1) + 1/(4K_2) = 0.28n² + 1.25. That is 1.53 (relevant, locked) for n = 1 and 2.36 (irrelevant) for n = 2, so R4 sits in the n ≥ 2 gap of Lemma L.

(c) Dynamics [derived near the U(1) class, assumed beyond]. In the repulsive U(1) parent the 2-string (bound pair) branch lies above the two-particle continuum, so pairs never become the lowest charge-2 excitation and their band never softens while the host liquid survives (C163 Fact 3). At the spin-wave level no second soft species appears at any period (§6), the free theory has c ≤ 1 at any period (§5), and exact decoupling gives only two Majorana chains (§7). A pair fluid in a repulsive host at order-one U(1) breaking is not excluded and is part of (A6).

Verdict on 9.3. No construction. The obstruction is that a repulsive host has no low-lying bound pairs (derived near U(1)), and beyond that it is (A6).

## 10. Summary of what is proved, derived and assumed

| construction | result | tag |
|---|---|---|
| on-site unitaries, any period | support invariant, no target reaches B_2 | proved |
| Clifford QCAs, KW, JW, cluster entangler | offset-diameter and component invariants exclude all targets; KW-stable B_2 ⊂ free locus; CZ acts as reflection | proved |
| blocking, any block size | qubit-interface invariant excludes all targets | proved |
| general QCA | equivalent to the original question; integrable targets excluded for symmetric J | proved / literature, assumed with DM |
| free fermions, any period | c ≤ 1 (c ≤ r at range r) | proved |
| free fermions in disguise, Pauli-aligned | none beyond JW-free | proved (enumeration) |
| spin-wave onsets, any period | at most one soft species, single magnon minimum | proved algebra, derived physics |
| exact commuting decoupling, any period | only XZ + ZX, c = 1 | proved |
| Ising × XX through XZ | needed term is four-site, induced sign attractive | proved / derived (C157) |
| Z_5 FZ at a melting tongue | chirality never symmetry-forbidden; Lagrange bound forces c = 1 floating phase in every controlled region | derived |
| Potts × half-floating boson | codimension 2 (T, K_I), ≥ 3 (N); needs a second species | derived |
| repulsive-host c = 2 | needs n ≥ 2 composites with K_2 < 1/4; pairs never soft near U(1) | derived, assumed beyond |
