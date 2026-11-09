# C172 derivations: closing A_S at order-one amplitude

This work was done by a Claude Opus 5.5 agent (C172, 2026-11-08/09), not by Codex. Astra was not used (out of credits).

Tags: [proved] complete argument or exact identity, [derived] physics-standard derivation with named inputs, [numerical] script and output named, [heuristic], [assumed]. Conventions: H = Σ_j [σ_jᵀ J σ_{j+1} + B·σ_j] with Pauli matrices, C165's family `fam` (copied verbatim into ed_kink.py), C170's corrections to C165 (codex/work_C170/RESULT.md) taken as given.

Target, in the corrected form of C170 (RESULT.md, "Consequences"): along the quantum switch set of the parent, namely the KT point Δ_eff = 1, h_eff = 0 and the Néel plateau edge h_eff = hz_c1(Δ_eff) for Δ_eff > 1, there is no second low-energy sector at order-one harmonic amplitude, in particular no bound magnon below the two-spinon threshold that condenses before the spinons.

## 1. The polar frame and the exact hidden symmetries

### 1.1 Frame [proved]

Write J = S O with O = R_m(θ) ∈ SO(3) and S symmetric (polar decomposition, S → −S and O → −O if needed). The on-site rotations σ_j = O^{−j} τ_j give exactly

    H = Σ_j [ τ_jᵀ S_j τ_{j+1} + B_j·τ_j ],   S_j = O^j S O^{−j},   B_j = O^j B.      (1.1)

Every bond is an XYZ-type exchange whose principal axes precess by θ about m from bond to bond, in a precessing field. The U(1)_m average of S_j is the parent of C154 §2, and its charge ±1, ±2 parts are the harmonics. For C165's family with only the charge-2 harmonic, the transverse block of J is M = −R(a) + g diag(1, −1), with singular values 1 ± g (checked: MMᵀ has eigenvalues (1 ± g)²). So the harmonic amplitude g is measured against the parent exchange by the polar decomposition itself, and g = 1 is where the transverse block becomes singular.

It is convenient to record the transverse block by two complex numbers. Writing the in-plane bond as A σ⁺σ⁻ + A* σ⁻σ⁺ + P σ⁺σ⁺ + P* σ⁻σ⁻,

    A = M_xx + M_yy + i(M_xy − M_yx),   P = M_xx − M_yy − i(M_xy + M_yx),   det M = (|A|² − |P|²)/4,      (1.2)

and the singular values of M are (|A| ± |P|)/2. C165's parent has |A| = 2 and its charge-2 harmonic gives |P| = 2(g² + g_xy²)^{1/2}. In the polar frame of the transverse block |P| ≤ |A| by construction, with equality exactly where M is singular.

### 1.2 Exact SU(2) exists only on a curve [proved]

Lead 2 of the spec asks for an exact identity that keeps a hidden SU(2) along a submanifold. There is none of positive dimension beyond the known curve.

Proposition 1. A bond (J, B) has an exact SU(2) symmetry with on-site generators if and only if J = λ O with O ∈ SO(3), λ ∈ ℝ, and B = 0.

Proof. An on-site SU(2) acting on a qubit is the spin-1/2 representation up to a site-dependent unitary, so its generators are Σ_j U_j σ^a_j U_j†, i.e. after the product unitary ⊗U_j the symmetry is the standard SU(2). The only SU(2)-invariant two-qubit Hermitian operators are 1 and σ·σ, and no single-site operator is invariant. Writing U_j σ U_j† = R_j σ, the transformed bond is σ_jᵀ R_jᵀ J R_{j+1} σ_{j+1} + (R_jᵀ B)·σ_j. Two-site operators σ^aσ^b are linearly independent of single-site ones, so invariance requires R_jᵀ J R_{j+1} = λ 1 and B = 0, hence J = λ R_j R_{j+1}ᵀ ∈ λ SO(3). Conversely J = λO is mapped to λ Σ σ·σ by R_j = O^{−j} (exact on open chains, a twist on rings). QED.

Modulo the global rotation, λO is labelled by the rotation angle, so the hidden-SU(2) set in B_2 is the twisted Heisenberg curve (AF for λ > 0, the h_F curve of C161 for λ < 0). The same argument with a Heisenberg chain in a uniform field gives B ∥ axis of O, which is the U(1) class at Δ = 1. Exact U(1) with on-site generators Σ n_j·σ_j reduces to XXZ, commuting Ising or the staggered-charge class (contrib/d2_gap_reduction results 3, 4). So at order-one harmonic amplitude no exact hidden symmetry is available, and lead 2 closes negatively. The SU(2)_1 structure at the KT point at O(1) amplitude is emergent (C170's selection rule), not exact.

Lead 1 (non-crossing of the SU(2)_1 towers) has the same fate. In N the only exact quantum number of a ring eigenstate is the lattice momentum. An Ising σ primary would enter at momentum π, where the SU(2)_1 tower already has N^z and the dimerization, so level repulsion is allowed and no counting argument protects the tower [proved, as a statement about available quantum numbers]. What survives of lead 1 is a different invariant, the kink parity of §3, which is exact on rings (odd versus even length).

## 2. Locking lemma: an independent momentum-π Ising sector against the liquid's own vertices

This is lead 3 (a bound on the momentum-π sector in terms of the liquid's vertex operators).

Setting [assumed for contradiction]: at a point of N the infrared is a floating boson (φ, ϑ; Luttinger parameter K, velocity v) plus an independent Ising sector (σ, ε, mass m_I). Translation acts on the Ising sector as σ → −σ (momentum π) and on vertices V_{n,m} = e^{i(2nφ + mϑ)} by their σ-frame momenta P_{n,m} = 2n k_F + mQ. Dimensions x_{n,m} = n²K + m²/4K. Define the mismatch δ_{n,m} = P_{n,m} − π mod 2π. The lattice generates every coupling

    H_c = Σ_{n,m} g_{n,m} ∫ dx σ(x) [V_{n,m}(x) e^{iδ_{n,m} x} + h.c.].      (2.1)

Lemma L [derived; inputs: standard commensurate-incommensurate RG, Giamarchi ch. 2 and 4, and C-IC soliton theory]. Let λ = 15/8 − x_{n,m} > 0 and ℓ_g = a (v/(a|g_{n,m}|))^{1/λ}, the length at which g_{n,m} reaches order one at δ = 0.
(i) If |δ| ℓ_g ≪ 1, the coupling reaches strong coupling before the oscillation cuts it off. σ and arg V lock (σ = sign cos arg V), the locked state is gapped with m_lock ~ v/ℓ_g, and δ acts on it as a chemical potential for its solitons.
(ii) For v|δ| above the soliton energy (∼ m_lock) the system enters by a Pokrovsky–Talapov transition an incommensurate phase of joint solitons, each a σ kink bound to a shift of arg V by π. This is a single boson, c = 1, and the Ising order parameter is its vertex (−1)^{N_sol}.
(iii) Only for |δ| ℓ_g ≫ 1 is (2.1) irrelevant before it grows, and the decoupled Ising × boson can survive. Second order in g then generates ε ∂_x(arg V) (dimension 2) and Ising-mass and boson renormalizations, none of which couples the sectors relevantly (consistent with C154 table 4.2, codimension 1).

So an Ising × boson point with K < 15/8 requires |δ_{n,m}| ℓ_g(n,m) ≫ 1 for every vertex with x_{n,m} < 15/8. For K ∈ [1/2, 1] (the hidden-Heisenberg region) these are (n,m) = (±1,0), (0,±1), (0,±2), (±1,±1).

Consequences.
- Stratum T [consistency check]. At B = 0 the liquid is half-floating with 2k_F = π (C159), so δ_{1,0} = 0 and (2.1) contains the commensurate σ cos 2φ, dimension 1/8 + K < 2 for K < 15/8. Lemma L then makes Ising × boson need the extra tuning g_{1,0} = 0. This is C159's count "two tunings unless K > 15/8", recovered with the second tuning identified.
- Neighbourhood of T in N [derived]. Off T, δ_{1,0} = 2π m_∥ with m_∥ the boson density induced by B (linear in the longitudinal component). A c = 3/2 hypersurface Σ ⊂ N must therefore lie in the wedge |g_{1,0}| ≲ v a^{−1}(|δ_{1,0}| a)^{λ}, which pinches to g_{1,0} = 0 as B → 0. The closure of Σ meets T only inside C159's two-tuning locus {ε = g_{1,0} = 0}. No hypersurface of c = 3/2 can approach T elsewhere.
- Dilute layers of Néel edges [derived]. If the liquid next to a Néel PT edge is a kink liquid (§3), the Néel order is its own vertex at momentum π + πn, n the kink density, and any extra Ising sector locks to it within |πn| ℓ_g ≲ 1.
- What it does not do [stated limit]. Deep in the liquid (all mismatches of order one) Lemma L allows decoupling unless the bare couplings are large, and the bare couplings are not computable at O(1) amplitude. That a lattice-scale momentum-π sector of a qubit chain couples at order one to the liquid's density vertex is the "one species per site" intuition, [heuristic].

## 3. The kink-parity dichotomy at the edge of a period-2 phase

This makes C170's "bound magnon" precise without any U(1), so that it can be posed at order-one amplitude in N.

Let G be a gapped phase in N with spontaneous period-2 order (Néel along any axis, possibly with uniform canting). Its two ground states are exchanged by translation T.

Definitions [proved, exact on rings]. On a ring of odd length the two ground states cannot be matched, and every state contains an odd number of domain walls (kinks). The kink parity is the parity of L, a topological quantum number. The kink gap is D_K = lim [E_0(L odd) − (E_0(L−1) + E_0(L+1))/2]. The kink-even gap D_2 is the gap above the quasi-degenerate doublet on even rings. Kink-even states include the two-kink continuum, which starts at 2D_K in the thermodynamic limit, so D_2 ≤ 2D_K, with D_2 < 2D_K if and only if a kink-even state lies below the continuum (a bound kink-antikink pair or another local mode, here called a bound magnon).

Proposition 2 [derived; inputs: a continuous transition out of a gapped phase is driven by a single softening excitation (C154 Lemma 3 setting), and the IR of a dilute gas of one species of quanta is one boson]. At a continuous edge of G, either
(a) kinks soften first (D_K → 0 with D_2 → 0 as 2D_K). Kink condensation restores translation, the period-2 order parameter becomes the string operator (−1)^{N_kink} = a vertex of the kink liquid's boson, the liquid is c = 1 (a single soft minimum in N gives one complex mode), and Lemma L locks any further Ising sector in the dilute layer; or
(b) a kink-even state softens first (D_2 → 0 with D_K > 0). Its quanta condense on top of an intact period-2 order: a floating + SSB phase, which is the ordered side of C154's P1 and the only route from G to a c = 3/2 point.
Case (b) requires a kink-even state below 2D_K at the softening point. Hence along the plateau-edge part of the switch set, A_S is equivalent to the two-body statement

    (B) no kink-even bound state below the two-kink continuum reaches zero energy while D_K > 0.

Any binding is fatal. If a bound state with binding E_b = 2D_K − D_2 > 0 persists into the edge, then D_2 reaches zero at D_K = E_b/2 > 0, before the kinks, and a floating + SSB window of width set by E_b opens. So (B) is equivalent to the absence of any kink-even state below the absolute two-kink bottom 2D_K in a neighbourhood of the edge. A numerical test with a resolution floor therefore cannot settle (B), and an argument has to come with a margin. §4.3 supplies one.

In the U(1) parent the kinks are the spinons (S^z = ±1/2) and the kink-even bound state is C170's bound magnon (S^z = 1). A longitudinal field shifts D_2 − 2D_K by zero in the U(1) case (both shift by the same S^z), so binding is field-independent there [proved]; numerically the P0 and P7 data below agree to all digits.

## 4. Proving (B) in the controlled regimes at order-one amplitude

### 4.1 The parent [derived, standard input]

The massive spin-1/2 XXZ chain (Δ > 1) has a spectrum generated by spinons with no bound states below the two-spinon continuum (Bethe ansatz; C170 R4.2 item 3). So (B) holds on the whole parent switch set.

### 4.2 Weak harmonics [derived, C170]

The parent's spinons scatter with S(0) = −1 and have no zero-energy resonance, so a weak attraction cannot bind them (1D binding by a weak potential needs a zero-energy resonance). (B) holds for harmonics below an O(1) threshold. §4.3 to 4.5 extend this to order-one amplitude in two regimes.

### 4.3 The Ising end of the switch set, any amplitude below the polar threshold [proved within first-order perturbation theory; derived for the chain]

Let the order axis be z with J_zz = Dz much larger than every other entry (Δ_eff → ∞, plateau edges at fields of order |A|). First-order degenerate perturbation theory is the projection of V = H − Dz ΣZZ onto the n-wall configurations (each wall costs 2Dz). In the Néel background:

- A flip-flop next to a wall moves it by two bonds, preserving its flavour (↑↑ or ↓↓): a ↓↓ wall hops right with A and left with A*, a ↑↑ wall with A* and A. Single-wall bands 2|A| cos(2k ∓ φ_A), minimum −2|A|, plus the field ∓hz. So E_1 = −2|A| − |hz|.
- A charge-2 term σ⁺σ⁺ on an isolated ↓↓ wall makes three walls (off-shell). On a contact pair (two walls at distance one, i.e. one reversed spin at site c) it moves the pair to c ± 1 and flips both flavours, amplitude P (P* backwards). This is the only first-order process that exists only at contact: the quantum version of a coherently moving magnon, and the quantum counterpart of C165's convex charge-2 quartic.
- Charge-1 terms and transverse fields move single walls by one bond with a flavour flip; at contact they annihilate the pair (off-shell, second order). The S^z = 0 two-wall sector (opposite flavours, odd separation) has no contact process at first order.

Lemma W. With charge-0 and charge-2 terms and a longitudinal field, the first-order two-wall spectrum has no state below the two-wall continuum bottom C = 2E_1 = −4|A| − 2|hz| whenever |P| ≤ |A|, that is det M ≥ 0. For |P| > |A| a bound state exists at suitable pitch (sharp at a = π).

Proof. Same-flavour sector f ∈ {↓, ↑}: walls at bonds b_1 < b_2 = b_1 + r, r odd ≥ 1. With ψ(X, r) = e^{iQX} u_f(r), X = b_1 + b_2, the four single-wall moves give

    (H u_f)(r) = t_f (u_f(r+2) + u_f(r−2)) + ε_f u_f(r),   u_f(−1) = 0,   t_{↓,↑} = 2|A| cos(2Q ∓ φ_A),   ε_{↓,↑} = ∓2hz,      (4.1)

and the contact process couples u_↓(1) and u_↑(1) with amplitude w = 2|P||cos 2Q| (phase removable). For a half-line chain with real hopping t one has the identity ⟨u, T u⟩ = |t| [Σ_r |u_r + s u_{r+2}|² − 2‖u‖² + |u_1|²] (s = sign t), so ⟨T⟩ ≥ −2|t|‖u‖² + |t||u_1|². With x = |u_↓(1)|, y = |u_↑(1)| and |t_f| ≤ 2|A|, ε_f + 2|hz| ≥ 0,

    ⟨H⟩ − C‖u‖² ≥ Σ_f [(4|A| − 2|t_f|)‖u_f‖² + |t_f| |u_f(1)|²] − 2wxy ≥ 2|A|(x² + y²) − w(x² + y²) ≥ 2(|A| − |P|)(x² + y²) ≥ 0,

using ‖u_f‖² ≥ |u_f(1)|². The S^z = 0 sector has no contact term and the same bound with w = 0. For the converse, at a = π (φ_A = 0) and Q = 0, t = 2|A|, w = 2|P| and the half-line with boundary potential −w binds for w > t with E = −(w + t²/w) < C exactly when |P| > |A|. QED.

[numerical, wall_pt.py, exact first-order projection on rings L = 28 to 41] Checks: E_1 = −3.999 (exact −4) at L = 41; at a = π the two-wall minimum is −8.000, −10.000, −13.333 at g = 1, 2, 3, matching −(w + t²/w) = −8, −10, −13.33; at a = 0.9 the threshold lies between g = 1.5 and 2, at a = π/2 at g = 2 (E = −4g, bound state at Q with t_Q = 0). Charge-1 terms, not covered by the proof: a random scan of 120 bonds (wall_scan.py, outputs wall_scan_1.out, wall_scan_2.out; pitch uniform, |P|/|A| ∈ [0, 1.4], charge-1 vector of norm up to 1.5, |hz| ≤ 1) finds no bound state (minimum of E_2 − 2E_1 is +0.047, a finite-size positive value), and a targeted scan at the sharp pitch (wall_target.py, wall_target.out; a = π ± 0.3, |P|/|A| = 0.8 to 1.2, charge-1 norm 0, 0.5, 1, 2, random hz) finds that charge-1 terms and fields only raise E_2 − 2E_1, with the only bound case at |P|/|A| = 1.2 and zero charge-1. So charge-1 terms do not lower the threshold [numerical].

Meaning. In the Ising regime the kinks win for every harmonic content as long as the transverse block of the bond is not singular. Since the polar parent of the transverse block has |P| ≤ |A| by construction (§1.1; for C165's family this is g² + g_xy² ≤ 1), the whole hidden-Heisenberg harmonic range of the Ising end of the switch set satisfies (B). The threshold |P| = |A| is the locus where the transverse block is singular (an XZ-type precessing chain); beyond it the natural parent is a different one, obtained by a π rotation of every other site about an in-plane axis, which exchanges A and P and turns the Néel order into a ferromagnet whose elementary excitation is the magnon. The bound magnon at |P| > |A| is then not a second sector but the first sector of the other parent.

Margin [proved within the first-order model]. The proof gives more than E ≥ C. It gives ⟨H⟩ − C‖u‖² ≥ 2(|A| − |P|)(x² + y²) plus the kinetic forms |t_f| Σ|u_f(r) + s u_f(r+2)|². An extra attractive potential −U_1 at contact is therefore harmless for U_1 ≤ 2(|A| − |P|). A potential −U_3 at separation 3 is controlled by splitting |u(3)|² between ‖u‖² (slack 4|A| − 2|t_f|, large where |t_f| is small) and the kinetic form (|u(3)|² ≤ 2(|u(1) + s u(3)|² + |u(1)|²), useful where |t_f| is large), so it is harmless below a threshold of order |A| − |P| [derived; the constant is not optimized]. This is the quantitative form of C170's "no zero-energy resonance": the hard-core (Dirichlet) boundary leaves a finite binding threshold.

For the chain at finite Dz, second order adds short-range terms of size O(V²/Dz), V the largest transverse coupling (virtual annihilation of a contact pair by charge-1 terms, flavour-changing charge-2 loops, longer hops). By the margin, Lemma W holds for the chain when V²/Dz ≪ |A| − |P| [derived]. It fails only in a thin layer next to the singular transverse block |P| = |A|.

### 4.4 The KT end of the switch set [derived, conditional]

Near Δ_eff → 1⁺ (plateau edge h_eff = hz_c1 → 0, exponentially small), D_K is exponentially small and the kinks have size ξ ≫ a. The kink-even channel at energies of order D_K is then fixed by the infrared theory. If the liquid at the KT point is in the SU(2)_1 basin (C170: exact selection rule, no momentum-π relevant operator, charged ones only in resonance tongues), this is the massive SU(2)-Thirring regime of sine-Gordon at β² → 8π, whose spinons have no bound states (breathers need β² < 4π). Irrelevant operators change the binding by (D_K/Λ)^y → 0, with Λ the lattice scale at which the harmonics are resolved. So (B) holds in a neighbourhood of the KT end at any amplitude for which the KT point is in the SU(2)_1 basin. The width of the neighbourhood is set by Λ and is not computed. An extra gapless sector exactly at the KT point would be a separate coincidence of codimension at least 3 in N, and it is not excluded here.

### 4.5 Intermediate anisotropy [numerical]

ed_kink.py: full ring ED, L = 10 to 18, b(L) = D_2(L) − [D_K(L−1) + D_K(L+1)]. The rings carry the pitch twist L(a + π), so b(L) has an oscillating finite-size part common to a bond and its parent. The test is the difference to the parent at the same a and Dz, which a bound state of binding E_b would shift by −E_b.

Parents (no harmonics; U(1), so unbound by §4.1): b(12, 14, 16) = +0.83, +0.76, +0.68 (a = 0.9, Dz = 2, both hz = 0 and hz = 0.6, identical as §3 requires), +0.63, +0.54, +0.53 (Dz = 1.6), +0.78, +0.49, +0.92 (a = 2.2, Dz = 2). These positive values are the finite-size kinetic energy of two hard-core kinks on L ≤ 18 rings with a twist, not physics.

| bond | Dz | ratio abs(P)/abs(A) | ‖harmonics‖_F | D_K(17) | D_2(16) | db(12, 14, 16) |
|---|---|---|---|---|---|---|
| P1 = C165 Q2 harmonics, hz 0.05 | 2.0 | 0.30 | 0.43 | 0.95 | 2.52 | −0.024, −0.018, −0.013 |
| P2 = same, hz 0.6 (control P7) | 2.0 | 0.30 | 0.43 | 0.41 | 1.44 | −0.020, −0.013, −0.008 |
| P5: g 0.6, G 0.3, hx 0.3, Dx 0.1 | 2.0 | 0.60 | 0.96 | 0.84 | 2.34 | −0.027, −0.029, −0.027 |
| P8: all nine harmonic types | 2.0 | 0.40 | 0.66 | 0.84 | 2.30 | −0.048, −0.040, −0.033 |
| Q1: g 0.9, G 0.4, hx 0.3, Dx 0.1 | 2.0 | 0.90 | 1.40 | 0.77 | 2.18 | +0.045, −0.038, −0.066 |
| Q2: g 0.8, gxy 0.4, G 0.5, G' −0.3, Dy 0.3, hx 0.4, hy 0.2, hz 0.3 | 2.0 | 0.89 | 1.28 | 0.53 | 1.66 | −0.011, −0.083, −0.105 |
| R1 = P5 harmonics | 1.6 | 0.60 | 0.96 | 0.44 | 1.38 | +0.055, −0.011, −0.049 |
| R2 = P8 harmonics | 1.6 | 0.40 | 0.66 | 0.43 | 1.36 | −0.048, −0.040, −0.033 |
| S1: pitch a = 2.2, g 0.7, G 0.3, hx 0.3, Dx 0.1 | 2.0 | 0.70 | 1.09 | 0.90 | 2.75 | −0.097, −0.012, −0.039 |

(ed_summary.txt; raw spectra and momenta in ed_*.out, ed_*.json.) Field path past the parent edge hz_c1 = 0.78 at Dz = 2: at hz = 1.0 with P1's harmonics D_K(13, 15, 17) = 0.03, −0.03, 0.02 while D_2(14, 16, 18) = 0.74, 0.66, 0.66, and with P5's harmonics D_K = 0.01, −0.06, −0.05 while D_2 = 0.70, 0.56, 0.61. At hz = 1.4 D_K is negative (the odd ring is the cheaper one, a magnetized liquid). So along both field paths the kinks soften first, case (a) [numerical].

Reading [numerical]. At harmonic amplitudes up to ‖·‖_F = 1.4 (charge-2 up to 90 percent of the polar threshold, charge-1 up to 0.5, all nine types) and Dz = 1.6 and 2, the harmonics shift the kink-even gap relative to twice the kink gap by at most 0.105, against 2D_K = 0.9 to 1.9. No strongly bound magnon exists at these points. A weak bound state (E_b ≲ 0.1) is not resolved, and by §3 even a weak one would be fatal, so these data support (B) at intermediate anisotropy without proving it. The trend of Q2 (db growing in magnitude with L) is the one to watch. It is also what a twist oscillation of the kink band does, as the parent S0 shows (b non-monotone in L).

## 5. Summary of what is derived, and the residual

- [proved] (§1.2) Exact on-site SU(2) exists only on the twisted Heisenberg curve, and exact U(1) only on the known classes. No exact identity reaches order-one amplitude. Ring momentum is the only exact quantum number in N, so SU(2)_1 towers have no non-crossing protection. Leads 1 and 2 of the spec close negatively, and their usable remnant is the kink parity.
- [derived] (§2, Lemma L) An Ising × boson point with K < 15/8 needs every relevant vertex of the liquid to be detuned from momentum π by more than its locking scale. On T this is C159's second tuning. Near T, any c = 3/2 hypersurface of N is pinched onto C159's two-tuning locus. In the dilute layer of a kink liquid every extra Ising sector locks.
- [derived] (§3, Proposition 2) At an edge of a period-2 phase in N, a c = 3/2 route exists only if a kink-even excitation softens before the kink, which needs a kink-even state below the two-kink continuum. This is C170's bound magnon in a form that needs no U(1). Any binding is fatal, so the statement must be exact.
- [proved in first-order perturbation theory, derived for the chain] (§4.3, Lemma W) In the Ising regime of the switch set, with charge-0 and charge-2 harmonics of any amplitude and any longitudinal field, no kink-even state lies below the continuum when |P| ≤ |A|, i.e. when the transverse block has non-negative determinant, which holds whenever C165's parent is the polar parent of the transverse block (g² + g_xy² ≤ 1). The bound comes with a finite margin 2(|A| − |P|), so second-order terms O(V²/Dz) do not change it away from |P| = |A|. The threshold is sharp (a bound state appears for |P| > |A| at a = π, energy −(w + t²/w) reproduced to four digits at L = 30). Charge-1 terms: numerical, 120 random and 72 targeted bonds, they only raise the two-wall energy.
- [derived, conditional on the KT point being in the SU(2)_1 basin] (§4.4) Near the KT end the kink-even channel is fixed by the massive SU(2)-Thirring spinons, which have no bound states, at any amplitude.
- [numerical] (§4.5) At Dz = 1.6 and 2 with amplitudes up to 1.4, no bound magnon with binding above about 0.1, and kinks soften first along field paths through the edge.

Residual A_S′ (strictly narrower than A_S as corrected by C170). On the plateau-edge part of the switch set at intermediate anisotropy (between the KT crossover of §4.4 and the Ising regime of §4.3), at order-one harmonic amplitude with non-singular transverse block (|P| < |A|), the two-kink channel has no bound state below the two-kink continuum. Equivalently, the hard-core repulsion of two spin-1/2 domain walls is not overcome there by the harmonics. It is a two-body statement about a gapped background, so small-system numerics with resolution well below 0.1 can test it directly (DMRG of the odd/even sectors or a two-kink variational ansatz on top of VUMPS), unlike the original A_S. The interior of the liquid away from the switch set remains under C154 (A6), with Lemma L's constraint that any independent Ising there be detuned from every relevant vertex by more than its locking scale.

What is not done. No proof at intermediate anisotropy. Lemma W with charge-1 terms is numerical. The KT-end statement assumes the SU(2)_1 basin at that amplitude and leaves the width of its neighbourhood uncomputed. The ED resolution (about 0.1) is not enough to exclude a weak bound state.
