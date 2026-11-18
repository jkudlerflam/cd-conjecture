# C201: proofs

2026-11-18. Manager: Claude Fable 5.2. Derivations delegated to Astra (gpt-6.2-astra via the shim, five calls: three derivations and two referee reports, each output ending with a trailer that names that model). Every Astra result below was read line by line by me and tested with my own scripts in this folder. Nothing here has had human review. Tags: [rigorous], [rigorous given cited result X], [numerical], [physics-level]. Each section names the model that produced it.

## 0. Conventions and the target

A bond is a Hermitian h on ℂ^d ⊗ ℂ^d. H_L = Σ_{j∈ℤ_L} h_{j,j+1} on the ring, with lowest eigenvalue E_L, and e₀ = lim E_L/L, which equals the infimum of tr(h ρ₂) over nearest-neighbour marginals of translation-invariant states on ℤ. M(h) is the minimum over Hermitian X of half the spectral spread of h + X ⊗ 1 − 1 ⊗ X. Ω_L and Ω_∞ are the sets of nearest-neighbour two-site marginals of translation-invariant states on the L-ring and on ℤ. A = Σ_a c_a v_a is defined through the infrared assumption E_L = L e₀ − πA/(6L) + o(1/L) on vacuum-compatible lengths.

Lemma P is the statement E_L ≥ L e₀ − K(d) M/L for all h and L. It gives A ≤ (6/π) K(d) M. It is open. Write

    g_d(L) = sup { L e₀(h) − E_L(h) : h on ℂ^d ⊗ ℂ^d, M(h) = 1 }.

Since E_L/L and e₀ are the minima of tr(hρ) over Ω_L and Ω_∞, and M is the norm dual to the trace norm on differences of marginals, g_d(L) = L·D_L with D_L the one-sided Hausdorff distance from Ω_L to Ω_∞ in trace norm. (My spec had L² here. Astra corrected it.) The cut-and-tile bound is g_d(L) ≤ 2, Lemma P is g_d(L) ≤ K(d)/L, and before this round nothing was known in between.

## 1. Targets, ranked

1. (B1) Lemma P for interacting classes by exact spectral reduction to classes already proved. Chosen first because the reduction is finite-dimensional algebra and can be tested by exact diagonalization. Outcome: Theorem B (exchange class, every d ≥ 2, unconditional at V = 0) and Theorem C (critical ferromagnetic Potts chains, d = 3, 4, given the project's XXZ theorem).
2. A fully general but weaker statement, g_d(L) → 0. Chosen because it is the weakest statement that bears on an extensive fast tower (it excludes one only at physics level, Corollary A3), it is purely kinematic, and no earlier round attempted it. Outcome: Theorem A, unconditional given a published recovery theorem, for the ring defect and for the doubling defect.
3. (B3) A fast tower. Outcome: none found. Two no-go lemmas (Section 6) and a correction to the finite-size signature used in earlier searches (Section 7.4).
4. (B2), (B4). Not attacked separately. Theorem A is the only statement here that holds for all chains, and it does not bound A (Section 2.10).

The full target, a finite K(d) for all chains, was not reached. Section 7 records why each attempted strengthening stops.

## 2. Theorem A: the ring defect and the positive part of the doubling defect vanish uniformly at fixed d

First proof by Astra (gpt-6.2-astra via the shim, `astra_A2_thirdlaw.out`), from a spec of mine that posed the uniform statement and listed the square-root loss of recovery maps as the known obstruction. A second, independent Astra call (gpt-6.2-astra, `astra_A4_referee_T1.out`) refereed that proof, found it correct as stated, and gave the shorter proof with a deterministic output length that is written out below, which also yields the doubling statement. I checked both versions line by line and tested the new lemmas numerically (`t1_lemma_checks.py`, `t1_chainrule_check.py`).

### 2.1 Statement

Fix d ≥ 2, an integer n ≥ 1 and 0 < δ ≤ 1. Define B₁ = 1 and, for i = 1, …, n,

    t_i = B_i + 1,   q_i = d^(t_i),   K_i = ⌈2 (2q_i − 1)² δ⁻² ln q_i⌉,   B_{i+1} = (K_i + 1) t_i,

and f = f_d(n, δ) = 3δ + 2√(2δ) + 4√(2 ln d / n).

Theorem A. [rigorous given cited theorem R] Let L ≥ 2 B_{n+1} and let h be any bond on ℂ^d ⊗ ℂ^d.

(A1) For every translation-invariant state ρ on the L-ring with nearest-neighbour marginal ρ₂, L (e₀(h) − tr(h ρ₂)) ≤ f M(h). In particular L e₀ − E_L ≤ f M, so g_d(L) ≤ f.

(A2) For every integer k ≥ 2, E_{kL} − k E_L ≤ k f M. In particular E_{2L} − 2E_L ≤ 2 f M.

Taking δ = n⁻² and n → ∞, g_d(L) → 0 and sup_h (E_{2L} − 2E_L)/M → 0 as L → ∞ at fixed d, and Ω_L lies in the (f/L)-neighbourhood of Ω_∞ in trace norm.

Theorem R (Junge, Renner, Sutter, Wilde, Winter, Ann. Henri Poincaré 19, 2955 (2018), arXiv:1509.07127, Corollary 4.1 for the inequality, together with Remark 2.4(2) for the exact recovery of the reference). For a state ω_ABC on a finite-dimensional space there is a completely positive map 𝓡 from B to BC, depending only on ω_BC, trace preserving on the support of ω_B, with 𝓡(ω_B) = ω_BC and I(A:C|B)_ω ≥ −2 ln F(ω_ABC, (id_A ⊗ 𝓡)(ω_AB)), where F(α, β) = ‖√α √β‖₁. It is an average of rotated Petz maps and is completed to a channel by a fixed preparation off the support of ω_B. The C203 red team (Claude Opus 5.5, `../work_C203/REDTEAM_A.md` §1) read the paper and confirmed these items: Corollary 4.1 gives the inequality for every finite-dimensional ω_ABC with no full-support hypothesis and with the unsquared fidelity, and the exact recovery 𝓡(ω_B) = ω_BC is Remark 2.4(2) applied to σ = 1_A ⊗ ω_BC and the partial trace over C. Universality of the map is never used below. The proof applies the theorem to the single state τ_ABC, so any rotated-Petz result with exact recovery of τ_BC would do.

The threshold B_{n+1} grows like a tower of exponentials of height n in n (each step is one exponential in the previous buffer size), so f decays in L more slowly than any iterated logarithm. The bound beats the trivial bound 2M only when 4√(2 ln d/n) + 3δ + 2√(2δ) < 2. At d = 2 this needs n ≥ 6 with δ of order 5·10⁻⁴ or smaller, or n ≥ 8 with δ = 10⁻², and in either case L ≥ 2B₇ or more, an iterated-exponential scale (already B₂ is of order 10⁶ to 10⁹ for these δ, and at n = 1, δ = 1 the threshold is L ≥ 548 with f ≈ 10.5). So Theorem A is qualitative and has no content at any length a computer or a physical system reaches (C203 red team, `../work_C203/REDTEAM_A.md` §5). The first proof had two exponentials per step and bounded g_d(L) by 3/n + 2√(2/n) + 8√(ln d/n), without (A2).

### 2.2 Idea

Tiling ℤ with cut copies of the ring state costs O(M) at one bond per L sites. To lower the cost per join to o(M), each join is made by a recovery channel that re-creates a short interval between two copies and reproduces the interior of that interval without error, so that only the two bonds at its ends are disturbed. The long arc of the ring is kept and the short interval is re-created, so there is one join per L sites whatever scale the argument selects. Two things have to be made small. The recovery error is controlled by a conditional mutual information, made small by choosing a scale on which the entropy increment has not moved (a pigeonhole over n scales, the only place where the factor ln d/n arises). And the two ends that meet at a join belong to independent copies, whereas the recovery channel expects the correlations that the two ends of an interval have in the ring. This is repaired by conditioning all copies on one shared classical measurement record, taken in the part of the ring that is cut out, with the length of that part chosen by a chain rule so that the record carries the correlations between the two ends.

### 2.3 A scale from the entropy

[rigorous] Let ρ be a translation-invariant state on the L-ring, S_j the entropy in nats of j consecutive sites (S₀ = 0) and u_j = S_{j+1} − S_j. Strong subadditivity gives S_{j+2} + S_j ≤ 2S_{j+1}, so u_{j+1} ≤ u_j. Weak monotonicity, S(a|D) + S(a|D′) ≥ 0 for a site a between two disjoint j-site intervals D, D′, gives u_j ≥ 0 whenever 2j + 1 ≤ L. And u₀ = S₁ ≤ ln d. The n nonnegative drops u_{B_i − 1} − u_{B_{i+1} − 1}, i = 1, …, n, sum to at most u₀ ≤ ln d, so one of them is at most ln d/n. Fix this i and write b = B_i, t = b + 1, q = d^t, K = K_i, B′ = B_{i+1} = (K + 1)t ≤ L/2.

For an interval of m consecutive sites with 2t ≤ m ≤ B′, let A = A_ℓ A_r be its two end sites, B = B_ℓ B_r the b sites next to each end, and C the remaining m − 2t sites. By subadditivity of conditional entropy and strong subadditivity, S(A|B) ≤ S(A_ℓ|B_ℓ) + S(A_r|B_r) = 2u_b. Removing the two end sites leaves an interval of m − 2 sites, so S(A|BC) = S_m − S_{m−2} = u_{m−1} + u_{m−2} ≥ 2u_{m−1}. Hence

    I(A:C|B) ≤ 2(u_b − u_{m−1}) ≤ 2(u_{B_i − 1} − u_{B_{i+1} − 1}) ≤ 2 ln d / n.          (2.1)

### 2.4 A measurement with a cheap inverse

[rigorous] For a Haar-random rank-one projector P on ℂ^q, ∫ P dP = 1/q and ∫ P ⊗ P dP = (1 + Swap)/(q(q+1)). By Carathéodory there are finitely many projectors P_y and weights w_y > 0 with the same two moments, so {w_y q P_y} is a measurement. With X_y = (q + 1)P_y − 1,

    Σ_y w_y q tr(γ P_y) X_y = q(q+1) tr₁[(γ ⊗ 1) Σ_y w_y P_y ⊗ P_y] − (tr γ) 1 = γ

for every operator γ, and ‖X_y‖₁ = q + (q − 1) = 2q − 1. Consequently, if Z = Σ_y Z_y ⊗ |y⟩⟨y| is the difference of two classical-quantum states obtained by measuring one tensor factor, the operator Σ_y Z_y ⊗ X_y it comes from has trace norm at most (2q − 1) ‖Z‖₁.

### 2.5 The gap: the two ends are nearly independent given a record

[rigorous] Take consecutive t-site windows U, W₁, …, W_K on the ring, (K + 1)t = B′ sites in all. Let Y_j be the outcome of the measurement of Section 2.4 on W_j. The Holevo bound and the chain rule give

    Σ_{j=1}^K I(U : Y_j | Y_{<j}) = I(U : Y₁ … Y_K) ≤ S(U) ≤ ln q,

so there is a j with I(U : Y_j | Y_{<j}) ≤ ln q / K. Fix it. Let v = Y_{<j} be the record of the windows W₁, …, W_{j−1}, with probability p_v, and ω^v the conditional state of the rest of the ring. For each v, Pinsker's inequality gives ‖ω^v_{U Y_j} − ω^v_U ⊗ ω^v_{Y_j}‖₁ ≤ √(2 I(U:Y_j)_{ω^v}), and averaging with Jensen's inequality gives Σ_v p_v ‖ω^v_{U Y_j} − ω^v_U ⊗ ω^v_{Y_j}‖₁ ≤ √(2 ln q/K). Undoing the measurement of W_j with Section 2.4, and using Σ_v p_v ω^v_{U W_j} = ρ_{U W_j} (the record is taken on other sites),

    ‖ρ_{U W_j} − Σ_v p_v ω^v_U ⊗ ω^v_{W_j}‖₁ ≤ (2q − 1) √(2 ln q / K) ≤ δ.              (2.2)

W_j is measured only inside this proof. In the construction only the windows W₁, …, W_{j−1} are measured. The two ends may keep order-one correlations (a GHZ-type state does). They are carried by the record.

### 2.6 The interval to be re-created and its recovery channel

[rigorous given theorem R] Let 𝕀 be the interval from the first site of U to the last site of W_j. It has m = (j + 1)t sites, with 2t ≤ m ≤ B′. Its first t sites are U = A_ℓ B_ℓ and its last t sites are W_j = B_r A_r, in the notation of Section 2.3, and C consists of the windows W₁, …, W_{j−1}. Let τ_ABC be the state of 𝕀 in ρ, so τ_AB = ρ_{U W_j}. By (2.1) and theorem R there is a channel 𝓡 from B to BC with 𝓡(τ_B) = τ_BC and, using ‖α − β‖₁ ≤ 2√(1 − F²) and 1 − e^(−x) ≤ x,

    ‖(id_A ⊗ 𝓡)(τ_AB) − τ_ABC‖₁ ≤ 2√(I(A:C|B)) ≤ 2√(2 ln d / n).                        (2.3)

### 2.7 Correcting the input marginal

At a join the input will be η_AB = Σ_v p_v ω^v_U ⊗ ω^v_{W_j}, not τ_AB. Feeding η_B to 𝓡 would leave an error on every bond inside the re-created interval, and there can be B′ of them. So the marginal on B is first corrected to τ_B without error.

Lemma (marginal correction). [rigorous] Let r, s be states on B and z = ‖r − s‖₁/2. There is a channel 𝒟 on B with 𝒟(r) = s such that every state ξ_AB with marginal r on B satisfies ‖(id_A ⊗ 𝒟)(ξ) − ξ‖₁ ≤ 2√z + z.

Proof. Let T = s + (r − s)₊. Then T − s = (r − s)₊ ≥ 0, T − r = (r − s)₋ ≥ 0 and tr(T − s) = z. Let K = √s T^(−1/2) on the support of T. Then K†K ≤ 1, K T K† = s and K r K† ≤ s. Define 𝒟(X) = K X K† + tr[(1 − K†K) X] ω with ω the normalized residual s − K r K† (any state if the residual vanishes). Then 𝒟(r) = s. The weight of the second branch on any extension of r is tr[(1 − K†K) r] ≤ tr[(1 − K†K) T] = z. Further ‖(1 − K)√r‖₂² = tr[(1 − K) r (1 − K)†] ≤ tr[(1 − K) T (1 − K)†] = ‖√T − √s‖₂² = tr T + tr s − 2 tr(√T √s) ≤ tr T − tr s = z, where the last inequality uses √T ≥ √s (operator monotonicity of the square root), so tr(√T √s) = tr(s^(1/4) √T s^(1/4)) ≥ tr s. On an extension ξ the same quantity equals ‖(1 ⊗ (1 − K))√ξ‖₂². Writing K ξ K† − ξ = (K − 1) ξ K† + ξ (K† − 1) and factoring each term through √ξ, the Hilbert–Schmidt Cauchy–Schwarz inequality bounds the trace norm by 2√z. Adding the second branch gives 2√z + z. ∎

[rigorous] Apply the lemma with r = η_B and s = τ_B. By (2.2), z ≤ δ/2. The join channel 𝒥 = 𝓡 ∘ 𝒟 satisfies

    𝒥(η_B) = τ_BC,   ‖(id_A ⊗ 𝒥)(η_AB) − τ_ABC‖₁ ≤ (√(2δ) + δ/2) + δ + 2√(2 ln d/n) = f/2,      (2.4)

by the triangle inequality through 𝓡(𝒟(η_AB)), 𝓡(η_AB), 𝓡(τ_AB), with contractivity of the trace norm, the lemma, (2.2) and (2.3). The equality controls every bond inside the re-created interval. The inequality is needed only for the two bonds at its ends.

### 2.8 Joining long arcs and counting the energy

[rigorous] Measure the windows W₁, …, W_{j−1} on the ring and keep the complementary arc, which starts at the first site of W_j, runs around the far side of the ring and ends at the last site of U. It has N = L − (j − 1)t sites. Let ξ_v be its conditional state given the record v. Its first t sites are in the state ω^v_{W_j} and its last t sites in the state ω^v_U, with the buffers B_r and B_ℓ at the two extreme ends and the sites A_r and A_ℓ next to them. The unconditional state Σ_v p_v ξ_v has marginal ρ₂ on every bond of the arc, since the measurement acts on other sites. Since N ≥ L − B′ ≥ L/2 ≥ B′ > 2t, the arc has at least two sites outside its buffers.

For k ≥ 2 place k arcs in cyclic order in the state Σ_v p_v ξ_v^{⊗k}. At each of the k joins the last t sites of one arc and the first t sites of the next belong to different copies, so their joint state is η_AB. Apply 𝒥 to the two buffers at every join. The k channels act on disjoint sites and preserve the trace, so the state on A, B, C at one join after all channels is (id_A ⊗ 𝒥)(η_AB). Each join adds m − 2t = (j − 1)t sites, so the result is a state on the ring of length k(N + (j − 1)t) = kL.

Choose the gauge and the additive constant of the bond so that ‖h‖ = M, and let e = tr(h ρ₂). Count the energy bond by bond.

- Each arc keeps N − 2b sites outside its buffers and N − 2b − 1 bonds between them. The channels act elsewhere, so each such bond has marginal ρ₂ and energy e.
- At each join the buffers and the re-created sites are in the state τ_BC, an interval of ρ with m − 2 sites, so each of its m − 3 bonds has energy e.
- At each join two bonds connect A_ℓ and A_r to the ends of that interval. In τ_ABC they have marginal ρ₂, so by (2.4) each has energy at most e + M f/2.

Per arc and join this is (N − 2b − 1) + (m − 3) + 2 = N + m − 2t = L bonds, so every bond of the kL-ring is counted once, and the energy of the trial state is at most kL e + k M f. Hence

    E_{kL} ≤ kL tr(h ρ₂) + k f M     for every translation-invariant state ρ on the L-ring and every k ≥ 2.

Taking ρ the translation average of a ground state of H_L gives (A2). Dividing by kL and letting k → ∞, with E_{kL}/(kL) → e₀, gives (A1). ∎

Remarks. The state ρ is arbitrary, pure or mixed, so (A1) is the kinematic inclusion of Ω_L in a neighbourhood of Ω_∞, and h is only a witness. On lengths where the ring is frustrated, E_L > L e₀, statement (A1) is empty and (A2) is still meaningful. For h = Z ⊗ Z on odd rings, L e₀ − E_L = −2M, so the absolute value of the defect does not tend to zero.

### 2.9 Checks and referee

[numerical, Claude Fable 5.2] `t1_lemma_checks.py`. (a) The marginal-correction lemma on 400 random extensions with z between 5·10⁻⁴ and 0.5: 𝒟(r) = s to 10⁻⁸, and the largest ratio of ‖(id ⊗ 𝒟)(ξ) − ξ‖₁ to 2√z + z is 0.755. (b) The covariance identities of the first proof on ring ground states at L = 12 (Heisenberg, XX, transverse Ising at g = 1/2, GHZ): they hold in all twelve cases. `t1_chainrule_check.py`: with single-site windows and the octahedral measurement on a qubit (q = 2, ‖X_y‖₁ = 3), on the same four ring states and j = 1, …, 4, the left side of (2.2) is below 3√(2 I(U:Y_j|Y_{<j})) in all sixteen cases (ratios between 0.19 and 0.49), and Σ_j I(U:Y_j|Y_{<j}) ≤ S(U). For the GHZ state the left side falls as (2/3)^(j−1), which shows the record absorbing the long-range order. These test the lemmas, not the theorem, whose lengths are out of reach.

[Astra gpt-6.2-astra, `astra_A4_referee_T1.out`] The referee call audited every step of the first proof and reported no gap, with three details to add (the completion of the recovery map off the support, the extraction of the connecting bonds from the direct sum, and the distinction between conditioning on one label and on all labels, which is harmless because the energy is linear and the lower bound affine in the length). It ran four sanity cases (GHZ order, product words of incommensurate period, the frustrated classical bond, and the half-filled hopping ring, whose marginal distance is of order L⁻²), all consistent. It also gave a fixed-d obstruction to shortening the scales: for d = 4 there are translation-invariant ring states (antisymmetric states of Q = 2^(t/4) encoded systems) whose two t-site end windows stay at constant trace distance from every shared-record mixture of products over spans polynomial in t. So the decoupling step (2.2) on full windows needs spans exponential in the buffer size. The obstruction concerns that step only, not every possible join. I read that construction and did not verify it independently.

[Claude Opus 5.5 red team, `../work_C203/REDTEAM_A.md`] An independent red team attacked the citation, the scale selection, the decoupling, the marginal correction, the join geometry, the bond count, the gauge handling of M, the k-fold length bookkeeping and the incommensurate-length cases, and found no gap. Its verdict is that (A1) and (A2) are correct given Corollary 4.1 with Remark 2.4(2), confirmed from the text of the paper. Its literature search (Fawzi, Fawzi and Scalet, arXiv:2305.06855, Blakaj and Wolf, arXiv:2308.14585, Kull, Schuch, Dive and Navascués, arXiv:2212.03014) found only O(1/n) rates, so the o(1/L) inclusion appears to be new as a statement, assembled from standard tools (universal recovery, the Brandão–Harrow chain rule, Markov-type gluing). The search was not exhaustive. Its corrections to the stated consequences are applied in Sections 2.1 and 2.10.

### 2.10 What Theorem A does and does not give

Corollary A1. [rigorous given theorem R] For every fixed bond, limsup_L (L e₀ − E_L) ≤ 0. No translation-invariant nearest-neighbour chain has a ring energy that stays a fixed amount below L e₀.

Corollary A1′ (third law, only for chains with a positive spatial transfer operator). [rigorous given theorem R and the positivity results of REPORT §9.88; the argument is from the C203 red team, §7(i), and I checked it] Let s̃_L = L E_{L+1} − (L+1) E_L, which under positivity is the entropy per unit Euclidean time of the crossed Gibbs state at crossed temperature 1/L. Then lim s̃_L = 0. Proof: positivity gives concavity of E_L in L and E_{2L} ≥ 2E_L. The second gives E_L/L ≤ E_{2^k L}/(2^k L) → e₀, so g_L = L e₀ − E_L ≥ 0, and the first makes g_L convex. Corollary A1 gives limsup g_L ≤ 0, so g_L → 0, and a convex nonnegative sequence that tends to zero is nonincreasing. A direct computation gives s̃_L = (L+1) g_L − L g_{L+1} = g_L + L (g_L − g_{L+1}) ≥ g_L ≥ 0, and s̃_L − s̃_{L−1} = L (2g_L − g_{L+1} − g_{L−1}) ≤ 0, so s̃_L decreases to a limit c ≥ 0. If c > 0 then g_L − g_{L+1} ≥ (c − g_L)/L ≥ c/(2L) for large L, and Σ (g_L − g_{L+1}) would diverge, which contradicts g_L ≥ 0. So c = 0. In the crossed channel this says that the ground state has no residual entropy. Without positivity the statement is false: for h = Z ⊗ Z, s̃_L equals 2L on even L and −2L − 2 on odd L.

Corollary A2. [rigorous given theorem R] The positive part of the doubling defect vanishes uniformly: limsup_L sup_h (E_{2L} − 2E_L)/M(h) ≤ 0 at fixed d, and the supremum is at least zero (frustration-free bonds), so it tends to zero. The same holds for E_{kL} − kE_L divided by k, uniformly in k. For an individual bond the defect can stay negative of order M (h = Z ⊗ Z on odd L has E_L = −L + 2 and E_{2L} = −2L, so E_{2L} − 2E_L = −4 with M = 1), so only the positive part is controlled. Before this the only general upper bound on the doubling defect was 2M.

Corollary A3. [physics-level] For a family of chains at fixed d whose ring energies have the form E_L = L e₀ − (π/6L) Σ_{a: ξ_a ≤ L} c_a v_a + r_L on compatible lengths, with ξ_a the crossover length of sector a and |r_L|/M → 0 uniformly in the family,

    Σ_{a: ξ_a ≤ L} c_a v_a ≤ (6/π) L (g_d(L) + |r_L|/M) M = o(L) M.

So the number of sectors with central charge at least c_min and velocity at least ηM that are visible at scale L is o(L)/(η c_min). An extensive tower is excluded, at physics level and only under the uniformity hypothesis on r_L, which does most of the work, and only asymptotically, since below the scale 2B₇ the bound is weaker than the trivial one. A fixed-chain asymptotic form with a remainder that is not uniform in the family does not justify this evaluation at the crossover scale.

What it does not give. A bound on A needs g_d(L) ≤ K/L. A family whose N-th member has deficits L e₀ − E_L = M N/(L + N²) has supremum over N equal to M/(2√L) → 0, so it obeys Theorem A, yet lim_L L (L e₀ − E_L) = M N is unbounded. Theorem A is strictly weaker than the target and is compatible with C(d) = ∞.

## 3. Theorem B: the exchange class at arbitrary d

Produced by Astra (gpt-6.2-astra via the shim, `astra_A3_adversary.out` Section 1), prompted by a weaker label-conserving version of mine (Claude Fable 5.2, `xxc_check.py`). I rechecked every step below by hand and numerically (`exchange_class_check.py`).

### 3.1 Setting

Split the site space as ℋ = ℋ₀ ⊕ ℋ₁ with d₀ = dim ℋ₀ ≥ 1, d₁ = dim ℋ₁ ≥ 1, d = d₀ + d₁, and let n be the projector on ℋ₁. Let W be any unitary from ℋ₁ ⊗ ℋ₀ onto ℋ₀ ⊗ ℋ₁, extended by zero to the rest of ℋ ⊗ ℋ. The bond is

    h_W = −t (W + W†) + V n ⊗ n + (μ/2)(n ⊗ 1 + 1 ⊗ n),   t > 0.

The scalar bond h_sc is the case d₀ = d₁ = 1, W = |01⟩⟨10|, which is the spin-1/2 XXZ chain in a field, with J = 2t and Δ = V/(2t) in the convention H = J Σ (SˣSˣ + SʸSʸ + Δ SᶻSᶻ) up to a flux and a field. For d ≥ 3 the unitary W can entangle the internal states of the two species at every exchange, so h_W is not a quasi-free bond in general, and it is not covered by the label-order argument, which needs W to be the flip.

### 3.2 Statement

Theorem B. For every d₀, d₁, W, V, μ, t and every L ≥ 2,

(B1) E_L(W) ≥ E_L(sc), e₀(W) = e₀(sc), and M(h_W) ≥ t. [rigorous]

(B2) If V = 0, then L (L e₀ − E_L) ≤ 2tL (csc(π/L) − L/π) ≤ (4 − 8/π) t ≤ (4 − 8/π) M, with equality at d = 2, μ = 0, L = 2, and limsup_L L (L e₀ − E_L) ≤ π t/3. [rigorous]

(B3) If V ≤ −2t, then E_L = L e₀ for every L. [rigorous]

(B4) If 0 ≤ V/(2t) ≤ 1, then E_L ≥ L e₀ − 1682 M / L. [rigorous given cited result X1, with the endpoint V = 2t by continuity, Section 3.6]

(B5) If −1 + δ ≤ V/(2t) ≤ 0 with 0 < δ ≤ 1, then E_L ≥ L e₀ − 2·10¹⁵ δ^(−25/3) M / L. [rigorous given cited result X2]

Cited results. (X1) REPORT §9.100, `codex/work_C142/PROOF.md` Theorem I.9, refereed in §9.104: for 0 ≤ Δ < 1, J = 1, every L ≥ 2 and 0 ≤ n ≤ L/2, |F_L(n) − L e_bulk(Δ, n/L)| ≤ 841/L, where F_L(n) is the minimum over flux of the sector ground energy of the XXZ ring with n reversed spins. (X2) REPORT §9.94, `codex/work_C136/PROOF.md`: the same with 841 replaced by 10¹⁵ δ^(−25/3) for −1 + δ ≤ Δ ≤ 0. Both are AI-produced and AI-refereed, with no human review.

Under the infrared assumption each of (B2), (B4), (B5) gives A ≤ (6/π) K M with the corresponding K, and (B2) gives A ≤ 2t ≤ 2M.

### 3.3 Proof of (B1)

Step 1, fibers. [rigorous] Fix the number m of sites in ℋ₁. For an occupation word s ∈ {0,1}^L with m ones let 𝒦_s = ⊗_j ℋ_{s_j}. The sector is ⊕_s 𝒦_s. The terms V n ⊗ n and (μ/2)(n ⊗ 1 + 1 ⊗ n) act on 𝒦_s as the scalar D(s) = V·#{bonds of type 11} + (μ/2) Σ_bonds (s_j + s_{j+1}). A hop on bond (j, j+1) takes a word s with (s_j, s_{j+1}) = (1, 0) to the word s′ with (0, 1) and acts as −t W_j, where W_j is a unitary from 𝒦_s onto 𝒦_{s′}. The reverse hop is −t W_j†.

Step 2, domination. [rigorous] For ψ = (ψ_s) put f_s = ‖ψ_s‖. Then

    ⟨ψ, H_L ψ⟩ = Σ_s D(s) f_s² − 2t Σ_{edges (s,s′)} Re⟨ψ_{s′}, W ψ_s⟩ ≥ Σ_s D(s) f_s² − 2t Σ_{edges} f_s f_{s′},

by Cauchy–Schwarz and unitarity of W on the fiber. The right side is the expectation in the vector f of the scalar ring with hopping amplitude −t on every bond, the seam included, and ‖f‖ = ‖ψ‖. Hence E_L(W) ≥ E_L(sc). At L = 2 the two oriented bonds give two edges between the words 10 and 01, and the argument is unchanged. The same inequality with W replaced by a phase shows that the scalar ring with hopping −t has the lowest sector energy among all scalar twists, so E_L(sc) restricted to a sector is the flux-minimized energy.

Step 3, the open chain is the scalar open chain tensored with an identity. [rigorous] On the open chain call a move a hop 10 → 01. It lowers the inversion number I(s) = #{i < j : s_i = 1, s_j = 0} by one, and every word reaches the sorted word 0^(L−m) 1^m. Claim: the product U_s of the unitaries W along a sorting path from s depends only on s. Induct on I(s). If two different first moves are available, they sit on bonds (j, j+1) and (k, k+1) that both carry the pattern 10. These bonds cannot share a site, since a shared site would have to be both 0 and 1. So W_j and W_k act on disjoint sites and commute, each move is still available after the other, and both orders reach the same word s″ with the same unitary W_k W_j. By the induction hypothesis the unitary from the word after the first move is path independent, so it may be computed along a path that starts with the other move. Both choices give U_{s″} W_k W_j. This proves the claim. In the basis identified through the U_s, a hop s → s′ at bond j is U_{s′}(−t W_j) U_s⁻¹ = −t, because U_s = U_{s′} W_j. The diagonal terms are scalars on fibers. So the open-chain Hamiltonian in sector m is the scalar open-chain Hamiltonian tensored with the identity on ℋ₀^{⊗(L−m)} ⊗ ℋ₁^{⊗m}.

Step 3′, on the ring the internal register enters only as a boundary twist. [rigorous; Astra gpt-6.2-astra, `astra_A3_adversary.out` Section 1, moved here at the request of the C203 referee, who verified it numerically in 36 cases; I checked the argument by hand] Fix a sector 1 ≤ m ≤ L − 1 and the orientation of the seam hop, which takes a word s with (s_L, s_1) = (1, 0) to the word s′ with (0, 1). In the basis of Step 3 it is represented by Θ_s = U_{s′} W_seam U_s⁻¹, a unitary on the sorted fiber. The source words of the seam hop have fixed occupations on sites L and 1 and m − 1 ones on the interior sites 2, …, L − 1, and any two of them are connected by open-chain hops on the bonds (2,3), …, (L−2, L−1). Such a hop W_j takes s to s₂ and s′ to s₂′, acts on sites disjoint from the seam bond, and so commutes with W_seam. By Step 3, U_s = U_{s₂} W_j and U_{s′} = U_{s₂′} W_j, hence Θ_s = U_{s₂′} W_j W_seam W_j⁻¹ U_{s₂}⁻¹ = Θ_{s₂}. So Θ_s is one unitary Θ = Θ_{L,m}, the same for every source word. In sector m the ring Hamiltonian is therefore the scalar ring with the seam hopping amplitude −t replaced by −t Θ, and diagonalizing Θ splits the sector into scalar rings with seam phases equal to the eigenvalues of Θ, each with the multiplicity of that eigenvalue. The sectors m = 0 and m = L have no hop. Consequently, on the open chain every W, entangling or not, is unitarily equivalent to the flip, and the whole content of W on the ring is the single unitary Θ.

Step 4, bulk densities. [rigorous] The ring Hamiltonian is the open-chain Hamiltonian plus one bond, so the two ground energies differ by at most ‖h‖, for h_W and for h_sc. Dividing by L and using Step 3, e₀(W) = e₀(sc).

Step 5, M ≥ t. [rigorous] Let D = d₀ d₁ and let P₁₀, P₀₁ be the projectors on ℋ₁ ⊗ ℋ₀ and ℋ₀ ⊗ ℋ₁. The states ρ_± = (P₁₀ + P₀₁ ± (W + W†))/(2D) are the normalized projectors on the graphs of ±W. They satisfy tr(h_W ρ_±) = μ/2 ∓ t, since tr[(W + W†)²] = 2D. The partial trace of W over either site vanishes, because the traced site would have to lie in ℋ₀ and in ℋ₁ at once. So the left and right one-site marginals of ρ_± are both (1_{ℋ₁}/d₁ + 1_{ℋ₀}/d₀)/2, and tr[(X ⊗ 1 − 1 ⊗ X) ρ_±] = 0 for every X. For every gauge representative h̃ of h_W, λ_max(h̃) − λ_min(h̃) ≥ tr h̃ (ρ₋ − ρ₊) = 2t. Hence M(h_W) ≥ t. (Astra also shows M(h_W) = M(h_sc); only M ≥ t is used below. Gauge minimization in `exchange_class_check.py` reproduces M(h_W) = M(h_sc) on nine random cases.)

### 3.4 Proof of (B2)

[rigorous] At V = 0 the scalar ring in sector m with hopping −t is a hard-core gas. In ordered particle coordinates, antisymmetrized, it is a free-fermion ring whose boundary phase is (−1)^(m−1), so the allowed momenta are k_j = (2j − m − 1)π/L, j = 1, …, m, in the lowest state, and

    E_L^sc(m) = μ m − 2t Σ_j cos k_j = μ m − 2t sin(π m/L)/sin(π/L).

Since the grids m/L become dense and L sin(π/L) → π, e₀ = min_{0≤x≤1} (μ x − (2t/π) sin π x). Let m minimize E_L^sc(m). Then L e₀ ≤ μ m − (2tL/π) sin(π m/L), and with (B1),

    L e₀ − E_L(W) ≤ L e₀ − E_L^sc(m) ≤ 2t sin(π m/L) (csc(π/L) − L/π) ≤ 2t (csc(π/L) − L/π).

Multiply by L and write x = π/L: the bound is 2π t (x/sin x − 1)/x². By Euler's product, x/sin x = Π_n (1 − x²/π²n²)⁻¹ has nonnegative Taylor coefficients in x², so (x/sin x − 1)/x² is increasing on (0, π/2]. Its value at x = π/2 (L = 2) gives (4 − 8/π) t and its limit at x → 0 gives π t/3. At d = 2, μ = 0, L = 2 one has M = t, e₀ = −2t/π and E₂ = −2t, so the constant is attained. The chemical potential dropped out before any limit, so the bound holds for a field that depends on L, including a field tuned to the z = 2 saturation point.

### 3.5 Proof of (B3)

[rigorous] For V ≤ −2t, h_W − ((V + μ)/2)(n ⊗ 1 + 1 ⊗ n) vanishes on the 00 and 11 sectors, and on the mixed sectors its eigenvalues are −V/2 ± t ≥ 0. Summing over bonds, H_L ≥ (V + μ) N̂ with N̂ the number of sites in ℋ₁, so E_L ≥ L min(0, V + μ). The empty and the full configurations attain the two values. Hence E_L = L e₀ on every ring.

### 3.6 Proof of (B4) and (B5)

[rigorous given X1 or X2] Write J = 2t and Δ = V/(2t). With n ⊗ n = (1/2 + Sᶻ) ⊗ (1/2 + Sᶻ), the scalar sector energy with hopping −t is

    G_L(m) = F_L(m) + (V + μ) m − V L/4,

where F_L(m) is the flux-minimized sector ground energy of J Σ (SˣSˣ + SʸSʸ + Δ SᶻSᶻ) with m spins up (Step 2 identifies hopping −t with the flux minimum, and F_L(m) = F_L(L − m) by a global spin flip). Let K_X be the cited constant (841 for X1, 10¹⁵ δ^(−25/3) for X2), and for rational ν ∈ [0, 1] put φ(ν) = J e_bulk(Δ, min(ν, 1 − ν)) + (V + μ) ν − V/4 and φ_min = inf_ν φ(ν). Then

    E_L(W) ≥ min_m G_L(m) ≥ min_m L φ(m/L) − K_X J/L ≥ L φ_min − K_X J/L.

For the reverse comparison of densities fix a rational ν = p/q and take L a multiple of q. Then E_L(sc) ≤ G_L(νL) ≤ L φ(ν) + K_X J/L, so e₀(sc) ≤ φ(ν) for every rational ν, hence e₀(sc) ≤ φ_min. No continuity of e_bulk in the density is used. With e₀(W) = e₀(sc) and J = 2t ≤ 2M,

    E_L(W) ≥ L e₀(W) − 2 K_X t / L ≥ L e₀(W) − 2 K_X M / L.

The endpoint V = 2t. [rigorous given X1] Theorem I.9 of C142 is stated for 0 ≤ Δ < 1 (C142 reaches Heisenberg only in its doubling corollary, by continuity at fixed L), so the display above is proved for 0 ≤ V < 2t. At fixed L, t, μ and W, E_L(W) is the lowest eigenvalue of a matrix that is affine in V, so it is continuous in V. The density e₀ is an infimum over translation-invariant states of functions affine in V, so it is concave and finite on the whole line, hence continuous. The constant 2 K_X t does not depend on V. Letting V increase to 2t gives the first inequality of the display at V = 2t, and t ≤ M gives the second.

### 3.7 Checks

[numerical, Claude Fable 5.2] `xxc_check.py`: for label-conserving chains with (d₀, d₁) = (1, 2), (2, 2), (1, 3), L = 3, 4, 5 and three parameter sets, the ring spectrum equals the union of twisted XXZ sector spectra to 5·10⁻¹⁰ in both directions. `exchange_class_check.py`: for Haar-random W with (d₀, d₁) = (2, 2), (1, 2), (2, 1), three parameter sets and L up to 6, E_L(W) − E_L(sc) is between 3·10⁻⁶ and 0.1 and never negative, the open-chain spectra agree with the scalar ones to rounding, and M(h_W) = M(h_sc) under gauge minimization.

### 3.8 What this class is and is not

Relation to the literature (added after the C203 referee report). The open-chain identification of Step 3 and the twist of Step 3′ are the squeezed-space, or spin-charge, factorization that is classical for the U = ∞ Hubbard ring (Ogata and Shiba, Phys. Rev. B 41, 2326 (1990), where the ring spectrum is that of spinless fermions with a twist given by the eigenvalues of the cyclic spin permutation) and for multispecies exclusion and XXC chains (Maassarani, Phys. Lett. A 244, 160 (1998), solv-int/9712008, and the multispecies models of Alcaraz and Bariev). The label-conserving case d₀ = 1, d₁ = 2, V = 0 is the U = ∞ Hubbard chain up to a fermion sign, which is itself a twist. The domination step is a diamagnetic inequality. The extension to an entangling W is a modest generalization, since on the open chain every W is equivalent to the flip. The class is therefore not new. What is new here, as far as the referee and I can tell, are the finite-size bounds in units of M: (B2) with its sharp constant, which is an elementary free-fermion computation once the reduction is in place, and (B4), (B5), which are inherited from the project's XXZ theorems. The Ogata–Shiba reference is cited from memory by the referee and was not opened.

The internal spaces can be arbitrarily large and the exchange can entangle them, yet the only propagating sector is the single charge mode of the scalar chain. The internal register enters ring spectra only through a boundary twist. So this is a general-d interacting class with a d-independent constant, and it shows in a controlled example that a large set of internal states does not produce fast sectors. It does not contain any chain with more than one fast sector, so it does not test the log₂ d growth of the conjectured bound.

## 4. Theorem C: the critical ferromagnetic Potts chains (d = 3, 4)

Produced by Astra (gpt-6.2-astra via the shim, `astra_A1_potts.out`), and refereed by a second, independent Astra call (gpt-6.2-astra, `astra_A5_referee_potts.out`), which found it correct as stated, verified the planar composition rule of Section 4.3 in detail, and supplied a second proof of the ground-state identity (C2) that uses neither that rule nor theorem M (Section 4.8). I checked the algebra of Sections 4.2, 4.4 to 4.7 and the coefficient tables of Section 4.8 line by line, read Section 4.3, and tested every conclusion numerically (Section 4.9).

### 4.1 Setting and statement

For Q ≥ 2 let u = Q^(−1/2) Σ_a |a⟩, E = √Q |u⟩⟨u| on one site and F = √Q Σ_a |aa⟩⟨aa| on two sites. Put N = 2L, e_{2j−1} = E_j, e_{2j} = F_{j,j+1} (indices mod N). The chain is H_P^(s)(L) = −s Σ_{i=1}^N e_i with bond h = −s[F + (E ⊗ 1 + 1 ⊗ E)/2]. For Q ≤ 4 write √Q = 2 cos γ, so γ = π/4, π/6, 0 for Q = 2, 3, 4. The sign s = +1 is the self-dual critical ferromagnetic Potts chain (c = 1/2, 4/5, 1).

The twisted XXZ ring on N spins is

    H_X(N; Δ, φ) = (1/2) Σ_{i=1}^{N−1} [σˣσˣ + σʸσʸ + Δ(σᶻσᶻ − 1)]_{i,i+1} + (1/2)[2e^{iφ} σ⁺_N σ⁻_1 + 2e^{−iφ} σ⁻_N σ⁺_1 + Δ(σᶻ_N σᶻ_1 − 1)],

and Spec(N; Δ, φ, m) is its spectrum in the sector Sᶻ = m.

Theorem C. Let Q ∈ {2, 3, 4} and a = cos γ.

(C1) For both signs and all L ≥ 2, every eigenvalue of H_P^(s)(L) lies in the union of Spec(2L; s a, φ, m) over the finite list: m = 0 with φ ∈ {2γ, −2γ, π}, and m = j ∈ {1, …, L} with φ = 2πk/j. All these twists are real. The sector m = 1 can be dropped (the two-through-line layer has zero image, by the referee call). [rigorous given theorem M]

(C2) E_P^(+)(L) = min Spec(2L; a, 2γ, 0) for every L ≥ 2. [rigorous, by Section 4.8, and independently rigorous given theorem M, by Sections 4.3 to 4.5]

(C3) e₀^P = 2 ε₀(a), with ε₀(a) the bulk ground energy per site of the untwisted XXZ ring at Δ = a. [rigorous]

(C4) M(h) = (√Q + 1)/2 for Q = 3, 4, and M(h) = 1 for Q = 2, for either sign. [rigorous]

(C5) E_P^(+)(L) ≥ L e₀^P − 841/L for all L ≥ 2, that is, L (L e₀ − E_L) ≤ 841, which is 841 M at Q = 2 (where M = 1), 615.7 M at Q = 3 and 560.7 M at Q = 4. [rigorous given cited result X1, with Q = 4 (Δ = 1) by continuity, Section 4.6]

(C6) For s = −1 and Q = 2, 3 the analogue of (C3) is false: e₀^(P,−) ≥ 0 while ε₀(−a) ≤ (a − 1)/2 < 0. The inclusion (C1) still holds, but the Potts representation keeps only excited levels of the XXZ sectors, so no finite-size bound centered at the Potts density follows. For s = −1 and Q = 4 the chain is frustration free with E_L = 0. [rigorous]

Theorem M (Morin-Duchesne and Saint-Aubin, J. Phys. A 46, 285207 (2013), arXiv:1203.4996). For even N and 2j through lines there is an intertwiner between the link module of the periodic Temperley–Lieb algebra with contractible loop weight κ² + κ⁻², essential loop weight v^N + v^(−N) and defect winding weight fixed by v, and the Sᶻ = j eigenspace of the twisted XXZ representation whose generator has antiparallel block ((κ², v²), (v⁻², κ⁻²)). It is an isomorphism outside an explicit set of critical parameters, which for fixed κ is a finite set of v. I confirmed from the arXiv abstract that the paper states an intertwiner that is an isomorphism for generic parameters with explicit critical curves. The two Astra calls quote the same items (Definitions 2.6 and 3.1, Proposition 3.2, Theorems 4.5 and 4.6, equations (8), (9)). I did not open the paper to confirm these numbers. Theorem M is needed only for the full inclusion (C1), not for the finite-size bound (C5).

Under the infrared assumption (C5) gives A ≤ (6/π)·615.7 M at Q = 3. The true value is A = c v = (4/5)(3/2) = 1.2 = 0.88 M, so the constant is far from sharp.

Relation to the literature and what is new (added after the C203 referee report). (C1), (C2) and (C3) are not new results. The correspondence between the periodic Potts chain and the XXZ chain with a seam, including the ground-state sector with twist 2γ, is due to Alcaraz, Barber and Batchelor (Phys. Rev. Lett. 58, 771 (1987), Ann. Phys. 182, 280 (1988)), at the level of exact level matching. The finite twist list of (C1) is the standard decomposition of the periodic Potts space into standard modules of the affine Temperley–Lieb algebra (Martin and Saleur, Graham and Lehrer, Richard and Jacobsen, Nucl. Phys. B 750, 250 (2006), Jacobsen, Ribault and Saleur, arXiv:2208.14298, and theorem M). What this section adds is a self-contained rigorous proof of (C2) (Section 4.8), the value of M in (C4), and the finite-size bound (C5). Only (C4) and (C5) are new, and (C5) is conditional on X1. The three-state Potts chain has a neutral, non-Gaussian infrared sector (c = 4/5, no on-site U(1)), and it now has a bound valid at every L, but that bound is inherited from the XXZ host: the Potts ground state is one level of a U(1) Gaussian ring, and it is the host that is bounded. This is not a new method for neutral sectors. It gives nothing for a neutral sector without such a host, which is the case the project needs (Section 7.6).

### 4.2 Temperley–Lieb relations

[rigorous] E² = √Q E and F² = √Q F. With P = |u⟩⟨u| and D = Σ_a |aa⟩⟨aa|, (P ⊗ 1) D (P ⊗ 1) = Q⁻¹ (P ⊗ 1), because P|a⟩⟨a|P = P/Q, and D (P ⊗ 1) D = Q⁻¹ D, because ⟨aa|(P ⊗ 1)|bb⟩ = δ_ab/Q. Multiplying by Q^(3/2) gives e_i e_{i±1} e_i = e_i. Generators at cyclic distance at least two commute: site terms on different sites, bond terms (all diagonal), and a site term with a bond not touching it. Summing the bond h over the ring counts each site term twice with weight 1/2, so H_P = −s Σ e_i.

### 4.3 Which twists occur

[rigorous; Astra, with the composition rule verified in `astra_A5_referee_potts.out` Section 1 by an Euler-characteristic count (gluing along g gap intervals lowers Σχ by g and the gap count by 2g, so the exponent χ − b/2 is additive), and attributed there to Jones's graph planar algebra (Theorem 3.1 of his paper on the planar algebra of a bipartite graph) and to Kodiyalam, Saini, Sruthymurali and Sunder, arXiv:1901.11191, Proposition 13. I read both accounts and checked item (i) by hand] Let G be the star graph with center o and Q leaves, with adjacency matrix A_G and Perron vector μ_o = √Q, μ_a = 1. Let V_n be spanned by closed walks of length n on G with a marked start. Annular Temperley–Lieb diagrams act between the V_n by summing over region colourings with the weights of equation (1) of the Astra output. Two facts are used.

(i) On V_N the generator acts as e_i|…, a, b, c, …⟩ = 1_{a=c} Σ_{b′∼a} (√(μ_b μ_{b′})/μ_a) |…, a, b′, a, …⟩. On walks with leaves at odd positions and the center at even positions this is the Potts chain: at an odd position a = c = o and every entry is 1/√Q, which is E, and at an even position the entry is √Q if the two neighbouring leaves agree and zero otherwise, which is F. [rigorous, checked by hand]

(ii) Filter V_N by the number r of through lines. The quotient at r > 0 is spanned by images of the standard modules W_{r/2,z}, where z is an eigenvalue of the cyclic rotation of walk labels on V_r, so z^r = 1. The part with r = 0 is spanned by images of W_{0,z} with z + z⁻¹ an eigenvalue of A_G, and the eigenvalues of the star are √Q, −√Q, 0. For Q ≤ 4 all of them lie in [−2, 2], so z is unimodular. The Perron value √Q gives z = e^{±iγ}.

Since characteristic polynomials multiply along a filtration and divide under surjections, every eigenvalue of B = Σ e_i on the Potts chain is an eigenvalue of B on one of these standard modules.

### 4.4 From standard modules to twisted XXZ

[rigorous given theorem M] Fix κ² = e^{iγ}. Theorem M gives similarity of B on the link module and on the spin module outside finitely many v. The coefficients of both characteristic polynomials are Laurent polynomials in v, so the characteristic polynomials agree for every v, including the exceptional twists. No invertible intertwiner at a root of unity is needed.

[rigorous, checked by hand] Conjugating the spin generators by v^{Σ_ℓ ℓ σᶻ_ℓ} makes every interior hopping coefficient 1 and the seam coefficient v^{2N} = z² = e^{iφ}. The diagonal part of each generator is (a/2)(1 − σᶻσᶻ) + ((κ² − κ⁻²)/4)(σᶻ_i − σᶻ_{i+1}), and the second term sums to zero on the ring. So Σ ē_i is similar to H_X(N; −a, φ). The product of σᶻ over odd sites flips the sign of every hopping and nothing else, so −H_X(N; −a, φ) is unitarily equivalent to H_X(N; a, φ), and N = 2L even makes this consistent at the seam. With H_P^(+) = −B this gives (C1) for s = +1 at Δ = a and for s = −1 at Δ = −a, with the dictionary m = r/2 and e^{iφ} = z².

### 4.5 The ground state (C2)

[rigorous, checked by hand] Take r = 0 and the Perron vector, so that every loop, contractible or not, has weight √Q. The diagram map C from W_{0,z} to the Potts space intertwines B and has nonnegative, nonzero columns. The Potts matrix B_P is symmetric, entrywise nonnegative and irreducible (single-site colour changes connect all configurations), so its top eigenvector x is strictly positive and its top eigenvalue is λ. Then y = Cᵀx > 0 and yᵀ B_W = λ yᵀ. Since B_W is entrywise nonnegative, Σ_i y_i |(B_W f)_i| ≤ λ Σ_j y_j |f_j| for every f, so every eigenvalue of B_W has modulus at most λ, and λ is an eigenvalue of B_W. So the top eigenvalue of B on the link module equals the top eigenvalue on the Potts chain, even if C has a kernel. With Section 4.4 this is (C2).

### 4.6 The bulk density (C3) and the finite-size bound (C5)

[rigorous] For any nearest-neighbour chain E_n/n converges to the infimum of the energy density over translation-invariant states (restrict to n sites, and tile with ring ground states, each at a cost of one bond). Let g_n be the unrestricted untwisted XXZ ring ground energy and g_n^φ(0) the minimum in the sector m = 0 with twist φ. A twist changes two matrix entries, so g_n^φ(0) ≥ g_n − 2. Conversely, place an untwisted ground vector of magnetization m and its global spin flip on consecutive blocks of a 2n-ring. The product has m = 0, and two bonds are deleted and two inserted, each of norm at most 1 + |a|, so g_{2n}^φ(0) ≤ 2 g_n + 4(1 + |a|). Hence g_{2n}^φ(0)/(2n) → ε₀(a), and with (C2), e₀^P = 2 ε₀(a).

[rigorous given X1] In the convention of X1 (J = 1, S = σ/2), H_X = 2 Σ(SˣSˣ + SʸSʸ + Δ SᶻSᶻ) − NΔ/2 up to the twist. For every even N, every sector and every real twist, min Spec(N; a, φ, m) ≥ 2 F_N(n) − N a/2 ≥ N [2 e_bulk(a, n/N) − a/2] − 1682/N with n = N/2 − |m|. As in Section 3.6, ε₀(a) ≤ 2 e_bulk(a, ν) − a/2 for every rational ν, by X1 applied along lengths with νN an integer, since a twist changes only the two seam entries and so costs at most 2 in energy. Hence min Spec(N; a, φ, m) ≥ N ε₀(a) − 1682/N. With N = 2L, (C2) and (C3) give E_P^(+)(L) ≥ L e₀^P − 841/L. X1 covers 0 ≤ Δ < 1, which contains a = cos γ for Q = 2, 3.

The case Q = 4, a = 1. [rigorous given X1] For Δ < 1 the argument above gives min Spec(N; Δ, φ, m) ≥ N ε₀(Δ) − 1682/N for every even N, sector and real twist. At fixed N, φ, m the left side is the lowest eigenvalue of a matrix affine in Δ, hence continuous, and ε₀(Δ) is an infimum of functions affine in Δ, hence concave, finite and continuous. Letting Δ increase to 1 gives the same inequality at Δ = 1, and (C2), (C3) at Q = 4 then give (C5).

### 4.7 The half-spread (C4) and the reversed sign (C6)

[rigorous, checked by hand] The bond is symmetric under exchange of the two sites and the gauge term X ⊗ 1 − 1 ⊗ X is antisymmetric. The spread is convex and unitarily invariant, so spr(h) ≤ [spr(h + G) + spr(h − G)]/2 = spr(h + G), and X = 0 is optimal. Let C₀ = F + (E ⊗ 1 + 1 ⊗ E)/2, U₀ = u ⊗ u and D₀ = Q^(−1/2) Σ_a |aa⟩. Then C₀ U₀ = √Q U₀ + D₀ and C₀ D₀ = U₀ + √Q D₀, so U₀ + D₀ is a positive eigenvector of the nonnegative matrix C₀ with eigenvalue √Q + 1, which is therefore the largest. C₀ ≥ 0, and for Q ≥ 3 the vector v ⊗ w − w ⊗ v with v, w independent in u^⊥ is annihilated by C₀. So the spectrum of C₀ spans [0, √Q + 1] and M = (√Q + 1)/2 for Q = 3, 4. At Q = 2 the space u^⊥ is one-dimensional and the antisymmetric vector does not exist. There E = (1 + X)/√2 and F = (1 + Z ⊗ Z)/√2, so C₀ = [2 + Z ⊗ Z + (X ⊗ 1 + 1 ⊗ X)/2]/√2. The operator Z ⊗ Z + (X ⊗ 1 + 1 ⊗ X)/2 has eigenvalues ±√2 on the sector X ⊗ X = 1 and ±1 on the sector X ⊗ X = −1, so C₀ has spectrum [√2 − 1, √2 + 1] and M = 1 (the C203 referee confirms M = 1 by semidefinite programming).

[rigorous] For s = −1, H_P = Σ e_i ≥ 0, so e₀ ≥ 0, while the alternating σˣ product state gives ε₀(−a) ≤ (a − 1)/2 < 0 for Q = 2, 3. At Q = 4 the vector ψ(s₁, …, s_L) = Π_j w_{s_j}(s_j − s_{j+1}) with w = (1, −3, 3, −1) is annihilated by every F (it vanishes unless neighbours differ) and by every E (Σ_t w_t t^k = 0 for k = 0, 1, 2), so E_L = 0.

### 4.8 A second proof of (C2), without the walk representation and without theorem M

[rigorous; Astra gpt-6.2-astra, `astra_A5_referee_potts.out` Section 5; I checked the coefficient tables by hand]

Step 1, from loop states to the Potts chain. Let W be the zero-through-line module of the periodic Temperley–Lieb algebra on N = 2L points with weight δ = √Q for contractible and for non-contractible loops. Its basis is the set of noncrossing pairings drawn on an annulus. Forgetting the hole maps it to noncrossing pairings of a disk and intertwines the action, since both kinds of loop have the same weight. Colour the regions of a disk pairing black and white so that the odd gaps, which are the Potts sites, are black. The pairing defines a partition π of the L sites, two sites being in one block when their gaps touch the same black region. Send the pairing to

    f_π(s₁, …, s_L) = δ^(−|π|) Π_{blocks A of π} 1{s is constant on A}.

On pairings, e_{2j} joins the black regions of sites j and j + 1, or removes a loop if they already coincide, and e_{2j−1} cuts site j off its black region as a singleton, or removes a loop if it already is one. The four cases match. F_{j,j+1} f_π = δ·1{s_j = s_{j+1}} f_π is f of the partition with the two blocks joined when they differ (one more constraint, one block fewer, δ^(−|π|)·δ = δ^(−(|π|−1))), and is δ f_π when they coincide. (E_j f_π)(s) = δ⁻¹ Σ_t f_π(s with s_j = t) is f of the partition with j detached when j shares its block (one value of t survives, and the prefactor becomes δ^(−(|π|+1))), and is (Q/δ) f_π = δ f_π when j is a singleton. So this map C from W to the Potts space intertwines B = Σ e_i, and each of its columns is entrywise nonnegative and positive on the constant configurations.

Step 2. The Perron–Frobenius argument of Section 4.5 applies to C verbatim and gives λ_max(B on W) = λ_max(B on the Potts space).

Step 3, from loop states to spins. Fix κ² = e^{iγ} and let v ≠ 0. The spin generator ē_i vanishes on parallel spins and has block ((κ², v²), (v⁻², κ⁻²)) on (↑↓, ↓↑), for every cyclic bond. It factors as the column (κv, (κv)⁻¹)ᵀ times the row (κ/v, v/κ). For an annular pairing, lift to the strip and orient each arc from p to q with 1 ≤ q − p < N, the disk side between them. To an arc of length ℓ = q − p assign A_ℓ |↑_p ↓_q⟩ + A_ℓ⁻¹ |↓_p ↑_q⟩ with A_ℓ = κ v^ℓ, and let J(v) send the pairing to the product over its arcs. Applying the row of ē_i at sites x, y = x + 1: on a short arc x → y it gives (κ/v)(κv) + (v/κ)(κv)⁻¹ = κ² + κ⁻² = δ. On the long arc y → x + N it gives (κ/v) A_{N−1}⁻¹ + (v/κ) A_{N−1} = v^(−N) + v^N. If x and y lie on different arcs, say x → k and y → j with x < y < j < k, the remaining endpoints j, k carry (κ/v) A_{k−x}/A_{j−y} = κ v^(k−j) = A_{k−j} for (↑_j ↓_k) and its reciprocal for the other assignment, which is the vector of the arc j → k, and the two other nestings (j → x, y → k and j → y, k → x) give A_{k−j} in the same way. The column of ē_i then creates the short arc x → y. So J(v) intertwines the loop action with non-contractible weight α = v^N + v^(−N) and the spin action in Sᶻ = 0, the seam included.

Annular pairings correspond bijectively to L-element subsets of the N positions (the closing ends of the arcs, by parenthesis matching), like the spin configurations with Sᶻ = 0. Dividing each column by Π A_ℓ and sending v → +∞ suppresses every reversed arc by A_ℓ⁻², so the rescaled J(v) tends to a permutation matrix and J(v) is invertible for large v. Hence B on W and Σ ē_i on Sᶻ = 0 have the same characteristic polynomial for large v, and since the coefficients are Laurent polynomials in v, for every v ≠ 0. Choose v on the unit circle with v^N = e^{iγ}. Then α = δ, and the gauge transformation and the odd-site rotation of Section 4.4 turn −Σ ē_i into H_X(2L; a, φ) with e^{iφ} = v^{2N} = e^{2iγ}.

Steps 1 to 3 give E_P^(+)(L) = −λ_max(B on Potts) = −λ_max(B on W) = min Spec(2L; a, 2γ, 0). The case Q = 4, γ = 0 is included.

### 4.9 Checks

[numerical, Claude Fable 5.2] `potts_xxz_inclusion.py`: for Q = 2, 3, 4, both signs, L up to 5, 4, 4, every Potts level matches a twisted XXZ level on a twist grid to 5·10⁻¹⁰ (rounding). `potts_table8_check.py`: the same with the twists restricted to the finite list of (C1), L up to 5 for Q = 2, 3 and 4 for Q = 4, same accuracy, and the equality (C2) for Q = 3 at L = 6, 7, Q = 4 at L = 5, 6, Q = 2 at L = 7, to 7·10⁻¹⁴. `potts_casimir.py`: M = 1.366025 and 1.5 for Q = 3, 4, and L (L e₀ − E_L)/M falls from 0.5005 at L = 2 to 0.4614 at L = 11 for Q = 3 (the limit π c v/(6M) with c = 4/5, v = 3/2 is 0.4599), from 0.606 to 0.553 for Q = 4 at L ≤ 9, and from 0.399 to 0.371 for Q = 2 at L ≤ 14. For s = −1 and Q = 2, 3 the same quantity is negative and grows in magnitude linearly in L, which is the failure (C6) and which I found numerically before reading Astra's proof of it.

## 5. Stacks

Proposition D. [rigorous, Claude Fable 5.2] Let h = h₁ ⊗ 1 + 1 ⊗ h₂ be a decoupled stack of two chains with Lemma P constants K₁ and K₂. Then E_L = E_L¹ + E_L², e₀ = e₀¹ + e₀², and L (L e₀ − E_L) ≤ K₁ M₁ + K₂ M₂ ≤ (K₁ + K₂) M(h), because M(h) ≥ max(M₁, M₂) (compress any gauge representative of h with a product vector of the other layer). By induction a stack of m layers drawn from class B, the three chains of Theorem C, quasi-free chains and the d = 2 U(1) class has L (L e₀ − E_L) ≤ K_max Σ_i M_i ≤ K_max m M(h) ≤ K_max log₂ d · M(h), since each layer has d_i ≥ 2. This is the conjectured log₂ d form, on a class where it carries no information beyond the single layers. When the extreme eigenvalues of every layer have eigenstates with equal left and right marginals (Heisenberg layers, for example), product witnesses give M(h) = Σ_i M_i and the constant is K_max, independent of d.

Proof of M(h) ≥ M₁. Let h̃ = h + X ⊗ 1 − 1 ⊗ X with X on the two-layer site. For a fixed unit vector χ of layer 2 on each site, the compression of h̃ to (ℂ^{d₁} ⊗ χ)^{⊗2} is h₁ + ⟨χχ|h₂|χχ⟩ + X_χ ⊗ 1 − 1 ⊗ X_χ, a gauge representative of h₁ up to a constant. A compression does not increase the spectral spread. ∎

## 6. Two no-go lemmas for the lane picture

Both are from Astra (gpt-6.2-astra, `astra_A3_adversary.out` Sections 2 and 3), and I checked both proofs by hand. They are statements about effective single-particle data, so their use for a microscopic chain is physics-level.

Lemma E (velocity tail from a bounded hopping moment). [rigorous] Let f be real, 2π-periodic and r times differentiable with r ≥ 3, and let p = 1/(r − 1). Then the sum of |f′(x)|^p over the zeros x of f in one period is at most 2π (‖f^(r)‖_∞/r)^p.

Proof. It suffices to bound the sum over every finite set of distinct zeros, since the bound is uniform in the set. Order such a set cyclically as x₁ < … < x_n and extend by x_{i+n} = x_i + 2π (the lifted nodes are distinct on the line even when n < r). Fix i and choose b so that g(x) = f(x) − b Π_{j=0}^{r−1}(x − x_{i+j}) has g′(x_i) = 0, that is, f′(x_i) = b Π_{j=1}^{r−1}(x_i − x_{i+j}). Then g has r + 1 zeros counted with multiplicity, so g^(r) vanishes somewhere, which gives |b| ≤ ‖f^(r)‖_∞/r!. Hence |f′(x_i)| ≤ (‖f^(r)‖_∞/r!) Π_{j=1}^{r−1}(x_{i+j} − x_i). Raise to the power p, sum over i and apply Hölder with r − 1 equal exponents. Since Σ_i (x_{i+j} − x_i) = 2π j, the sum is at most (‖f^(r)‖_∞/r!)^p Π_{j=1}^{r−1}(2π j)^p = 2π (‖f^(r)‖_∞/r)^p. ∎

The same proof works for r = 2, where p = 1 and no Hölder step is needed. It gives Σ|f′| ≤ π ‖f″‖_∞ over the zeros, which is a bound on A with no gain in the tail exponent, so the statement is of interest for r ≥ 3.

Consequence. [physics-level] If the gapless modes of a family of chains are the Fermi points of at most s(d) effective bands f_j(k) = Σ_ℓ t_{j,ℓ} e^{iℓk}, each carrying central charge at most γ per crossing, and Σ_ℓ |ℓ|^r |t_{j,ℓ}| ≤ D(d) M with the range ℓ measured in lattice sites, then Σ_{a: v_a ≥ u} c_a ≤ 2π γ s(d) (D(d)/r)^p (M/u)^p with p = 1/(r − 1) < 1. With the Lieb–Robinson bound v_a ≤ 4eM this is the tail estimate that REPORT §9.93 names as the missing lemma, and it gives A ≤ 2π γ s (D/r)^p (4e)^(1−p) M. The sawtooth of N teeth and slope M has third moment of order M N², so the cost of a fast tower in this picture is a hopping moment that grows with the number of lanes. The exponent is sharp: f_n(k) = M n^(−r) sin(nk) satisfies the hypothesis with D = 1 and has 2n zeros of speed M n^(1−r).

Lemma F (rank across a cell boundary). [rigorous] Let B(k) = B₀ + T e^{ik} + T† e^{−ik} be a Bloch matrix of any size with rank T ≤ R. Then for every energy there are at most 2R linearly dispersing crossings, counted with multiplicity.

Proof. In the column expansion of det(B(k) − E) a term that takes more than R columns from zT, z = e^{ik}, vanishes, and likewise for z⁻¹T†. So z^R det(B − E) is a polynomial in z of degree at most 2R, not identically zero if E is not in the spectrum of B(0). A momentum with b crossing branches is a zero of multiplicity at least b. Flat bands at E are avoided by shifting E slightly, which does not remove any branch of nonzero slope. ∎

Consequence. [physics-level] A superlattice of any period whose inter-cell hopping has rank at most R(d) carries at most 2R(d) fast crossings, and A ≤ 8eγ R(d) M. This extends the free-fermion channel count of REPORT §9.71 to arbitrary cell size. It does not apply to an interacting superlattice, where a single physical bond does not bound the rank of a dressed hopping matrix.

## 7. Where each strengthening stops

### 7.1 Permutation-invariant ring states are harmless

Proposition G. [rigorous; Astra gpt-6.2-astra, `astra_A2_thirdlaw.out` Section 5 and `astra_A4_referee_T1.out` Section 3; checked by hand] Every permutation-invariant state on (ℂ^d)^{⊗L}, pure or mixed, has a two-site marginal within 30 d⁵/L² of Ω_∞ in trace norm.

Proof. For L ≥ 3 let R = ⌊L/3⌋ and take three disjoint groups U, V, W of R sites. On each site use the measurement with effects (1 ± G_α)/(2d²) and outcome matrices X = ± d² G_α, where G_α runs over a Hilbert–Schmidt orthonormal basis of Hermitian d × d matrices. The outcome matrix has mean equal to the state of the site, also conditionally on outcomes elsewhere, and ‖X‖₂² = d⁴. By permutation invariance a = 𝔼‖X_i‖₂² = d⁴ and c = 𝔼 tr(X_i X_j) is the same for all i ≠ j, with |c| ≤ a. With X̄ the mean over a group, expanding the square gives 𝔼‖X̄_U − X̄_V‖₂² = 2(a − c)/R ≤ 4d⁴/R, and the same for W and V. Let v be the record on V and γ(v) the conditional state of any unmeasured site. Then ρ₂ − Σ_v p_v γ(v)^{⊗2} = 𝔼[(X̄_U − γ(v)) ⊗ (X̄_W − γ(v))], the conditional expectation γ(v) = 𝔼[X̄_U | v] has mean-square error at most that of the estimator X̄_V, and ‖Y‖₁ ≤ √d ‖Y‖₂ with Cauchy–Schwarz give ‖ρ₂ − Σ_v p_v γ(v)^{⊗2}‖₁ ≤ d · 4d⁴/R = 4d⁵/R. Keep the L − R sites outside V, take tensor powers of their conditional state with a common record v, mix over v and average over translations. Internal bonds have marginal ρ₂, and one bond in L − R has marginal Σ_v p_v γ(v)^{⊗2}. So the distance is at most 4d⁵/(R(L − R)) ≤ 30d⁵/L², using R ≥ L/5 and L − R ≥ 2L/3. At L = 2 the diameter 2 is below 30d⁵/4. ∎

The W state of ring momentum p = 2πj/L is covered as follows. Take the infinite-chain state that approximates the zero-momentum W state, average it over global number rotations (this fixes the target marginal and contracts the error, and makes the state invariant under global number rotations), and then apply the site-dependent rotation e^{ipx n̂}. A translation changes this rotation pattern by a global number rotation, under which the twirled state is invariant, so the rotated state is still translation invariant, and its bond marginal and the target are conjugated by the same two-site unitary.

So Dicke states, W states of any ring momentum and ground states of mean-field type are not candidates for a ring marginal at distance of order 1/L from Ω_∞. This answers ingredient (I4) of my spec in the negative.

### 7.2 Inclusion of spectra is not enough: the reversed Potts sign

For H = +Σ e_i at Q = 2, 3 the spectral inclusion (C1) holds, and Lemma P holds for the host XXZ rings at Δ = −a, yet no bound follows for the Potts chain, because the Potts representation contains only excited levels of the host sectors and its bulk density differs from the host's ((C6), and L (L e₀^host − E_L) grows linearly in L in `potts_casimir.py`). The transfer of Lemma P by spectral reduction needs two things: every level of the reduced chain lies in a host sector that obeys the bound, and the reduced chain has the host's bulk density. In Theorems B and C the second is supplied by a Perron–Frobenius or domination argument, which is available only for one sign.

### 7.3 The square-root barrier, sharpened

[rigorous for the general statement, physics-level for the conformal estimate; Claude Fable 5.2] Every construction of the type of Section 2 pays, per join and per L sites, 2M times a trace-norm error, and the trace norm is controlled only by the square root of an information quantity. Two such quantities enter: the conditional mutual information I(A:C|B) = 2(u_b − u_{b′}) between the sites next to a buffer of b sites and the inserted interval, and the error z in decorrelating the two ends that meet at a join. The general theorem pays a tower of exponentials because both are made small by pigeonhole over scales and the second needs a window dimension q = d^(b+1).

For a conformal ring state one can estimate both directly. The entropy increments are u_j ≈ C/(3j), so the recovery term is of order √(C/b). Two buffers of b sites at separation s have mutual information of order (b/s)^(4x), with x the smallest scaling dimension that couples them, so the decorrelation term 2√z is of order (b/s)^x. With s at most L, the best buffer is b ~ L^(2x/(2x+1)) and the join error is of order L^(−x/(2x+1)), below L^(−1/2) for every x. Even if energy errors were second order in the amplitude errors (linear in the information, which is true for local channels that nearly fix a ground state, by the identity ΔE = Σ_i ⟨χ_i|(H − E)|χ_i⟩ with K_i ψ = c_i ψ + χ_i, but is not available across a change of ring length), the exponent would be 2x/(2x+1) < 1. Reaching 1/L needs buffers of order L, where the two ends are correlated at order one and must be reconnected, not decorrelated. That is the original problem. And a linear-response version with buffers of order L would give at best L e₀ − E_L ≲ M C/L with C the total central charge, that is, A ≲ M Σ_a c_a, which the Lieb–Robinson bound already gives and which is not a bound in terms of d.

A separate loss is exhibited in the referee call for Theorem A (Section 2.9, read and not verified independently): for certain ring states at d = 4 the decoupling of two full end windows fails on every span polynomial in the buffer size, so the full-window decoupling step (2.2) cannot be carried out on polynomial scales. This does not exclude a join that decouples only the end sites, or controls only the observables to which the two connecting bonds pull back, or reconnects the ends, and the example gives no lower bound on the energy defect.

Conclusion: recovery-map periodization, including the shared-record version, cannot reach Lemma P. It does reach the uniform o(1) statement, which the earlier version of this obstruction (REPORT §9.90, and (O3) in my spec, where I had asserted that pigeonhole gives only the trivial order) had missed. The step I had missed is to retain the long interval and insert a short one, so that the period stays of order L whatever scale the pigeonhole selects.

### 7.4 The finite-size signature of a tower

[rigorous as a statement about sequences; Astra, checked] A family with deficits L e₀ − E_L = M N/(L + N²) has deficit at most M/(2√L) at every length for every N, and lim_L L (L e₀ − E_L) = M N. So a tower need not show an order-M plateau at any length. My spec's necessary condition (N4), and the heuristic of REPORT §9.93 that a tower with N sectors of crossover length N shows L e₀ − E_L of order M up to L ~ N, describe only the extreme scenario. Adversarial searches that look for a visible order-M deficit at accessible sizes can miss a tower. Theorem A excludes the extreme scenario and nothing more.

### 7.5 Site insertion

[rigorous; Astra, checked] If E_{j+1} ≤ (j+1) E_j/j + K M/j² for all j ≥ L₀, then dividing by j + 1 and summing from L to infinity, with E_j/j → e₀, gives e₀ − E_L/L ≤ K M Σ_{j≥L} 1/(j²(j+1)) ≤ K M/L², which is Lemma P on L ≥ L₀. (I had claimed in my spec that an additive constant stays undetermined. That was wrong.) The CFT value of the insertion excess is πA(2L+1)/(6L²(L+1)), so the hypothesis is of the right order. As stated it is false: for the classical bond that penalizes equal neighbouring bits, E_L = 0 on even L and E_{L+1} = 1. A correct form has to be restricted to compatible lengths, or to insertion of a full period.

### 7.6 Charged and neutral sectors

[rigorous given the stiffness bound of REPORT §9.75; Astra] In a normalization with Σ_α D_α = Σ_a w_a v_a for unit-central-charge modes with Luttinger weights w_a, the stiffness bound and w_a ≥ w_min > 0 give A ≤ 4M/w_min. The response data w_a = 2^(−a), v_a = M satisfy the stiffness bound for every N with A = N M. No bond realizing them is known. The stiffness route therefore still needs a lower bound on Luttinger weights, as recorded in §9.77.

[physics-level; Astra] Neutral interacting factors, simultaneously occupied dressed strings and higher-order Schrieffer–Wolff constructions remain open. For each, finite perturbative order for one member of a family gives neither a uniform number of dressed branches nor the hopping-moment hypothesis of Lemma E in lattice units.

### 7.7 The remaining estimate

The general target is unchanged: g_d(L) ≤ K(d)/L. In the language of Section 2 it is a join whose trace-norm error at the two disturbed bonds is O(1/L) per L sites, for the ring ground state. No route in this round, and none in REPORT §9.68 to §9.98, supplies it.

## 8. Checks run (all scripts in this folder, Claude Fable 5.2)

| Script | What it tests | Result |
|---|---|---|
| `potts_xxz_inclusion.py` | Potts spectrum inside twisted XXZ sectors on a twist grid, Q = 2, 3, 4, both signs, L ≤ 5, 4, 4 | all levels matched to 5·10⁻¹⁰ |
| `potts_table8_check.py` | the same with only the predicted finite twist list, and E_P = min Spec(2L; a, 2γ, 0) to L = 7 | matched, 7·10⁻¹⁴ for the ground state |
| `potts_casimir.py` | M of the Potts bond, e₀ from the Bethe integral, L(L e₀ − E_L)/M | 0.5005 → 0.4614 (Q = 3, L = 2 to 11), 0.606 → 0.553 (Q = 4, L ≤ 9), 0.399 → 0.371 (Q = 2, L ≤ 14), far below the proved 615.7 and 560.7, and the reversed sign fails as in (C6) |
| `xxc_check.py` | label-conserving chains, d = 3, 4, against twisted XXZ sectors | equal spectra to 5·10⁻¹⁰ |
| `exchange_class_check.py` | Haar-random exchange unitaries, d = 3, 4, L ≤ 6 | E_L(W) ≥ E_L(sc) in all 24 cases, open spectra equal, M equal |
| `t1_lemma_checks.py` | marginal-correction lemma, covariance identities of the first proof of Theorem A | bound holds (largest ratio 0.755), identities hold in 12 cases |
| `t1_chainrule_check.py` | decoupling estimate (2.2) with q = 2 on four ring states | holds in 16 cases, ratios 0.19 to 0.49 |
| `b2_bound_check.py` | the sharp V = 0 bound (B2) with random W and eight fields, and the SU(3) swap chain as a control | largest ratio 1.000000 (saturated at d = 2, L = 2), SU(3) gives 2.67, 2.32, 2.25 M at L = 3, 6, 9, above the class constant 1.45 M |

Not tested numerically: Theorem A itself (its lengths are out of reach), the cited XXZ theorems X1 and X2 (taken from REPORT §9.94 and §9.100 with their own checks), Lemmas E and F beyond the worked example in the proof of Lemma E. The chains named in my brief are covered as follows: XX, Heisenberg and XXZ are the scalar members of the exchange class and the hosts of the Potts chains, the transverse Ising chain is the Q = 2 Potts chain, and the SU(3) swap chain is outside every class treated here (a three-species exchange skeleton, for which no finite-size theorem exists), so Theorems B and C say nothing about it.
