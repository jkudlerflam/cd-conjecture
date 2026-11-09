This was written by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C178 derivations: how to read the ring spectrum of the cusp model

Model (3.2) of C173 on a ring of L sites, hard-core bosons (the periodic spin chain), no two adjacent, units |t_2| = 1, w = |w| e^{iφ}, uniform flux θ (each hop by d carries e^{iθd/L}). N is conserved, translation T conserves the momentum K = 2πm/L.

## 1. Symmetries and the flux line

(1a) [proved] At φ = π/2, L even, S = I · (−1)^{Σ_{j odd} n_j} (I the site inversion j → −j) commutes with H at θ = 0 and sends K → Nπ − K. Proof: I maps w b†_{j+1}b_j to w b†_{j−1}b_j, the staggered sign flips the sign of every nearest-neighbour hop, and w̄ → −w̄ = w when w = i|w|. The NNN hop, V_2, V_3 and the blockade are invariant. Translation by one site changes the staggered sign by (−1)^N, which gives K → −K + Nπ. Checked numerically (validate.out: [H,S] = 0 to 3e-14, S² = 1).

(1b) [proved] With flux, S sends θ → −θ. The large gauge transformation U = exp(2πi Σ_j j n_j/L) satisfies U H(θ) U† = H(θ + 2π) and T U = e^{−2πiN/L} U T, so the spectrum F(m, θ) in sector m obeys

    F(m, θ + 2π) = F(m + N, θ).                                                        (1.1)

Combining, F(m, θ) = F(NL/2 − m − N, 2π − θ) for 0 < θ < 2π and F(m, 0) = F(NL/2 − m, 0). Both checked numerically (validate.out, validate2.out). Under S the kinetic valley label C1 = Re e^{iθ/L}⟨Σ_j b†_{j+1}b_j⟩ = Σ_k cos(k) n_k changes sign, since S sends k → π − k.

(1c) Consequence of (1.1) [proved]. Following a level continuously in θ over the whole real line visits the sectors m, m + N, m + 2N, ... So the minimum of a branch over all fluxes depends only on κ = m mod N (L is a multiple of N in every run here). This "twist-relaxed" energy removes the charge current zero mode exactly, whatever the interactions are.

## 2. Spin-charge factorisation in the impenetrable regime

[derived, leading order in the dilute impenetrable limit of C173 PROOF 4] When inter-valley pairs are impenetrable, the order of the particles along the ring is fixed up to pseudospin exchange processes. The eigenstates factorise in the Ogata-Shiba form

    Ψ(x_1 < ... < x_N; σ_1 ... σ_N) = det[e^{i k_α x_l}] · Π_l e^{i k_{σ_l} x_l} · χ(σ_1, ..., σ_N),

with k_a ≈ −s, k_b ≈ π + s the valley centres, a Slater determinant of small relative momenta (hard rods, Girardeau), and χ an eigenstate of the pseudospin chain on the N particles. Moving the last particle once around the ring cyclically shifts χ by one pseudospin site and multiplies by e^{i k_σ L}. The uniform part e^{i k_a L}, the cyclic shift eigenvalue e^{iq} (q = 2πκ/N the pseudospin chain momentum) and the flux θ together fix the boundary condition of the charge determinant. The rest, e^{i(k_b − k_a)L} per flipped pseudospin, is a residual pseudospin twist Φ_s = (k_b − k_a)L mod 2π in the transverse coupling. Hence:

- the energy is E_c(θ + q + k_a L + π(N−1)) + E_s(q, S^z; Φ_s) + O(charge-pseudospin coupling), with S^z = M/2;
- minimising over θ removes E_c's zero mode, so the twist-relaxed levels are the pseudospin chain levels E_s(q, S^z) up to a constant;
- the total momentum of an M = 0 level relative to the ground state is ΔP = 2πκ/L mod 2πN/L, so κ = (m − m_0) mod N is the pseudospin chain momentum q = 2πκ/N. The Néel (valley-crystal) partner of the ground state sits at κ = N/2.

This is the "compactification lattice" of the two bosons in the form that a ring can test: charge quantum numbers (θ-branch, ΔN) on one side, pseudospin (κ, M) on the other, with the residual cross-term Φ_s set by the valley momenta.

Consistency check [numerical, local/n6_pi2.jsonl]. At L = 36, N = 6, θ = 0 the ground state is a doublet at ΔP = ±πn (κ = 3, Marshall sign of an antiferromagnetic pseudospin chain with N/2 odd), 0.011 below the κ = 0 state, while at θ = π the κ = 3 state is lower by 0.115 and alone. The factorisation predicts E(0, π) − E(π, 0) = E(π, π) − E(0, 0) = E_s(π) − E_s(0). The data give −0.135 and −0.127. At L = 24, N = 4 they give +0.184 and +0.175. So the factorisation holds to about 5 percent of the pseudospin splitting at these sizes, and E_s(π) − E_s(0) alternates in sign with N/2, as in an antiferromagnetic (or XY) pseudospin ring.

## 3. What each phase predicts for the relaxed spectrum

Write the pseudospin chain as an XXZ ring, J > 0, anisotropy Δ, twist Φ_s (the sign of J_⊥ is a gauge choice on an even ring up to Φ_s → Φ_s + πN).

- Easy-plane, −1 < Δ < 1 (c = 2). Levels scale as 1/N at fixed n: D_ps (S^z = ±1) = (2π v/N)/(4K_s), D_4 (S^z = ±2) = 4 D_ps, D_κ=1 = 2π v/N, D_Néel (κ = N/2, S^z = 0) = (2π v/N) K_s with K_s > 1/2. So R = D_Néel/D_ps = 4K_s² > 1.
- Heisenberg, Δ = 1: R = 1, with logarithmic corrections, and the BKT transition to Néel order is exponentially slow.
- Easy-axis, Δ > 1 (valley crystal, c = 1): D_Néel → 0 exponentially in N, D_ps → a finite gap, R → 0.
- Ferro, Δ < −1 (one valley, c = 1): the fully polarised state |M| = N becomes the ground state.

At N = 4 to 10 the CFT ratios are not reached, so the comparison is made with the N-site XXZ ring itself (xxz_fit.py), fitting J, Δ and Φ_s to the relaxed levels.

## 4. Velocities without a normalisation convention

[derived, standard Luttinger-liquid zero modes] Let the uniform flux θ couple to the total charge, and let the twist-relaxed ground-state energy be E(θ) near its minimum and E(N) at fixed L. Define

    v_J = π L ∂²E/∂θ²,      v_N = (L/π) ∂²E/∂N²   (from [E(N+2) + E(N−2) − 2E(N)]/4).

For one species of velocity v and Luttinger parameter K, v_J = vK and v_N = v/K. For two independent species with the same v, v_J = 2v and v_N = v/2. For impenetrable two-component particles (charge = Girardeau fermions of the total density), v_J = v_N = v_F. In all three cases v_c = (v_J v_N)^{1/2} is the charge velocity, independent of how many components carry the charge. The pseudospin velocity is read off the relaxed pseudospin levels (D_κ=1 ≈ 2π v_s/L at large N).

Calibration [numerical, local/v2_rods.jsonl, local/v2_legs.jsonl]. Hard rods (c = 1, n = 1/6): v_J = 0.98 (K = v_J/v = 0.70, against the hard-rod value (1 − n)² = 0.69), the κ = 1 phonon gives v = 1.34, 1.36, 1.38 at L = 24, 30, 36 (hard-rod value 2 sin(π n/(1 − n))/(1 − n) = 1.41), and the Casimir fit of the relaxed ground state over L = 24, 36 gives (c v)_eff = 1.42. Two XX legs (c = 2, n = 1/6, v = 2 per leg): v_J = 4.02 to 4.05 (2v), the phonon gives 1.77 and 1.90 at L = 24, 36, and (c v)_eff = 4.05 against 4. So the relaxed Casimir term reproduces Σ c_i v_i to 1 percent in both controls at these sizes.

## 5. Pair liquid against two-component liquid

C175 found an inter-valley bound pair at K = π (E_b = 0.0195 at |w| = 1, 0.1225 at |w| = 2, size about 21 and 8 sites). If the pairs survive at finite density the infrared is one charge-2 boson with a pseudospin gap Δ_s (Luther-Emery). The ring signatures [derived, standard]:

- Parity gap. Δ_pb = E(N+1) + E(N−1) − 2E(N) at even N tends to 2Δ_s, so L Δ_pb grows linearly. In a gapless two-component liquid L Δ_pb tends to a constant (charge and pseudospin zero modes).
- Pseudospin gap. D_ps (lowest |M| = 2) and every M = 0 level at κ ≠ 0 saturate at a scale of order Δ_s, instead of falling as 1/L.
- Flux period. A liquid of charge-2 bosons has ground-state energy periodic in θ with period π, up to corrections exponentially small in L/ξ_s, so δ_π = min_K E(θ = π) − min_K E(θ = 0) falls faster than 1/L. In a two-component liquid L δ_π tends to a nonzero constant (it is the q = π pseudospin zero mode of section 2).
- One boson. The relaxed spectrum below the gap is a single tower (κ = ±1 phonon at 2π v_c/L, with the next level at the same κ at three times the first), with no dense pseudospin tower below the charge phonon.

All four require L ≫ ξ_s. With ξ_s ≳ pair size, |w| = 1 needs L well above 21 sites per pair at low density, and C175's BCS estimate at n ≈ 0.17 (Δ_s of order 1e-4) puts it out of reach of any ring. The |w| = 2 and n = 1/12 runs are the calibration: there the pair gap should be visible if these diagnostics work.

## 6. Reading the R drift

[heuristic, standard bosonization] For two valleys exchanged by S, the pseudospin boson has one marginal operator at the SU(2)-like line, the inter-valley backscattering cos(2√2 φ_−) of C180's table. Its Kosterlitz-Thouless flow is the only way a 1/L-gapless pseudospin at small sizes can turn into a gapped one at large sizes. On the gapless side the flow stops on the fixed line, and the effective anisotropy (here Δ_R, or equivalently R) approaches a constant, with corrections that decelerate in ln L. On the gapped side the effective coupling grows and crosses the SU(2) value at a scale ξ ~ exp(const/|g|). The data (RESULT claim 6) show R falling by an amount per step in ln N that does not decrease (n = 1/4: −0.74 then −1.05 per unit ln N; n = 1/5: −0.40 then −0.56), at every density, with R already below 1 at n = 1/12 and the crossing scale moving to smaller N as n falls or |w| grows. This is the pattern of the gapped side with a crossover scale beyond L ≈ 40 at n ≥ 1/5. It is consistent with C175's bound inter-valley pair (a Luther-Emery pseudospin gap). At these sizes it is suggestive, not decisive: three sizes cannot separate a slow approach to an easy-plane fixed point from a crossing at N ≈ 10 to 15.

The calibration at |w| = 2 (pair binding 0.12, C175) shows what a pair liquid looks like in the same diagnostics: R collapses (0.35 → 0.003 at n = 1/6, N = 4 → 6), the flux period halves (L·δ_π → 0.1), while D_ps still scales as 1/L at these sizes. So R and δ_π are the sensitive diagnostics of pairing on rings, and D_ps alone is not.
