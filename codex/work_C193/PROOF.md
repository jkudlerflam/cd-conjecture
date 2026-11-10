This run was done by a Claude Opus 5.5 agent (C193, 2026-11-10), not by Codex. Astra was not used.

# C193 derivations

Tags [derived], [numerical], [heuristic] open the paragraph they apply to.

## 1. What the stored correlators are

[derived, from `codex/work_C190/z2dmrg.py` measure()] Sites j = 0, ..., L − 1, particle n_j = (1 − Z_j)/2, SP, SM the spin raising and lowering matrices (rotated back by the gauge phases, so everything is in the frame of the Pauli Hamiltonian). On a central window of W sites starting at i_lo (W = 48 at L = 128 and 256, 32 at L = 64):

- Gb[i, j] = ⟨σ⁻_i σ⁺_j⟩ for i < j, Hermitian, with diagonal (1 − ⟨Z_i⟩)/2 = n_i. No string. Up to particle–hole relabelling (irrelevant for an exponent) this is the hard-core boson one-body matrix ⟨b†_i b_j⟩ of the particles n_j.
- Gf[i, j] = ⟨σ⁻_i (Π_{i<k<j} Z_k) σ⁺_j⟩ (the string is applied on the sites strictly between, via `_apply_site(E, B, ZZ_)`), the Jordan–Wigner fermion one-body matrix ⟨c†_i c_j⟩ up to a sign convention.
- CP(r) = ⟨P_{i₀} P_{i₀+r}⟩ − ⟨P_{i₀}⟩⟨P_{i₀+r}⟩ for r = 2, ..., L/2 − 1, i₀ = L/4, with P_j = 2(σ⁺_jσ⁻_{j+1} + σ⁻_jσ⁺_{j+1}) = X_jX_{j+1} + Y_jY_{j+1}.

So CP is a bond–bond (hop–hop) correlator, not a four-boson pair function. This is the main limitation of the stored data (§2).

Map to the effective model [derived, C190 PROOF §2]. With c_j = (Π_{k<j} Z_k)σ⁺_j, the bond is C152's fermion form, and inside the blockaded space at second order in Δ/V it is C173 (3.2). The Schrieffer–Wolff rotation e^S, S = O(Δ/V) Σ (c_j c_{j+1} − h.c.), dresses σ^∓_j by terms O(Δ/V) that change the particle number by the opposite sign (e.g. σ⁻_j → σ⁻_j + O(Δ/V) σ⁺_{j±1}·(local)). In ⟨σ⁻_iσ⁺_j⟩ the cross terms are anomalous (⟨b b⟩), which vanish in the N-conserving effective model, and the diagonal correction is O((Δ/V)²) = 6e-4 (V/Δ = 40) in amplitude with the same exponents. So at leading order Gb is the effective-model boson, Gf the effective-model Jordan–Wigner fermion, and P_j = b†_{j+1}b_j + h.c. (the w hop), each up to an amplitude renormalization 1 + O((Δ/V)²).

## 2. Their components in C191's frame

[derived] Frame (C191, §9.160): two valleys A, B at boson momenta k_A, k_B, weakly coupled as mutually bosonic fields, b_j ≈ e^{ik_A j} e^{iθ_A} + e^{ik_B j} e^{iθ_B} (times density factors), θ_c,s = (θ_A ± θ_B)/√2, φ likewise. Normalization: e^{i√2θ_c} has two-point exponent 1/K_c, e^{i√2θ_s} has 1/K_s, ∂φ_ν has exponent 2 with amplitude ∝ K_ν.

- Boson. e^{iθ_A} = e^{i(θ_c + θ_s)/√2}, so η_b = 1/(4K_c) + 1/(4K_s), at k_A and at k_B.
- Inter-valley pair b_A b_B = e^{i√2θ_c}, exponent 1/K_c, momentum Q = k_A + k_B (π on the symmetric line). Not stored.
- Inter-valley exchange b†_A b_B = e^{−i√2θ_s}, exponent 1/K_s, momentum k_B − k_A. In terms of the measured numbers, 1/K_s = 4η_b − 1/K_c.
- P_j. Inserting the two-valley form,
  b†_{j+1}b_j + b†_j b_{j+1} ⊃ e^{i(k_B − k_A)j} ψ†_Aψ_B [e^{−ik_A} + e^{ik_B}] + h.c.,
  and |e^{−ik_A} + e^{ik_B}| = 2|cos((k_A + k_B)/2)|. On the symmetric line k_A + k_B = π (C192: Q = π) and the exchange amplitude vanishes at leading order. Off the line it is 2|sin((Q − π)/2)| ≈ |Q − π| ≈ 0.005 to 0.015 for the φ = 1.53 bonds. The next term in the gradient expansion, ∂θ e^{i√2θ_s}, has exponent 1/K_s + 2. The intra-valley parts give 2cos k_A (n_A − n_B) to leading order, since cos k_B = −cos k_A on the line: the valley polarization, whose smooth part is −(2cos k_A)²·K_s/(π² r²) in the bulk (n_A − n_B = −√2∂φ_s/π, ⟨∂φ_s∂φ_s⟩ = −K_s/(2r²)), exponent 2. Density harmonics e^{2iφ_ν} have exponents of order K_c + K_s > 2.
- Symmetry check. On the line the valley exchange k → π − k is the antiunitary K·(−1)^j, and P_j → −P_j under it, so P is valley-odd, consistent with the polarization reading.
- Gf. The Jordan–Wigner string over the other valley adds to the boson exponent (C192 item 9, fermion exponent η_b + min(K_c, K_s) > 1). Used only as a frame check.

Conclusion. From the stored data, η_b is measurable, 1/K_c is not, and CP carries 1/K_s only through an exchange piece that is suppressed by |Q − π| (or by a gradient) and is expected to be invisible under the exponent-2 polarization. On the two-leg control (legs = even and odd sites, each a hard-core XX chain), P_j is the inter-leg hop b†_A b_B and CP decays as r^{−1/2}·r^{−1/2} = r^{−1} at q = 0, which is what the control calibrates.

## 3. Open-chain forms used in the fits

[derived] For a field θ obeying Neumann conditions at both ends (φ pinned by the open end, θ free), on sites x = j + 1 ∈ [1, L] with walls at 0 and L + 1, the image construction gives ⟨θ(x)θ(y)⟩ = −g[ln f(x − y) + ln f(x + y)] + const with f(u) = sin(πu/(2(L + 1))). Then

⟨e^{iθ(x)} e^{−iθ(y)}⟩ = exp(−½⟨(θ(x) − θ(y))²⟩) ∝ [ √(f(2x) f(2y)) / (f(x − y) f(x + y)) ]^η = F(x, y)^{−η},

with η = g the bulk exponent. For pairs symmetric about the centre, x + y = L + 1 gives f(x + y) = 1 and f(2x) = f(2y) = cos(πr/(2(L + 1))), so F = tan(πr/(2(L + 1))), C192's form. The same F with η = 2 is used for the smooth part of a ∂φ correlator only as an amplitude normalization (§6), where the Dirichlet image enters with the opposite sign. At r ≪ L, F ∝ r and the forms reduce to bulk power laws.

Estimators (c193.py). E1: symmetric pairs about the window centre, averaged over 9 centres, demodulated at each of the two largest Fourier peaks (refined by the phase slope), four binomial low-pass passes (response cos⁸(Δk/2) at the other valley, 3.5e-3 at Δk = 2.66), fit of ln|h| against ln F. E2: same on CP(r) with F(i₀ + 1, i₀ + 1 + r). Frozen ranges are in RESULT.md.

## 4. K_s from η_b

[derived] 1/K_s = 4η_b − 1/K_c, so K_s = 1/(4η_b − 1/K_c) and δK_s = 4K_s² δη_b. With 1/K_c = 0.62 and η_b ≈ 0.31, K_s ≈ 1.6, an error of 0.005 in η_b is 0.05 in K_s and a shift of 0.03 in 1/K_c is 0.08 in K_s. K_s > 1 needs 1/K_c > 4η_b − 1, which at η_b = 0.31 to 0.34 is 1/K_c > 0.24 to 0.36, far below any measured value.

## 5. Kosterlitz–Thouless flow and the drift bound

[derived, Giamarchi normalization as in C191 ktflow.py] Pseudospin sine-Gordon with g₁ cos(√8φ_s) of dimension 2K_s, y₁ = g₁/(πu_s):

dK_s/dl = −K_s² y₁²/2,   dy₁/dl = (2 − 2K_s) y₁.

Near K_s = 1 put K_s = 1 + y_∥/2. Then dy_∥/dl = 2dK_s/dl = −K_s²y₁² ≈ −y₁² and dy₁/dl = −y_∥y₁, the BKT equations with the invariant y_∥² − y₁². The flow ends on the gapless line iff y_∥ ≥ |y₁| (with y_∥ > 0), at y_∥* = √(y_∥² − y₁²).

Bound from a drift. If K_s(l) is measured at two scales l₁, l₂ = l₁ + Δl, integrating dK_s/dl with y₁ and K_s nearly constant over the interval gives

y₁² ≈ −2ΔK_s / (K̄_s² Δl),   K̄_s = (K_s(l₁) + K_s(l₂))/2.

A drift ΔK_s ≥ 0 bounds |y₁| by the error alone, y₁² ≤ 2σ_ΔK/(K̄²Δl). The upper value quoted uses ΔK − σ. With Δl = ln 2 for L = 128 → 256. The linearization is not accurate at K_s ≈ 1.6 (y_∥ ≈ 1.2), so the criterion y_∥ > |y₁| is applied as registered and the full flow is used only for the extrapolation below.

Self-consistency of a KT reading [derived]. If a drift dK/dl is caused by y₁, then y₁ ∝ e^{−(2K−2)l} and the drift still to come is ∫_l^∞ (dK/dl) dl' = (dK/dl)/(2(2K − 2)). So a drift of −0.17 per doubling at K_s = 1.6 implies only about −0.1 more, K_s(∞) ≈ 1.5. A KT-driven drift that ends in a gap needs K_s to reach 1 while y₁ grows, which requires the drift to accelerate, not to decay.

What the 48-site window does to the drift [derived]. E1 at L = 128 and 256 uses the same separations, r ∈ [8, 34]. In the RG the scale is set by r, not by L, once r ≪ L, and the open-chain form removes the r/L dependence. So K_s(256) − K_s(128) at fixed window is not the flow over Δl = ln 2. In a pure conformal liquid it must vanish. A nonzero value measures an L dependence of the state at fixed r (boundary, density, N mixing, χ), not the KT flow. The secondary drift, between R3 = [8, 28] and R4 = [12, 34] at L = 256, is a scale change Δl = ln(√(12·34)/√(8·28)) = 0.30, and it is the one that maps onto the flow (with the caveat that subleading corrections to scaling also make the local exponent run).

## 6. The smooth part of C_P as a K_s ratio (unregistered)

[derived] CP_smooth(x, y) ≈ −Z_P² (K_s/π²) × [bulk 1/r² with images] and on the open chain F ∝ (π/(2(L + 1)))·r at r ≪ L, so fitting |CP| = A F^{−2} gives A(L) ∝ Z_P² K_s/(L + 1)². Z_P (≈ 2cos k_A) is a lattice constant independent of L, so K_s(256)/K_s(128) = [A(256)/A(128)] (257/129)². This needs only the exponent-2 form, which E2 reads (2.04 to 2.07), not the unknown Z_P. It was not calibrated on a control (on the XX chain q = 0 and 2k_F compete), so it is a secondary reading only.
