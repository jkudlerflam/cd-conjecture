# C194 PROOF: derivations behind the gap test

Model and conventions as in RESULT.md (C173 PROOF eq. 3.2, |w| = 1, φ = π/2, |t₂| = 1). Two valleys a, b at the band minima ε_min = −9/4. C191's frame: K_s is defined by x(cos √8φ_s) = 2K_s, two decoupled XX legs have K_c = K_s = 1, an attractive two-component gas has K_s < 1.

## 1. The dilute limit of the pair-breaking gap

Take N = 2M particles, M ≥ 2, at density n → 0 on a chain of length L → ∞ with n fixed first. The two-body problem has one bound state of the a–b pair, with energy ε_p = 2ε_min − E_b and E_b = 0.0195 at |w| = 1 (C175, checked by two independent codes and by the coordinator's ring ED). The trimer (pair plus atom) is bound by E_3 = 6.0e-4 below pair + atom (C175 claim 7). The tetramer is not resolved and is taken unbound at |w| = 1 (C175 open item 2), so the dilute ground state of even N is a gas of M pairs.

Interactions between clusters vanish as n → 0 (their energy per cluster is O(n²) for impenetrable composites in 1D, kinetic energy of the relative motion, and O(e^{−r/a₊}) from the overlap tails). So, to leading order,

- E(2M) = M ε_p,
- E(2M ± 1) = (M − 1) ε_p + ε_t with ε_t = ε_p + ε_min − E_3 for 2M + 1 (trimer plus M − 1 pairs), and for 2M − 1 the same with M − 1 replaced by M − 2, i.e. E(2M − 1) = (M − 2)ε_p + ε_t.

If the trimer is ignored (E_3 = 0) the odd states are pairs plus one free atom, E(2M + 1) = Mε_p + ε_min and E(2M − 1) = (M − 1)ε_p + ε_min. Then

Δ_pb = E(2M + 1) + E(2M − 1) − 2E(2M) = 2ε_min − ε_p = E_b.

With the trimer, Δ_pb = 2ε_t − 3ε_p = 2(ε_min − E_3) − ε_p = E_b − 2E_3. So

E_b − 2E_3 = 0.0183 ≤ Δ_pb(n → 0) ≤ E_b = 0.0195.

In the language of a two-component liquid, Δ_pb = 2Δ_s with Δ_s the energy of one unpaired particle (one pseudospin soliton, with the charge part separated off, §2), so Δ_s → E_b/2 up to the trimer correction. That is the "half the two-body binding energy" of the spec. For M = 1 (N = 2, so N − 1 = 1) the odd state N − 1 is a free atom and N + 1 = 3 is a trimer, which gives E_b − E_3 instead; this matters only at L = 64 for n = 1/32 and L = 48 for n = 1/24.

Finite size. In an open box each cluster carries a confinement energy π²/(2m_X L²) (m_X its band mass), and the gas of clusters has the 1/L charge term of §2. So at fixed small n, Δ_pb(L) = Δ_∞ + a/L + b/L², the registered form.

## 2. The finite-size form of Δ_pb on open chains

Gapless two-component liquid (c = 2). Bosonize each valley, φ_c = (φ_a + φ_b)/√2, φ_s = (φ_a − φ_b)/√2. With open ends the zero-mode energy of a sector with δN = N − N̄ particles and valley imbalance m = N_a − N_b is (C191 PROOF 4 normalization, free fermions per valley give E = (πv/2L) Σ_α δN_α²)

E(δN, m) = (π/(4L)) [ (v_c/K_c) δN² + (v_s/K_s) m² ] + O(1/L²).

At even N the ground state has m = 0; N ± 1 has m = ±1. So

Δ_pb(L) = (π/(2L)) (v_c/K_c + v_s/K_s) + O(1/L²),

which is C191's L·Δ_P formula (boundary with φ_s pinned at short scales, C191 RESULT, comparison item 2). Check on the two-leg XX control: v_c = v_s = v, K = 1 gives Δ_pb = πv/L. For two legs of ℓ = L/2 sites at leg filling ν, v (in chain sites) = 2 · 2 sin(πν), and at n = 1/4 (ν = 1/4) L·Δ_pb = 2π·√2 = 8.886, against the exact 8.870 at L = 256 (controls_exact.py). Hard rods: one component, Δ_pb = πv/(K L) → 0 as 1/L, exact values in controls_exact.py.

Gapped pseudospin (c = 1, Luther–Emery). The charge sector keeps its zero mode, (π/(4L))(v_c/K_c)δN², and the odd sector costs one soliton, Δ_s, plus its confinement energy π²/(2m_s L²) and an L-independent interaction with the walls which falls off as e^{−L/ξ_s}. So

Δ_pb(L) = 2Δ_s + (π v_c/(2K_c))/L + O(1/L²) for L ≫ ξ_s,

and the same a/L + b/L² form with Δ_∞ = 2Δ_s. For L ≲ ξ_s the gapless form applies and the extrapolation cannot see the gap. So a density with ξ_s of several hundred sites is outside the resolution of L ≤ 256 by construction, and Δ_∞ there will read zero within its error. That is the main limitation of the direct test, and why the infinite chain is used too.

## 3. The Kosterlitz–Thouless form of the gap near the boundary

C191's flow (ktflow.py) in Giamarchi normalization: dK/dl = −K²y²/2, dy/dl = (2 − 2K)y, with y = g₁/(πu). Near K = 1 put x = 2(K − 1), so to leading order

dx/dl = −y², dy/dl = −xy.

Then d(x² − y²)/dl = 2x(−y²) − 2y(−xy) = 0, so C = x² − y² is an invariant. The separatrix is C = 0 with x = |y|. On the gapped side write C = −μ² < 0. Along the flow dl = −dx/y² = −dx/(x² + μ²), so the scale at which the flow leaves the weak-coupling region (x → −O(1)) is

l* = ∫_{−O(1)}^{x₀} dx/(x² + μ²) = (1/μ)[arctan(x₀/μ) + arctan(O(1)/μ)] → π/μ as μ → 0.

The bare couplings depend smoothly on n, so near the crossing μ² ∝ (n_c − n) at first order, and

ξ_s ∝ e^{l*} ∝ exp(B/√(n_c − n)), Δ_s ∝ u_s/ξ_s ∝ exp(−B/√(n_c − n)).

That is the registered KT fit form, with A, B, n_c free. On the gapless side C = +μ² > 0, the flow ends on the fixed line at x* = μ, so the infrared K_s* = 1 + μ/2 ∝ 1 + √(n − n_c) · const. The measured K_s at L = 128 to 256 is the running K at l = ln(L/a₀), which lies between the bare K_s⁰ and K_s*, so close to n_c it drifts down with L (logarithmically, d K/d ln L = −K² y²/2 with y² = x² − μ²). Away from n_c (x₀ ≫ |y₀|) the drift is negligible, as C191 found (3e-5 at n = 1/8).

Two consequences for the analysis.

(i) A naive scan places n_c too low. At n_c − n = 0.01 with B of order 1, ξ_s ∝ e^{10}, far beyond any size here, so the last density at which a gap is seen lies well below n_c. The fit to exp(−B/√(n_c − n)) extrapolates the visible gaps to the essential singularity. With few points it constrains n_c weakly, and its error is reported from the fit covariance and the leave-one-out spread.

(ii) From the gapless side. At n = n_c the bare value is K_s⁰ = 1 + |y₀|/2 (x₀ = |y₀|), and along the separatrix K runs down toward 1 without reaching it at any finite scale. So the measured K_s(L; n_c) lies in (1, 1 + |y₀|/2). Since K_s(L; n) rises with n, the density where the measured K_s crosses 1 lies below n_c (a lower estimate), and the density where it crosses 1 + |y₁|/2 lies above n_c (an upper estimate), provided |y₁| is the right bare coupling. Below n_c the running K_s still exceeds 1 at scales short of the crossover, so the lower estimate can sit well below n_c when ξ_s ≫ L. Since (K_s* − 1)² ∝ (n − n_c), a fit of (K_s − 1)² linear in n is a third, unregistered, estimate.
