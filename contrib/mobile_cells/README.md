# Mobile cells: spin–charge separation removes the ring-length obstruction

## Claim

Class: translation-invariant, nearest-neighbour, Hermitian chains on rings of every length L whose low-energy states are cells between conserved markers (C7's selector with C8's marker hop), each carrying a payload. c is the Virasoro central charge.

1. If a marker hop acts on payloads as a scalar, the Hilbert space at fixed cell number M is (markers and cells) ⊗ (payload)^⊗M modulo one cyclic identification, an Ogata–Shiba twist linking the factors. [checked]
2. If the markers form a length fluid and the payload chain is critical and stable alone with central charge c_p, no relevant coupling between the sectors survives at generic long-cell density ν, and every ring carries the spectrum of one c = c_p + 1 theory in L-dependent zero-mode sectors. A finite set of commensurate densities is excluded. [heuristic]
3. Toy with c_p = 1 (d = 6). Frozen markers give ground-state degeneracy 3 when 3 | L and between L/2 and L otherwise. Mobile markers give a unique ground state on every ring tested (L = 8 to 18, primes included), gaps at q = 2π/L scaling as 1/L, and c = 2 from a Casimir test with measured velocities, ring entropy, and finite-entanglement scaling (FES). [numerical]
4. Toy with c_p = 2 (d = 10, two incommensurate flavours). Unique ground state on L = 8 to 12; c = 3 is not yet established. [numerical, partial]

## Status

As tagged. The next upgrade is the same test on a cell carrying the C8 §6 multi-flavour payload.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a sub-session of the coordinating session, 2026-10-27.
- AI assistance: one closed-book gpt-6-astra (high) run, M1, on kinematics, couplings and finite-size spectra (`astra_M1/`, transcript 20261027T200718-63360), which agrees with my derivation on every point used. Astra C13 (`codex/work_C13/RESULT.md`) independently reaches the same conditional result. Numerics are mine.
- Checks: the squeezed-coordinate ED against brute force on the d = 6 and d = 10 chains (spectra, translation, entanglement, 10⁻¹⁴); the twist against the complete spectrum at u = t, J = 0 (4 × 10⁻¹⁴); VUMPS calibrated on XX (c = 0.99), two XX chains with velocity ratio 0.4 (2.02) and t–XX (1.96).

## Argument

**Kinematics.** A state is a reference marker x, cell data d and payloads z, with (x; d; z) ~ (x + ℓ₀; Sd; Sz), S the cyclic shift. A momentum-P state with payload cell momentum Q obeys f(Sd) = e^{i(Pℓ₀ − Q)} f(d): the length fluid sees a flux set by P and Q. This is Ogata–Shiba's U → ∞ Hubbard factorization, cell lengths as holes and payloads as spins. In the toy at u = t, J = 0 the complete spectrum is free fermions with k_a = (2πI_a + Q + π(W − 1))/L summed over Q (`check_twist.py`). C8's `joint_fluid.md` has the same gluing at first order.

**Couplings.** Bosonize the length fluid on the cell lattice (density ν, parameter K_ℓ, field φ), with payload primaries O_α of dimension Δ_α at cell momentum q_α. Every coupling is e^{i(2πmν + q_α)j} e^{−2imφ} ∂ⁿO_α, of dimension m²K_ℓ + Δ_α + n, and survives only if 2πmν + q_α ≡ 0 mod 2π. At generic ν only m = 0 survives: O_α renormalizes the payload, ∂φ O_α has dimension 1 + Δ_α, and ∂φ times a payload current is marginal and Gaussian, leaving c unchanged. So a payload bond that depends on cell length is harmless. Resonances are isolated. mν ∈ ℤ locks the fluid, relevant only near its edges. mν ∈ ℤ + ½ couples the length density wave to a staggered payload operator; with one payload site per cell this is a spin-Peierls point at ν = ½ of dimension K_ℓ + Δ_ε, relevant for a Heisenberg payload. C8's natural density ν* = ½ − 1/[π²(b + ½)] is near ½, so μ must be detuned off that window.

**Conditions.** (i) A marker hop maps the two payload codes to themselves as a scalar (C8 §16, §17). (ii) μ lies in the fluid window between the b- and (b+1)-crystals, of width set by J_b = uA_b², with a positive density Hessian (C13 §5); outside it an incompatible ring carries a massive defect with a quadratic band (C8 §12). (iii) ν and the payload fillings avoid the resonant set. (iv) The payload is stable alone.

**Rings.** Then E = e_∞L + (2π/L)[v_ℓ(x_ℓ − 1/12) + v_p(x_p − c_p/12)], with zero modes set by M − ρL, the currents and the twist. They move with L, so each ground state is the vacuum of an L-dependent sector of one theory, as in any incommensurate Luttinger liquid, and along subsequences of L a neighbouring charge sector comes below the 1/L scale (C13 §7). A payload needing M ∈ Nℤ takes the nearest allowed M at cost O(N²/L).

## The toy

Letters L, R, #↑, #↓, X↑, X↓ (C8 §15 arrows) and

h = KQ + (μ/2)(n_# ⊗ 1 + 1 ⊗ n_#) − t Σ_σ(|L Xσ⟩⟨Xσ R| + h.c.) − u Σ_σ(|#σ L⟩⟨R #σ| + h.c.) + (J/2)(S⁺S⁻ + h.c.) on adjacent (#, X), (X, #).

At u = 0 this is C7's selector with a spin on each marker and head. The hop carries the marker's spin, so transport is coherent (L and R are the vacancies). The payload is an XX chain of 2M squeezed spins (c_p = 1) whose bonds act only when a head touches a marker, with first-order weight (2/ℓ) sin²(π/ℓ), the strongest possible length dependence. The reflection i → −i with L ↔ R is site-centred on the squeezed chain and forbids its dimerization for all u, t. Rotating the # spins flips J, so the chain is stoquastic and each (W, S^z) sector has a unique ground state for u, J > 0. Throughout t = 1, J = 0.5.

| Point | u | μ | Regime |
|---|---|---|---|
| S | 0 | −1 | frozen period-3 crystal |
| A | 1 | −1 | equal hopping (C8 §15) |
| B | 0.3 | −0.4 | melted selector, mean cell length 3.3 |
| C | 0.1 | 0 | (3, 4) fluid at the C8 crossing, ν ≈ 0.41 |
| D | 0.3 | −0.4 | B with exchange 0.75 on (#, X) and 0.25 on (X, #) |

## Results

**Frozen versus mobile.** Ground-state degeneracy over all sectors:

| L | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| S | | 3 | 5 | 11 | 3 | 13 | 14 | 3 | 8 | | |
| A | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| B | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| C | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | | |

In S a second defect arrangement lies within 0.005 at L = 10, 13 and 16. Every mobile ground state has P = 0. Up to L = 16 all sectors near the optimum were compared; at L = 17, W = 10 (S^z = 0, 1) and 12; at L = 18, W = 10, 12, 14 at S^z = 0. Through L = 17 the first excitation of A and B sits at L·ΔE/2π = 0.07 to 0.10. At L = 18, B's ground state has W = 10 with the W = 12 sector at 0.009, a charge-sector near-crossing of the kind C13 §7 predicts. At C, rings L = 9 and 12 admit only all-3 tilings and the first gap is 0.016 and 0.008: the small-u fluid needs a few cells of each length before it is asymptotic.

**z = 1.** The wall density ρ_q couples to the length boson and S^z_q to the payload. Lanczos from O_q|0⟩ at q = 2π/L gives one dominant level (weight 0.75 to 1.0) and v = Lω/2π (`run_modes.py`). At W/L = 2/3:

| L | v_ρ (A) | v_S (A) | v_ρ (B) | v_S (B) |
|---|---|---|---|---|
| 9 | 1.630 | 0.430 | 0.874 | 0.386 |
| 12 | 1.641 | 0.416 | 0.909 | 0.384 |
| 15 | 1.647 | 0.410 | 0.927 | 0.383 |

The velocities move by at most 6% from L = 9 to 15, in shrinking steps, so ω ∝ 1/L for both modes. At other lengths they follow the density smoothly.

**Central charge.**

| Point | Casimir C/(v_ρ + v_S) | c_half, ring entropy | c, FES |
|---|---|---|---|
| S | | 0.98, 0.99 | 1.00 |
| A | 1.02 | 2.11, 1.94, 2.06 | 2.12 |
| B | 0.98 | 2.09, 1.90, 2.01 | 2.08 |
| C | | | 2.04 |
| D | | | 1.00 |

Casimir: E₀/L = e_∞ − πC/(6L²) fitted at W/L = 2/3 on L = 9 to 18 (C = 2.097 for A, 1.283 for B), divided by the velocities measured at L = 15, so 1 means one c = 1 mode per velocity. Consecutive pairs give C within 0.02 of the fit for A and 0.03 for B. c_half = 3[S(L₂/2) − S(L₁/2)]/ln(L₂/L₁) on the same sequence; for S, cuts at cell walls of one marker pattern, L = 12→18→24, where only the payload is gapless. FES is 6 dS/d log ξ over the last four bond dimensions (χ ≤ 128 for A and C, 160 for B, 96 for S, 64 for D), illegal-pair density below 10⁻¹². Single χ pairs scatter widely in every two-velocity run, the two-XX calibration included (1.4 to 2.7), so only fits are quoted. In D the broken reflection dimerizes the squeezed chain, the payload gaps, and only the length fluid's c = 1 remains, as the symmetry analysis requires. With Heisenberg exchange instead (SU(2)₁ payload, marginal logarithms) point A gives FES 2.00 at χ ≤ 96.

**First-order reduction.** At C, `joint_fluid.md`'s exact first-order spectrum with J_b = uA₃² (A₃ = 0.854) predicts the length velocity at q = 2π/L; the microscopic ρ_q peak at finite J agrees within 16%, even/odd pattern included:

| L (long cells) | 10 (1) | 11 (2) | 13 (1) | 14 (2) |
|---|---|---|---|---|
| first order | 0.304 | 0.422 | 0.265 | 0.459 |
| ED | 0.316 | 0.385 | 0.299 | 0.387 |

**Two flavours.** Each wall carries qubits a and b, exchanged only within a flavour (J_a = 0.5, J_b = 0.35), coupled on each wall by 0.3 S^z_a S^z_b, with a field 0.06 on b so the two fillings differ (`ed2.py`, u = 0.3, μ = −0.4). The ground state is unique with P = 0 on L = 8 to 12, with a neighbouring magnetization sector at L·ΔE/2π = 0.003 to 0.026, the near-crossings of C13 §7 now in three charges. The c = 3 test is not done: up to χ = 48 the b magnetization still drifts (0.069 to 0.046 per site) and ξ stalls, so FES has no usable slope yet. Job 17513224 continues to χ = 128 and L = 13; a second estimator still has to be added.

## Verdict, and what it adds to d₀

The ring-length obstruction belongs to frozen markers, not to cells. With the C8 hop the cell lengths form a c = 1 fluid, the payload lives on the squeezed lattice, and every ring carries the spectrum of a c = c_p + 1 theory. The toy shows this at finite coupling with a strongly length-dependent payload bond, from u = 0.1 to u = t. It fails at commensurate densities where a length density wave locks to a payload operator, and when a payload operator is forbidden by a symmetry that the motion breaks (control D).

The selector costs nothing in d₀: the hop uses C7's four letters. Transporting a payload spread over a cell costs one vacancy symbol, since the transferred site must be empty (C8 §16, §17). A whole-site vacancy takes the programmable-cell d₀ = 4k_c q to 4(k_c q + 1) while c gains one, so c − (d₀ − 1) moves by a constant. Two items of C8 §19 stay open: a frustration-free core, so that §17 applies, and pinning the number of occupied sites per cell while its length fluctuates.

## Reproducing

numpy and scipy only. ED: `ed.py`, `ed2.py` (each runs its brute-force check), `check_ed.py`, `check_twist.py`, `run_ed.py`, `run_ed2.py`, `ed_one.py` (one cluster sector), `static_entropy.py`, `run_modes.py`. VUMPS: `run_fes.py` (imports `../d2_followup_220_1282/vumps.py`). `analyze.py` writes `data/summary.json`. L ≥ 16 and large χ ran on Typhon (`slurm_*.sh`, `tasks_*.txt`; jobs 17512777, 17512786, 17513038, 17513163, 17513224).

## Relation to the rest of the repository

This settles the ring-length item of REPORT §9 at the physics standard for constructions meeting condition (i); for the programmable cell, (i) reduces to the two open items of C8 §19. It tests `codex/work_C8/joint_fluid.md` at finite coupling, agrees with C13, and adds the length mode and its resonances to `contrib/ir_control/`'s checklist.
