# Mobile cells: spin–charge separation removes the ring-length obstruction

## Claim

Class: one-site translation-invariant, nearest-neighbour, Hermitian chains on rings of every length L, whose low-energy states consist of spontaneously formed cells between conserved markers (C7's selector with C8's marker hop), each cell carrying a payload. Central charge means Virasoro c.

1. Squeezed space is exact. If a marker hop acts on the payloads as a scalar, the legal Hilbert space at fixed cell number M is (marker positions and cell states) ⊗ (payload)^⊗M modulo one cyclic identification. At fixed momentum this identification is an Ogata–Shiba twist that links the two factors, nothing more. [checked]
2. If, in addition, the markers form a length fluid and the payload chain is critical and stable on its own with central charge c_p, then no relevant coupling between the two sectors survives at generic long-cell density ν, and the spectrum on every ring is that of one c = c_p + 1 theory, with L-dependent zero modes and bounded ground-state degeneracy. The exceptions are a finite set of commensurate ν. [heuristic]
3. Toy (d = 6, each cell carries two qubits, c_p = 1). The frozen chain (u = 0) has ground-state degeneracy 3 on rings with 3 | L and up to L otherwise. The mobile chain has a unique ground state on every ring tested (L = 8 to [[18]], primes included), gaps at q = 2π/L that scale as 1/L with two measured velocities, and c = 2 from three estimators. [numerical]

## Status

Item 1 is exact and checked numerically. Item 2 is standard bosonization and power counting. Item 3 is numerical. What would upgrade it: a mobile cell whose payload has c_p > 1 (the C8 §6 multi-flavour target), at finite coupling.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a sub-session of the coordinating session, 2026-10-27.
- AI assistance: one closed-book gpt-6-astra (high) run, M1, on the kinematics, the coupling classification and the finite-size spectrum (`astra_M1/`, transcript 20261027T200718-63360). It agrees with my derivation on every point used below and supplied two qualifications (marginal current–current mixing, and the meaning of "finite union of towers"). All numerics are mine.
- Checks: the squeezed-coordinate ED against brute-force diagonalization of the d = 6 chain (spectra, translation and entanglement agree to 10⁻¹⁴, `check_ed.py`); the Ogata–Shiba twist against the complete spectrum at u = t, J = 0 (`check_twist.py`, 4 × 10⁻¹⁴); momentum-resolved Lanczos against full diagonalization; VUMPS calibrated on XX (c = 0.99), two decoupled XX chains with velocity ratio 0.4 (2.02), and the t–XX chain (1.96).

## Argument

**Kinematics.** Fix M. A legal state is a reference marker x, cell data d = (ℓ₀, h₀; …; ℓ_{M−1}, h_{M−1}) (lengths and internal states), and payloads z = (z₀, …, z_{M−1}), with (x; d; z) ~ (x + ℓ₀; Sd; Sz) and S the cyclic shift. The cell Hamiltonian never touches z, and payload bonds act on neighbouring z's with coefficients that may depend on d. A momentum-P state with payload cell-momentum Q obeys f(Sd) = e^{i(Pℓ₀ − Q)} f(d), so the length fluid sees a flux set by P and Q. This is Ogata–Shiba's factorization of the U → ∞ Hubbard chain, with cell lengths playing the holes and payloads the spins. In the toy at u = t, J = 0 the complete ring spectrum is free fermions with k_a = (2πI_a + Q + π(W − 1))/L, summed over squeezed momenta Q with their multiplicities. Astra C8's `joint_fluid.md` derives the same gluing for the first-order two-length model with an arbitrary payload.

**Allowed couplings.** Bosonize the length fluid on the cell lattice (density ν, parameter K_ℓ, field φ) and let the payload have primaries O_α of dimension Δ_α at cell momentum q_α. The effective Hamiltonian is a sum over cells of one local function of nearby lengths and payloads, so every coupling has the form e^{i(2πmν + q_α)j} e^{−2imφ} ∂ⁿO_α, of dimension m²K_ℓ + Δ_α + n, and survives only if 2πmν + q_α ≡ 0 mod 2π.

- At generic ν only m = 0 survives. O_α itself renormalizes the payload, which must be stable anyway. ∂φ O_α has dimension 1 + Δ_α and is harmless unless the payload already has a relevant symmetric operator. ∂φ times a payload current is marginal and Gaussian, so it mixes the two bosons without changing c. A payload bond that depends on the cell length therefore does no harm.
- Resonances are isolated. mν ∈ ℤ locks the length fluid (relevant only for m²K_ℓ < 2, the edges of the fluid). mν ∈ ℤ + ½ couples the length density wave to a staggered payload operator. For one payload site per cell this is a spin-Peierls point at ν = ½ of dimension K_ℓ + Δ_ε, relevant for a Heisenberg payload and marginal for XX. C8's natural density at the crossing, ν* = ½ − 1/[π²(b + ½)], sits near ½, so μ has to be detuned off that window.

**Conditions.** (i) Coherent transport: a marker hop maps the payload codes of the two cells to themselves as a scalar (C8 §16 padding, §17 vacancy parent). Otherwise the first-order hop carries a payload operator at the hop rate. (ii) A fluid: μ inside the window between the b- and (b+1)-crystals, whose width is set by J_b = uA_b², and length-dependent payload energies weaker than the length bandwidth (a strong effective attraction between equal lengths phase-separates them). Outside the window an incompatible ring carries a massive mobile defect with a quadratic band, C8 §12. (iii) ν off the resonant set. (iv) A stable payload.

**Rings.** Under (i) to (iv), E = e_∞L + (2π/L)[v_ℓ(x_ℓ − 1/12) + v_p(x_p − c_p/12)], where the zero-mode dimensions depend on M − ρL, the current quantum numbers and the twist. They vary with L (quasi-periodically at irrational density), so the ground state on each ring is the vacuum of an L-dependent sector of the same theory, as in any Luttinger liquid at incommensurate filling, and degeneracies are finite multiplets. A payload that needs M ∈ Nℤ (SU(N)₁) takes the nearest allowed M at cost O(N²/L).

## The toy

Letters L, R, #↑, #↓, X↑, X↓ with C8 §15's arrows, and

h = KQ + (μ/2)(n_# ⊗ 1 + 1 ⊗ n_#) − t Σ_σ(|L Xσ⟩⟨Xσ R| + h.c.) − u Σ_σ(|#σ L⟩⟨R #σ| + h.c.) + (J/2)(S⁺S⁻ + h.c.) on adjacent (#, X) and (X, #) pairs.

At u = 0 this is C7's selector with a spin on each marker and head. The C8 hop moves a marker with its spin, so transport is coherent by construction (L and R are the vacancies). The squeezed payload is an XX chain on 2M spins (c_p = 1), and a bond is active only when a head touches a marker, with first-order weight (2/ℓ) sin²(π/ℓ), so the payload coupling depends on the cell length as strongly as it can. The exact reflection i → −i with L ↔ R maps (#, X) bonds to (X, #) bonds and acts on the squeezed chain as a site-centred reflection, which forbids its dimerization for every u, t. Rotating the # spins by π flips the sign of J, so the chain is stoquastic. Parameters: t = 1, J = 0.5.

| Point | u | μ | Regime |
|---|---|---|---|
| S | 0 | −1 | frozen period-3 crystal |
| A | 1 | −1 | equal hopping (C8 §15), wall density 0.62 |
| B | 0.3 | −0.4 | melted selector, mean cell length 3.3 |
| C | 0.1 | 0 | two-length (3, 4) fluid at the C8 crossing, ν ≈ 0.41 |
| D | 0.3 | −0.4 | B with (#, X) and (X, #) exchange 0.75 and 0.25, reflection broken |

**Ground-state degeneracy.** Frozen chain S, exact over all marker patterns (`run_ed.py static`):

| L | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|
| degeneracy | 3 | 5 | 11 | 3 | 13 | 14 | 3 | 8 |

At L = 10, 13 and 16 a second defect arrangement lies within 0.005 above, so the nearly degenerate manifold is larger still. Mobile chains A, B and C: the ground state is unique and at P = 0 on every L from 8 to [[16]], and for A and B also at [[17, 18]] (`run_ed.py mobile`, Typhon). Perron–Frobenius makes it unique within each (W, S^z) sector, and the sectors never cross at the minimum.

**z = 1 and velocities.** The two modes are probed by their own operators at q = 2π/L: the wall density ρ_q couples to the length boson, and S^z_q to the payload. A Lanczos run from O_q|0⟩ gives the dominant level (weight 0.75 to 1.0) and v = Lω/2π (`run_modes.py`). Velocities depend on density, so the clean z test is the fixed-density sequence W/L = 2/3:

| L | v_ρ (A) | v_S (A) | v_ρ (B) | v_S (B) |
|---|---|---|---|---|
| 9 | 1.630 | 0.430 | 0.874 | 0.386 |
| 12 | 1.641 | 0.416 | 0.909 | 0.384 |
| 15 | [[ ]] | [[ ]] | [[ ]] | [[ ]] |
| 18 | [[ ]] | [[ ]] | [[ ]] | [[ ]] |

[[z statement]] At the other lengths (density 0.55 to 0.62) the same operators give gaps that are linear in 1/L with velocities that follow the density smoothly.

**Central charge.**

| Point | Casimir C / (v_ρ + v_S) | c_half, ring entropy | c from FES (χ ≤ [[ ]]) |
|---|---|---|---|
| S | | 0.98, 0.99 (cuts at cell walls) | [[1.0]] |
| A | [[ ]] | 2.11, 1.94 (9→12→15) [[18]] | [[ ]] |
| B | [[ ]] | 2.09, 1.90 [[18]] | [[ ]] |
| C | | | [[ ]] |
| D | | | [[1.00]] |

The Casimir column fits E₀/L = e_∞ − πC/(6L²) on W/L = 2/3 and divides by the sum of the measured velocities, so a value of 1 means one c = 1 mode per velocity and c = 2 in total. The two pairs (9, 12) and (12, 15) give the same C within [[ ]], so there is no parity alternation at this density. The ring entropy is the half-chain estimator of `contrib/counterexample_search/`, here on the same sequence. FES is S against log ξ for a one-site VUMPS (three-site for S), with the illegal-pair density below 10⁻¹². Control D shows that the symmetry analysis is doing real work: once the reflection is broken the squeezed chain dimerizes, the payload gaps out and FES gives c = 1, the length fluid alone.

**Comparison with the first-order reduction.** At point C the two-length picture predicts a length velocity 2J_b sin(πν)(b + ν) ≈ 0.48, with J_b = uA₃² and A₃ = 0.854. [[measured]]

## Verdict, and what it adds to d₀

The ring-length obstruction is a property of frozen markers, not of cells. With the C8 hop the cell lengths form a c = 1 fluid, the payload lives on the squeezed lattice of cells, and on every ring the spectrum is that of a c = c_p + 1 theory. In the toy this holds at finite coupling with a strongly length-dependent payload bond, from u = 0.1 to u = t. The two places where it can fail are worth keeping in view: commensurate densities where the length density wave locks to a payload operator, and payloads whose relevant operators are allowed once the cells move (control D).

For d₀ the selector costs nothing, since the hop uses C7's four letters. Coherent transport of a payload spread over the cell costs one vacancy symbol in the tracks that carry logical data, because the transferred site must be empty (C8 §16 and §17). With a vacancy that empties the whole site this takes the programmable-cell count d₀ = 4k_c q to 4(k_c q + 1), while c gains one, so the margin c − (d₀ − 1) changes by a constant. Two items from C8 §19 are not addressed here and remain open: a frustration-free core, so that the §17 extension applies, and a mechanism that pins the number of occupied sites per cell while its length fluctuates.

## Reproducing

numpy and scipy only. `ed.py` is the squeezed-coordinate ED (`python3 ed.py` runs the brute-force check), `check_ed.py` and `check_twist.py` the checks, `run_ed.py static|mobile` the ring spectra, `static_entropy.py` the frozen-chain entropies, `run_modes.py` the velocities, `run_fes.py` the VUMPS runs (it imports `../d2_followup_220_1282/vumps.py`), and `analyze.py` writes `data/summary.json`. Laptop runs take seconds to an hour each. L ≥ 16 ran on Typhon (`slurm_*.sh`, `tasks_*.txt`; jobs 17512777, 17512786, 17513038). Data are in `data/`.

## Relation to the rest of the repository

This closes the ring-length item of REPORT §9 and README "Where things stand" at the physics standard, conditional on coherent transport (condition i), which C8 §16 and §17 supply for frustration-free cores. It is the finite-coupling test of `codex/work_C8/joint_fluid.md`, and it adds one item to `contrib/ir_control/`: the length mode and its commensurate resonances with the payload.
