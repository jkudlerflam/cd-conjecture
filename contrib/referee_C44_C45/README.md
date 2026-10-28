# Referee report on C44 and C45, the mobile wrapper and the infrared step at d0 = 4464

## Claim

Under review are `codex/work_C44/` and `codex/work_C45/`, which complete the C39/C42 route. The claim is that for every n there is one Hermitian nearest-neighbour bond h_n on ℂ^4464 ⊗ ℂ^4464 (fixed alphabet, entries depending on n), one-site translation invariant on periodic rings, whose infrared theory is n payload compact bosons plus one length boson, so c = n + 1 in Latorre and Sierra's sense c_T = Σ_a c_a, with every coupling fixed before L → ∞. C44 audits the C8/C32 mobile wrapper on C39's 46-letter core. C45 adds a finite-order local Schrieffer–Wolff (SW) rotation with a volume-independent remainder, a ground-relative auxiliary threshold and an ordered schedule, and concludes c = n + 1 from one stated continuum input S. Then no bound c ≤ f(d) exists, and c ≤ d − 1 fails at d0 = 4464 once n ≥ 4463.

## Status

Refereed. No fatal error. Every finite bound I could rebuild holds, the local SW bound is correctly stated and used, and the schedule is not circular. The infrared conclusion is a correct reduction to S, which is standard physics in content. It is not yet an established counterexample, for reasons of verification rather than logic. The upgrades are listed in the verdict.

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-28.
- AI assistance: none beyond this session. No Codex, since each check is a short script whose specification would be as long as the code.
- Checked by: me. I read C44 and C45 in full, the C40 and C42 results, `referee_C39/`, `referee_d24/`, `ir_control/` and parts of STATUS and `routeA_assembly/`. I rederived C45 (7) to (16) and (22) and C44 (1) to (3), (7), (8), (13), (14) by hand. The scripts here are written from the prose, import nothing from C44 or C45, and I reran no C44 or C45 script.

## Report

### 1. The wrapper on this core (C44)

The wrapper's hypotheses hold for C39's core. Positivity, the biased initialized code (C8 needs termwise annihilation, not equal amplitudes), endpoint sources, the core gap and exact U(1) survive, and the exclusion bounds (H_wrap ≥ −4U, the secant inequality (1), the cut bound (2)) do not depend on D, J_f or any positive core norm. That independence is what the lacunary schedule needs, since D grows without bound with n. The selector tuning gives E_b = eb and E_{b+1} = e(b+1), with e < 0 < μ_sel.

| Rebuilt check (C44 equation) | Result |
| :-- | :-- |
| Charge gaps at n = 1, ρ = 0.2, against C44's cell ED | 0.0147580634663 and 3.36309542068e−5, agreeing in every quoted digit |
| Lower bound (5), n = 1 to 4, p = 1, 2, ρ = 10⁻¹ to 10⁻³ | Holds in all 32 cases, with slack 10⁸ to 10⁶⁰ |
| Asymptotic (6) at ρ = 10⁻³ | Coefficient 2.000 for p = 1, and 0.500, 0.250, 0.125, 0.0625 for p = 2 at n = 1 to 4 |
| Adjacency and endpoint weights (7) | 1 − d ≤ ζ² and r_p ≥ (1 − ζ²)/2 in every case |
| Occupation filter, k = 3 to 8, r = 2 to k + 2 | Left-pinned gap ≥ 1.76 κ_r, F_r ≥ 3.85 times its bound, J_f F_r ≥ 10.5 × 2U, zero energy only at r = k |
| Selector, b = 4 to 4470, l ≤ 6b | μ_c ≥ 0.599, η_b b⁴ ≥ 7.66, (1) holds (equality at its defining point, to 2e−10 roundoff). Reproduces t_sel η_b − U = 1.69357336, 1.39840043, 1.30289697, 1.25125262 |
| Equalization (13), 12000 random instances in the wrapper's block form | Left side at most 0.54 of the right side |

The equalization bound is a correct block inverse and is used correctly. Errors are measured against J0 = 2g²V0 D^(−2n−4), not g², and d_eq = J0 ε_n/(16 C_b) with ε_n ≤ v_* δ/(1024 A_n) puts them below the slowest speed times the smallest margin. The scalar density interaction c_22 − c_21 − c_12 + c_11 in (14), which no chemical potential removes, is budgeted, and one-port mixed terms vanish by bus conservation.

I did not rebuild the source formula (11), the transport overlap (9) or the vanishing first-order compression, which need the 46-letter rule table. C44 checks them only at n = 1, 2, at diagnostic parameters. The bound (5) carries the exponent 2kp where the splitting scales as ρ^(kp). It is valid, and only squares an already absurd scale.

### 2. The local SW step and the schedule (C45)

The construction (9) is the standard finite-order local SW (Datta et al., Bravyi–DiVincenzo–Loss) with one-cell reduced resolvents. The commutator bound (7), the BCH tail constants and the threshold (15) are correct, and the completely positive extension proof of (15) allows defect motion.

| Rebuilt check (C45 equation) | Result |
| :-- | :-- |
| Identity (10), N = 2 to 5 cells, auxiliary gaps 1 and 9 | Error below 7e−16, for a generic and a port-like V |
| Remainder (13), λ = 0.04 down to 0.002 | ‖R‖/λ³ converges, to 2.37, 4.97, 7.68, 10.52 at N = 2 to 5, so linear in N |
| Constants (11) | Actual ‖R‖ is 2.5e−8 to 3.4e−7 of 2(N − 1) C_SW λ³, also at λ = 0.002, inside the range (12) |
| Threshold (15) | Smallest eigenvalue above −7e−16 |
| Three-cell part of P V Q A⁻¹ Q V P | 0.34 relative for a generic V, 1e−16 for a port-like V |
| Gaussian metric (22), 2000 random stiffness matrices, speed ratios to 10⁴ | Smallest ratio 0.77496, against √(3/5) = 0.77460 |

Identity (10) holds for any V with PVP = 0 once A is on-cell, since a surviving cross term has all its excitations in X ∩ Y. Reading it as the sum of the four length-pair bonds needs the port rule that one V excites both cells it touches, which C44 proves. A is on-cell only inside the invariant selector-legal sector (the helix pins at a cut are conditioned on the neighbouring selector letter), so the construction should be stated there.

The bound (13) is used correctly, as a per-cell remainder C_SW λ³ against W = λ² J0/s, and C45 rightly rejects the false criterion of an extensive norm below 1/L. The order is C42's witness, ν and ε_n, then D, g, ρ, s, then Δ, β_n, t_m = W/β_n and the readout, and λ last. Each correction is bounded by a parameter chosen after the quantity it must beat, and nothing later feeds back. Carrying out the prescription at n = 1, 2, 3 confirms this (log10 values below), with the remainder condition binding λ.

| n | ε_n | D | J0 | ρ | s | Δ | β_n | λ | W |
| --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| 1 | −37 | 114 | −1472 | −826 | −12525 | −23277 | −14881 | −38263 | −65472 |
| 2 | −132 | 239 | −4057 | −2234 | −42715 | −80715 | −53640 | −134583 | −230507 |
| 3 | −304 | 441 | −9259 | −5029 | −116127 | −221749 | −150880 | −373056 | −639244 |

S's constants are set to 1, with K_c = D, J_h = U = 1 and the root-patch data replaced by the margins. At n = 4463, D ≥ 16384^(2n+5) alone forces log10 J0 < −6.7 × 10^8 and charge splittings below 10^(−10^12). Circularity could enter only through S's constants A_n in ε_n, which must depend on the normalized reference alone, not on Δ, s or ρ. C45 asserts this. It is physically right, but it is part of S.

### 3. Is S standard?

Stated precisely, S says this. After the rotation the chain is D + R, with D block diagonal, an auxiliary threshold Δ/2 above the reference ground, and a retained block equal to W times a non-resonant quadratic reference plus relative corrections ε, while R is symmetric, quasi-local and ‖R‖_a ≤ εW. Then (a) integrating out the auxiliary sector shifts the connected low-frequency retained vertices by at most A_dec ‖R‖_a, uniformly as W/Δ → 0, and (b) the multi-component compact-boson theory obeys weak-coupling RG (marginal majorant, irrelevant decay) with constants fixed by the normalized reference.

This is how physicists get the Heisenberg chain from the large-U Hubbard chain. One writes H_eff = (4t²/U) Σ (S_i·S_j − 1/4) + O(t⁴/U³), treats higher orders as small local SU(2)-symmetric corrections, ignores that a doublon can decay into about U/J spinons, and invokes the stability of SU(2)₁. Part (a) is the first two steps and part (b) the last, as in the Ogata–Shiba picture of spin–charge separation, which this construction resembles (payload on the squeezed cell lattice, lengths as charge). So S is standard in content, and its uniform form states carefully what is normally left implicit.

The application is unusual in three ways. The protection is a set of momentum margins (10⁻¹⁶, 10⁻²⁸, 10⁻⁴¹ at n = 1, 2, 3, never computed beyond), not a symmetry. Many small parameters are chosen in sequence, not one. And nothing can be checked, whereas the Hubbard reduction is confirmed by the Bethe ansatz. Rigorous work covers a related multi-channel Gaussian problem (Mastropietro–Porta) but not massive decoupling under a gapless extensive band (`ir_control/`, H6). Lack of rigor is not my objection. S carries the entire infrared step, as C45's own tags say.

### 4. Order of limits and the meaning of c

Fixing n and every coupling, taking L → ∞, and only then comparing n is the right order for refuting c ≤ f(d), and no uniformity in n is needed. Latorre and Sierra count c_T = Σ_a c_a for "a collection of unitary CFTs", which matches the convention for n + 1 modes of different speeds. The Casimir energy is velocity weighted, and at irrational ν the zero-mode sectors drift with L, so L(E_0 − eL) has no limit. The central charge should be read from the towers or from the interval entanglement of the infinite chain, ((n + 1)/3) log ℓ, which is speed independent because the normal-mode transformation acts pointwise. The conformal regime begins only beyond about 10^35 cells at n = 1 (10^303 at n = 3).

### 5. Ranked problems

| # | Problem | Rank |
| :-- | :-- | :-- |
| 1 | The infrared step is the assumed input S, untested on this family at any n | Serious |
| 2 | No microscopic c for any member of the 4464-letter family. C44's numerics are projected models or single cells at parameters that violate the budget | Serious as evidence |
| 3 | C44 (9), (11) and PVP = 0 were checked only by the pipeline that wrote them, at n ≤ 2. C39 to C45 come from one model family, and C45 was written with C44's conclusions in hand | Serious for acceptance, fixable |
| 4 | Violating c ≤ d − 1 needs n ≥ 4463. C42's search has run only at n ≤ 3 and has no depth bound (referee_C39, item 3) | Serious as evidence, not a logical gap |
| 5 | Independence of A_n from Δ, s and ρ is asserted. Without it the choice of ε_n would be circular | Fixable, state inside S |
| 6 | C45's RESULT and commit da52641 say the infrared step "closes", while its tags keep S assumed and C44 calls these inputs "substantial" | Fixable, wording |
| 7 | c should be defined through towers or infinite-chain entanglement, not the Casimir coefficient | Fixable |
| 8 | C44 (5) uses the exponent 2kp, twice the true kp | Cosmetic |
| 9 | C44 (17)–(18) and C45 (21) give different mixed-margin lists; C45's is the right minimal one | Cosmetic |
| 10 | A is on-cell only inside the selector-legal sector, where the local SW should be stated | Cosmetic |

### Verdict

No fatal error. C44's wrapper is sound where I could test it, and C45 correctly reduces the infrared question to one standard input. A careful theoretical physicist should accept the argument as a correct derivation conditional on S, and S as the kind of input the field uses routinely. I would not yet call unbounded c at d0 = 4464 established, for the reason the d24 referee gave. The input carrying all of the infrared physics has never been tested on this family, the scheduled couplings (W near 10^(−65000) at n = 1) are beyond any direct test, and the core-dependent parts of the wrapper have only been checked by the pipeline that built them. I would accept it after (i) a calibrated c = 2 from two estimators for the n = 1 member, at couplings where the payload hopping is resolved, (ii) an independent rebuild of C44 (9) and (11) from the 46-letter rule table, and (iii) a human read of C39 to C45.

## Reproducing

numpy and scipy only, from this folder with `python3 -B`. Each `ref_*.py` writes the matching `.json`. `ref_charge.py`, `ref_filter_selector.py` and `ref_equalize_gauss.py` take seconds, `ref_localsw.py` about 70 s, and `ref_schedule.py` (which reads `../../codex/work_C42/budgets_n*.json`) under a second.

## Relation to the rest of the repository

This reviews C44 and C45 in the manner of `contrib/referee_C39/` and `contrib/referee_d24/`. STATUS's route-A all-n line could say that the C44 wrapper audit is refereed with no fatal error, C45's local reduction and schedule are correct, and c = n + 1 at fixed d0 = 4464 is a credible physics-level result conditional on the standard input S, untested on this family. It does not support changing "Where things stand" to "counterexample found".
