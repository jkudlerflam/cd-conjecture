# The c = 16/5 target at d = 4: no tricritical singlet on the boundary of the spin-3/2 c = 3 phase

## Claim

Class: one-site translation-invariant, nearest-neighbour, Hermitian chains on (ℂ⁴)^⊗L with periodic boundary conditions and on-site SU(2) acting as spin 3/2, h = Σ_{J=0}^{3} a_J P_J. Removing the identity and the scale leaves a two-sphere, charted here as (a₁, a₂, a₃) = (sin θ + η, cos θ, sin θ − η) with a₀ = 0. The circle η = 0 is the Sp(4) family cos θ P₅ + sin θ P₁₀, (θ, η) = (π/2, 0) is SU(4), and the chart misses only the two points ±(P₁ − P₃). The central charge is the Virasoro c of the infrared theory.

1. [Astra, checked] The glued Sp(4)₁ × tricritical Ising point (c = 16/5), approached from the c = 3 phase without breaking translation, has exactly two relevant scalar operators invariant under SU(2) (or Sp(4)), translation, both inversions and time reversal. Both have Δ = 6/5. One is the tricritical perturbation ε′, the other is the product M ε of the Sp(4)-singlet mass M = Σ_a iψ^a_L ψ^a_R of the five Majoranas and the tricritical energy ε. The point therefore has codimension two. It is not reached generically on the Sp(4) circle, and in the SU(2) family it can only occur at isolated points. The SU(2) commutant of (spin 3/2)^⊗2 is four-dimensional, so this two-parameter family is already the largest nearest-neighbour family with this symmetry, and no further extension inside the symmetry class exists.
2. [numerical] At L = 8 and 12 the c = 3 phase occupies roughly 0.43π < θ < 1.15π, −0.5 < η < 0.9. It ends in a dimerized phase on the left through a continuous, BKT-like transition, in ferromagnetic or partially polarized phases at the top and right through first-order level crossings, and at the bottom in a region where c_half falls smoothly to 2.0 to 2.5.
3. [numerical] No point of the family shows the tricritical-singlet fingerprint. Wherever the momentum-π singlet is soft and scale invariant, its descendant and zero-momentum partners sit at 2 to 4 times its energy, the free-Majorana pattern, never at the factor 6 of tricritical Ising. Every value of c_half above 3.15 sits next to a ferromagnetic boundary or a forming dimer doublet, with all velocities collapsing or z_eff far from 1.

So no c = 16/5 point was found, and no counterexample to c ≤ d − 1 at d = 4 in this family.

## Status

Numerical and negative, with the codimension count proved at the level of the operator algebra. What would upgrade it is ED at L = 16 or DMRG along the left boundary, which would push the resolution below the current grid. The finite-size argument in the search section says that a tricritical point with O(1) metric factors would already have been visible.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a sub-session of the coordinating session, 2026-10-27.
- AI assistance: one closed-book gpt-6-astra (high) run, F1, on the field theory (`specs/F1_glued_point.md`, output `specs/F1_glued_point_astra_out.md`, transcript 20261027T181712-15382). I derived the operator list and the anomaly assignment myself before reading its answer. The two agree on every item I checked, namely the eight extended primaries and their Z₂ × Z₂ × Fibonacci fusion, ν = 7 for the tricritical chiral parity that reverses the left supercurrent (I recomputed the twisted-sector spins, all ≡ −1/16 mod 1/2), the inherited translation as five-Majorana left parity times tricritical right parity with ν = 6 and order four, the Ramond fields at momentum ±π/2, and the two Δ = 6/5 relevant singlets. The full table of inversion and time-reversal signs I checked only for those two singlets, by continuity with staggered bond operators.
- Numerics: mine. ED reuses `ced.py` and `models.py` from `../counterexample_search` (validated there against brute force). `ed_lean.py` (a lower-memory version for L = 16, not used for results) agrees with `ced.py` to 1e-12 at L = 8. A tricritical control, the O'Brien–Fendley chain, calibrates the spectroscopic ratios.

## Field theory (task 1)

The glued theory is the bosonization of five free Majoranas times the fermionic tricritical Ising model (the N = 1 minimal model). Its relevant scalars, with their quantum numbers, are as follows (Δ, momentum relative to the ground state, SU(2) spins).

| Operator | Δ | momentum | spins | symmetric? |
|---|---|---|---|---|
| ε (tricritical energy) | 1/5 | π | 0 | no, translation odd |
| Ramond pair (s, σ) | 7/10 | ±π/2 | 0, 1, 2, 3 | no, translation ±i |
| M_ab = iψ^a_L ψ^b_R | 1 | π | 0 to 4 | no |
| ε′ | 6/5 | 0 | 0 | yes |
| M ε (trace part) | 6/5 | 0 | 0 | yes |
| traceless M_ab ε, and K^a_± | 6/5 | 0 or π | 1 to 4, and 2 | no SU(2) singlet |
| Ramond pair (s, σ′) | 3/2 | ±π/2 | 0 to 3 | no |

The marginal symmetric scalars are the current-current terms, one with Sp(4) and two (spin-1 and spin-3 channels) with SU(2) only. The momentum assignments follow from the anomaly. In the c = 3 phase translation acts on the six Majoranas as a chiral fermion parity with ν = 6 mod 8. The singlet sector carries ν = 1, and since the left-supercurrent parity of the fermionic tricritical model has ν = 7, translation must act there as the right-supercurrent parity. The alternative assignment, translation through the X± simple current of the extended algebra as in A1 (`../counterexample_search/specs`), has ν = 4. It would put the Ramond scalars at momentum 0 and π and make the Δ = 7/10 singlet translation invariant, a third relevant direction. Continuity with SU(4)₁ excludes it.

The two relevant couplings have the same RG eigenvalue 4/5. With one sign, ε′ alone gives the massless flow tricritical Ising → Ising, back into the c = 3 phase. With the other sign it orders the singlet sector, breaking translation, and then M ε gives the five Majoranas a mass ∝ ⟨ε⟩. Away from the multicritical point the exit from the c = 3 phase is controlled at weak coupling by the Gross–Neveu flow ġ₅ = 3g₅² + g×², ġ× = 4g₅g×, whose massless basin is g₅ < 0, |g×| ≤ −g₅. Its edge is a marginal, BKT-like transition.

The finite-size fingerprint of the glued point on a ring with L ∈ 4ℤ is a spin-0 level φ at momentum π with x = 1/5, its two descendants at π ± 2π/L, and ε′ at momentum 0, with E(ε′)/E(φ) = E(desc)/E(φ) = 6 independently of the two velocities. In the c = 3 phase the same levels give the Ising ratio 2.

## Map of the c = 3 phase (task 2)

The grid is θ ∈ [0.30π, 1.40π] in steps of 0.025π and η ∈ [−1.2, 1.4] in steps of 0.1 (1215 points), plus finer cuts across the left boundary, all at L = 8 and 12 with momentum and spin labels. Sizes 8 and 12 are the commensurate ones for the period-four c = 3 phase, whose ground state sits at momentum 0 for L = 8 and π for L = 12. `map.png` shows the classes and c_half(8 → 12). The calibrated controls of `../counterexample_search` give c_half(8 → 12) = 3.05 at SU(4)₁ and 2.98 at the Sp(4) point θ = π.

- Left, dimerized. The momentum-π singlet descends smoothly as θ decreases, with its descendant at π ± 2π/L near twice its energy, until a doublet collapses (π singlet below 0.005 × 2π/L at L = 12 by θ ≈ 0.40π for −0.4 ≤ η ≤ 0.4). At L ≤ 12 the apparently critical line sits near θ ≈ 0.425π for η ∈ [−0.4, 0.2], inside the true dimerized phase, which on the Sp(4) circle extends up to θ = π/2. This is the drift expected from a marginal transition. The dimerized region reaches θ ≈ 0.45π around η ≈ 0.4 to 0.7.
- Top and right, ferromagnetic and partially polarized. The ground-state spin jumps from 0 to 3L/2 (or to L), so these transitions are first order. Approaching them every velocity falls together. Along θ = π/2 the current level falls from 0.52 to 0.13 (in units 2π/L) between η = 0.5 and 0.9, the singlet-to-current ratio stays at 0.81 to 0.86, and c_half rises from 3.10 to 3.29. The rise is a proximity effect of the z = 2 ferromagnet, not a new sector.
- Bottom. Below η ≈ −0.4, c_half falls to 2.0 to 2.5 by η ≈ −0.8, the lowest level becomes a spin-2 multiplet at momentum π, and the singlet stays stiff (singlet-to-current ratio near 1.5). I did not settle whether this is a crossover towards c = 1 or a separate c ≈ 2 phase. The bottom-left corner meets the period-two c ≈ 1 phase that contains the Heisenberg point.
- Top-left corner, near (0.45π to 0.48π, 0.8). The dimerized tongue, the c = 3 phase, the ferromagnet and a period-two region with c_half 2.0 to 2.6 meet. The spectrum there is a dense set of multiplets with spins up to 6, not a CFT tower.

## Search for a tricritical singlet

The O'Brien–Fendley chain at its tricritical point (`controls/tci_obf.py`) gives ε′/ε = 8.08, 6.86, 6.49 and descendant ratios 6.78, 6.51, 6.34 at L = 8, 12, 16, with L·Δ_ε/2π = 0.325, 0.327, 0.327. So at these sizes a tricritical singlet shows ratios of 6.5 to 8 that decrease with L, and a scale-invariant φ.

`tci_score.py` scans every computed point for a soft (below 0.6 of the current level), scale-invariant (|φ₁₂/φ₈ − 1| < 0.2) momentum-π singlet. Of the 46 such points none has both ratios above 4.5. The largest values are F/φ = 4.3 and descendant ratio 3.8, and along the c = 3/dimer boundary they are typically 3.0 to 3.3 and 2.0 to 2.6, close to the Ising value 2. The two boundary points where the ratios reach about 5 at L = 12, (0.425π, 0.3) and (0.425π, −0.5), have φ falling by 20 to 25 per cent from L = 8 to 12 and ratios growing with L, the opposite of the control, so they are dimerizing. Across the whole c = 3 phase the median descendant ratio is 1.6.

The region in which a tricritical fixed point would control the spectrum at L = 12 has size of order L^{−4/5} ≈ 0.14 in units of the bare couplings. That is larger than the grid spacing, so with O(1) metric factors the grid could not have missed it. This is the main caveat: a point with an unusually small basin, or one reached only at L ≫ 12, is not excluded.

The two residual spin-3/2 points of the earlier search, 188 and 235 (c_half 3.40 and 3.20), lie in a period-two region next to the SU(4) ferromagnet (θ ≈ 1.3π to 1.4π, η < 0), outside the c = 3 phase. There the lowest level is a spin-1 multiplet at momentum π with z_eff ≈ 1.8 to 1.95, so the excess there is not a central charge either.

VUMPS with a four-site cell (`vumps_run.py`) is less useful than hoped. On the dimerized side of the boundary, at (0.43π, 0), the bond alternation decreases only slowly with χ over the accessible range (D = 0.113, 0.109, 0.089 at ξ = 2.5, 5.6, 6.1), which over this range is indistinguishable from the slow decay ξ^{−1/5} a tricritical energy operator would give, so I do not use it as a discriminator.

## Reproducing

numpy, scipy and matplotlib only. Set `OPENBLAS_NUM_THREADS=1` for ED.

- `chart.py --sizes 8 12 --points P.json --chunk i n --nlev 8 --out f.jsonl` runs the chart points of P.json (pairs [θ, η], θ in radians). About 70 s per point at L = 12 on one laptop core, about four times that on Typhon under load.
- `analyze.py`, `detail.py`, `spectrum.py θ/π η`, `mapplot.py` and `tci_score.py` produce the tables, spectra, `map.png` and the fingerprint search from `data/`.
- `controls/tci_obf.py 0.428` is the tricritical control (a few minutes).
- `vumps_run.py --theta T --eta E --n 4 --chis ...` (uses `../d2_followup_220_1282/vumps.py`), and `ed_lean.py` for L = 16 momentum blocks (about 25 to 40 GB per block, untested at that size).
- `cluster/` holds the point lists and the pack scripts. Typhon jobs: 17511595 to 17511598 and 17511668 to 17511679 (grid), 17511707 to 17511710 (left cuts), 17511685 and 17511878 to 17511881 (VUMPS). 17511599 to 17511602 were cancelled and resubmitted as the second set.
- Data: `data/grid/` (grid, JSON lines per point and size), `data/fine/` (left cuts), `data/vumps/`, `data/map_table.json`.

## Relation to the rest of the repository

This closes the "sharpest remaining target" of `../counterexample_search` (REPORT §9.4) in the negative within the spin-3/2 family, and sharpens it. The glued point needs two tunings, not one, so the Sp(4) line is the wrong place to look, and the full SU(2) sphere, the natural two-parameter home, does not seem to contain it. A remaining route would break SU(2) to U(1) ⋊ Z₂ (for example single-ion anisotropy, which keeps the LSM anomaly). That raises the codimension to four but adds many more parameters. It does not change "Where things stand" in the top-level README beyond recording that the 16/5 target was searched for and not found.
