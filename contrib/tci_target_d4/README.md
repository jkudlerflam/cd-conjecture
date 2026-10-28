# The c = 16/5 target at d = 4: no tricritical singlet on the boundary of the spin-3/2 c = 3 phase

## Claim

Class: one-site translation-invariant, nearest-neighbour, Hermitian chains on (ℂ⁴)^⊗L with periodic boundary conditions and on-site SU(2) acting as spin 3/2, h = Σ_{J=0}^{3} a_J P_J. Up to shift and scale this is a two-sphere, charted as (a₁, a₂, a₃) = (sin θ + η, cos θ, sin θ − η), a₀ = 0. The circle η = 0 is the Sp(4) family cos θ P₅ + sin θ P₁₀, (π/2, 0) is SU(4), and the chart misses only ±(P₁ − P₃). The central charge is the Virasoro c of the infrared theory.

1. [Astra, checked] The glued Sp(4)₁ × tricritical Ising point (c = 16/5), approached from the c = 3 phase without breaking translation, has exactly two relevant scalars invariant under SU(2) (or Sp(4)), translation, both inversions and time reversal, both with Δ = 6/5. One is the tricritical ε′, the other the product Mε of the Sp(4)-singlet mass M = Σ_a iψ^a_L ψ^a_R of the five Majoranas with the tricritical energy ε. So the point has codimension two. It is not generic on the Sp(4) circle, and in the SU(2) family it can only be an isolated point. That family is already the largest nearest-neighbour family with spin-3/2 SU(2) symmetry (the commutant of (spin 3/2)^⊗2 is four-dimensional), so there is no further extension inside the symmetry class.
2. [numerical] At L = 8, 12 the c = 3 phase (c_half within 0.1 of 3) occupies roughly 0.43π < θ < 1.1π, −0.4 < η < 0.8. It ends in a dimerized phase on the left, apparently continuously, in ferromagnetic or partially polarized phases at the top and right through first-order level crossings, and at the bottom in a region where c_half falls smoothly to 2.0 to 2.5.
3. [numerical] No point shows the tricritical-singlet fingerprint. All along the c = 3/dimer boundary the soft momentum-π singlet has an effective dimension 0.5 to 0.7 in units of its own velocity, against 1/5 for the tricritical energy (0.18 in a control at the same sizes). Every c_half above 3.15 sits next to a ferromagnetic or partially polarized boundary, next to a forming dimer doublet, or at an incommensurate level crossing.

So no c = 16/5 point was found, and no counterexample to c ≤ d − 1 at d = 4 in this family.

## Status

Numerical and negative, with the codimension count established at the level of the operator algebra. L = 16 data on the left boundary would upgrade it. ED at L = 16 for the most extreme boundary point is running (see Reproducing).

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a sub-session of the coordinating session, 2026-10-27.
- AI assistance: one closed-book gpt-6-astra (high) run, F1, on the field theory (`specs/F1_glued_point.md`, output `specs/F1_glued_point_astra_out.md`, transcript 20261027T181712-15382). I derived the operator list and the anomaly assignment before reading its answer, and we agree on everything I checked: the eight extended primaries with Z₂ × Z₂ × Fibonacci fusion, ν = 7 for the tricritical parity that reverses the left supercurrent (I recomputed the twisted-sector spins, all −1/16 mod 1/2), translation as five-Majorana left parity times tricritical right parity (ν = 6, order four), the Ramond fields at ±π/2, and the two Δ = 6/5 singlets. Inversion and time-reversal signs I checked only for those two, by continuity with staggered bond operators.
- Numerics: mine. ED reuses `ced.py` and `models.py` from `../counterexample_search`. `ed_lean.py` (lower memory, for L = 16) agrees with `ced.py` to 1e-12 at L = 8. The O'Brien–Fendley chain at its tricritical point calibrates the spectroscopy.

## Field theory (task 1)

The glued theory is the bosonization of five free Majoranas times the fermionic tricritical Ising model. Its relevant scalars, with momentum relative to the ground state and SU(2) spins, are these.

| Operator | Δ | momentum | spins | lattice-symmetric? |
|---|---|---|---|---|
| ε | 1/5 | π | 0 | no, translation odd |
| Ramond (s, σ) | 7/10 | ±π/2 | 0 to 3 | no, translation ±i |
| M_ab = iψ^a_L ψ^b_R | 1 | π | 0 to 4 | no |
| ε′ | 6/5 | 0 | 0 | yes |
| Mε | 6/5 | 0 | 0 | yes |
| traceless M_ab ε, K^a_± | 6/5 | 0 or π | 1 to 4, 2 | no SU(2) singlet |
| Ramond (s, σ′) | 3/2 | ±π/2 | 0 to 3 | no |

The marginal symmetric scalars are current-current terms, one with Sp(4) and two (spin-1 and spin-3 channels) with SU(2). The momenta follow from the anomaly. In the c = 3 phase translation is a chiral fermion parity with ν = 6 on the six Majoranas, and the singlet carries ν = 1. Since the left-supercurrent parity of the tricritical model has ν = 7, translation must act there as the right-supercurrent parity. The alternative, translation through the X± simple current of the extended algebra (as in A1, `../counterexample_search/specs`), has ν = 4, puts the Ramond scalars at 0 and π, and leaves a Δ = 7/10 singlet at momentum 0, a likely third relevant direction. Continuity with SU(4)₁ excludes it.

Both relevant couplings have RG eigenvalue 4/5. With one sign ε′ gives the flow tricritical Ising → Ising, back into the c = 3 phase. With the other it orders the singlet, breaking translation, and Mε then gaps the five Majoranas with mass ∝ ⟨ε⟩. The earlier note missed Mε when it said one tuning would suffice. Away from the multicritical point, the exit from c = 3 is governed at weak coupling by the Gross–Neveu flow ġ₅ = 3g₅² + g×², ġ× = 4g₅g×, whose massless basin g₅ < 0, |g×| ≤ −g₅ ends in a marginal, BKT-like transition.

On a ring with L ∈ 4ℤ the glued point has a spin-0 level φ at momentum π (x = 1/5), two descendants at π ± 2π/L, and ε′ at momentum 0, with E(ε′)/E(φ) = E(desc)/E(φ) = 6 for any velocities. In general E(desc)/E(φ) = 1 + 1/x_φ, so x_eff = 1/(E(desc)/E(φ) − 1) is the dimension of φ in units of its own velocity. A free singlet Majorana (the c = 3 phase) gives 1.

## Map of the c = 3 phase (task 2)

ED at L = 8 and 12 with momentum and spin labels on a grid θ ∈ [0.30π, 1.40π] (steps 0.025π) by η ∈ [−1.2, 1.4] (steps 0.1), 988 of 1215 points so far, plus 292 points on cuts across the left boundary (steps 0.01π and 0.005π). These sizes are the commensurate ones for the period-four c = 3 phase, whose ground state sits at momentum 0 for L = 8 and π for L = 12. `map.png` shows the classes and c_half(8 → 12). The controls of `../counterexample_search` give c_half(8 → 12) = 3.05 at SU(4)₁ and 2.98 at Sp(4) θ = π.

- Left, dimerized. As θ decreases the π singlet descends smoothly until a doublet collapses (it falls by factors of 2 to 12 between L = 8 and 12 at θ ≈ 0.40π to 0.44π). There is no first-order jump at L ≤ 12, consistent with the marginal mechanism, though a weak first-order transition cannot be excluded. The last scale-invariant singlet sits near θ ≈ 0.42π for −0.3 ≤ η ≤ 0.1 and moves to 0.455π by η = 0.6. On the Sp(4) circle the true transition is at π/2, so this drift is the finite-size shadow of a marginal transition.
- Top and right, ferromagnetic and partially polarized. The ground-state spin jumps from 0 to 3L/2 or L, so these are first order. On approach, multiplets of spin 4 to 13 from the polarized phases invade the low spectrum. Along θ = 0.7π, c_half climbs from 3.05 (η = 0.5) to 3.54 (η = 1.1) while the current level falls from 0.64 to 0.13 (units 2π/L) and a spin-6 level at zero momentum drops below x = 1, a state neither the c = 3 nor the glued theory has below x = 2.
- Bottom. Below η ≈ −0.4, c_half falls to 2.0 to 2.5, the lowest level becomes a spin-2 multiplet at π, and the singlet stays stiff (about 1.5 times the current level). I did not settle whether this is a crossover towards c = 1 or a separate c ≈ 2 phase. The bottom-left corner meets the period-two c ≈ 1 phase containing the Heisenberg point.
- Top-left corner, near (0.45π, 0.8). The dimerized tongue, the c = 3 phase, the ferromagnet and a period-two region with c_half 2.0 to 2.6 meet, in a dense set of multiplets with spins up to 6.

## Search for a tricritical singlet

The O'Brien–Fendley chain at its tricritical point (`controls/tci_obf.py`) gives ε′/ε = 8.08, 6.86, 6.49 and descendant ratios 6.78, 6.51, 6.34 at L = 8, 12, 16 (x_eff = 0.17 to 0.19), with L·Δ_ε/2π = 0.325, 0.327, 0.327. So a tricritical singlet is recognizable at these sizes, with ratios decreasing towards 6.

`boundary.py` picks, on each η cut, the softest π singlet still scale invariant between L = 8 and 12 (`boundary.png`, `data/boundary_table.json`). For −0.3 ≤ η ≤ 0.6 the zero-momentum singlet sits at 3.2 to 3.6 times φ and the descendant at 2.5 to 3.1 times φ at L = 12, so x_eff = 0.49 to 0.68 (0.67 to 0.86 at L = 8). The descendant ratio stays a factor of two below the tricritical value everywhere. `tci_score.py` applies the same test to all 157 points with a soft, scale-invariant π singlet. Only one passes both ratios above 4.5, at (1.325π, −0.5) in the ferromagnet-adjacent period-two region, where c_half(8 → 12) ≈ 5 and the L = 8 and 12 ratios disagree.

At those boundary points c_half is 3.08 to 3.15 for −0.3 ≤ η ≤ 0.3, inside the c = 3 band, and 3.21, 3.31, 3.33 at η = 0.4, 0.5, 0.6. The excess near the top-left corner comes with unchanged singlet ratios (x_eff ≈ 0.5), and c_half keeps rising (to 3.6 to 3.8) as one moves into the forming doublet, so I read it as corner soft modes plus doublet formation. The most extreme point, (0.445π, 0.5), with φ a quarter of the current level and c_half = 3.31, is the one sent to L = 16.

A tricritical fixed point would control the L = 12 spectrum within a region of size about L^{−4/5} ≈ 0.14 in bare-coupling units, larger than the grid spacing. With O(1) metric factors the scan could not have missed it. A point with an unusually small basin, or one visible only at L ≫ 12, is not excluded, and that is the main caveat.

The earlier residual points 188 and 235 (c_half 3.40, 3.20) lie in the same ferromagnet-adjacent period-two region (θ ≈ 1.3π to 1.4π, η < 0), outside the c = 3 phase. There c_half ranges from 1.5 to 5.8 with z_eff jumping between −5.9 and 1.95 from one grid point to the next. Next to the SU(4) ferromagnet the ground state lies close to the permutation-symmetric manifold, whose Dicke-type states have S ≈ (3/2) ln L at d = 4, enough by itself to fake c_half ≈ 4.5. So that excess is not a central charge either.

VUMPS with a four-site cell (`vumps_run.py`) separates the c = 3 interior from the boundary but not the boundary from a tricritical point. At the Sp(4) point θ = π the finite-χ bond alternation obeys D·ξ = 0.14 ± 0.01 for χ = 16 to 64, the x = 1 decay of a Majorana mass. At (0.43π, 0), on the dimerized side, D falls only from 0.113 to 0.079 as ξ grows from 2.5 to 9.8, which over this range mimics the ξ^{−1/5} of a tricritical energy.

## Reproducing

numpy, scipy and matplotlib only, `OPENBLAS_NUM_THREADS=1` for ED.

- `chart.py --sizes 8 12 --points P.json --chunk i n --nlev 8 --out f.jsonl` (P.json lists [θ, η], θ in radians), about 70 s per L = 12 point on a free laptop core.
- `analyze.py`, `detail.py`, `spectrum.py θ/π η`, `mapplot.py`, `boundary.py`, `tci_score.py` build the tables, spectra and figures from `data/`. `controls/tci_obf.py 0.428` is the control.
- `vumps_run.py` uses `../d2_followup_220_1282/vumps.py`. `ed_lean.py --L 16 --momenta k` builds one L = 16 block (2.4 × 10⁷ states, 40 to 110 GB).
- `cluster/` has the point lists and pack scripts. Typhon jobs, under `/scratch/lustre/jkudlerflam/tci_d4/tci_target_d4/`: grid 17511595 to 17511598 and 17511668 to 17511679, with incomplete chunks rerun as 17512610 to 17512612 (writing `data/grid1r_c*.jsonl`); left cuts 17511707 to 17511710; L = 16 at (0.445π, 0.5), (0.435π, 0.3) and (0.6π, 0) as 17512248 to 17512250 (momenta 0 with entropy, 8 and 7, writing `data/l16_k*.jsonl`). 17511599 to 17511602 and the VUMPS jobs 17511685, 17511878 to 17511881 were cancelled.
- Data: `data/grid/`, `data/left/`, `data/fine/` (JSON lines per point and size), `data/vumps/`, `data/map_table.json`, `data/boundary_table.json`.

## Relation to the rest of the repository

This answers the "sharpest remaining target" of `../counterexample_search` (REPORT §9.4) in the negative within the spin-3/2 family, and corrects it: the glued point needs two tunings, not one. The Sp(4) line is therefore the wrong place to look, and the full SU(2) sphere, its natural two-parameter home, does not seem to contain it. A remaining route breaks SU(2) to U(1) ⋊ Z₂ (single-ion anisotropy keeps the LSM anomaly), which raises the codimension to about four but adds many parameters. Nothing changes in "Where things stand" beyond recording that the 16/5 target was searched for and not found.
