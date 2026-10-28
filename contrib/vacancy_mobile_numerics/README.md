# One vacancy on the route-A demonstrator: numerics of the mobile completion

## Claim

Class: one real nearest-neighbour bond on ℂ¹¹ ⊗ ℂ¹¹, one-site translation invariant, rings of every length, Virasoro c. It is the ten-letter demonstrator bond h of `../routeA_demonstrator/` (λ = 0.15, other coefficients as there) lifted by the collaborator's one-vacancy construction, proposed for the 23-letter bond in `../counterexample_campaign/doping/`:

h₁₁ = −t Σ_a (|0 a⟩⟨a 0| + h.c.) + J h_occ − (μ_v/2)(n ⊗ 1 + 1 ⊗ n).

Unless stated, t = 1, J = 0.25, μ_v = 0.6 (occupied density ρ ≈ 0.6). The frozen chain (no vacancy) has c = 1 on top of a threefold broken translation.

1. [numerical, exact ED] On every ring L = 5 to 11 (L = 11 in its best first-order sector), primes included, the ground state has momentum zero and is unique, except for the exact Q = ±½ doublet of the logical U(1) when the occupied number M is odd (as in any odd XXZ ring). The frozen chain is 3- to 22-fold degenerate there.
2. [numerical] The frozen triplet's partners become charge-twist levels (flux ±2π/3) with x within 1.5% of the free value, which tends to 1/9.
3. [numerical] First-order spin–charge separation holds level by level, with errors ∝ J².
4. [numerical] The charge boson has z = 1 and calibrated Casimir c = 0.998.
5. [numerical, negative] The payload is about 360 times slower than the charge, and neither Casimir nor finite-entanglement scaling (FES) on the full chain resolves its c. So c = c_payload + 1 = 2 rests on item 3 plus the demonstrator's c_payload = 1, not on a direct measurement.
6. [numerical plus estimate] Rings shorter than L* ≈ 700 sites (t = 4J) carry a pattern defect, which keeps the ground state unique but adds a slow rotor band.
7. [numerical] For the 23-letter bond at n = 1 plus a vacancy (d = 24), the frozen 2- and 5-fold ground states become unique on L = 5 to 8.

## Status

Numerical, as tagged. A direct c = 2 needs FES at bond dimension of order 10³, or a payload much faster than the demonstrator's second-order hopping.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28, testing the collaborator's proposal. No Codex call.
- `brute.py` builds H from the literal 121 × 121 bond on the physical basis and matches the squeezed solver to 2 × 10⁻¹³ in five sectors (frozen, mobile, defect words). `d24.py check` does the same for the literal 576 × 576 lift of `bond23.py`. `vacbond.py` equals `bond24.vacancy_lift` entry by entry.
- FES is calibrated on XX and on two decoupled XX chains. Nobody else has rerun anything.

## Method

Vacancy swaps carry a letter with its type, so the cyclic P/B word along the occupied sites is conserved. In squeezed coordinates (`sqed.py`) a state is positions ⊗ labels with a word offset, and exactly H = −t(T_bulk + B ⊗ R + h.c.) + J Σ_j A_j h_{j,j+1} − μ_v M, with R the cyclic label shift. `scan.py` screens every word of runs B Pᵏ (k ≤ 4, plus the all-port word) by first-order energy and diagonalizes all competitive (M, word, 2Q) sectors in all momenta. `separated.py` builds the first-order spectrum, free charge states at flux 2πr/M plus Jα times frozen-ring levels with R eigenvalue e^{2πir/M}. `modes.py` reads velocities off the density structure factor. `run_fes.py` runs VUMPS (`../d2_followup_220_1282/vumps.py`) started in the P B P sector, since a random start falls into a metastable all-port sector.

## Results

**Ground-state degeneracy.**

| L | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|
| frozen | 10 | 3 | 14 | 8 | 6 | 10 | 22 | 3* |
| mobile | 2 | 1 | 1 | 2 | 2 | 1 | 2* | running |
| word (M) | BPP (3) | BPPP (4) | BPPP (4) | BPBPP (5) | BPBPP (5) | BPPBPP (6) | BPPBPPP (7) | |

Every 2 is the Q = ±½ pair. Entries marked * come from the best first-order sector in all momenta, with competitors (0.007 and 0.16 higher at first order) still running. J = 1, 2 and 4 give the same pattern on L = 6 to 10. The mechanism is the Ogata–Shiba gluing. Translates of a word are exactly degenerate inside the payload, each R eigenvalue e^{iθ} acts as a charge flux θ, and the hard-core charge fluid picks θ = 0.

**What becomes of the period.** The frozen triplet's partners sit at flux ±2π/3. Their gaps match free hard-core bosons to 4 × 10⁻⁴ (0.13446 against 0.13451 at L = 10). With x = L·gap/(2πv_c) we find 0.1263, 0.1165, 0.1160, 0.1155, 0.1151, 0.1147 at L = 5 and 9 to 13, against free values 0.1251, 0.1152, 0.1144, 0.1138, 0.1134, 0.1131. These are 1/L conformal levels of a charge boson compactified by the three-site cell.

**Spin–charge separation** (L = 10, M = 6). Largest deviation of the six lowest levels in each momentum except π, where degenerate charge states mix:

| J | 0.05 | 0.1 | 0.25 | 0.5 | 1.0 |
|---|---|---|---|---|---|
| max \|E_ED − E_pred\| | 9.8e-5 | 3.9e-4 | 2.3e-3 | 9.0e-3 | 3.4e-2 |
| max gap error | 3.8e-5 | 1.5e-4 | 9.0e-4 | 3.5e-3 | 1.3e-2 |

Doubling J multiplies the error by 3.7 to 4.0.

**Velocities and z.** The density probe at q = 2π/L excites one level (weight at least 0.94), giving v_c = 1.762, 1.887, 1.677, 1.846, 1.925, 1.948, 1.937 at (L, M) = (5, 3), (6, 3) and (9, 6) to (13, 6). Each is 1.0–1.5% below the free value, so the gaps scale like the free ones and z = 1. The logical-charge probe peaks at 0.035 (L = 10), mostly the charge twist (2πv_c/L)/M² dragged along by the payload's lowest current state. The payload's own share is about 1.6 × 10⁻³.

**c.** Casimir pairs E₀/L = e∞ − πC/(6L²) at fixed density give C/v_c = 1.075 (L = 5, 10) and 1.051 (L = 6, 12). Free hard-core bosons at the same sizes give 1.078 and 1.053, so the calibrated charge c is 0.9975 and 0.9983. The payload's expected share, c_payload v_p/v_c ≈ 0.003, is below this resolution. FES, with c from S = (c/6) ln ξ over the last four χ:

| chain | χ | ξ (sites) | c |
|---|---|---|---|
| XX | 8 to 96 | 11 to 345 | 0.99 |
| two XX, velocity ratio 0.3 | 8 to 128 | 3.8 to 46 | 1.95 |
| two XX, velocity ratio 0.05 | 8 to 128 | 5.3 to 65 | 2.14 (pairwise 1.0 to 4.6) |
| d = 11, ρ = 0.61, J = 0.25 | 18 to 32 | 12 to 33 | 1.14, gradient still 10⁻⁴ |
| d = 11, ρ = 0.92, J = 0.5 | 18 to 48 | 25 to 87 | 1.02 (pairwise 0.9 to 1.4) |

VUMPS stays in the P B P sector (n_B/n = 0.333) but converges slowly, as the payload's energy scale is below 10⁻³. Balancing truncation costs v/ξ² suggests ξ_payload ≈ √r ξ_charge for velocity ratio r [estimate], a few sites at r = 1/360 and a few port-dressed cells at the dilute point (r ≈ 1/36). So the full-chain FES sees only the charge boson.

**Short rings carry a defect.** One defect cell costs Δ_def = Jα × 0.070 (a B P run) or Jα × 0.086 (a B P³ run), from the demonstrator's run energies, about 0.008 here. Moving M to a multiple of 3 costs up to πv_c/L, and ED confirms that M follows ρL whatever its residue. The M offsets of a defect word are exactly degenerate, so their twist levels (2πv_c/L)(r/M)² form a slow band. Defects disappear only for L ≫ L* ≈ πv_c/Δ_def ≈ 700, which grows like t/J.

**d = 24** (`d24.py`, bond23 at n = 1, cell length 2, λ = 0.2, J = 0.25, μ_v = 0). Arrow legality of every squeezed bond is conserved, so the legal sector is exact. The checked working point of `../explicit_candidate/` has logical hopping 2 × 10⁻⁷λ², too small to resolve, so I used D_x = D_M = 0.5, g = 0.3, Δ = 0.4 (code still the cell ground space, gap 0.2, |T| = 1.0 × 10⁻³). The frozen chain is 2-fold at L = 4, 6 and 5-fold at L = 5. The vacancy chain is unique on L = 5 to 8 (words DDM, DDM, DMDM, DMDM). The would-be translates sit at charge-twist energies, 0.279 and 0.241 at L = 5, 6 (x ≈ 1/9) and 0.444 at L = 7 (x ≈ 1/4, two two-site cells), with only payload levels near 10⁻⁶ below. The payload is about 10⁵ times slower than the charge, so c is out of reach, and at n = 23 the logical bandwidth is about 10⁻¹⁵².

## Window

- J/t: the ground-state structure is unchanged from J = 0.25 to 4 on L ≤ 10, and separation holds to O(J²/t) from J = 0.05 to 1. Small J costs a payload velocity ∝ J and a defect length L* ∝ t/J.
- μ (demonstrator): to first order the payload is the frozen chain times Jα, so its selection window 3λ ≲ μ < 3ε is unchanged. μ_v only sets ρ.
- Density: α(ρ) → 0 as ρ → 0 decouples the payload, and v_c → 0 as ρ → 1. The pattern harmonic e^{−2isφ/3} at lattice momentum 2πsρ/3 has dimension s²K/9. It is resonant only for s = 4 at ρ = 3/4, with dimension 16K/9 ≈ 1.8, so it is relevant and locks a period-4 crystal P B P 0. Harmonics dressed by the staggered payload operator (ν = s + ½, dimension ν²K/9 + 0.91) are all irrelevant, the closest being ρ = 6/7 at 2.27. So c = 2 should hold for 0 < ρ < 1 away from 3/4, at the physics standard and conditional on bosonization.

## Verdict

The frozen payload's period becomes harmless. Vacancies glue every translate and defect position to a charge flux, so each ring has one ground state (up to the payload's Q = ±½ doublet), and the frozen triplet survives only as x ≈ 1/9 twist levels. What this cheapest completion cannot deliver numerically is the payload's half of c = 2, since its velocity is about 1/360 of the charge velocity and rings below L* ≈ 700 can keep a pattern defect. At d = 24 both hierarchies are astronomically worse, so c = 24 has to rest on the infrared argument rather than on simulation.

## Reproducing

Only numpy and scipy are used. Scripts are `brute.py` (2 min), `scan.py`, `gs_quick.py` (no momentum projection, L ≥ 13), `separated.py`, `modes.py`, `run_fes.py`, `d24.py`, and `analyze.py` (writes `data/summary.json`). Typhon arrays 17517900 (scans), 17517901 (L = 15, 16) and 17517984 (L = 13 to 15), with 17517902 cancelled. FES, L ≤ 10, the J scans and d = 24 ran on the login node at nice 19 (`tasks_login.txt`). Mobile L = 12 to 16, the frozen L = 12 competitors and the ρ = 0.92 FES beyond χ = 48 were still running when this was written. d = 24 at L = 9 needs a legality-aware enumeration (product space 8.5 × 10⁸).

## Relation to the rest of the repository

This addresses item 3 of STATUS "Where help is most useful" for the cheapest mobile completion and extends `../mobile_cells/` from a toy to the route-A demonstrator. It supports the gluing assumed in `../counterexample_campaign/doping/README.md` §3, and adds the defect length L* and the payload-velocity hierarchy to the infrared transfer still open there.
