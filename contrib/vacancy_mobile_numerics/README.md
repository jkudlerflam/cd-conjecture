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

Update (2026-10-28, evening, section "Late results and calibrated c" at the end). On the logical vacancy model, an effective model of the vacancy lift with a resolvable payload, c = 2.01 ± 0.02 (ring entanglement) and 1.9 to 2.0 ± 0.2 (FES) at n = 1, each calibrated and each separated from a c = 1 control. At n = 2 neither estimator delivers a calibrated c = 3 at the sizes reached.

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
| frozen | 10 | 3 | 14 | 8 | 6 | 10 | 22 | 3 |
| mobile | 2 | 1 | 1 | 2 | 2 | 1 | 2 | 2† |
| word (M) | BPP (3) | BPPP (4) | BPPP (4) | BPBPP (5) | BPBPP (5) | BPPBPP (6) | BPPBPPP (7) | BPPBPPP (7) |

The L = 11 and 12 entries were completed from the late Typhon results (section "Late results and calibrated c", below), and † marks a partial sector search. Every 2 is the Q = ±½ pair. When this was first written, the L = 11 and 12 entries came from the best first-order sector in all momenta, with competitors (0.007 and 0.16 higher at first order) still running. J = 1, 2 and 4 give the same pattern on L = 6 to 10. The mechanism is the Ogata–Shiba gluing. Translates of a word are exactly degenerate inside the payload, each R eigenvalue e^{iθ} acts as a charge flux θ, and the hard-core charge fluid picks θ = 0.

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

## Late results and calibrated c (2026-10-28, evening)

Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, same day, answering the referee's request in `../referee_d24/` for a calibrated c = n + 1 for the vacancy lift at n = 1 and 2. No Codex call (the work is repository-bound numerics). Nobody else has rerun it.

### 1. The Typhon jobs of the first pass

- [numerical, ED] Mobile L = 11. The competitor word BPBPBPP finished, 0.0081 above the BPPBPPP ground state, so the entry is 2 (the Q = ±½ doublet at zero momentum). The second first-order competitor (0.16 above at first order) was not rerun.
- [numerical, ED, partial] Mobile L = 12. Only the best first-order word BPPBPPP (M = 7) finished, in all twelve momenta. Its ground state is the Q = ±½ doublet at zero momentum, with the next level 0.0033 above. The other words hit the 6 h limit, hence 2†.
- [numerical, ED without momentum projection, two candidate words] Mobile L = 13: M = 8, word BPBPPBPP, 2Q = 0, unique with zero momentum, next level 1.1 × 10⁻³ above. The winning word carries a B P run, a defect as in item 6. Mobile L = 14: M = 9, word BPPBPPBPP, Q = ±½ doublet at zero momentum, next pair (k = ±2π/14) 0.0115 above. The second candidate ran out of memory (18 GB Lanczos basis). L = 15 is still running (job 17517901_0), and L = 16 ran out of memory.
- [numerical, ED] Frozen L = 12 is 3-fold (table updated). The competing words lie 0.156 and 0.160 above. Frozen L = 13, partial: 26-fold inside BPPBPPBPPBPPP (all 13 momenta times the Q doublet), with 2Q = 3 0.0073 above, other words cut by the time limit. Frozen L = 15, partial and still running: 6-fold inside (BPP)⁵ at k = 0, ±2π/3 times the Q doublet, the expected threefold breaking. Frozen L = 14 finished only one sector.
- [numerical] Charge velocities from the density probe: 0.839, 1.367 and 1.905 at (L, M) = (7, 6), (8, 6), (14, 6), 1.611 at (13, 9) and 1.763 at (14, 9) (the last from the log, job still running). They lie 0.0 to 1.5% below the free values (0.839, 1.378, 1.934, 1.630, 1.787), as before. Twist-level x = 0.1179, 0.1172, 0.1144, 0.1143.
- [numerical] FES calibration, two XX chains at velocity ratio 0.01, now to χ = 160: last-four fit 2.03, pairwise slopes 1.4 to 3.6.
- [numerical, negative] ρ = 0.92 FES (`vac_P2`) at χ = 64 and 96: ξ = 88.5 and 122, S = 2.153 and 2.317, but the gradient is 3 × 10⁻⁷ and 7.8 × 10⁻⁶ with the iteration cap reached, and the pairwise slopes are 34.6 and 3.1. This run is not converged, so the 1.02 quoted above from χ ≤ 48 should not be read as a measurement. χ = 128 is still running on the login node.
- [numerical, negative] d = 24 beyond L = 8 did not run. Frozen L = 7 lost ARPACK convergence in one sector, and frozen L = 8, 9 and vacancy L = 9, 10 need 51 GB to 285 GB index arrays in `sqed.Internal`. The legality-aware enumeration is still needed.

### 2. Why an effective model

The real lift cannot be measured at any n. At n = 1 the d = 24 chain's payload is about 10⁵ times slower than the charge (above), and at n = 2 the bond's logical hopping is 1.1 × 10⁻⁸λ² with a band ratio λ_min/λ_max = 1.4 × 10⁻² (`../explicit_candidate/README.md` §2.6). So c is measured on an effective model in which the payload is resolvable, as the referee asks.

**Model LV_n** (`logical_vac.py`). A site is a vacancy or a logical cell carrying n Jordan–Wigner fermion flavours (d = 1 + 2ⁿ, so d = 3 and 5). The payload bond is the Gaussian part of the logical chain of `../explicit_candidate/README.md` §2.4,

h_log = −Σ_f τ_f (c†_{j,f} c_{j+1,f} + h.c.) + ½ Σ_f ε_f (n_{j,f} + n_{j+1,f}),

ordered by cell and then flavour, so the bond carries the intra-cell Z strings. It is lifted exactly as in `../counterexample_campaign/doping/bond24.py`: −t vacancy swaps, J h_log on adjacent occupied pairs only, −μ_v per occupied site. Vacancies carry no modes, so the swap needs no fermion sign. M and every flavour number are conserved.

[heuristic] Relation to the real lift. To first order in J/t the lift of any frozen chain is free hard-core charge plus Jα H_frozen (Method, and item 3), and the low-energy H_frozen of bond23 at member n is the logical chain. LV_n keeps the n flavours with their strings, the vacancy dynamics and the bond-only payload coupling. It differs in three ways. Each b = n + 1 site cell is one site (this changes the charge compactification, not c). T is diagonal, τ = (1, 0.5) instead of the bond's T, which makes the flavour numbers separately conserved and so forbids inter-band Cooper transfer by symmetry, where the bond uses distinct centres instead. The O(g⁴) non-Gaussian part is dropped.

Couplings: t = 1, J = 0.5 (and 0.25 at n = 1), so the payload is slower than the charge, the referee's ordering. [estimate, first order] With α = ρ − sin²(πρ)/(π²ρ) = 0.553 at ρ = 2/3, the payload-to-charge velocity ratios are 2Jατ_f sin(πν_f)/(2ρ sin πρ), that is 0.48 for τ = 1 at half filling and 0.24 for τ = 0.5. [numerical] The ED adjacency gives α = 0.555 to 0.559 for L = 9 to 18, within 1% of the first-order value.

[numerical] Checks. `logical_vac.py check` compares the sector Hamiltonian with the literal d^L Hamiltonian on rings (n = 1, 2, with fields and Ising terms on), difference 0 and no leakage. At n = 2 (L = 12, M = 8, N₀ = 4) the energy steps for adding pairs of the slow flavour, N₁ = 0 → 2 → 4 → 6, are −0.5323, −0.2138, +0.2138 in ED against −0.5171, −0.2142, +0.2142 for a free band times Jα, so the slow flavour separates as a free gapless band at this size (agreement 0.2% to 3%).

**Estimators.** (a) Ring interval entanglement [numerical, ED]: S(ℓ) of the ring ground state in a fixed (M, N_f) sector. The main estimator is the half-ring pair c = 3[S(L₂/2) − S(L₁/2)]/ln(L₂/L₁) at equal density and fillings. It is calibrated by the same quantity for the exactly solvable reference, a sum of free-fermion (XX) rings with the same L and particle numbers, one per gapless mode (M particles for the charge, N_f per flavour, which reproduces the physical Friedel wavevectors), through c_cal = m c_raw/c_ref with m the number of reference modes. A per-L fit of S(ℓ) to (c/3) ln[(L/π) sin(πℓ/L)] plus Friedel terms is also recorded in `data/lv_summary.json`, but the control below shows it is biased upward by up to 0.9 at these sizes, so I do not use it. (b) FES [numerical, VUMPS, `run_fes.py --model lv`]: c from S = (c/6) ln ξ fitted over χ = 16 to 96, calibrated the same way by decoupled XX chains with the first-order velocity ratios (`run_fes.py --model kxx`) over the same χ, with the fit's standard error propagated. The two-chain references themselves scatter between 1.92 and 2.23 across four-χ windows, so I add 0.1 per mode as a systematic to the FES numbers.

**Controls.** One flavour is gapped by an Ising term, an XXZ payload with Δ = vz/(2τ) = 3 (Néel). Expected c = 1 at n = 1 and c = 2 at n = 2. This is a stringent control, since the gapped payload still entangles with the charge through the parity of the number of cells in an interval [heuristic].

**Results, n = 1** (target 2).

| model (LV_1) | estimator | sizes | c |
|---|---|---|---|
| J = 0.5, ρ = 2/3 | half-ring pairs | L = (6, 9), (9, 12), (12, 15), (15, 18) | 2.024, 2.021, 2.009, 2.007 |
| J = 0.5, ρ = 1/2 | half-ring pairs | L = (8, 12), (12, 16) | 2.037, 2.016 |
| J = 0.25, ρ = 2/3 | half-ring pairs | L = (12, 15), (15, 18) | 2.013, 2.010 |
| J = 0.5, ρ = 0.668 | FES | χ = 16 to 96, ξ ≤ 40 | 1.92 ± 0.06 (fit), ± 0.2 (systematic) |
| J = 0.25, ρ = 0.629 | FES | χ = 16 to 96, ξ ≤ 47 | 2.04 ± 0.09, ± 0.2 |
| control, Néel payload, ρ = 2/3 | half-ring pairs | L = (9, 12), (12, 15), (15, 18) | 1.77, 1.64, 1.53 (falling) |
| control, Néel payload, ρ = 0.663 | FES | χ = 16 to 128 | 1.05 ± 0.01 |

[numerical] So c = 2.01 ± 0.02 from rings (the spread of the pairs, at two densities and two J, all from L ≥ 6 and flat in L) and 1.9 to 2.0 ± 0.2 from FES, while the c = 1 control is separated by both estimators, cleanly by FES and by a downward drift in the rings. The per-L ring fits give 2.06 to 2.14 for the model and 1.81 to 2.00 for the control, which is the bias mentioned above.

**Results, n = 2** (target 3). [numerical, negative] Not established.

| model (LV_2) | estimator | sizes | c |
|---|---|---|---|
| both flavours half filled, ρ = 2/3 | half-ring pair | L = (6, 12) | 3.06 |
| control, Néel second flavour | half-ring pair | L = (6, 12) | 2.78 |
| fillings 1/2 and 0.37, ρ = 0.648 | FES | χ = 16 to 64 | 2.14 ± 0.17 (fit), ± 0.3 |
| both flavours half filled, ρ = 0.668 | FES | χ = 16 to 96 | 2.06 ± 0.16, ± 0.3 |
| control, Néel second flavour, ρ = 0.79 | FES | χ = 16 to 96 | 1.91 ± 0.12, ± 0.2 |

Two obstructions. Closed shells on rings need even flavour numbers, so at fixed fillings only L = 6 and 12 (ρ = 2/3) are reachable, a single pair, and the same pair puts the control 0.78 above its target. Both flavours at half filling are also a poor choice, since 2k_F0 + 2k_F1 = 2π makes an inter-flavour umklapp marginal [heuristic, bosonization], so the generic fillings (1/2, 1/4) were run as well, but only per L (L = 10 to 14, per-L c_cal = 3.13 to 3.18, against 2.86 to 2.96 for the control). FES with three velocities calibrates (the 1 : 0.48 : 0.22 and 1 : 0.48 : 0.24 references give raw 3.7 to 4.0 over the same χ), but on the model it resolves only about two of the three modes up to χ = 96, the same failure as for the n = 2 chain of `../routeA_n2_chains/`. The separation check in §2 indicates that the slow flavour is a gapless free band at L = 12, so I read the FES deficit as under-resolution of the slowest mode rather than a gap [heuristic].

### 3. What remains

(i) n = 2: FES at χ of several hundred (the generic-filling run and its reference are still climbing on the Typhon login node, `tasks_lv_fes.txt` lines 9 and 10), or finite DMRG on open chains of 100 to 200 sites for the interval estimator, which is velocity independent. Ring ED cannot go further (L = 16 at fixed fillings is 6 × 10⁷ states and was queued, job 17522036_7). (ii) A version of LV_n with the bond's own T and on-site μ = iA (non-commuting, U(1) only), which needs complex sector ED, and with the band ratio of the real bond. (iii) The step from LV_n back to the real lift rests on first-order separation (item 3 and the n = 2 check above) and is heuristic. (iv) The scans of §1 marked as running were still running when this was written.

### Reproducing (this section)

`logical_vac.py ring ...` (task lists `tasks_lv_ring.txt`, local, and `tasks_lv_ring_big.txt`, Typhon login node at nice 19, up to 20 min and 30 GB each; the batch script now changes to the submission directory), `run_fes.py NAME --model lv|kxx ...` (`tasks_lv_fes.txt`, `tasks_lv_fes_ctrl.txt`; Typhon array 17522020 and login-node copies with suffix L, which `analyze_lv.py` merges by taking the longer run), then `analyze_lv.py` (writes `data/lv_summary.json`). The FES calibration table above is the χ ≥ 16 common-window fit printed by `analyze_lv.py`. Compute-node jobs were slowed by the per-user job limit, which is why most runs moved to the login node.
