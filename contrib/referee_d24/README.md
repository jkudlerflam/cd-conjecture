# Referee report on the d = 24 vacancy counterexample

## Claim

Under review is `contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md`. It adds one vacancy state to the explicit 23-letter route-A bond h_23 of `contrib/explicit_candidate/`, with species-blind vacancy swaps of amplitude −t, a copy J h_23 acting only on physically adjacent occupied pairs (on-site halves included), and an occupied-number chemical potential μ. The chain is Hermitian, nearest-neighbour and one-site translation invariant on (C^24)^⊗L with periodic boundary conditions. At C33's certified n = 23 point, with nonzero g, then λ, then J/t chosen small in that order and fixed before L → ∞, it argues that the low-energy theory is 23 payload compact bosons plus one charge boson. The total Virasoro c would then be 24 > d − 1 = 23, contradicting Latorre and Sierra (arXiv:2403.17242).

## Status

Refereed. I found no fatal error, and the vacancy step works at the physics standard, but the claim is not yet an established counterexample. Two things would upgrade it, a calibrated measurement of c = n + 1 for the vacancy lift of this bond family at n = 1 or 2, and a check of the C33 certificate and of the matching inputs by someone outside the pipeline that produced them.

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-28. No Codex call, since every check is a short script.
- Read: the claim with `doping/`, `doping_check/`, `infrared/` and `route_a/`, then `explicit_candidate/README.md`, `codex/work_C33/`, `codex/work_C37/`, `redteam_routeA_cells/`, `mobile_cells/`, `codex/work_C20/RESULT.md` and `ir_control/`. I formed the verdict before opening `codex/work_C41/` and `contrib/vacancy_mobile_numerics/`, and compared afterwards.
- Independent checks: the bond rebuilt from prose, the quoted constants recomputed, the resonance rule rederived, the C33 certificate rerun, and a toy analogue studied by exact diagonalization.

## Report

### 1. Is the bond as stated?

Yes. `rebuild_bond.py` writes h_23 from the term list of `explicit_candidate/README.md` §1.2 and the lift from §1 of the claim, importing nothing from `bond23.py` or `bond24.py`. At the fixture's parameters it agrees with `doping/bond24_fixture.npz` entry by entry (largest difference 0, 923 nonzero entries in both). It is Hermitian and commutes with occupied number, internal charge and marker number, also at μ ≠ 0. Every off-diagonal entry of h_23 preserves each site's marker bit, so marker positions in squeezed coordinates are conserved. The n = 23 parameter list agrees with C33's certified point. The quoted constants reproduce (γ_24 = 7.5667157e−5, α(1/2) = 1/2 − 2/π², log10 d_24 = 72.0416, 2tδ_4/(πb) = 8.371e−9 t, δ_4/(2b) = 6.5746e−9). `rerun_c33.py` regenerates C33's 512-bit loop matrices bit for bit, and the root certificate passes again (46 roots, δ_2 > 0.0140194, δ_4 > 3.1559e−7). This is a rerun of the collaborators' code, not an independent derivation.

Two caveats. The fixture uses moderate parameters, so it checks structure, not the n = 23 point. More importantly t, J, μ and the final λ are never given. What is defined is a family "for sufficiently small couplings" with uncomputed constants, not one Hamiltonian.

### 2. Does the vacancy remove the frozen-cell objection?

In the letter of the conjecture, yes. The squeezed-coordinate identity is exact (I checked the same identity for a toy against brute force). The cyclic marker word is conserved, and in the perfect period-24 sector the cut hop permutes the 24 offsets cyclically. They become charge flux branches K = (P + 2πℓ)/24 whose sea energies differ at O(t/L), not 24 disconnected copies. The ground state is translation invariant, marker order is algebraic with exponent 2K_c/b² ≈ 1/288, and every sufficiently long ring admits an occupied number divisible by 24. Uniqueness on every ring is not proved (the campaign's own `infrared/README.md` §7 and C41 say so), but finite accidental degeneracies do not change c.

The squeezed period survives intact in material coordinates, as a superselection label, and the payload remains a chain on M = N/24 cells. The vacancy does not wash the period out. The charge fluid re-imprints any period-24 structure of the bond onto the payload at relative order J/t (toy below). That is harmless here only because the payload's unit cell is the material cell, so the imprint is a small cell-periodic correction. For a payload with a shorter period it would be a relevant perturbation. The claim should say so.

Rounding N to a multiple of 24 costs O(t b²/L) and beats a wrong-length cell, which costs J α Δ_def. Two consequences go unstated. The charge sector is a "charge-24" liquid, since adding one particle forces a defect at finite cost and only 24-particle composites are gapless. This is consistent with c = 1 (it is the compactification change the gluing implies). And the crossover length is L* ~ t b²/(J α γ_24). Below it incompatible rings carry frozen defects. With the hierarchy of §3, L* exceeds 10^517 sites.

A fixed μ does select ρ = 1/2 in the thermodynamic limit. The relaxed compressibility is positive, so ρ(μ) is continuous and increasing and μ* exists, though its value is unknown. Exact half filling is not needed anyway. The dimension argument holds on an open window (ρ < 1/√2 in `infrared/README.md`, 2/5 < ρ < 3/5 with stiffness margins in C41). Only the free-charge threshold theorem uses ρ = 1/2, and the conclusion does not need it.

### 3. Is the infrared argument acceptable?

The operator algebra is right, and I rederived it. A payload vertex of cell momentum q acquires exp[iρ(q + 2πℓ)x/b] exp[−i(q + 2πℓ)φ_c/(πb)] in physical coordinates, with charge dimension [(q + 2πℓ)/(2πb)]² K_c. A uniform term at ρ = 1/2 needs q + 2πℓ = 4πbh. Either h = 0, which is the payload's own zero-momentum selection, or the charge part alone has dimension 4h²K_c = 4 + O(J/t). With C33's margins every uniform non-Gaussian vertex has dimension at least three, and the detuned marginal ones are cut off at δ_4/(2b). K_gap = πv/e'' = 1/ρ² = 4, α, e_gap and v_gap are correct.

The order of limits is consistent, since each parameter is bounded above by earlier ones and none below by a later one, and the Bravyi–DiVincenzo–Loss local Schrieffer–Wolff step applies to the blocked frozen cells. What the hierarchy demands is extreme. C33's slowest payload speed is 10^−431.7 λ² in internal units. Controlling BDL's generic second-order remainder against it, which is what the claim invokes, needs λ near 10^−435. Even if every data-dependent higher-order term keeps the express suppression of T, as the red team found at n = 2, 3, λ ≲ 10^−31 is needed. Inequality (10) of `route_a/DEFECT_TRANSFER.md`, with M ~ D_x = 1000, then puts J/t below a bound somewhere between 10^−510 and 10^−1316 depending on λ, and the slowest physical speed is 10^−1000 t or smaller. None of this is inconsistent, but no numerical check of the n = 23 chain will ever be possible. The conclusion rests entirely on asymptotic arguments.

Those arguments use two inputs that are assumed, not derived. P is that the frozen chain at the certified point has c = 23 (bosonization of 46 nonresonant Fermi points plus local elimination, C40's "M2"). M is the Wilsonian matching of the gapless charge fluid onto J α H_23 plus one Gaussian boson, with small normalized remainders and positive dressed defect costs. The campaign's own `infrared/README.md` tags M "[assumed, sufficient physics input]" and the microscopic conclusion "[heuristic]", C41 agrees, and the synthesis relabels this a "controlled physics argument". Both are standard effective-field-theory inputs, and I found no mechanism that defeats either. Neither has been tested on this bond family at any n. `explicit_candidate/README.md` records that no central charge has been measured, and `routeA_n2_chains/` reports that finite-entanglement scaling fails its calibration for the n = 2 "spatial" chain.

### 4. Anything fatal?

Nothing fatal. What they missed is items 4 and 6 to 9 of the table. The review chain is also not independent, since the three "independent" workers are GPT-6-family agents under one coordinator and C41 was unfinished when cited.

### Ranked problems

| # | Problem | Rank |
|:--|:--|:--|
| 1 | Input P (c = 23 for the frozen chain at the certified point) is inherited and untested numerically for this bond at any n | Serious |
| 2 | Vacancy matching M is assumed; the synthesis upgrades the campaign's own "[assumed]" and "[heuristic]" tags to "controlled" | Serious |
| 3 | No explicit λ, J/t or μ, and uncomputed matching constants (C41: no justified interval for J/t) | Serious, fixable in principle |
| 4 | λ ≲ 10^−31 at least, J/t ≲ 10^−510, payload speeds ≲ 10^−1000 t, L* > 10^517, so nothing is checkable | Serious as evidence, not a logical gap |
| 5 | Reviews not independent; no human check of the C33 and C37 certificates | Serious for acceptance |
| 6 | Unique ground state on every large ring not proved | Fixable, not needed for c |
| 7 | Charge-24 liquid and crossover length L* not stated | Fixable |
| 8 | Exact half filling assumed but unnecessary; restate on a density window | Fixable |
| 9 | The vacancy re-imprints the cell period on the payload at O(J/t), harmless only because the payload period is the cell | Cosmetic, should be stated |
| 10 | Fixture at moderate parameters; "[proved]" tags on matching steps | Cosmetic |

### A small analogue

[numerical] `toy_ed.py` treats a d = 1 + 2b chain in exact squeezed coordinates, checked against brute force at L ≤ 6 with zero matrix error and zero leakage. A letter is a colour c ∈ Z_b and a spin-1/2. On adjacent occupied pairs only, the bond has a penalty K when colours fail to advance, an XXZ payload, and an oriented colour-pair energy η. Colours never change, so this is a frozen period-b crystal with a c = 1 payload, lifted as in the claim.

At b = 2 (d = 5) with t = 1, J = 0.5, K = 10, Jxy = 1, Jz = 0.5, η = ±0.3 and one μ fixed from L = 12, the grand-canonical ground state on every ring L = 7 to 16, primes included, lies in an even-N sector and is non-degenerate with zero momentum. The best odd-N (defect) sector lies 0.25 to 0.69 higher. The two-particle curvature E(N+2) + E(N−2) − 2E(N) falls like 22/L (2.74 and 1.84 at L = 8 and 12). At b = 3 (d = 7, N = 6, L = 11 to 14) the legal sector again wins, but every ground state is a ±2π/L time-reversal doublet. That is the flux degeneracy C41 warns about. Uniqueness depends on the payload's cell momentum and the number parity, not on the mechanism.

The imprint test uses b = 2 with a colour-blind payload. In the frozen chain the squeezed ⟨S^zS^z⟩ alternation between AB and BA bonds vanishes identically for any η, which is then a constant. With vacancies it equals (0.51 to 0.52) Jη/t at L = 12 for Jη ≤ 0.1, and it grows with L (0.017, 0.026 and 0.033 at L = 8, 12 and 16 for J = 0.1, η = 0.5), as a relevant dimerization of a period-1 payload should. All numbers are in `toy_summary.json`.

I did not measure c. A faithful analogue needs a payload whose period equals the crystal's at incommensurate filling, which puts the rings needed by two calibrated estimators beyond exact diagonalization. `contrib/vacancy_mobile_numerics/` (in progress, read after drafting) attempts this with VUMPS for the d = 11 demonstrator lift. Its small-ring scans already show defect-carrying words winning at L ≤ 9, as the L* argument predicts when the defect cost is small.

### Verdict

[assessment] The vacancy completion is a sound and economical way to remove the frozen-cell objection, and its operator analysis is correct. The d = 24 claim therefore reduces to one question, whether the frozen 23-letter chain at C33's point has c = 23 at fixed small coupling. That is answered only by standard but untested effective-theory reasoning spanning a hierarchy of at least a thousand orders of magnitude, assembled and reviewed entirely by AI agents. A careful theoretical physicist should treat this as a credible, well-argued candidate, not as an established counterexample to c ≤ d − 1. I would accept it at the physics standard after (i) a calibrated c = n + 1 for the vacancy lift of the same bond family at n = 1 or 2, at couplings where the payload is resolvable and with the same ordering, and (ii) a human check of the C33 certificate and of inputs P and M.

## Reproducing

numpy and scipy only, run from this folder with `python3 -B`.

- `rebuild_bond.py` (one second) writes `rebuild_bond.json`.
- `rerun_c33.py` (under a minute) copies `codex/work_C33` scripts to a temporary directory, reruns `full_loop_integer.py` and `certify_roots.py`, and writes `c33_rerun.json`.
- `toy_ed.py` runs the brute-force self-test. `toy_runs.py sectors L...`, `sectors_b b lo hi L...` and `dimer L J eta` append to `results.jsonl` (minutes locally up to L = 16), and `analyze_toy.py` writes `toy_summary.json`.

## Relation to the rest of the repository

This is the referee that `contrib/routeA_final/` lists as pending for Form V. It does not support changing "Where things stand" to "counterexample found". It supports recording the d = 24 vacancy chain as the leading candidate, conditional on inputs P and M, with the two upgrades above as the next steps.
