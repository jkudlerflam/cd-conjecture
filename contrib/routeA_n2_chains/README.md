# The n = 2 logical chains of four route-A cells at finite coupling

## Claim

Class: the n = 2 logical chains of C28, C29, `contrib/routeA_spatial/` and `contrib/routeA_compile/`. Each is a translation-invariant nearest-neighbour chain on C⁴ (two logical fermion modes per cell) with exact total U(1), built from the complete order-λ² bond at a fixed conversion coupling g or τ, every generated term kept. The question is the Virasoro c of the ground state.

1. [numerical] Wherever both runs converge, VUMPS on the interacting chain reproduces VUMPS on its own self-consistent Hartree–Fock (HF) chain. The bond entropy agrees to 2 × 10⁻³ (0.02 for spatial) at every χ up to 128 to 160, (ξ up to 30 to 70 cells), as do densities, measured Fermi points and even finite-χ artifacts. The quartic couplings at the Fermi points are 10⁻⁴ to 0.06, and 0.14 at most (C29 at g = 1). So at reachable scales, c is the number of Fermi seas of the HF chain.
2. [numerical, HF level] Every departure from c = 2 comes from the generated quadratic terms. With the sources' fixed fields, C28's narrow band fills completely for g ≳ 0.15 (c = 0), and at g = 1 C28 has a single Fermi sea at ν = 0.45 and 0.6 (c = 1). C29's source fields put ν = 1.03, so its two flavour bands cross just below the Fermi level. The finite-g inter-flavour hopping hybridizes them, and at g = 0.5 and 1 HF opens a gap at ν = 1 (c = 0), which VUMPS and open-chain DMRG confirm at g = 0.5.
3. [numerical, exact to 10⁻¹³] The compile bond has a unitary particle–hole symmetry, staggered between cells. So μ_c = 0 is exactly ν = 1, where the two-band umklapp R₀R₁ ← L₀L₁ is exactly resonant. Its coupling is 10⁻⁴, which is invisible numerically, but it is the wrong operating point.
4. [heuristic, one loop plus margins] With one uniform field retuned (C28 to ν = 0.45, C29 to ν = 1.2, compile to ν = 0.9), every chain tested has four Fermi points except C28 at g = 1. Every 2 → 2 resonance is detuned by at least 0.01, and the couplings are at most 0.14. Standard weak-coupling RG then gives c = 2.

## Status

Numerical. FES passes its calibration for C28 and C29 and fails it for spatial and compile. The ring/Casimir estimator fails its calibration everywhere, and open-chain DMRG finished for only one chain. Three things would upgrade this: U(1)-symmetric iMPS at rational fillings, fixed-N open-chain DMRG, and χ well above 200 for spatial and compile.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28. No Codex call.
- Checks: the bonds were rebuilt from source code copied unchanged into `src/`. The C28 and C29 leading Fermi points and velocities match the sources to all printed digits, and compile's Fermi-sea centres match `switch_break_m2.json`. Every estimator was run on XX, on the sources' free chains and on the HF chain of every coupling. Not rerun by anyone else.

## Method

**Chains** (`bonds.py`). C28 g = 0.2, 0.5, 1 (field retuned by HF to the design density ν = 0.45). C29 g = 0.2, 0.5 at the source fields, and 0.5, 1 at ν = 1.2. Spatial g = 0.05 (working point), 0.1, 0.2 at the working field. Compile τ = 0.003 (source), 0.01 at μ_c = 0, and 0.01 at ν = 0.9. Each bond is divided by its largest hopping. Fermi points, resonance mismatches and Fermi-point couplings (C28 eq. 13) come from HF on the full bond (`bands.py`, `resonances.py`).

**Estimator 1, FES.** Damped VUMPS (`vumps.py`, copied from `d2_followup_220_1282`), χ = 16 to 192, c = 6 dS/d ln ξ over the last three or four converged χ. A two-cell unit cell is used where the one-cell uMPS is non-injective, and it also serves as the translation check. Interacting Fermi points come from the jumps of n(k), built from ⟨c†₀c_r⟩ with strings.

**Estimator 2.** The ring/Casimir estimator fails on the exactly solvable free parts at every reachable size (`ring_casimir_free.py`). At L = 20 to 40, minimized over flux, it gives 0.1 to 0.9 (C28, C29), 0.5 to 7 (spatial) and 25 to 330 (compile). XX gives 1.0 only at even L at half filling. Close Fermi-point pairs (δ₂ = 0.08 to 0.6) and slow velocities keep the two seas' zero modes from aligning. We substituted open-chain DMRG at L = 64 (`run_chain.py obc`) with a Calabrese–Cardy fit, compared with the exact entropies of the HF chain on the same open chain (`exact_open_table.py`). This exact calibration gives 1.90 to 2.25 for all four-point chains, except compile at 4.1 to 4.4 (2.7 at L = 256).

**z.** The energy error obeys e(χ) − e∞ ∝ ξ^(−p) with p = 1.93 for XX and for the C28 free chain (exact e∞), and p = 2.0 from free fits, so z = 1 there. The interacting chains track their HF chains, whose crossings are linear.

## Results

The second column gives the HF chain's number of Fermi seas, which is its exact c. FES lists the interacting chain (χ range, largest ξ in cells) and then the HF chain. OBC lists DMRG (L = 64, χ = 48) and then the exact HF value. "fails" means the estimator also fails on the solvable chain. Δ is the largest density difference between the cells of a two-cell unit cell at χ = 64 (interacting / HF).

| cell, coupling, field | HF seas (= exact c of HF chain), ν | FES full / HF | OBC full / HF exact | Δ |
|---|---|---|---|---|
| calibration XX | 1 | 0.99 (48–128, 574) | – / 1.03 | – |
| calibration C28 free | 2 | 2.02 (64–160, 70) | – / 2.11 | – |
| calibration C29 free | 2 | 1.98 (48–128, 43) | – / 1.92 | – |
| C28 g = 0.2, source field | 0 (band full) | – | – | – |
| C28 g = 0.2, ν = 0.45 | 2 | 1.89 (24–64, 31) / S equal to 4 × 10⁻³ | – / 2.06 | 10⁻⁶ / 10⁻⁶ |
| C28 g = 0.5, ν = 0.45 | 2 | 1.97 / 2.00 (two-cell, 16–64, 33) | – / 1.99 | 10⁻⁹ / 10⁻⁹ |
| C28 g = 1, ν = 0.45 | 1 | fails (non-injective, both) | – / 1.02 | – |
| C29 g = 0.2, source | 2 (pockets, δ₂ = 0.08), 1.03 | 1.73 / 1.73 (48–128, 33) | no fixed N, discarded / 1.93 | 10⁻⁷ / 10⁻⁹ |
| C29 g = 0.5, source | 0 (gap 0.08), 1.00 | 0.6 / 0.6, S saturates at ξ = 15 | −0.04 / 0.01 | 10⁻⁶ / 10⁻¹¹ |
| C29 g = 0.5, ν = 1.2 | 2 | 2.06 / 2.02 (24–48, 51) | – / 1.97 | – |
| C29 g = 1, ν = 1.2 | 2 | 2.01 / 2.01 (two-cell, 16–96, 66) | – / 1.90 | 10⁻⁷ / 10⁻⁷ |
| spatial g = 0.05, 0.1, 0.2 | 2, 0.92 to 0.95 | fails (ξ stuck at 11 to 19), S equal to 0.02 | – / 2.14, 2.12, 2.25 | ≤ 10⁻⁵, except 0.2 at g = 0.2 (both) |
| compile τ = 0.003, 0.01, μ_c = 0 | 2, 1 exactly | fails, S equal to 10⁻⁴ | – / 4.36 | 0.2 / 0.2 |
| compile τ = 0.01, ν = 0.9 | 2 | fails (narrow band left empty, both) | – / 4.13 | – |

The n(k) Fermi points at χ = 64 agree with HF to 0.03 (spatial g = 0.05 gives ±1.18, ±1.70, against ±1.18, ±1.72). For C29 at g = 0.2, ξ = 20 does not resolve its 0.08-wide pockets. The period-2 modulations (spatial g = 0.2 and compile) appear equally in the exactly uniform HF chains, so they are finite-χ artifacts.

**Resonances** (momentum mismatch, 0 = resonant; g_max is the largest dimensionless Fermi-point coupling):

| chain | pair transfer | umklapp 2π(ν−1) | other inter-pair | δ₂ | g_max |
|---|---|---|---|---|---|
| C28 g = 0.2 / 0.5 | 0.10 / 0.10 | 2.83 | 0.19 / 1.17 | 1.36 / 1.14 | 0.002 / 0.015 |
| C29 g = 0.2, source | 0.034 | 0.19 | 1.22 | 0.08 | 0.006 |
| C29 g = 0.5 / 1, ν = 1.2 | 0.051 / 0.15 | 1.25 | 0.61 / 0.16 | 0.60 / 0.55 | 0.04 / 0.14 |
| spatial g = 0.05 / 0.1 / 0.2 | 0.011 / 0.026 / 0.091 | 0.49 / 0.46 / 0.35 | ≥ 0.86 | ≥ 0.38 | 0.004 / 0.015 / 0.06 |
| compile, μ_c = 0 / ν = 0.9 | 0.055 / 0.033 | 0 / 0.63 | 0.45 / 0.78 | 0.20 / 0.37 | 10⁻⁴ / 5 × 10⁻⁴ |

Velocities, in units of the largest hopping, are 0.008 (C28), 0.6 and 1.8 (C29), 0.05 and 0.12 (spatial) and 0.0015 to 0.02 (compile). At g = 0.2 we reproduce C28's one-loop numbers (x = −0.0038, y = 0.0034, irrelevant). For spatial, x > −|y| with y = 0.004 to 0.06, so its nearly resonant pair transfer would grow if it were resonant. The detuning stops the flow after 10 to 90 cells, far below exp(1/y).

## Verdict

- **C28.** It keeps c = 2 at g = 0.2 and 0.5 once the chemical potential is retuned with g, by FES and the exact HF reference, with no translation breaking. With fixed fields it becomes a band insulator beyond g ≈ 0.15, and at g = 1 the generated hopping leaves one Fermi sea (c = 1).
- **C29.** At its source fields it is a nearly compensated semimetal 2π × 0.03 from the ν = 1 umklapp, and at g ≥ 0.5 it flows to a hybridization insulator (c = 0, both estimators). At ν = 1.2 it keeps c = 2 at g = 0.5 and 1.
- **routeA_spatial.** Four Fermi points, couplings ≤ 0.06, and a closest resonance in the pair transfer (0.011 at the working point). Its numerics track the HF chain, but FES fails its own calibration, so c = 2 rests on the HF chain plus the margins.
- **routeA_compile.** Exactly at ν = 1 when μ_c = 0. With μ_c ≠ 0 it is a detuned c = 2 chain at the HF level, with couplings ≤ 10⁻³. Neither estimator calibrates, because the velocities are 10⁻³ of the bandwidth.

No cell shows an interaction-driven flow away from c = 2 at reachable scales. The flows we see are band-level effects of the generated quadratic terms, cured by retuning one uniform field, except C28 at g = 1.

## Reproducing

numpy and scipy only. `resonances.py`, `exact_open_table.py`, `ring_casimir_free.py` and `table.py` take minutes locally. `run_chain.py SPEC fes|fes2|obc` runs one chain (task lists in `slurm/`, fits in `analyze.py`). Typhon jobs: 17517487 and 17517544 (two-cell), 17517501 and 17517662 (one-cell, interacting and HF, plus two-cell extensions), 17517663 (n(k)) and 17517502 (open chains). Outputs are in `/scratch/lustre/jkudlerflam/routeA_n2_chains/data/`, with the copies used here in `data/`. Of the open-chain runs, only C29 at g = 0.5 is usable. C29 at g = 0.2 ended in a superposition of particle numbers, ⟨N⟩ = 65.87, and the others had not finished χ = 48 when this was written.

## Relation to the rest of the repository

This addresses STATUS help item 2 for the n = 2 cells of REPORT §9.19 and §9.20. Suggested edits there: move C29's and compile's fields off ν = 1, retune C28's with g, and note that at n = 2 the threat to c = 2 is band-level, not an infrared flow.
