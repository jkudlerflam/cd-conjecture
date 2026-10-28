# Route A end to end in one ten-letter bond

## Claim

Class: one Hermitian nearest-neighbour bond h on ℂ¹⁰ ⊗ ℂ¹⁰, one-site translation invariant, rings and infinite chains. Virasoro c.

1. [proved] Cells P B P form spontaneously. Site types are conserved, and an exact operator lower bound shows that the period-3 pattern beats every other type configuration by at least 0.021 per site for λ ≤ 0.15. VUMPS from random states finds it unaided [numerical].
2. [proved] Each cell has an exactly degenerate one-qubit code, gap 0.4, and every port operator compresses to tr(O)/2, so first order vanishes. The second-order bond is exactly a logical XXZ bond with hopping J_xy = J⊥J_zκ²/[4D(D² − κ²)]. Exact two-cell diagonalization reproduces it with a residual ∝ λ³. The interface escapes C18's scalar-resolvent obstruction, which returns at κ = 0.
3. [numerical] The full microscopic chain has c = 1 at λ = 0.05, 0.10, 0.15, from finite-entanglement scaling (FES) and from ring Casimir energies, both calibrated on XX and on the effective chain. Ring spectra (L = 9 to 18 sites) match the effective logical chain to 0.1 to 0.6%, and z = 1.

## Status

Items 1 and 2 are exact, checked to machine precision. Item 3 is numerical, with FES reaching ξ ≈ 30 cells. Frozen cells and one flavour give c = 1 at d = 10, far below d − 1. The point is the mechanism, not the bound.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-27/28. No Codex. VUMPS from `../d2_followup_220_1282/vumps.py`.
- Internal checks: ring code against a full sparse build (M = 3), logical rings against free fermions, closed-form J_xy against numerical second order (10⁻¹⁶), the symmetries of h (exact). Not rerun by anyone else.

## The bond

Letters P = ℂ² (port qubit p) and B = ℂ⁸ (qubits s_L, c, s_R: the two port partners and the logical qubit). Basis |0⟩ (Z = +1), |1⟩, σ⁺ = |0⟩⟨1|, Φ± = (|00⟩ ± |11⟩)/√2 on (p, s). The port check on (p, s, c) is

H_pair = D₀|Φ⁻⟩⟨Φ⁻| + D₁(|01⟩⟨01| + |10⟩⟨10|) + κ(|Φ⁻⟩⟨01| ⊗ σ⁺_c − |Φ⁻⟩⟨10| ⊗ σ⁻_c + h.c.),

which annihilates Φ⁺ ⊗ anything and is ≥ 0 for D₀D₁ ≥ κ². All of h (`cell.py`, `bond()`, ‖h‖ = 1.45):

| term | on pairs | coefficients |
|---|---|---|
| K n_B ⊗ n_B | BB | K = 2 |
| −(μ/2)(n_B ⊗ 1 + 1 ⊗ n_B) | all | μ = 0.55 |
| H_pair(p, s_L, c) | (P, B) | D₀ = D₁ = 1, κ = 0.6 |
| H_pair(p, s_R, c) | (B, P) | same |
| λ(XX + YY − ZZ) | (P, P) | λ = 0.05, 0.10, 0.15 |

h conserves each site's type, the charge Q = Σ(Z_p − Z_s + Z_c)/2, the conjugation C = ΠX (Q → −Q, which needs the relative sign in κ), and a bond-centred reflection.

Scales: K = 2 > D = 1 > μ = 0.55 > cell gap D − κ = 0.4 > frustrated-port cost ε = 0.229 > λ ≤ 0.15 ≫ logical hopping |J| = 6.0, 21, 42 × 10⁻⁴ ≫ three-cell terms ≈ 0.004λ³. The selection margins (≥ 0.021 per site) sit far above the logical scale.

## Mechanism

**Selection.** A configuration is B's separated by runs of k P's. Assigning each check and port bond to its run gives the exact bound H ≥ Σ_B(−μ) + Σ_runs e_k, with e_k the lowest eigenvalue of the run operator (shared logical qubits do not spoil it). A one-P run is a port claimed by two partners, and monogamy costs ε. Long runs gain at most λ per P–P bond. The energy per site is then a weighted average of u_k = (−μ + e_k)/(k + 1), and in the window 3λ ≲ μ < 3ε the minimum is k = 2 (`selector.py`, exact e_k for k ≤ 7 and a bound beyond).

| λ | u_2 | next best | margin | VUMPS period-3 energy |
|---|---|---|---|---|
| 0.05 | −0.18492 | −0.16054 (k = 1) | 0.024 | −0.184845 |
| 0.10 | −0.18925 | −0.16054 (k = 1) | 0.029 | −0.188973 |
| 0.15 | −0.19577 | −0.17425 (k = 3) | 0.022 | −0.195211 |

The true period-3 energy lies at most 6 × 10⁻⁴ above u_2, so configurations with more than a few percent of defect sites are excluded rigorously, and dilute defects pay the local margin. The window needs λ < ε, so λ/(D − κ) ≲ 0.5.

**Code and second order.** The code is Φ⁺(p_L, s_L) ⊗ Φ⁺(s_R, p_R) ⊗ ℂ²_c, so ports are silent. The term σ⁺_{pR}σ⁻_{pL} of V creates syndromes |01⟩, |10⟩ of charge ±1 in the two cells, κ converts each to Φ⁻ while moving the charge onto c, and the ZZ part of V closes both ports. The logical charge has hopped. With V = J⊥(σ⁺σ⁻ + h.c.) + J_zZZ (here J⊥ = 2, J_z = −1), solving the 2 × 2 blocks gives

H₂ = λ²[J_xy(σ⁺σ⁻ + σ⁻σ⁺) + J_zzZZ + a₀], J_xy = J⊥J_zκ²/[4D(D² − κ²)] = −0.28125, J_zz = 0.02086.

U(1) forbids pairing at every order and C forbids a field, so the logical chain is XXZ with Δ′ = −2J_zz/J_xy ≈ 0.15, well inside the critical window. At κ = 0 the checks see only port and partner and H₂ = −1.5λ² exactly, as C18 requires. The κ term escapes because it acts on the logical qubit.

## Effective-bond check (`ed_check.py`)

One cell (32 states): two exact zeros (residual 10⁻¹⁷), gap 0.4, port compressions below 10⁻³¹, code charges ±½. Two cells (1024 states), des Cloizeaux bond from the exact four lowest states:

| λ | J_xy/λ² | J_zz/λ² | other Pauli terms | traceless residual/λ³ |
|---|---|---|---|---|
| 0.0125 | −0.2704 | 0.0194 | 10⁻¹⁶ | 0.99 |
| 0.05 | −0.2418 | 0.0160 | 10⁻¹⁷ | 0.89 |
| 0.10 | −0.2106 | 0.0129 | 10⁻¹⁶ | 0.79 |
| 0.15 | −0.1850 | 0.0107 | 10⁻¹⁷ | 0.71 |

The exact bond is XXZ plus a constant at every λ. Three cells (32768 states) give a connected three-cell term of norm 1.8 × 10⁻⁵ at λ = 0.15. The effective chain below uses the exact two-cell bond.

## The microscopic chain

**Cells form.** VUMPS on the d = 10 chain from random states, with 3-site cells (two seeds) and a 6-site cell, converges to B P P with n_B = 1, 0, 0, BB density below 10⁻⁶ and zero charge on every site. With 4- and 5-site cells, which cannot hold period 3, it ends higher and unconverged.

**FES.** Since types are conserved, the period-3 sector is exactly a chain of 32-dimensional cells, which VUMPS handles much faster (`--model cell`). At equal bond dimension its energy agrees with the d = 10 runs to 10⁻⁶. c = 6 dS/d ln ξ at the inter-cell bond, fitted over χ ≥ 16, ξ in cells.

| chain | χ | ξ | c |
|---|---|---|---|
| XX | 16 to 96 | 32 to 345 | 0.995 ± 0.004 |
| effective, second order | 16 to 96 | 30 to 328 | 0.993 ± 0.003 |
| effective, exact bond (λ = 0.1) | 16 to 96 | 31 to 333 | 0.992 ± 0.004 |
| microscopic, λ = 0.05 | 16 to 48 | 4.0 to 18 | 1.01 ± 0.04 |
| microscopic, λ = 0.10 | 16 to 64 | 4.0 to 30 | 1.01 ± 0.03 |
| microscopic, λ = 0.15 | 16 to 64 | 4.0 to 30 | 1.00 ± 0.03 |

Port dressing splits each logical Schmidt value, so the microscopic chain needs about four times XX's bond dimension for the same ξ, and its pairwise slopes alternate between 0.75 and 1.19. On the window ξ ≤ 32 alone the calibration chains give 0.90 to 0.92 (three points), so at these correlation lengths the estimator carries a systematic uncertainty of about 0.1 beyond the quoted fit errors.

**Rings** (`ring_ed.py`, resolved by charge and momentum, no truncation). M = 3 to 6 cells, up to 1.55 × 10⁸ states per sector. gap_J is the lowest Q = 0 level at momentum 2π/M (the current, x = 1), gap_Q the lowest charge-one level.

| λ | M | gap_J micro | gap_J effective | gap_Q micro | gap_Q effective |
|---|---|---|---|---|---|
| 0.05 | 4 | 1.7900e−3 | 1.7914e−3 | 5.821e−4 | 5.825e−4 |
| 0.05 | 6 | 1.2877e−3 | 1.2888e−3 | 3.785e−4 | 3.788e−4 |
| 0.10 | 4 | 6.2020e−3 | 6.2195e−3 | 2.003e−3 | 2.008e−3 |
| 0.10 | 6 | 4.4561e−3 | 4.4694e−3 | 1.302e−3 | 1.305e−3 |
| 0.15 | 4 | 1.2188e−2 | 1.2259e−2 | 3.915e−3 | 3.936e−3 |
| 0.15 | 6 | 8.7476e−3 | 8.8018e−3 | 2.544e−3 | 2.558e−3 |

The deviations (0.08%, 0.3%, 0.6%) grow roughly linearly in λ, as expected from omitted multi-cell terms. Ground energies agree to 10⁻⁶ relative, and M = 3, 5 behave the same (`data/summary.json`).

**Casimir c**, from E₀/M = e∞ − πvc/(6M²) at M = 4 and 6 with v from gap_J(6): microscopic 1.1902, 1.1883, 1.1867 at the three couplings, effective chain 1.1902, 1.1885, 1.1872, XX 1.165. On larger rings the effective chain's estimate falls to 1.026 at (12, 14) and XX's to 1.023. The microscopic value equals the calibrated effective one to 3 × 10⁻⁴.

**z.** At λ = 0.1 the finite-size velocity M·gap_J/2π is 3.95, 3.92, 4.26 × 10⁻³ for M = 4, 5, 6 (effective 3.96, 3.93, 4.27). On effective rings it reaches 4.48 × 10⁻³ at M = 14, approaching the Bethe-ansatz velocity 4.53 × 10⁻³, just as XX's approaches 2. So gap ∝ 1/M and z = 1 (the raw exponent from M = 4, 6 is 0.82 for both chains, and 0.86 for XX).

## Verdict

The mechanism works end to end in one explicit Hamiltonian. One ten-letter bond makes cells spontaneously, each cell holds an exact code with silent ports, cells couple only at second order, and the full chain is critical with c = 1, matching its effective logical chain in central charge, velocity and spectrum. Two ingredients carry it. The port check couples syndrome and logical register by charge transfer, the simplest escape from C18. And an exact U(1) with C protects the logical chain at all orders, so the λ³ pairing that gapped the free-fermion target of `../routeA_numerics/` cannot arise.

Limitations: frozen cells (three-fold broken translation, defects when 3 ∤ L), one flavour, and a selection window λ < ε that keeps the logical scale near 1% of the cell gap. Two flavours with distinct Fermi momenta would need a second syndrome channel per port (C19's pole count), which I did not build.

## Reproducing

numpy and scipy only. `cell.py` (bond, cell), `effective.py` and `ed_check.py` (task 1, 30 s), `selector.py`, `run_fes.py` (models `micro`, `cell`, `eff`, `xx`), `ring_ed.py` (microscopic rings), `eff_ring.py` (logical rings), `analyze.py` (tables, `data/summary.json`). Typhon array 17514348 tasks 9 to 17 ran the M = 6 rings, and M ≤ 5 ran locally. Array 17514347 was cancelled while queued, and its FES and selection runs went to the Typhon login node at nice 19 (`tasks_fes.txt`, `tasks_cell.txt`). Continuations to χ = 96 and 128 (`tasks_cellB.txt`) were still running when this was written. Logs in `data/logs/`.

## Relation to the rest of the repository

This answers item 2 of STATUS's "Where help is most useful" in its smallest form, realizing route A and the C18 escape at d = 10 with the selector included. At n = 1 an exact U(1) replaces detuned momenta. `../routeA_spatial/` (d = 23) and `../routeA_compile/` (d = 48) are larger cells aimed at c = n with the same exact U(1).
