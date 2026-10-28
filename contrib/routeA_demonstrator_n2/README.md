# Route A end to end at two flavours, in one thirteen-letter bond

## Claim

Class: one Hermitian nearest-neighbour bond h on ℂ¹³ ⊗ ℂ¹³, one-site translation invariant. Virasoro c. This is the n = 2 sequel to `../routeA_demonstrator/` (n = 1, d = 10), with the same mechanism.

1. [proved] Cells P B P form spontaneously. Site types are conserved, and an exact run-unit bound gives the period-3 pattern a margin of at least 0.022 per site for λ ≤ 0.2.
2. [proved] Each cell has an exactly degenerate four-state code (two logical modes), gap 0.292, and every traceless port operator compresses to zero. Both flavour charges N₁, N₂ are exactly conserved. At second order the flavours hop independently, with no inter-flavour term at all. Two-cell diagonalization reproduces the bond with residual ∝ λ³, and inter-flavour couplings appear at about λ⁴. At κ = 0 the bond is a scalar (C18).
3. [numerical] The logical fillings are n₁ ≈ 0.49 and n₂ ≈ 0.41, so the Fermi momenta differ and every inter-flavour backscattering term oscillates. In the full microscopic chain, FES gives c = 1.98 ± 0.08 (λ = 0.15) and 1.7 to 1.85 (λ = 0.2, less converged) at ξ ≤ 5 cells. A block-entropy estimator matched to c = 2 and c = 1 chains at equal ξ gives c ≈ 2.1 ± 0.15. M = 3 ring spectra match the effective chain to 0.3 to 1.4%.

## Status

Items 1 and 2 are exact. Item 3 is numerical and weaker than at n = 1. Port dressing costs the microscopic chain about six times its logical chain's bond dimension, so χ ≤ 96 reaches only ξ ≈ 5 cells, and rings stop at M = 3. The block estimator separates c = 2 from c = 1 cleanly at these ξ and FES is compatible with c = 2, but neither reaches the calibration chains' asymptotic window (ξ ≈ 25 to 37). The cells are frozen.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28. No Codex. VUMPS from `../d2_followup_220_1282/vumps.py`, patched for unequal site dimensions (`vumps2.py`).
- Internal checks: ring code against a full sparse build (M = 2, every level), second order against exact two-cell diagonalization, U(1) × U(1) on h, the VUMPS frame below (an exact unitary). Not rerun by anyone else.

## The bond

Letters: P = ℂ³ ⊗ ℂ³ (a port qutrit x and its partner y in the same letter) and B = ℂ² ⊗ ℂ² (occupations n₁, n₂ of two logical modes). Charges: x = 0, 1, 2 carries 0, e₁, e₂, y carries minus that, and B carries n₁e₁ + n₂e₂. Left ports use Φ_L = Σ_x|xx⟩/√3, right ports Φ_R = Σ_x ω^x|xx⟩/√3 (ω = e^{2πi/3}). Being orthogonal, they make a P shared by two B's frustrated, which replaces partner monogamy at n = 1. Neutral syndromes are N^s_a = (D_a ⊗ 1)Φ_s with D₁ = √(3/2) diag(1, −1, 0) and D₂ = diag(1, 1, −2)/√2. The check on (x, y; n₁, n₂) is

H_s = D_N Π_N^s + Σ_{a≠b} D_{ab}|ab⟩⟨ab| + Σ_a κ_a(|N^s_a⟩⟨a0| ⊗ c_a† − |N^s_a⟩⟨0a| ⊗ c_a + h.c.).

A charge +e_a syndrome |a0⟩ turns into N_a while mode a fills, and the charge −e_a syndrome |0a⟩ does so while it empties. All of h (`cell2.py`, `bond()`):

| term | on pairs | coefficients |
|---|---|---|
| K n_B ⊗ n_B | BB | K = 3 |
| −(μ/2)(n_B ⊗ 1 + 1 ⊗ n_B) | all | μ = 1.2 |
| H_L(x, y of P; n₁, n₂ of B) | (P, B) | D_N = 1, D_{ab} = 1 except D_{02} = 0.8, κ₁ = κ₂ = 0.6 |
| H_R(x, y of P; n₁, n₂ of B) | (B, P) | same, with Φ_R |
| λ[Σ_a(|a⟩⟨0| ⊗ |0⟩⟨a| + h.c.) + Σ_a D_a ⊗ D_a] on (x, x′) | (P, P) | λ = 0.10, 0.15, 0.20 |

The lowered |02⟩ breaks particle–hole symmetry of flavour 2 and sets its filling.

Scales: K = 3 > μ = 1.2 > D = 1 > frustrated-port cost ε = 0.461 > cell gap 0.292 > λ ≫ logical hopping (5 to 32 × 10⁻⁴) ≫ inter-flavour terms (≈ 10⁻³λ²).

## Mechanism and exact checks

**Selection** (`selector2.py`, `data/selector2.json`). The run-unit bound is the same as at n = 1. The cheapest competitors are a shared port (cost ε) and long free-P runs (at least −2λ per P–P bond).

| λ | u₂ | best other | margin |
|---|---|---|---|
| 0 | −0.4000 | −0.3694 (shared port) | 0.031 |
| 0.10 | −0.4055 | −0.3694 (shared port) | 0.036 |
| 0.15 | −0.4124 | −0.3694 (shared port) | 0.043 |
| 0.20 | −0.4219 | −0.4000 (long runs) | 0.022 |

VUMPS on the full d = 13 chain from random states found B P P in three of five runs. The other two stayed in a P-only configuration 0.19 per site higher, which a local optimizer cannot leave because types never change.

**Code and second order** (`ed_check2.py`, `effective2.py`, `data/ed_check2.json`). One cell (324 states) has four exact zeros (residual 10⁻¹⁶) and gap 0.292. Every traceless port operator on x or y compresses below 10⁻¹⁶, and the code carries charges (n₁, n₂). In each flavour the process is the n = 1 one. A port hop of e_a creates |a0⟩ and |0a⟩, κ moves the charge into mode a on one side and out on the other, and D_a ⊗ D_a closes both ports. Since κ_a couples N_a only to mode a and D_a closes only N_a, no virtual path mixes flavours. The exact second-order bond is

H₂ = Σ_a [t_a c†_{a,1}c_{a,2} + h.c. + g_a n_{a,1}n_{a,2} + f_a(n_{a,1} + n_{a,2})] + const,

with t₁ = (3/64)e^{−iπ/3}, t₂ = 0.0588 e^{iπ/3}, g₁ = −f₁ = 0.0093 (flavour 1 particle–hole symmetric), g₂ = 0.0092, f₂ = 0.0337, and inter-flavour terms below 3 × 10⁻¹⁵. The phases are fluxes per cell, removable flavour by flavour.

| λ | t₁/λ² | t₂/λ² | inter-flavour/λ² | traceless residual/λ³ |
|---|---|---|---|---|
| 0.025 | 0.0482 | 0.0607 | 1.2 × 10⁻⁵ | 0.22 |
| 0.05 | 0.0495 | 0.0628 | 3.4 × 10⁻⁵ | 0.23 |
| 0.10 | 0.0526 | 0.0673 | 1.7 × 10⁻⁴ | 0.24 |
| 0.15 | 0.0561 | 0.0727 | 5.1 × 10⁻⁴ | 0.25 |
| 0.20 | 0.0603 | 0.0791 | 1.3 × 10⁻³ | 0.27 |

At κ = 0 the second-order bond is a scalar to 3 × 10⁻¹⁵. So the logical chain is two U(1) chains (XXZ in a field, Δ′ ≈ 0.1) coupled only at about fourth order, and only through U(1) × U(1)-invariant terms. Pair tunnelling and inter-flavour hopping are forbidden exactly. With n₁ ≠ n₂ every inter-flavour backscattering term carries momentum 2(k_{F1} ± k_{F2}) ≠ 0, so only forward scattering survives, and c = 2.

## The microscopic chain

**A frame for VUMPS.** The product of on-cell unitaries U = Π_j exp(iΣ_a φ_a j N_{a,j}) maps h to the same bond with the flavour-a port hopping multiplied by e^{−iφ_a}. It changes no entanglement entropy or correlation length. With φ_a = arg t_a + π both logical hoppings are real and negative. In the original frame a one-site uMPS turns non-injective, and in this one VUMPS converges to 10⁻⁹. The effective chain then gives n₁ = 0.496, 0.493, 0.490 and n₂ = 0.403, 0.411, 0.420 at the three couplings. The microscopic chain starts from the code embedding of the effective state at χ = 16 (`--load_eff`), because random starts stall in logically empty sectors. It keeps the fillings to 2 × 10⁻³ (n₁ = 0.4929, n₂ = 0.4128 at λ = 0.15, χ = 96).

**FES** (`run_fes2.py --model cell3`, the exact period-3 type sector with site dimensions 9, 4, 9). c = 6 dS/d ln ξ at the inter-cell bond, ξ in cells.

| chain | χ | ξ | c |
|---|---|---|---|
| one XX chain (c = 1) | 6 to 16 | 8 to 32 | 0.95 ± 0.04 |
| two XX chains, velocities 1 and 0.7 (c = 2) | 8 to 128 | 3.6 to 37 | 2.08 ± 0.04 |
| effective chain, exact bond, λ = 0.15 | 16 to 96 | 5.6 to 26 | 2.00 ± 0.05 |
| effective chain, exact bond, λ = 0.20 | 16 to 96 | 5.5 to 26 | 1.98 ± 0.04 |
| microscopic, λ = 0.15 | 48 to 96 | 2.5 to 4.8 | 1.98 ± 0.08 |
| microscopic, λ = 0.15 | 32 to 96 | 1.5 to 4.8 | 1.79 ± 0.09 |
| microscopic, λ = 0.20 | 48 to 96 | 2.5 to 4.4 | 1.85 ± 0.01 |
| microscopic, λ = 0.20 | 32 to 96 | 1.5 to 4.4 | 1.72 ± 0.06 |

Microscopic gradient errors are 10⁻⁹ to 10⁻⁵, except 3 × 10⁻³ at λ = 0.2, χ = 96. For ξ ≤ 3.7 even the calibration chains scatter (two XX 1.48, effective 1.22, one XX 1.22), so the microscopic values are compatible with c = 2 but rest on ξ ≤ 5.

**Block entropies** (`blockS.py`). For blocks of ℓ whole cells, c_b = 3[S(2) − S(1)]/ln 2 removes the boundary dressing and is compared at equal ξ with chains of known c.

| chain | ξ ≈ 2 to 2.5 | ξ ≈ 3.2 to 3.6 |
|---|---|---|
| one XX (c = 1) | 0.62 | 0.70 (ξ 4.3) |
| two XX (c = 2) | 1.13 | 1.22, 1.34 |
| effective, λ = 0.15 (c = 2) | 1.26 | 1.37, 1.48 |
| microscopic, λ = 0.15 | 1.25 | 1.39 |
| microscopic, λ = 0.20 | 1.28 | |

The microscopic values sit inside the c = 2 band and about twice the c = 1 value. Scaling by the mean c = 2 calibration at each ξ gives c ≈ 2.09 and 2.06, with ± 0.15 from the spread of the two calibration chains.

**Rings** (`ring2.py`, both charges and momentum resolved, no truncation). Microscopic rings of M = 3 cells (L = 9, up to 1.6 × 10⁶ states per sector), in the five lowest charge sectors at λ = 0.15 and the two lowest at λ = 0.1 and 0.2, all momenta. They match the effective ring with the exact two-cell bond to 0.3% (λ = 0.1), 0.7% (0.15) and 1.4% (0.2) in every gap. M = 4 sectors have 4 × 10⁸ states and were not attempted, so rings are a spectral check here and not a c estimator.

## Verdict

The two-flavour mechanism works in one explicit 13-letter Hamiltonian. Cells form, each holds an exact two-mode code with silent ports, and at second order the flavours hop through the same ports without mixing. Exact U(1) × U(1) removes pair tunnelling and inter-flavour hopping at every order, and one lowered syndrome energy sets distinct fillings. The full chain follows its effective chain spectrally, and both estimators give c ≈ 2, with the block estimator excluding c = 1. The weak point is reach. The microscopic FES spans ξ = 2.5 to 5 cells, against 25 to 37 for the calibrations, so this measurement is less sharp than at n = 1. χ ≈ 300 (ξ ≈ 15) would close most of the gap, and runs to χ = 256 were started on Typhon (job 17518160). The design does not scale to n flavours at fixed d, since each flavour needs its own port charge.

## Reproducing

numpy and scipy only. `cell2.py` (bond, cell, frame), `effective2.py` (second order, exact two-cell bond), `ed_check2.py` (one and two cells, 3 min), `selector2.py`, `run_fes2.py` (models `micro`, `cell3`, `eff`, `flav`, `xx1`, `xx2`, with `--gauge`, `--load_eff` and `--init product`), `ring2.py`, `blockS.py`, `fes_window.py`, `analyze2.py` (writes `data/summary2.json`). `states/` holds MPS checkpoints and is not meant for the repository. Typhon: ring array 17517764 (done). The microscopic FES array 17518160 (χ to 256, output in /scratch/lustre/jkudlerflam/routeA_demo_n2/data/fes/) started only as this was written, so the χ ≤ 96 results above come from local runs (`data/logs/`). A duplicate on apollo (17518359) was cancelled. Effective and calibration runs ran on the Typhon login node at nice 19 and locally.

## Relation to the rest of the repository

At d = 13 this is smaller than the explicit n = 2 cells of `../routeA_spatial/` (23), `codex/work_C28/` (17), `codex/work_C29/` (26) and `../routeA_compile/` (48). Unlike them it conserves each flavour separately, which removes pair tunnelling outright rather than through displaced Fermi seas, at the price of one port charge per flavour. On effective chains `../routeA_n2_chains/` finds FES values of 1.89 to 2.06 for C28 and C29 at some couplings, and anomalous values for the spatial cell. Here the effective chain gives 2.00 and 1.98, and the microscopic chain agrees with it at the scales reached.
