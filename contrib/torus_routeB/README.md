# Species torus against the route-B instability

## Claim

The proposal under test is the species-torus selection rule of [`../chat_followup/torus_encoding.md`](../chat_followup/torus_encoding.md), imported from the collaborator's chat "Prove c d conjecture" (the idea is theirs). A register with states |0,β⟩, |1,α⟩ (b = 2r) carries N = r² channels F_{αβ}. The 2r − 1 independent species phases are an exact symmetry, so inter-flavour pair tunnelling is forbidden while flavour exchange must be detuned by distinct Fermi momenta k_a. Here r = 2, N = 4.

1. Effective theory, Fermi points ±k_a. With only total charge and parity, pair tunnelling stays resonant and the one-loop flow runs away for both signs of the pair-hopping amplitude. With the torus and distinct k_a only forward couplings stay resonant and the flow is bounded, umklapp included. With the torus and equal k_a, exchange and rectangle couplings run away for both signs. VUMPS and open-chain entropies at g = 0.3 agree in the cases where the densities stay near their design values (§1).
2. Route B. The graded port C^{2|2} of a q = 4 route-B cell is the torus register. A cell with flavour-dependent dressing keeps an exactly 16-fold degenerate code and the exact torus, and every second-order term (two- and three-cell) commutes with the torus to 10⁻¹⁴. The allowed marginal couplings are umklapp, exchange, rectangle and forward ones. In the untuned cell all are resonant and the flow runs away at RG time l* ≈ 8.3/λ (Δ = 1). In the tuned cell with a common chemical potential only forward couplings remain and the flow is bounded. At physical λ ≤ 0.4Δ the effective chains track the free chain at accessible scales, as the long l* predicts (§2).
3. Verdict. Torus plus distinct k_a plus doping off half filling gives a stable c = N at r = 2 at the physics standard. The torus alone does not, since the untuned (4,2) cell already has it (§3).

## Status

Numerical and one-loop. Upgrades would be open-chain entropies at χ ≥ 256 (or a symmetry-resolved DMRG), and a check at r = 3, where rectangles first have more than one shape.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27/28. No Codex run, since every piece is a short script whose specification would be as long as the code.
- Checks. The one-loop sign is calibrated on the Hubbard chain (U > 0 bounded and U < 0 runaway at k = 0.3π, U > 0 runaway at half filling), in both the u(4) and the so(8) codes. The normalization dG/dl = (1/4πv) f f G G follows from the OPE and reproduces dg/dl = g²/2πv for g J_R·J_L with J = ψ†(σ/2)ψ. The second-order code reproduces `../encoded_chain_numerics/effective.py` at q = 3 to 10⁻¹⁴. Every estimator is run on decoupled free chains. Not rerun by anyone else.

## Method and results

### 1. Effective theory (`fermions.py`, `tcases.py`, `oneloop_cases.py`)

Four flavours a = 2α + β on a chain, local dimension 16, hopping t_a = 1/(2 sin k_a) and μ_a = −2t_a cos k_a, so every flavour has Fermi points ±k_a and velocity 1. Distinct momenta k/π = (0.30, 0.375, 0.41, 0.555) maximize the margins (0.035π for every two- and four-field momentum sum, 0.07π for torus-allowed ones). Equal momenta are k = 0.4π. Interactions are random nearest-neighbour couplings of strength g in classes U, V (densities), J (exchange), R (rectangles f†₀₀f†₁₁f₁₀f₀₁), P (pair hopping) and X (other number-conserving flavour-changing terms). Case (i) has all classes and only total charge, with s = ±1 on P. Cases (ii) (distinct k) and (iii) (equal k) keep only U, V, J, R, which are torus-symmetric (checked on the lattice operator, |[H, T_s]| = 0), with s = ±1 on J and R. Case (ip) adds pairing, which leaves only parity.

One loop. Lattice terms are projected on currents with f_{x,a} → e^{ik_a x}ψ_{aR} + e^{−ik_a x}ψ_{aL}, both on u(4) (resonant terms only) and on so(8), which includes umklapp ψ†_Rψ†_Rψ_Lψ_L. In the so(8) flow each coupling is dropped once its momentum exceeds the running cutoff e^{−l}.

| case | resonant non-forward couplings | so(8) flow, seed 0, g = 0.1 / 0.2 / 0.3 | runaway, 100 seeds, g = 0.1 and 0.2 |
|---|---|---|---|
| (i) s = ±1 | pair tunnelling J^{ab}_R J^{ab}_L | l* = 7.0, 3.7, 2.5 (+); 5.4, 2.8, 1.9 (−) | 100 of 100, both signs |
| (ii) s = ±1 | none | bounded to l = 60 | 0 of 100, both signs |
| (iii) s = ±1 | exchange, rectangles | l* = 3.8, 1.9, 1.3 (+); 6.8, 3.4, 2.3 (−) | 100 of 100, both signs |
| (ip) | pairing ψ_{aR}ψ_{aL} at zero momentum | relevant at tree level | |

The window of case (ii) closes at strong coupling: at g = 0.5 with s = +1 the u(4) flow runs away at l = 1.1, before the exchange couplings are detuned.

Numerics (Typhon). VUMPS on a two-site unit cell (`run_vumps.py`, damped mode, which converges at incommensurate filling where plain VUMPS oscillates), χ = 48 to 256, c_FES = 6 dS/d ln ξ over the last three χ. Open chains of L = 24 by two-site DMRG (`run_obc.py`), fitted to (c/6) ln chord plus Friedel terms at 2k_a. Calibration: free chains give c_FES = 4.56 (distinct k, χ ≤ 256, ξ ≤ 8.9) and 4.21 (equal k, χ ≤ 192). The open-chain fit gives 3.67 (distinct) and 4.02 (equal) on the exact free chain at L = 24, but DMRG at χ = 64 only reaches S(L/2) = 2.32 of the exact 3.31, so at that χ the open chains are compared with the free chain at equal χ.

RESULTS-TASK1

### 2. Microscopic cell (`torus_cell.py`, `routeB_torus.py`)

The cell is that of `../encoded_chain_numerics/cells.py` with q = 4, q₊ = 2, n = 4 (dimension 4096). Its eight odd channels are the real and imaginary parts of F_{αβ} = |1,α⟩⟨0,β|, flavour a = 2α + β is qubit a of the register, and the dressing is g_a(|Ω_A⟩⟨Ω| + h.c.) ⊗ γ_A with A ∈ {2a, 2a + 1}. The torus acts as T_s = Π_s(p) − Π_s(s) on both port pairs and −Σ_a W_{as} n_a on the register.

| check | untuned (g_a = 0.5) | tuned (g = 0.35, 0.5, 0.65, 0.8) |
|---|---|---|
| code splitting, gap | 8 × 10⁻¹⁵, 1.786 | 8 × 10⁻¹⁵, 1.842 |
| [T_s, H_cell], [T_s, bond] | 0, 0 | 0, 0 |
| first-order hopping / λ per flavour | 0.104 each | 0.043, 0.082, 0.126, 0.173 |
| relative [T, T1], [T, T2], [T, T3] | 3 × 10⁻¹⁵, 3 × 10⁻¹⁵, 7 × 10⁻¹⁴ | 4 × 10⁻¹⁵, 4 × 10⁻¹⁵, 3 × 10⁻¹⁴ |
| largest torus-charged current coupling | 1.6 × 10⁻¹⁶ | 1.4 × 10⁻¹⁶ |
| allowed couplings per λ²: umklapp, exchange, rectangle, forward | 0.38, 0.15, 0.09, 0.17 | 0.40, 0.14, 0.07, 0.19 |

So the torus is exact at second order. Since it is a symmetry of the microscopic chain, it holds at every order. The physical upper-species charge restricts to −Q on the code (residual 2 × 10⁻¹⁴), so μ₀Q is a torus-symmetric physical term and is added exactly after Schrieffer–Wolff. We take μ₀ = 0.3 × 4λ⟨t⟩. The quadratic part then puts flavour a's Fermi points at δk_a and π − δk_a with δk/π = 0.27, 0.13, 0.08, 0.06 (free values), and the smallest mismatch of an allowed non-forward coupling is 0.14.

| cell | one loop, λ = 0.1 / 0.2 / 0.4 |
|---|---|
| untuned (δk = 0) | everything resonant, runaway at l* = 83 (beyond l = 80), 41, 21 |
| uniform g plus μ₀ (equal δk) | umklapp detuned, exchange and rectangles resonant, runaway at 2l* |
| tuned plus μ₀ | only forward couplings survive, bounded to l = 80 |

The dimensionless coupling is max|G|/2πv ≈ 0.03λ/Δ per flavour hopping, so l* grows like 1/λ and the untuned instability sits at ξ ~ e^{8/λ}, far beyond any simulation. Amplifying the λ² terms tenfold brings l* down to 1.4 to 2, but at that strength the on-site second-order terms exceed the slow flavour's bandwidth and every chain becomes a trivial product state (S = 0.32, ξ < 1), so this stress test discriminates nothing.

RESULTS-TASK2

The full chain of 4096-dimensional cells was not simulated. With c = 4 it would need χ ≳ 100 on a site 12 times larger than the q = 3 cell of `../encoded_chain_numerics/`, which already took an hour per χ step at χ ≈ 32.

### 3. Verdict and what remains

VERDICT

## Reproducing

numpy and scipy only. `python3 oneloop_cases.py` (a few minutes, `oneloop_cases_stdout.txt`). `python3 routeB_torus.py --cache DIR` (about 25 minutes, three 4096 × 4096 diagonalizations, `routeB_torus_stdout.txt`, `data/eff_*.npz`, `data/routeB_torus_summary.json`). `run_vumps.py` and `run_obc.py` with the task files `fes*_tasks.txt`, `obc*_tasks.txt`, `eff_tasks.txt`, `effobc_tasks.txt` and the slurm scripts beside them, then `python3 analyze.py fes|obc data/*.json` and `analyze.py exact`. Typhon jobs: JOBIDS.

## Relation to the rest of the repository

This tests STATUS step 4's repair of route B and the red team's §2 (`../redteam_construction/`). STATUS could say that the species torus of `../chat_followup/` is exact in a route-B cell with a (2|2) port, that with distinct hoppings and a common chemical potential it leaves only forward couplings at second order, and that the untuned (4,2) cell already had the torus, so the runaway the red team found there came from equal Fermi points (umklapp, exchange, rectangles), not from pair tunnelling.
