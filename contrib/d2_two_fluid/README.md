# Two fluids at d = 2: unpaired particles plus adjacent pairs

## Claim

For H = Σ_j h_{j,j+1} on (ℂ²)^⊗L with one fixed Hermitian h and one-site translation invariance, the two-fluid mechanism (unpaired up spins plus a fluid of adjacent up-spin pairs, c = 2) gives no c = 2 phase anywhere we looked. With U(1) the pair minimum is pinned to twice the single-particle minimum and binding implies clustering [exact]. Charge-one "conditional creation" terms unpin it and give two Fermi seas in the two-body problem on a small open set [numerical], but in every case tested a dense canted state, built by the same terms, undercuts the vacuum before the second sea fills, and the surviving dilute liquid has c = 1 [numerical]. Pairing terms shift the pair minimum only at second order [heuristic]. No counterexample.

## Status

Numerical, with exact two-body statements (section 2). A proof that the canted state always preempts the second sea would upgrade it.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), subagent of the coordinating session, 2026-10-27/28. No Codex/GPT was used.
- Checked by nobody else. Ring ED (`ring_ed.py`) matches dense ED to 2·10⁻¹³ (L ≤ 8, d = 2 and 3), the two-body solver matches cutoff diagonalization to 10⁻⁷ and the U(1) closed form to 2·10⁻¹⁵. `vumps.py`, `fdmrg.py`, `idmrg.py` are copied unchanged from `../d2_followup_220_1282/` and `../counterexample_search/`.

## Method and results

### 1. Setup

With |1⟩ = up = particle (n = (1 − σ^z)/2 in the convention of `specs/_d2_background.md`),

h = −τ|10⟩⟨01| + V|11⟩⟨11| − μ(n⊗1 + 1⊗n)/2 + g_L|11⟩⟨01| + g_R|11⟩⟨10| + Δ|11⟩⟨00| + f_L|10⟩⟨00| + f_R|01⟩⟨00| + h.c.,

τ = e^{iφ}. Every bond has this form up to constants and telescoping terms (`twobody.to_pauli` gives (J, B)). The g terms create a particle next to an existing one (g_L = −g_R real is a Dzyaloshinskii–Moriya vector along y). With Δ = f = 0 the all-down vacuum is exact, a six-parameter slice. There the two-up-spin problem at total momentum K lives on {one particle} ∪ {two particles at separation r ≥ 1}, a tridiagonal chain with single energy −2cos(K + φ) − μ, coupling A(K) = g_L e^{iK/2} + g_R e^{−iK/2} to r = 1, V at r = 1 and relative hopping −2cos(K/2 + φ). `twobody.py` solves it on the half-infinite chain. With U(1) this is exact, otherwise it is a truncation, checked below on rings.

Seas centred at k_a (unpaired) and k_b (pairs) stay decoupled only if m k_a + n k_b ≠ 0 mod 2π for every vertex e^{i(mθ_a + nθ_b)} of dimension below two near K = 1, (m, n) = (1,0), (0,1), (1,±1), (2,0), (0,2), (2,±1), (1,±2). With U(1) only the conversion (2,−1) is allowed.

### 2. U(1) bonds [exact]

These are twisted XXZ. The bound pair E(K) = V + 4cos²(K/2 + φ)/V − 2μ has its minimum at K = −2φ = 2k_a, so the conversion momentum vanishes identically. Binding below the global threshold needs V < −2, and then at the pair onset the up-polarized state is already below the vacuum by (V² − 4)/(2|V|) per site. `check_u1.py` confirms this on 200 random (φ, V) (levels to 1.8·10⁻¹⁵, Q = −2φ, all 65 points with V < −2 preempted by at least 0.18), as in C15 §3.2 and C17 §4.2.

### 3. Charge-one terms: the two-body screen

`screen_tda.py` finds the threshold μ* where the lowest excitation reaches zero, then fills the isolated bands rigidly (μ = μ* + δμ ≤ μ* + 2) until a second sea opens.

| | random | slice g_R = −g_L e^{iδ}, abs(δ) < 0.3 |
|---|---|---|
| bonds | 10,100 | 1,382 |
| two seas | 2,523 | 401 |
| second sea bound by > 0.05 | 429 | 85 |
| and abs(k_b − 2k_a) > 0.2 | 324 | 71 |
| and up-polarized state above vacuum at μ* | 145 | 54 |
| and first-sea density < 0.25 at onset | 99 | 19 |
| and all listed vertex momenta > 0.2 | 57 | 8 |

We sampled φ ∈ [0, π], V ∈ [−4, 4] ([−3, 3] on the slice), abs(g_L), abs(g_R) ∈ [0, 3] and uniform phases. Conversion does at d = 2 what C12's construction needed at d = 3, moving the pair minimum by O(1) and opening a second sea at low density. The random survivors mostly have φ between 1.4 and 2.8 and V < 0 (median −1.6, 11 of 57 repulsive). On the slice they also occur at V ≈ 1.7, g ≈ 2, bound by conversion alone.

On rings of L = 12 (`screen_ed.py`, μ scanned through μ*), 23 of 53 random survivors already have a non-vacuum ground state at μ* − 0.6, 16 jump into a dense state (⟨N⟩/L > 0.3) and 14 have a dilute onset. Their low states carry ⟨N⟩ ≈ 1.3 to 2 and 2.5 to 4, so the "pair" is a short string of fluctuating length.

### 4. Why it fails: a canted dense state

A uniform product state of density n has energy per site, relative to the vacuum,

e(n) = ε(0) n + (V + 2cos φ) n² − 2 abs(g_L + g_R) n^{3/2}(1 − n)^{1/2},  ε(0) = −2cos φ − μ.

The last term, conditional creation evaluated classically, is non-analytic and negative, so a canted state with n ≈ 0.6 to 0.8 drops below the vacuum unless ε(0) is large compared with the conversion. When g_L + g_R = 0, period-three states with density gradients do the same through g_L n_{j+1} + g_R n_{j−1}. `classical.py` minimizes over product states of period up to six and returns the μ_d at which the first one undercuts the vacuum. A dilute two-fluid regime needs μ_d − μ* to exceed the two-sea onset δμ.

| set | bonds checked | μ_d > μ* (largest window, its onset) | window > onset |
|---|---|---|---|
| random survivors | 57 | 6 (0.40, 0.68) | 0 |
| random, bound, abs(k_b − 2k_a) > 0.2 | 168 of 324 | 10 (0.79, 1.24) | 0 |
| slice survivors | 8 | 3 (0.17, 1.0) | 0 |
| slice, bound, abs(k_b − 2k_a) > 0.2 | 37 of 71 | 18 (0.31, 1.68) | 0 |

The product bound is variational, so the true windows are smaller. The terms that unpin the pair are the terms that build the canted state.

### 5. Finite density, calibrated

The main estimator is finite-entanglement scaling with damped VUMPS (`run_vumps.py`, n-site cell, c(χ) = 6ΔS/Δln ξ for consecutive χ, fit over the last three), with ring ED to L = 20 (`ring_fd.py`) second. Two-site iDMRG (`run_idmrg.py`) does not converge its entropy for incommensurate liquids, even on the U(1) control, so it gives only densities, and open-chain scaling (`run_fdmrg.py`) was too slow.

| model (exact c) | VUMPS cell, χ | c, VUMPS | c_S1, rings |
|---|---|---|---|
| XX (1) | | | 1.02 to 1.03 (L = 14 to 20), Casimir 1.01 |
| XXZ, Δ = 0.5 (1) | | | 1.02 to 1.03, Casimir 1.02 |
| XX in field 0.8 (1) | 4, 32 to 64 | 0.99 | |
| XXZ in a field, V = −1.03, ρ = 0.18 (1) | 4, 24 to 48 | 1.00 | 1.12 to 1.15 |
| two XX chains as d = 4, unequal velocities and fields (2) | 2, 32 to 64 | 2.02 | 2.18, 2.12 (L = 8, 10) |
| two identical XX chains as d = 4 (2) | 2, 32 to 64 | 2.14 | 2.11, 2.08 |

Candidate C is φ = 1.40, V = −1.03, g_L = 1.43, g_R = 0.38 e^{1.40i} (at μ = −2.30, J = [[−0.085, 0.493, −0.715], [−0.493, −0.085, 0], [−0.032, −0.187, −0.258]], B = (0.747, 0.187, −0.635)). Its exact L = 12 spectrum has two valleys, ⟨N⟩ ≈ 1.5 at K = −π/2 and ⟨N⟩ ≈ 2.2 to 2.8 near K = −2π/3 to π, 0.11 to 0.15 higher. An iDMRG scan (χ ≤ 48) finds the vacuum at μ = −2.50, a liquid with density rising from 0.13 (μ = −2.44) to about 0.46 (μ = −2.14), and a canted, nearly product state (ρ = 0.66, S = 0.012) at μ = −2.08, against μ_d = −2.26 from product states. Candidate A (φ = 2.98, V = −1.50, g_L = 1.40, g_R = 1.35 e^{−0.32i}) is canted and gapped (ρ ≈ 0.8, ξ ≈ 1) for all μ ∈ [−2.5, −1.3]. S1 to S3 are the random survivors with the widest product windows and AS1, AS2 the slice survivors with positive windows, each run inside its window.

| bond (φ, V, g_L, g_R) | μ | ρ | cell, χ | c(χ), consecutive | c, fit |
|---|---|---|---|---|---|
| C | −2.40 | 0.18 | 9, 8 to 24 | 1.10, 1.01, 1.04 | 1.03 |
| C | −2.30 | 0.28 | 7, 8 to 48 | 1.06, 1.03, 1.00, 1.02, 1.00 | 1.01 |
| C | −2.22 | 0.34 | 7, 8 to 24 | 1.08, 1.03, 1.03 | 1.03 |
| S3 (2.615, −1.835, 0.857, 0.575 e^{1.682i}) | −2.163 | 0.51 | 7, 8 to 48 | 1.16, 1.12, 1.12, 1.04, 1.02 | 1.03 |
| S2 (2.823, −2.318, 0.923, 1.068 e^{−1.758i}) | −2.742 | 0.44 | 7, 8 to 16 | 1.42, 1.30 | unconverged |
| AS1 (0.551, 1.733, 1.814, 1.814 e^{2.938i}) | −2.146 | 0.19 | 7, 8 to 16 | 1.00, 1.13 | 1.04 |
| AS2 (0.819, 1.693, 1.874, 1.874 e^{2.871i}) | −2.324 | 0.20 | 7, 8 to 16 | 0.91, 1.05 | 0.96 |

The growth ξ ∝ χ^κ also fits one conformal mode, κ = 1.33 (C, μ = −2.30), 1.35 (C, −2.40) and 1.27 (S3), against 6/(c(√(12/c) + 1)) = 1.34 for c = 1 and 0.87 for c = 2 (control 1.34). S1 (2.844, −1.768, 0.105, 1.567 e^{2.179i}) gave no usable estimate (unconverged in a four-site cell, non-injective in a seven-site cell). A wrong cell gives early slopes that fall with χ (C in a four-site cell gives 1.25, 1.23, 1.11, 1.03), as S2 does so far. Ring ED for C has ground-state momenta that wander with L, so its entropy fits (1.2 to 5) and gaps (L·Δ/2π between 0.29 and 0.68 at L = 14 to 20, without trend) fix neither c nor z, while the U(1) control gives L·Δ/2π = 0.73 to 0.78, z = 1.

### 6. Pairing terms

With Δ ≠ 0 (Z₂ = (−1)^N) pair creation from the vacuum is first order but carries momentum zero, so at order Δ it only mixes the K = 0 pair with the vacuum. The pair minimum moves off 2k_a only through virtual four-particle states, k_b − 2k_a = O(Δ²/E), so the conversion stays nearly resonant, and binding at small Δ still needs V < −2 with the clustering of section 2. Strong pairing is C15's problem, whose singular corner gives ferromagnetic degeneracies rather than two fluids.

### Verdict

No c = 2 phase was found, so this mechanism gives no counterexample at d = 2. With U(1) it is excluded by an exact two-body argument. Without U(1) the two-body ingredients exist on a small open set, but unpinning the pair needs conversion terms that also pull a canted dense state below the vacuum, which preempts the second sea, and the leftover dilute liquid has c = 1. Open: the generic bond with a dressed vacuum (f, Δ ≠ 0) at finite density, and converged runs for S1 and S2.

## Reproducing

numpy and scipy only. Typhon jobs 17513249 and 17513262 (screen), 17513988 and 17514104 (slice), 17513491 to 17513493 (iDMRG scans), 17513495 (ring calibrations), 17513538 (ring check), 17513857 (product states). With the queue saturated, VUMPS, ring ED for C and the control, and the bound-set product states ran niced on the login node.

- `check_u1.py`. `screen_tda.py SEED 100 OUT.json [as]` (4 min). Raw screens are in `data/tda_*.json.gz`, and `analyze_tda.py FILE` reproduces section 3.
- `screen_ed.py data/tda_random_survivors.json OUT.json 12`, then `analyze_ed.py`. `classical.py IN.json OUT.json` (a minute per bond).
- `run_idmrg.py SPEC --chis 16 32 48 --out OUT.jsonl`. `run_vumps.py SPEC --n 7 --out data/vu_NAME.json`, then `analyze_vu.py NAME`. `ring_fd.py SPEC L`, then `analyze_ring.py NAME`.

Specs are in `models.py` ('p:t,φ,V,μ,g_L,abs(g_R),arg(g_R)' is the particle bond) and batch scripts in `slurm/`.

## Relation to the rest of the repository

This addresses the two-fluid loophole at d = 2 (STATUS, "Small d", item 5 of "Where help is most useful"). It reproduces the U(1) statements of C15 and C17, transplants C12's d = 3 idea to d = 2 with U(1)-breaking conversion, and locates where it fails. "Where things stand" is unchanged.
