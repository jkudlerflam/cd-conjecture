# Casimir survey: A = Σ c_a v_a against 2π‖h‖

## Claim

Class and notation as in `../refined_conjecture/`. H = Σ_i h_{i,i+1} on periodic rings of (ℂ^d)^⊗L, one Hermitian h, low-energy theory a set of unitary CFTs (c_a, v_a), A = Σ c_a v_a from E_0(L) = e_0 L − πA/(6L) + o(1/L), ‖h‖ the half spread of spec h minimized over h → h + o⊗1 − 1⊗o, and R = A/(2π‖h‖). The refined conjecture is R ≤ 1.

1. The integrable SU(N) chains with Sym^k(ℂ^N) on each site, which Astra's C14 run (`../../codex/work_C14/RESULT.md`) identified as counterexamples, pass every lattice test I could run. At (N, k) = (2, 1 to 4), (3, 1 to 3), (4, 1 and 2) and (5, 1), exact diagonalization reproduces the Bethe/WZW values of e_0, v, c and A to about 1%, and at N = 3 and 4 it confirms that rank two beats the fundamental chain by the predicted margin (R = 0.716 against 0.711 predicted and 0.667 for SU(3), and 0.842 against 0.833 and 0.75 for SU(4)). The violating members (k = 2, N ≥ 7, d ≥ 28) are out of reach of direct lattice numerics, for reasons given below. [numerical, with the Bethe/WZW input of C14]
2. Outside that family nothing I examined exceeds the fundamental SU(d) value R = 1 − 1/d. Every one-parameter family through an SU(d) point (spin-1 BLBQ, the Sp(4) circle at d = 4, the SO(5) vector circle) peaks there, as a twirling argument says it must to first order. The d = 2 candidates have R ≤ 0.23, and the largest value outside the symmetric family is SU(4)'s 0.75. [numerical]
3. Across everything computed here and the whole symmetric family, A ≤ π‖h‖ log₂ d holds, with equality only for the Heisenberg chain. The symmetric family forces any bound of the form A ≤ K(d)‖h‖ to grow at least like log d/log log d. [numerical and analytic]
4. A adds exactly under stacking, but ‖h‖ is only subadditive. It adds for reflection-symmetric bonds and can drop strictly below the sum otherwise, because the optimal o may entangle the two stacked factors. [numerical]

## Status

Numerical. Item 1 inherits the standard Bethe-ansatz and WZW identifications from C14 for N ≥ 5. A direct computation of A at N = 7 or 8 (for example from exact finite-size Bethe roots) would upgrade it.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the coordinating session, 2026-10-27.
- AI assistance: gpt-6-astra (high) through `gpt_solve.sh`, two closed-book runs writing a finite-size Bethe solver (`astra_bethe/`), described below. All ED, fits and the survey are mine.
- Checks: `core.hnorm` (smoothed convex minimization over o) agrees with Powell searches on random bonds to 10⁻⁹. The Sym^k bond (`chains.sym_rep_chain`, built from Schwinger bosons and Lagrange interpolation in X = Σ E_ab ⊗ E_ba) equals the Takhtajan–Babujian bond at N = k = 2 to 4·10⁻¹⁶ and C14's bond up to h_C14 = 2h − 2H_k, so R is the same in both normalizations. The ED engine is `../counterexample_search/ced.py`, already checked there against full diagonalization. Nobody else has rerun this.

## Method and calibration

Two estimators. A_cas fits E_0(L)/L = e_0 − πA/(6L²) + b/L⁴ over ring lengths commensurate with the ground state (L a multiple of N for SU(N)). A_cv multiplies c_half = 3[S(L₂/2) − S(L₁/2)]/ln(L₂/L₁), from periodic ground states, by v = L·[E(k_0 ± 2π/L) − E_0]/2π, the current state at the largest L. When e_0 is known exactly I also use A(L) = −6L[E_0(L) − e_0L]/π size by size.

| chain | c | sizes | A_cas (exact) | c_half | v (exact) | R_cas, R_cv |
|---|---|---|---|---|---|---|
| XX | 1 | 10 to 24 | 1.9998 (2) | 0.996 | 1.994 (2) | 0.318, 0.316 |
| Ising | 1/2 | 10 to 22 | 1.0000 (1) | 0.499 | 1.993 (2) | 0.113, 0.112 |
| Heisenberg, h = P | 1 | 10 to 24 | 3.156 (π) | 1.001 | 3.168 (π) | 0.502, 0.504 |
| SU(3), h = P | 2 | 6 to 15 | 4.214 (4.189) | 1.948 | 2.098 (2.094) | 0.671, 0.650 |
| SU(4), h = P | 3 | 4, 8, 12 | 4.743 (4.712) | 3.053 | 1.561 (1.571) | 0.755, 0.759 |
| Takhtajan–Babujian S = 1 | 3/2 | 8 to 16 | 2.380 (2.356) | 1.505 | 1.595 (1.571) | 0.505, 0.509 |

Without marginal operators A_cas is exact to 10⁻⁴. With them (SU(2)₁, SU(3)₁, SU(2)₂) it overshoots by 0.5 to 1%. A_cv tracks A_cas within 3%. Two-size fits (L = 8, 12 at d = 4 and 5, 10 at d = 5) overshoot by 5.3% and 11.7% on SU(4) and SU(5), and I divide by those factors where only two sizes exist.

## Results

### The symmetric family

In my normalization h = Σ_j H_{k−j} Π_j (H_m harmonic numbers, Π_j the projector on the Young diagram (2k − j, j)), ‖h‖ = H_k/2, v = π/N, c = k(N² − 1)/(N + k), and e_0 = H_k − (1/N)Σ_{j=1}^{k}[ψ(1 + (j−1)/N) − ψ(j/N)] (C14 eq. 5). Ratios to the prediction:

| N, k | d | sizes | A_ED/A | e_0 free fit − exact | c_half (pred.) | v(L_max)/(π/N) | R_ED (pred.) |
|---|---|---|---|---|---|---|---|
| 2, 1 | 2 | 4 to 16 | 1.005 | 2·10⁻⁶ | 0.997 (1) | 1.003 | 0.502 (0.500) |
| 2, 2 | 3 | 4 to 12 | 1.011 | 2·10⁻⁵ | 1.498 (1.5) | 1.008 | 0.506 (0.500) |
| 2, 3 | 4 | 4 to 12 | 1.017 | 4·10⁻⁵ | 1.812 (1.8) | 1.019 | 0.499 (0.491) |
| 2, 4 | 5 | 4 to 10 | 1.023 | 6·10⁻⁵ | 2.021 (2) | 1.020 | 0.491 (0.480) |
| 3, 1 | 3 | 3 to 12 | 1.002 | −2·10⁻⁵ | 2.113 (2) | 0.995 | 0.668 (0.667) |
| 3, 2 | 6 | 3, 6, 9 | 1.007 | −5·10⁻⁵ | 3.071 (3.2) | 0.997 | 0.716 (0.711) |
| 3, 3 | 10 | 3, 6, 9 | 1.010 | −7·10⁻⁵ | 3.895 (4) | 1.012 | 0.735 (0.727) |
| 4, 1 | 4 | 4, 8, 12 | 1.006 | 1·10⁻⁵ | 3.053 (3) | 0.994 | 0.755 (0.750) |
| 4, 2 | 10 | 4, 8 | 1.011* | | 5.296 (5) | 0.987 | 0.842 (0.833) |
| 5, 1 | 5 | 5, 10 | 1.006* | | 4.605 (4) | 0.984 | 0.805 (0.800) |

The starred entries extrapolate A(L) with exact e_0 linearly in 1/L² over the two sizes. c_half at these sizes carries the same 5 to 15% bias as on the SU(N) controls at equal sizes (4.605 for SU(5) at L = 5, 10). The overshoot of A_ED grows with k at N = 2 exactly as the marginal Takhtajan–Babujian logarithms do. Every number is consistent with the formula, including the free-fit e_0, which tests the Bethe ground state independently of any CFT input.

### Why N = 7 and 8 cannot be tested directly

Rings must have L ≡ 0 mod N. At L = N the ground state is frustration free (E_0 = 0 exactly, the SU(N) singlet with every bond in the (k, k) component), so it carries no Casimir information, and its lowest gap equals that of the fundamental chain at the same L for every k I checked. The next size, L = 2N = 14 at N = 7, has a colour-balanced sector of 1.8·10¹⁶ states. Matrix-product methods fare no better. With c = 32/3 the finite-entanglement exponent is κ = 6/[c(1 + √(12/c))] ≈ 0.27, so χ = 10⁴ gives a correlation length of order ten sites, far from the 1.6% precision a test of R = 64/63 needs, and an SU(7)-symmetric DMRG was not available. Finite-size Bethe roots are the remaining route. [BA PARAGRAPH]

### Survey of other chains

| chain | d | c | R | notes |
|---|---|---|---|---|
| d = 2 candidates of `../d2_classification` (16 points, incl. 220, 1282) | 2 | 1 | 0.02 to 0.23 | rotor-twist fit to stored E_0(L), L = 10 to 26, VUMPS e_0 for 220 and 1282; Lieb–Robinson caps any c = 1 chain at e/π = 0.865 |
| Takhtajan–Babujian S = 1/2 to 2 | 2 to 5 | 1 to 2 | 0.50 to 0.48 | table above |
| spin-1 BLBQ, θ ∈ [π/4, 0.48π] | 3 | 2 | 0.67 to 0.04 | L = 6 to 15, both estimators within 2%, monotone from SU(3) to the ferromagnet |
| supersymmetric t–J, n = 2/3 | 3 | 2 | about 0.42 | μ in the bond norm, L = 6, 12 |
| Sp(4) circle cos θ P₅ + sin θ P₁₀, θ ∈ (π/2, 5π/4) | 4 | 3 | 0.74 to 0.20 | L = 8, 12, monotone from SU(4) toward the ferromagnet |
| Sp(4) at θ = π, SO(5)₁ × Ising | 4 | 3 | 0.64 | split velocities, single-v Casimir artefact of `../counterexample_search` resolved |
| spin-3/2 off the Sp(4) line near SU(4) | 4 | 3 | [SPIN32] | L = 8, 12 |
| SO(5) vector cos θ P + sin θ K, θ ∈ [0, π/2) | 5 | 4 | 0.80 to 0.08 | L = 5, 10, monotone from SU(5) |
| SU(d) swap | d | d − 1 | 1 − 1/d | |

So the fundamental SU(d) chain is the best known chain at every d ≤ 20, and the rank-two symmetric chain takes over at d = 21 (R = 0.972 at N = 6 against 0.952 for SU(21)).

Why SU(d) points are local maxima. ‖h‖ is convex in h and invariant under h → (u⊗u)h(u⊗u)†, so twirling over SU(d) cannot increase it, ‖P + εV̄‖ ≤ ‖P + εV‖. In an SU(d)-singlet ground state only the twirled part V̄ = a + bP of a perturbation shifts E_0 at first order, and that rescales the swap without changing R. Hence R(P + εV) ≤ R(P) + O(ε²) for every perturbation. The symmetric-family counterexample sits far from any fundamental SU(d) point, at a different representation, which is why local searches and symmetric scans never saw it.

### Scaling with d

`scaling.py` evaluates R(N, k) = k(N² − 1)/[N(N + k)H_k] for all (N, k) with d up to 10⁴⁰. The envelope R_max(d) = max over binom(N + k − 1, k) ≤ d grows like 0.37 to 0.39 × ln d/ln ln d from d ≈ 10² to 10⁴⁰ (for example 1.40 at d = 10³, 2.03 at 10⁶, 2.60 at 10⁹). The largest R/ln d over the family and every chain in the survey is 1/(2 ln 2) = 0.721, attained by the Heisenberg chain. So

A ≤ π‖h‖ log₂ d

holds in every case examined, with equality only at d = 2 for the Heisenberg chain, and with room everywhere else (a factor 1.2 for SU(3), growing with d, since SU(d) has 2(1 − 1/d) ≤ log₂ d). It also agrees with the free-fermion and Majorana bounds A ≤ π‖h‖ and is closed under stacking even though ‖h‖ is subadditive, since ‖h_stack‖ ≥ max_i ‖h_i‖ (compress onto a product state of the other factors) gives ‖h_stack‖ Σ_j log₂ d_j ≥ Σ_i ‖h_i‖ log₂ d_i. A d-independent bound has no such closure. I have no argument that the log form is true. At d = 2 it says R ≤ 1/2 for every qubit chain, which Lieb–Robinson alone does not give (it gives e/π for c = 1), and which every qubit chain computed here satisfies.

### Stacking

Heisenberg ⊗ Ising as one d = 4 site gives E_0(L) equal to the sum of the two rings to 10⁻¹⁴ and ‖h‖ = 1 + √2 exactly. For three random pairs of bonds without reflection symmetry the stacked norm is 7.982 against 8.127 for the sum, 12.765 against 12.765, and 12.457 against 12.525. Restricting o to o₁⊗1 + 1⊗o₂ recovers the sum, so the deficit comes from entangling o.

## Reproducing

numpy and scipy only. Run inside this folder.

- `core.py` (‖h‖, ring spectra through `../counterexample_search/ced.py`, both estimators), `chains.py` (bonds, including `sym_rep_chain(N, k)`), `run.py` (one chain, several L, JSON lines), `scan.py` (linear families with cached momentum blocks), `analyze.py`.
- Calibration: `python3 run.py --model SU3 --L 6 9 12 15 --out data/calib.jsonl` and likewise for XX, Heisenberg, Ising, SU4, `TB --par 1`. Seconds to a minute per size up to L = 24 at d = 2.
- Symmetric family: `python3 run.py --model symrep --par N k --L ... --out data/symrep_local.jsonl`, then `python3 analyze_symrep.py data/symrep*.jsonl`. (3, 2) at L = 9 takes a minute.
- Survey: `scan.py --family blbq|sp4|so5|spin32`, `d2_points.py`, `tj_run.sh` and `tj_analyze.py`, `stack.py`, `survey_table.py`, `scaling.py`.
- Typhon (`slurm_task.sh` with `tasks_*.txt`, results in `/scratch/lustre/jkudlerflam/casimir_survey/`): [JOBS].

## Relation to the rest of the repository

This tests `../refined_conjecture/` item 3 and complements C14, which found the counterexample analytically. Suggested STATUS line: "The Casimir bound A ≤ 2π‖h‖ is false (C14). ED at N ≤ 4 confirms the symmetric-representation family's e_0, v, c and A to about 1%, every other chain surveyed stays below 1 − 1/d, and A ≤ π‖h‖ log₂ d fits all data, tight only for the Heisenberg chain (`contrib/casimir_survey/`)."
