# Casimir survey: A = Σ c_a v_a against 2π‖h‖

## Claim

Class and notation as in `../refined_conjecture/`: periodic rings of (ℂ^d)^⊗L with H = Σ h_{i,i+1}, unitary CFTs (c_a, v_a) in the infrared, A = Σ c_a v_a from E_0(L) = e_0 L − πA/(6L) + o(1/L), ‖h‖ the half spread of spec h minimized over h → h + o⊗1 − 1⊗o, and R = A/(2π‖h‖). The refined conjecture is R ≤ 1.

1. The integrable SU(N) chains with Sym^k(ℂ^N) on each site, found by Astra's C14 run (`../../codex/work_C14/RESULT.md`), violate it. ED at (N, k) = (2, 1 to 4), (3, 1 to 3), (4, 1 and 2) and (5, 1) reproduces the Bethe/WZW e_0, v, c and A to about 1%, and at N = 3 and 4 shows rank two beating the fundamental chain by the predicted margin. For the violating members, exact finite-size Bethe roots validated against this ED give A directly from E_0(L) on rings up to L = 105 and 120: R = 1.0156 at N = 7, d = 28 (1.0139 to 1.0177 over fit forms, predicted 64/63) and R = 1.0497 at N = 8, d = 36 (1.0481 to 1.0516, predicted 1.05). No WZW input enters this last step. [numerical]
2. Outside that family nothing examined exceeds the SU(d) value 1 − 1/d. Every family through an SU(d) point peaks there, as a twirling argument requires to first order. [numerical]
3. A ≤ π‖h‖ log₂ d holds for every chain computed and the whole symmetric family, with equality only for the Heisenberg chain. The symmetric family forces any bound A ≤ K(d)‖h‖ to grow at least like log d/log log d. [numerical and analytic]
4. A adds exactly under stacking, but ‖h‖ is only subadditive. It adds for reflection-symmetric bonds and can fall strictly below the sum otherwise. [numerical]

## Status

Numerical. At N ≥ 5 item 1 assumes that the continued Bethe root configuration stays the ground state, which ED confirms at every size it reaches. An SU(N)-symmetric DMRG at N = 7 would remove that assumption.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the coordinating session, 2026-10-27.
- AI assistance: gpt-6-astra (high) through `gpt_solve.sh`, two closed-book runs that wrote the Bethe solvers in `astra_bethe/` (specs included). The second ended with a backend capacity error after writing its script and results, so it left no report. I reran its validation (11 ED energies, to 4·10⁻¹¹) and its N = 7 production (identical energies, root residuals below 2·10⁻¹³). The fits (`bethe_refit.py`), all ED and the survey are mine.
- Checks: `core.hnorm` agrees with Powell searches on random bonds to 10⁻⁹. `chains.sym_rep_chain` (Schwinger bosons, Lagrange interpolation in X = Σ E_ab ⊗ E_ba) equals the Takhtajan–Babujian bond at N = k = 2 to 4·10⁻¹⁶, and C14's bond is h_C14 = 2h − 2H_k, so R agrees between normalizations. The ED engine is `../counterexample_search/ced.py`. Nobody else has rerun this.

## Method and calibration

A_cas fits E_0(L)/L = e_0 − πA/(6L²) + b/L⁴ over commensurate ring lengths. A_cv is c_half = 3[S(L₂/2) − S(L₁/2)]/ln(L₂/L₁) times v = L[E(k_0 ± 2π/L) − E_0]/2π at the largest L. With exact e_0 I also use A(L) = −6L[E_0(L) − e_0L]/π.

| chain | c | sizes | A_cas (exact) | c_half | v (exact) | R_cas, R_cv |
|---|---|---|---|---|---|---|
| XX | 1 | 10 to 24 | 1.9998 (2) | 0.996 | 1.994 (2) | 0.318, 0.316 |
| Ising | 1/2 | 10 to 22 | 1.0000 (1) | 0.499 | 1.993 (2) | 0.113, 0.112 |
| Heisenberg, h = P | 1 | 10 to 24 | 3.156 (π) | 1.001 | 3.168 (π) | 0.502, 0.504 |
| SU(3), h = P | 2 | 6 to 15 | 4.214 (4.189) | 1.948 | 2.098 (2.094) | 0.671, 0.650 |
| SU(4), h = P | 3 | 4, 8, 12 | 4.743 (4.712) | 3.053 | 1.561 (1.571) | 0.755, 0.759 |
| Takhtajan–Babujian S = 1 | 3/2 | 8 to 16 | 2.380 (2.356) | 1.505 | 1.595 (1.571) | 0.505, 0.509 |

A_cas is exact to 10⁻⁴ without marginal operators and 0.5 to 1% high with them. Two-size fits (L = 8, 12 at d = 4, L = 5, 10 at d = 5) overshoot SU(4) and SU(5) by 5.3% and 11.7%, and I divide by those factors where only two sizes exist.

## Results

### The symmetric family by ED

In my normalization h = Σ_j H_{k−j} Π_j (H_m harmonic numbers, Π_j the projector on the diagram (2k − j, j)), ‖h‖ = H_k/2, v = π/N, c = k(N² − 1)/(N + k), and e_0 is C14's eq. 5 shifted.

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

Starred entries extrapolate A(L) with exact e_0 linearly in 1/L². c_half carries the bias of the SU(N) controls at equal sizes. The A_ED overshoot grows with k as the marginal Takhtajan–Babujian logarithms do. The free-fit e_0 tests the Bethe ground state without any CFT input.

### The violating members, N = 7 and 8

Rings need L ≡ 0 mod N. At L = N the ground state is frustration free (E_0 = 0), so it carries no Casimir information, and at L = 14 the colour-balanced sector of N = 7 has 1.8·10¹⁶ states. Matrix-product states are no better: at c = 32/3 the finite-entanglement exponent is about 0.27, so χ = 10⁴ reaches a correlation length of order ten sites, far from the 1.6% a test of 64/63 needs.

So I used Bethe roots. Astra's first solver, with ideal k-strings, misses ED by 0.4/L to 0.7/L in E_0 at k ≥ 2, which corrupts the Casimir term, and was discarded. The second (`astra_bethe/bethe_exact.py`) solves the full nested equations for complex roots, with strings only as seeds, and matches all eleven ED energies at k = 2, 3, 4 to 4·10⁻¹¹. Continuing that configuration in L, I fit A(L) with the Bethe e_0 (which free fits confirm to 10⁻⁸):

| N, k | d | L | R(L_max) | R, 1/ln³L + 1/L² fit | range over forms | predicted |
|---|---|---|---|---|---|---|
| 2, 2 | 3 | 4 to 60 | 0.5023 | 0.5003 | 0.4998 to 0.5019 | 0.5 |
| 3, 2 | 6 | 6 to 90 | 0.7133 | 0.7112 | 0.7106 to 0.7130 | 0.7111 |
| 4, 2 | 10 | 8 to 120 | 0.8353 | 0.8333 | 0.8327 to 0.8350 | 0.8333 |
| 5, 2 | 15 | 10 to 105 | 0.9166 | 0.9142 | 0.9130 to 0.9161 | 0.9143 |
| 7, 2 | 28 | 14 to 105 | 1.0183 | 1.0156 | 1.0139 to 1.0177 | 1.0159 |
| 8, 2 | 36 | 16 to 120 | 1.0522 | 1.0497 | 1.0481 to 1.0516 | 1.05 |

The forms are A_∞ + a/L² (+ b/L⁴), A_∞ + a/ln³L and A_∞ + a/ln³L + b/L², each over the upper half and upper two thirds of the sizes. The mixed form hits the known answer to 3·10⁻⁴ in R for N ≤ 5, and the others bracket it. Every fit at N = 7 and 8 lies above 1.

### Other chains

| chain | d | c | R | notes |
|---|---|---|---|---|
| 16 candidates of `../d2_classification` (incl. 220, 1282) | 2 | 1 | 0.02 to 0.23 | rough rotor-twist fits to stored E_0(L), L = 10 to 26 |
| Takhtajan–Babujian S = 1/2 to 2 | 2 to 5 | 1 to 2 | 0.50 to 0.48 | first table |
| spin-1 BLBQ, θ from π/4 to 0.48π | 3 | 2 | 0.67 to 0.04 | L = 6 to 15, estimators within 2% |
| supersymmetric t–J, n = 2/3, μ in the norm | 3 | 2 | about 0.42 | L = 6, 12 only |
| Sp(4) circle cos θ P₅ + sin θ P₁₀, θ from 0.52π to 1.17π | 4 | 3 | 0.74 to 0.20 | L = 8, 12 |
| Sp(4) at θ = π (SO(5)₁ × Ising) | 4 | 3 | 0.64 | split velocities |
| spin 3/2, a₁ − a₃ = −0.1, −0.3, −0.6 off SU(4) | 4 | 3 | 0.73, 0.67, 0.58 | L = 8, 12 |
| SO(5) vector cos θ P + sin θ K, θ from 0 to 0.42π | 5 | 4 | 0.80 to 0.08 | L = 5, 10 |
| SU(d) swap | d | d − 1 | 1 − 1/d | |

Each family falls monotonically away from its SU(d) point. Lieb–Robinson caps every c = 1 chain at e/π = 0.865. The fundamental SU(d) chain is the best known chain for d ≤ 20, and the rank-two chain overtakes it at d = 21 (0.972 against 0.952).

Why SU(d) points are local maxima. ‖h‖ is convex and invariant under h → (u⊗u)h(u⊗u)†, so twirling cannot increase it, ‖P + εV̄‖ ≤ ‖P + εV‖. In an SU(d)-singlet ground state only V̄ = a + bP shifts E_0 at first order, and that leaves R unchanged. So R(P + εV) ≤ R(P) + O(ε²). The counterexample lives in another representation, which is why local scans never saw it.

### Scaling with d

`scaling.py` evaluates R(N, k) = k(N² − 1)/[N(N + k)H_k] for d up to 10⁴⁰. The envelope max over binom(N + k − 1, k) ≤ d grows like 0.37 to 0.39 × ln d/ln ln d (1.40 at d = 10³, 2.03 at 10⁶, 2.60 at 10⁹). The largest R/ln d anywhere, family or survey, is 1/(2 ln 2), from the Heisenberg chain. So A ≤ π‖h‖ log₂ d holds in every case, tight only at d = 2, with a factor 1.2 to spare for SU(3) and more beyond. It contains the free-fermion and Majorana bounds A ≤ π‖h‖, and it is closed under stacking although ‖h‖ is subadditive, because compressing onto a product state of the other factors gives ‖h_stack‖ ≥ max_i ‖h_i‖. A d-independent bound has no such closure. I have no argument that the log form is true. At d = 2 it says R ≤ 1/2 for every qubit chain, stronger than Lieb–Robinson, and every qubit chain here satisfies it.

### Stacking

Heisenberg ⊗ Ising as one d = 4 site has E_0(L) equal to the sum of the two rings to 10⁻¹⁴ and ‖h‖ = 1 + √2 exactly. For three random pairs of bonds without reflection symmetry the stacked norm is 7.982 against a sum of 8.127, 12.765 against 12.765, and 12.457 against 12.525. Restricting o to o₁⊗1 + 1⊗o₂ restores the sum, so entangling o causes the deficit.

## Reproducing

numpy and scipy only, run inside this folder. `core.py` (‖h‖, ring spectra through `../counterexample_search/ced.py`, estimators), `chains.py` (bonds), `run.py` (one chain, several L), `scan.py` (linear families), `analyze.py`, `analyze_symrep.py`, `survey_table.py`, `d2_points.py`, `tj_run.sh` and `tj_analyze.py`, `stack.py`, `scaling.py`, `bethe_refit.py`. Examples: `python3 run.py --model symrep --par 3 2 --L 3 6 9 --out data/symrep_local.jsonl` (a minute), `python3 astra_bethe/bethe_exact.py --validate-only`, then `--families 7 --max-multiple 15` (seconds). Typhon job 17513209 (task 2) ran (3, 3) at L = 9 (19 minutes, 7.7 GB). Larger calibration sizes and (3, 2) at L = 12 sat in the queue and were cancelled, and nothing here depends on them.

## Relation to the rest of the repository

This tests `../refined_conjecture/` item 3 alongside C14. Suggested STATUS line: "A ≤ 2π‖h‖ is false. Exact finite-size Bethe roots give R = 1.016 at d = 28 and 1.050 at d = 36 for the rank-two symmetric SU(N) chains, ED confirms that family to about 1% at N ≤ 4, no other chain surveyed exceeds 1 − 1/d, and A ≤ π‖h‖ log₂ d fits everything, tight only for the Heisenberg chain (`contrib/casimir_survey/`)."
