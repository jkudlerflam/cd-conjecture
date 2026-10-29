# Attacking A <= pi M log2 d: stacking, grading, spin-charge separation and bond optimization

## Claim

Class. One-site translation-invariant nearest-neighbour chains H = sum_i h_{i,i+1} on (C^d)^L, h Hermitian, whose low-energy theory is a set of CFT sectors (c_a, v_a). A = sum_a c_a v_a is the Casimir coefficient, E_0(L) = e_0 L - pi A/(6L) + o(1/L) on a ring, and M = min_o (1/2)[lambda_max - lambda_min](h + o (x) 1 - 1 (x) o). The conjecture under attack is

    R := A/(pi M log2 d) <= 1.

Findings.

1. [proved] The "many commuting fast modes" mechanism cannot work. For a stack h = h_1 (x) 1 + 1 (x) h_2 on d = d_1 d_2, A adds while M(h) >= max(M_1, M_2) (compression onto a product state of the other factor), so a stack of chains obeying the bound obeys it with strict slack. Twists h -> (1 (x) w) h (1 (x) w^dag) keep the bulk A and never lower M for bonds whose extremal eigenspaces are the symmetric and antisymmetric subspaces (Heisenberg, every SU(N) permutation chain). Blocking b sites divides R by b^2 M'/M >= b^2 M_2(h)/M. Lemmas and proofs in `PROOFS.md`, numerical checks in `lemmas_check.py` and `data/lemmas_check.json`.
2. [proved] Equality at d = 2 is not isolated. Every twisted Heisenberg bond -(1 (x) w)|s><s|(1 (x) w^dag), w in SU(2), none of which is SU(2)-invariant as a bond, has A = pi/2 and M = 1/2, hence R = 1. The d = 2 optimizer rediscovers this family from random rank-one projector starts.
3. [numerical] Graded, spin-charge-separated and spin-orbital chains all stay below the SU(d) permutation chain at their d. Supersymmetric t-J (d = 3): R = 0.55 at electron density 2/3, once the chemical potential that selects the density is charged to M (0.64 if it is not). The SU(3) t-J analogue gl(3|1) at d = 4: 0.54 to 0.56. gl(2|2): 0.36. Hubbard at d = 4: at most 0.34 (U = 0, half filling), falling with U at every filling. Kugel-Khomskii spin-orbital bonds and coupled Heisenberg pairs peak at the SU(4) point, R = 0.755 at L = 4, 8, 12 (0.75 exact). The Hermitian OSp(1|2)-invariant bonds are only span{1, P_g} (the third invariant is not Hermitian), so there is no OSp(1|2) circle to scan; the d = 3 "pairing" circle cos(theta) P_g + sin(theta) K that was scanned instead is critical nowhere with a CFT-shaped E_0(L).
4. [numerical] Direct maximization of R over bonds at d = 2 (12 real or 15 complex parameters after removing the identity and the gauge, 80 Nelder-Mead starts from random bonds, random projectors, XXZ points and noisy Heisenberg points, rings L = 8 to 16) never exceeds the Heisenberg value under the same finite-size estimator: best 1.0138 against the Heisenberg reference 1.0138 (L = 8 to 14), reached from a random reflection-symmetric start and from a random projector start, and verified on rings to L = 18 to be a twisted Heisenberg chain (c_half -> 1, ratio 1.002 by the full fit, identical two-point values). At d = 3 (36 real parameters, L = 6 to 12) nothing exceeds the SU(3) value (reference score 0.860 at L = 6, 9, 12; best found 0.66 from noisy SU(3) starts, 0.23 from spin-1 XXZ starts). Every inflated candidate of an earlier run that used the plain least-squares ratio failed the drift test (two-point estimates A(L, L+2) not constant in L).
5. [heuristic] The largest R over every chain known to the repository is 1, at d = 2, and the envelope of known chains falls with d: 0.84 (d = 3, SU(3)), 0.75 (d = 4, SU(4)), 2(1 - 1/d)/log2 d along the SU(d) chains, then about 0.76 ln d/(ln ln d . log2 d) along the rank-two symmetric SU(N) chains (0.42 at d = 28). The data also admit the tighter candidate A <= pi M H_{d-1}, H_{d-1} = sum_{j<d} 1/j: equal to log2 d at d = 2, below it for all d >= 3, superadditive (so stacking-closed by Lemma 1), satisfied with slack 0.89 at d = 3 and 0.82 at d = 4 by the SU(d) chains, and generous by a factor ln ln d against the symmetric family at large d (`candidate_f.py`, `data/candidate_f.json`). The d-independent candidate 2(1 - 1/d) fails from d = 21 on.

The conjecture survives this attack. Largest ratio found: R = 1.000, the Heisenberg chain and its twisted relatives at d = 2. Largest at d >= 3: the SU(3) permutation chain, R = 4/(3 log2 3) = 0.841.

## Status

Items 1 and 2 proved. Items 3 and 4 numerical: exact diagonalization on rings up to L = 18 (d = 2), 12 (d = 3) and 12 (d = 4), A from fits of E_0(L); the bond search is a finite-size local optimization and does not exclude narrow maxima it never visited. Item 5 heuristic. Two cluster jobs with a larger search (d = 2 with L up to 16, d = 3, more seeds, and the remaining d = 4 points at L = 12; task lists in `typhon/`) were submitted but had not finished when this note was written; their outputs are not included, and nothing here depends on them. What would upgrade the result: a proof of A <= pi M at d = 2 (where the bound is sharp), or a small-d analogue of the rank-k symmetric chains with velocities above the SU(d) ones.

## How it was established and checked

- Author: Claude (Anthropic, Fable 5.1), an adversarial subagent of the coordinating session, 2026-10-29.
- AI assistance: none beyond the author. No Codex run.
- Checked: the ED engines reproduce the Heisenberg chain (R = 1.005 at L <= 16 by the full fit, 1.000 exact), the ULS chain (0.846 at L = 6, 9, 12 against 0.841), SU(4) (0.755 at L = 4, 8, 12 against 0.750), free fermions at U = 0 (ring energies exact with the parity-dependent boundary condition) and the survey's t-J value. The lemmas were checked numerically on random bonds. The symmetric engine is the repository's `../counterexample_search/ced.py`. Nobody else has rerun this.

## Method

Bond scale. `core.hnorm` minimizes the log-sum-exp smoothed spread over o by BFGS at six annealed inverse temperatures (the spread is convex in o) and returns the exact spread at the optimum; this is the survey's algorithm, re-implemented. For reflection-symmetric bonds with simple extremes the gauge is provably trivial (Lemma 2), which the code confirms; for random bonds without reflection symmetry it lowers M by up to 40 percent.

Casimir coefficient. Symmetric families use the charge- and momentum-resolved engine through `symring.py`; generic bonds use the matrix-free full-space Lanczos in `core.Ring`. A is fitted from E_0(L)/L = e_0 - pi A/(6 L^2) + b/L^4 on commensurate rings, with the two-size fit as a second estimate and the half-ring entanglement c_half and the gap times L as criticality diagnostics. For fermionic chains at fixed density the canonical sector is used, ring fluxes 0 and pi are compared and the lower energy taken (closed shells), and the chemical potentials mu_a = [E(N_a + 1) - E(N_a - 1)]/2 at the largest ring are added to the bond before computing M, so that the density studied is the ground state of the bond that is charged (`canonical.py`). Calibration: the two-size d = 4 fits at L = 4, 8 overshoot by about 14 percent (0.856 against 0.75 for SU(4)); L = 4, 8, 12 is accurate to 1 percent.

Bond optimization. `optimize_random.py` parametrizes h in an orthonormal Hermitian product basis with the identity removed and the one-site part symmetrized (the antisymmetric one-site part is the gauge and changes neither the ring energies nor M). The score is the smallest two-point Casimir ratio A(L_i, L_{i+1})/(pi M log2 d) over consecutive ring sizes minus half the spread of those estimates. This penalizes level crossings, incommensurate ground states and gapped chains with inflated least-squares fits; an earlier run with the plain fitted ratio was driven to such artifacts (fitted "R" up to 3.7 with two-point values 1.1, 4.2, 7.2), and its logged points were re-scored (`rescore.py`, `data/old_score/`). Nelder-Mead from starts of several kinds (random bond; minus a Haar-random projector of random rank; spin-S XXZ with random anisotropy; Heisenberg or SU(d) permutation with 30 to 60 percent relative noise), then a local polish. Reference scores under the same estimator: Heisenberg 1.0138 (L = 8 to 14) and 1.0062 (L = 8 to 16); SU(3) 0.860 (L = 6, 9, 12).

## Results

### Largest ratio per mechanism

| mechanism | best chain | d | R | tag |
|---|---|---|---|---|
| stacking of saturating chains | Heisenberg (x) Heisenberg | 4 | 0.500 (= 1/2 exactly, M = 1) | proved |
| twists of the Heisenberg bond | -(1 (x) w)|s><s|(1 (x) w^dag), any w in SU(2) | 2 | 1.000 | proved |
| blocking, b = 2, Heisenberg | blocked bond, M' = 3/4 | 4 | 0.168 (fit) = 1/6 (Lemma 4) | proved + numerical |
| graded gl(N|M) permutation chains with chemical potentials | t-J at n = 2/3 | 3 | 0.55 | numerical |
| | gl(3|1), one fermion per four sites | 4 | 0.54 to 0.56 | numerical |
| | gl(2|2), quarter filling each | 4 | 0.36 | numerical |
| spin-charge separation | Hubbard U = 0, half filling | 4 | 0.34 (exact 0.318) | numerical |
| spin-orbital, Kugel-Khomskii (S.S + x)(T.T + y) | SU(4) point x = y = 1/4 | 4 | 0.755 (exact 0.750) | numerical |
| coupled Heisenberg pairs (zz, rung, zigzag cross) | SU(4) point again (zz, lambda = 4) | 4 | 0.75 | numerical |
| spin-1 bilinear-biquadratic circle | ULS point theta = pi/4 | 3 | 0.846 (exact 0.841) | numerical |
| spin-1 XXZ with single-ion D (XY phases) | delta = 0.5, D = 0.5 or 1.0 | 3 | 0.31 | numerical |
| bond optimization, d = 2 (real, complex, reflection-symmetric) | twisted Heisenberg | 2 | 1.000 | numerical |
| bond optimization, d = 3 (real) | below SU(3) | 3 | < 0.84 | numerical |

### Selected family points

R by the full fit and by the two largest sizes; c_half from the two largest rings; gap times L at each ring. Full table from `summarize.py` (`data/summary.json`).

| chain | d | rings | M | R (fit, two-size) | c_half | gap x L |
|---|---|---|---|---|---|---|
| ULS, theta = pi/4 | 3 | 6, 9, 12 | 0.707 | 0.846, 0.876 | 2.11 | 5.94, 5.93, 5.94 |
| BLBQ theta = 0.9 | 3 | 6, 9, 12 | 0.864 | 0.639, 0.661 | 2.11 | 5.16, 5.11, 5.10 |
| BLBQ theta = 1.2 | 3 | 6, 9, 12 | 1.217 | 0.291, 0.300 | 2.14 | 3.00, 2.99, 3.00 |
| t-J, n = 1/3 | 3 | 6, 12 | 1.073 | 0.383, 0.384 | 2.12 | 6.00, 5.59 |
| t-J, n = 1/2 | 3 | 4, 8, 12 | 1.145 | 0.475, 0.490 | 2.02 | 8.00, 7.07, 6.84 |
| t-J, n = 2/3 | 3 | 6, 12 | 1.229 | 0.546, 0.549 | 1.98 | 7.34, 6.66 |
| t-J, opposite sign, n = 2/3 | 3 | 6, 12 | 1.250 | 0.288, 0.290 | 2.31 | 5.46, 2.91 |
| spin-1 XXZ+D, delta = 0.5, D = 0.5 | 3 | 6 to 12 | 1.414 | 0.305, 0.311 | 0.98 | 2.15, 2.14, 2.14, 2.14 |
| SU(4) = KK x = y = 1/4 | 4 | 4, 8, 12 | 0.25 | 0.755, 0.789 | 3.05 | 2.00, 1.85, 1.83 |
| KK x = y = 1/2 | 4 | 4, 8, 12 | 0.375 | 0.379, 0.469 | 1.42 | 1.39, 0.70, 0.31 |
| KK x = 0, y = 1 | 4 | 4, 8, 12 | 0.625 | 0.503, 0.521 | 2.32 | 1.00, 0.64, 0.40 |
| two Heisenberg chains, zz coupling lambda = 2 | 4 | 4, 8 | 0.75 | 0.675, 0.690 | 2.94 | 2.79, 1.39 |
| zigzag cross coupling lambda = 0.25 | 4 | 4, 8, 12 | 1.076 | 0.293, 0.313 | 1.05 | 3.18, 3.08, 3.27 |
| gl(3|1), n_f = 1/4 | 4 | 4, 8, 12 | 1.142 | 0.541, 0.559 | 3.09 | 8.00, 5.38, 4.72 |
| gl(2|2), n = 1/4 each | 4 | 4, 8 | 1.309 | 0.357, 0.365 | 3.44 | 8.00, 4.21 |
| Hubbard U = 0, n = 1 | 4 | 4, 8 | 2.0 | 0.341, 0.349 | 1.82 | 11.3, 12.3 |
| Hubbard U = 4, n = 1 | 4 | 4, 8, 12 | 2.236 | 0.153, 0.183 | 1.16 | 5.03, 4.26, 3.81 |
| Hubbard U = 4, n = 1/2 | 4 | 4, 8, 12 | 3.029 | 0.334, 0.237 | 2.03 | 4.49, 5.62, 4.14 |

The KK points away from SU(4) have gap times L falling with L (gapped or dimerizing), so their R values are upper estimates. For the graded chains the quoted M includes the chemical potentials; without them (M_mu0) the t-J value would be 0.64 and gl(3|1) 0.67, still below SU(3) and SU(4).

### Bond optimization, d = 2

71 finals in `data/opt_d2_*.jsonl` (real, complex, reflection-symmetric, L up to 14 or 16). Scores (drift-penalized, same estimator as the Heisenberg reference 1.0138 at L = 8 to 14):

| start type | runs | best scores |
|---|---|---|
| noisy Heisenberg, real | 14 | 1.013, 1.013, 1.012, 1.012, 1.012 |
| noisy Heisenberg, reflection-symmetric | 1 | 1.0138 |
| random projector, real | 24 | 1.011, 1.010, 0.46, 0.39, 0.37 |
| random projector, reflection-symmetric | 3 | 1.0138, 0, 0 |
| random projector, complex | 6 | 0.91, 0.87, 0, 0, 0 |
| random bond, real | 8 | 0.47, 0.41, 0.02, 0, 0 |
| random bond, reflection-symmetric | 1 | 1.0138 |
| random bond, complex | 2 | 0.46, 0 |
| spin-1/2 XXZ with random anisotropy | 8 | 1.012, 1.011, 0.94, 0.59, 0.33 |

Verification on rings L = 8 to 18 (`verify_candidates.py`, `data/verify_local_d2.json`) of the three best distinct bonds: full-fit ratios 1.0029, 1.0017, 1.0014, two-point sequences 1.044, 1.031, 1.023, 1.017, 1.014 (Heisenberg: 1.047, 1.033, 1.025, 1.020, 1.014), c_half 0.981 to 0.998, i.e. all three are the Heisenberg chain up to a twist and a local unitary. Scores of exactly 0 are gapped chains (random bonds and generic projectors are gapped), where Nelder-Mead has no gradient; the search is therefore a search among the critical points it can reach from its starts.

### Bond optimization, d = 3

Three runs (`data/opt_d3_*.jsonl`): noisy SU(3) starts at L = 6, 9, 12 reach 0.66 (SU(3) itself scores 0.860 there); spin-1 XXZ starts at L = 6 to 12 reach 0.23. Nothing approaches the SU(3) value, let alone 1.

## Reproducing

numpy and scipy only, run inside this folder (the symmetric engine is imported from `../counterexample_search/ced.py`).

- `lemmas_check.py` (two minutes): the checks of `PROOFS.md`, into `data/lemmas_check.json`.
- `scan.py --family <name> --par k=v ... --Ls ...` and `canonical.py --family <name> --dens ... --Ls ...` for the families in `families.py` (`tj`, `osp12` pairing circle, `blbq`, `spin1xxzd`, `hubbard`, `kk`, `twochain`, `gperm`); `tasks_local_d3.sh`, `tasks_local_d4.sh`, `tasks_local_d4_L12.sh` reproduce `data/families_*.jsonl` (minutes per point, up to an hour per L = 12 point at d = 4). `summarize.py` prints the table.
- `optimize_random.py --d 2 --Ls 8 10 12 14 --starts 12 --evals 400 --init mix --seed <s> --out ...` (three to four minutes per start at d = 2, L <= 14; `--complex`, `--reflect`, `--init heis|su|xxz|proj|random|mix`). `rescore.py 'data/opt_d2_*.jsonl'` re-scores every logged point; `verify_candidates.py <rescore.json> --d 2 --Ls 8 10 12 14 16 18` recomputes the best bonds with gaps and c_half. `compact_logs.py` trims logs to finals plus the last improvements.
- `candidate_f.py`: the envelope of the known chains against log2 d, H_{d-1} and 2(1 - 1/d).
- `typhon/`: task lists and the one-node SLURM wrapper for the larger cluster runs (not completed at writing).

## Relation to the rest of the repository

Tests the velocity-weighted bound of `../refined_conjecture/` (refined in `../casimir_survey/` item 3 and `../../specs/_bound_background.md`) from the adversarial side, complementing C101's survey of integrable classes with structural lemmas (stacking, twisting and blocking cannot raise R) and a direct optimization over bonds at d = 2 and 3. Suggested STATUS line: "Adversarial search for A > pi M log2 d (stacks, twists, blockings, graded and spin-charge chains, bond optimization at d = 2, 3): nothing above the Heisenberg value R = 1, which is attained on the whole family of twisted Heisenberg bonds; the envelope of known chains falls with d and also fits the tighter A <= pi M H_{d-1} (`contrib/fable_bound_attack/`)."
