# Reproducing the finite claims behind `note.pdf`

This is an ordered checklist of the scripts that regenerate every load-bearing finite claim in the note. It follows the note's structure, Claim A first, then Claim B, then the small-d mechanism checks. The infrared input S, and inputs P and M of the d = 24 chain, are not finite claims and no script tests them.

Everything was produced with AI assistance and has not been refereed by a human. The scripts are the authors' own, and where a referee rebuilt something independently that script is listed separately.

## Before you start

Requirements are Python 3 with numpy and scipy (the exact certificates use the standard library only). The commands below were tested on 2026-10-28 with Python 3.12.3, numpy 2.4.6 and scipy 1.17.1 on an Apple M1 Pro, one BLAS thread per process, with up to four test batches running at once. Times marked "tested" are wall-clock times from that run. Times marked "quoted" come from the source folder's README and were not rerun for this checklist.

Most scripts write their output files (`.json`, `.npz`) into their own folder and overwrite the committed copies. To leave the repository untouched, work in a copy.

```bash
REPO=/path/to/cd-conjecture/notes/cd_conjecture      # the folder containing STATUS.md
WORK=$(mktemp -d)/cd && rsync -a --exclude .git "$REPO/" "$WORK/" && cd "$WORK"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
```

Afterwards `diff -rq "$REPO" "$WORK"` lists the regenerated files. In the test run every regenerated output that was compared (25 files) agreed with the committed copy except for timing fields and last-digit roundoff (at most 4e-15).

Paths. Run every command from the folder named in its step (`cd "$WORK/<folder>"`). This is required for all scripts under `codex/work_C*/`, `contrib/referee_C39/` and `contrib/referee_C44_C45/`, which read and write relative to the current directory. Tested failure modes from the repository root were `python3 -B contrib/referee_C39/ref_bond.py` (FileNotFoundError, it reads `../../codex/work_C39/bond_n*.npz` relative to the current directory), and `codex/work_C39/certify.py 1`, `codex/work_C45/checks.py` and `contrib/referee_C44_C45/ref_schedule.py`, which ran but wrote their output into the root. Scripts under `contrib/counterexample_campaign/`, `contrib/referee_d24/`, `contrib/routeA_final/`, `contrib/routeA_demonstrator/`, `contrib/vacancy_mobile_numerics/` and `contrib/margins_audit/` locate themselves through `__file__` and run from any directory (tested for `rebuild_bond.py`, `infrared/checks.py`, `doping/bond24.py`, `doping_check/check.py` and `mobile_modular.py`), but they read sibling folders, so the repository layout must be intact.

## A. Claim A, unbounded c at d0 = 4464

### A1. The all-n band lemma, exact certificates (C39)

Folder `codex/work_C39/`.

- [ ] `python3 -B certify.py 1`, `python3 -B certify.py 2`, `python3 -B certify.py 3`. Tested 1 s, 2 s, 6 s. Exact Fraction arithmetic, Sturm isolation of all 2n roots and rational comparison of every pair sum at w = 1/4, v = 1/1024, u = 1/64, eta = 1/100, mu = 10^-(8n+10). Expected printed `delta4_display` 5.1197e-16, 1.2364e-28, 2.3631e-41, and `delta2_display` 3.1416, 0.016609, 1.4260e-14. Writes `certificate_n*.json`.

### A2. Independent high-precision margins (referee of C39)

Folder `contrib/referee_C39/`.

- [ ] `python3 -B ref_pencil.py fixed 1,2,3,4 500`. Tested 8 s. Expected one JSON line per n with `"delta4"` 5.119744e-16, 1.236376e-28, 2.36e-41, 4.2e-54, `"lacunary_ok": true`, `"delta4_ge_bound17": true`, and `r_over_kappa_v4j` equal to 1 within 2e-6. Writes `fixed_1_2_3_4.json`.
- [ ] `python3 -B ref_pencil.py fixed 8 900`. Tested 114 s. Expected `"delta4": "3.087374e-105"`, `"lacunary_ok": true` and `"delta4_ge_bound17": true`, the certified tuple carried to n = 8. Writes `fixed_8.json`.
- [ ] `python3 -B ref_pencil.py search 8 700`. Tested 54 s. Expected the search to stop at `"eta": "1/32"`, `"u": "1/32"`, `"v": "1/4"` with `"delta4": "1.948375e-21"`. Writes `search_n8.json`.

### A3. The 46-letter bond and its second-order coefficient

Folder `codex/work_C39/`.

- [ ] `python3 -B literal.py`. Tested 1 s. Expected three lines with `'alphabet': 46`, `'bond_dimension': 2116`, `'core_dimension'` 54, 160, 424, `'entry_error'` at most 6e-16, `'charge_commutator': 0.0` and `'code_residual'` at most 4e-16.
- [ ] `python3 -B checks.py`. Tested 3 s. Expected JSON with Hermiticity and charge-commutator errors near 1e-16 and the Gram normalization check `relative_error` 0.979, 0.333, 0.038 at D = 10, 100, 1000 (it converges like 1/D). Writes `bond_checks.json`.

Folder `contrib/referee_C39/` (rebuilt from the rule table, independent of C39's code).

- [ ] `python3 -B ref_bond.py`. Tested 3 s. Expected for n = 3 `"max_abs_diff"` 6.6e-17 against C39's `bond_n3.npz` and `"hop_rel_err"` 9.1e-10. Writes `ref_bond_compare.json`.
- [ ] `python3 -B ref_largeD.py`. Tested 5 s. Expected `rel_err` between 1.4e-9 and 3.5e-9 for the g^2 hopping formula at D = 0.2 and 2, and `times_D` approaching 22.5 (n = 1) and 38 (n = 2), so the Gram limit is approached like 1/D.
- [ ] `python3 -B ref_gauss.py`. Tested 3 s. Expected, for n = 1, 2, 3, an order-g^0 part that is a scalar only and an order-g^2 part that is quadratic with `nonGaussian` at the 1e-21 level.
- [ ] `python3 -B ref_g4scale.py`. Tested 10 s. Expected `nonGauss4_over_hop2` near 0.07 to 0.1 at n = 1, 2 and growing (3.3, 13.7, 79.6 at D = 2, 8, 32) at n = 3, the scale problem the referee flags.

### A4. Finite D and g for every n (C42)

Folder `codex/work_C42/`.

- [ ] `python3 -B budgets.py 1`, `python3 -B budgets.py 2`, `python3 -B budgets.py 3`. Tested under 1 s each (n = 1, 3). Expected `'g_dyadic_exponent'` 1138, 1738, 2761, `'D_log10'` 102.15, 120.51, 155.36 and `'finite_D_error_over_tau': 0.25`. Writes `budgets_n*.json`.
- [ ] `python3 -B leading_checks.py`. Tested 2 s. Expected `quadratic_residual` below 1e-19 and `relative_bus_phase_error` near 1e-24 for each case. Writes `leading_checks.json`.
- [ ] `python3 -B ed.py 2 3 1 --seconds 180` and `python3 -B ed.py 3 3 1 --seconds 180` (arguments n, cells, Q). Tested 85 s and 78 s. Expected `error_over_lam4` 0.033706 and 0.023737 at lambda = 0.1, matching the table in C42 RESULT. Other (n, cells, Q) combinations, n in {1, 2, 3}, cells in {2, 3}, Q in {0, 1}, regenerate the rest of that table (quoted up to about 2 min each). As the C39 referee notes, these cluster runs resolve the scalar bond and readout, not the small hopping, which A3 checks directly.

### A5. The mobile wrapper on the 46-letter core (C44)

Folder `codex/work_C44/`.

- [ ] `python3 -B bounds_checks.py`. Tested 1 s. Expected one record per n = 1 to 4 with positive `illegal_selector_margin`, `wrong_occupation_margin` and `min_secant_ratio` (0.0047 at n = 4). Writes `bounds_results.json`.
- [ ] `python3 -B consistency_checks.py`. Tested 1 s. Expected `smallest_local_eigenvalue` near -2e-15 (positivity), `termwise_code_residual` and `bare_source_comparison_error` below 1e-16, `readout_compression_error` 0. Writes `consistency_results.json`.
- [ ] `python3 -B transport_ring.py`. Tested 3 s. Expected projected marker-graph spectra of dimensions 104, 432, 1024 with `spectrum_error` at most 2.3e-14. Writes `transport_ring_results.json`.
- [ ] `python3 -B audit_numeric.py cell`. Tested 1 s. Expected n = 1 legal cells of dimension 1728 and 4320, each with a two-fold zero eigenvalue (the code) and next levels 0.0147581 and 3.3631e-5, eigenpair residuals below 2e-13.
- [ ] `python3 -B audit_numeric.py`. Tested 36 s. Source-bond calculation at n = 1, 2 for four wrapper parameter choices. Writes `source_bonds.npz` and `source_results.json`.

Folder `contrib/referee_C44_C45/` (rebuilt from the prose, imports nothing from C44 or C45).

- [ ] `python3 -B ref_charge.py`. Tested 5 s. Expected `C44_cellED_n1_p1` 0.0147580634663 and `_p2` 3.36309542068e-5 agreeing with C44 to 3e-10 relative, `"all_bound5": true`, and asymptotic coefficients 2.000 (p = 1) and 0.5, 0.25, 0.125, 0.0625 (p = 2).
- [ ] `python3 -B ref_filter_selector.py`. Tested 6 s. Expected `selector_b_range` [4, 4470], `min_mu_c` 0.599, `min_eta_b4` 7.66 and `min_secant_slack` near -2e-10 (equality at the defining point, roundoff).
- [ ] `python3 -B ref_equalize_gauss.py`. Tested 3 s. Expected `equalize_max_ratio` 0.537 over 12000 trials (C44 eq. 13 holds), `gauss_min_ratio` 0.77496 against `sqrt_3_5` 0.77460, and `min_dimension_of_D3_vertex` 2.325.

Not yet independently reproduced. C44 eqs. 9 and 11 (transport overlap and source formula) are checked only by C44's own scripts. A rebuild is in progress in `contrib/verify_C44/`, which has no report yet.

### A6. Local Schrieffer–Wolff, threshold and schedule (C45)

Folder `codex/work_C45/`.

- [ ] `python3 -B checks.py`. Tested 0.5 s. Expected `second_order_identity_error` below 7e-17, `remainder_over_lambda3` approaching 0.343 (two cells) and 0.621 (three cells) as lambda halves, `defect_inequality_min_eigenvalue` above -4e-15, and `"exact_budget_checks": "passed"`. Writes `checks.json`.

Folder `contrib/referee_C44_C45/`.

- [ ] `python3 -B ref_localsw.py`. Tested 167 s (quoted 70 s). Expected, for a generic and a port-like V on N = 2 to 5 cells, identity error `id10` below 7e-16, `R/l3` converging as lambda decreases and growing linearly in N (about 2.4, 4.9, 7.6, 10.4 for the generic V), the actual remainder between 2.5e-8 and 3.4e-7 of the C_SW bound (`bnd`), and threshold eigenvalues `thr` above -7e-16. Writes `ref_localsw.json`.
- [ ] `python3 -B ref_schedule.py`. Tested 0.3 s. Reads `../../codex/work_C42/budgets_n*.json`. Expected a table of log10 scales, at n = 1 `D=114.1 g=-384.3 J0=-1472 lam=-3.826e+04 W=-6.547e+04 crossover_cells_length_mode=35.4`, with W -2.305e+05 and -6.392e+05 at n = 2, 3, and `'all_positive': True` with every schedule inequality satisfied. Writes `ref_schedule.json`.

### A7. Partial results toward S (C46)

Folder `codex/work_C46/`.

- [ ] `python3 -B checks.py`. Tested 1.3 s. Expected `"status": "numerical, independent toy controls only"`, nine partition-function comparisons with `log_ratio` below `log_upper`, the raw-coefficient growth of the one-cell example, and the chiral-speed thermal coefficient converging to 1.11112 (error 1.7e-5 at the smallest T). Writes `checks.json`.

## B. Claim B, the d = 24 chain

Most of this folder is the work of a collaborator (GitHub: salpace)'s campaign with GPT-6-family Codex workers.

### B1. The bond

- [ ] `contrib/counterexample_campaign/doping/`: `python3 -B bond24.py`. Tested 1 s. Expected `"bond_shape": [576, 576]`, `"nonzero_entries": 923`, zero Hermiticity, occupied-block, particle-number, charge and marker errors, and `"alpha_half": 0.29735763271532445` (that is 1/2 - 2/pi^2). Writes `bond24_fixture.npz` and `bond24_checks.json`. The fixture uses moderate parameters, so it checks structure, not the n = 23 point.
- [ ] `contrib/referee_d24/`: `python3 -B rebuild_bond.py`. Tested 0.5 s. Independent rebuild from the prose. Expected `"fixture_max_abs_difference": 0.0`, 923 nonzero entries in both, and every commutator 0.0.

### B2. Squeezing, gluing and the charge threshold

Folder `contrib/counterexample_campaign/doping_check/`.

- [ ] `python3 -B check.py`. Tested 0.3 s. Random complex internal bonds, the literal ring Hamiltonian against its squeezed form. Expected six lines (L, q, M, dimension, error) with dimensions 24 to 270 and errors at most 2.3e-15.
- [ ] `python3 -B color_gluing.py`. Tested 0.2 s. Expected errors at most 4.5e-15 on prime and composite rings.
- [ ] `python3 -B charge_threshold.py`. Tested 1.6 s. Prints L, theta/pi, the fixed-momentum threshold and the proved lower bound. Expected for L = 16 to 256 the third column above the fourth in every row.

### B3. Cell-length selection and illegal words at fixed nonzero lambda

- [ ] `contrib/counterexample_campaign/route_a/`: `python3 -B check_global_sectors.py`. Tested 2 s. Rational n = 23 selector checks and finite diagnostics. Expected `excited_min` above `proved_excited_lower` 0.9 at each n and `port_norm` 1.7320508 (sqrt 3). The all-length statement itself is proved in that folder's README.

### B4. The n = 23 Fermi-point certificate (C33), rerun

- [ ] `contrib/counterexample_campaign/route_a/`: `python3 -B audit_c33.py`. Tested 88 s. Copies `codex/work_C33/` to a temporary folder, regenerates the 512-bit loop matrices and reruns the integer interval certificate. Expected `"generated_arrays_equal"` all true, `"root_count": 46`, `"delta2_lower": 0.01401938...`, `"delta4_lower": 3.15587...e-07` and an inertia sequence rising 0 to 23 and falling back.
- [ ] `contrib/referee_d24/`: `python3 -B rerun_c33.py`. Tested 99 s. The referee's own rerun of the same code. Same expected numbers. Writes `c33_rerun.json`.

### B5. Harmonic selection and the vacancy margins

- [ ] `contrib/counterexample_campaign/infrared/`: `python3 -B checks.py`. Tested 0.6 s. Expected `"physical_mismatch_lower": 6.5747e-09`, `"free_charge_restoration_threshold_over_t": 8.3712e-09`, `"minimum_uniform_nonGaussian_dimension": 3`, `"minimum_nonzero_h_charge_dimension": 4`.

### B6. Margins of the 23-letter chain for n up to 23

Folder `contrib/margins_audit/`.

- [ ] `python3 -B literal_check.py 4`. Tested 80 s. Structured second-order coefficient against a literal two-cell computation, n = 1 to 4. Expected `T_rel_err` near 1e-15 and `quadratic_residual_over_T` at most 4e-11.
- [ ] `python3 -B precision_check.py 23`. Tested 742 s (quoted 12 min). Reruns n = 23 at 313 and 373 digits, both at the audit's schedule and at C33's exact point. Expected two lines with `"max_root_diff": 0.0`, `"max_centre_diff": 0.0`, delta4 `1.81935444773e-7` (audit schedule) and `3.99638534799e-7` (C33's point, above its certified 3.1558e-7), and `log10_lmin_F` -63.6148. Appends to `precision_check.json`.
- [ ] `python3 -B run_schedule.py`. Quoted 1 h on 9 cores, not rerun. Regenerates `results/n*.json` and the margins table of the README (delta4 at least 4.2e-8, delta2 at least 0.014 for every n from 2 to 23).

### B7. Modular non-identity certificates (C37 method)

Folder `contrib/routeA_final/`.

- [ ] `python3 -B mobile_modular.py check`. Tested 1 s. Expected `"matches_C37": true` in every regression case against C37's saved frozen-core polynomials and `"collision detected"` in every no-loop and no-field control at n = 2, 3, 4.
- [ ] `python3 -B mobile_modular.py 23`. Tested 15 s. Expected `"gcd": 1` and `"status": "no_pair_product_collision"`. `python3 -B mobile_modular.py frozen 23 24 ... 40` extends the frozen-core certificate (quoted a few minutes at the top end).
- [ ] `python3 -B mobile_leading.py`. Tested 2 s. Structured order-g^2 symbol against literal extraction, n = 1 to 3. Expected `T_rel_err` near 1e-16 and `M_err_over_T` at most 3e-13.

## C. The mechanism at small d

### C1. Ten-letter frozen demonstrator (c = 1)

Folder `contrib/routeA_demonstrator/`.

- [ ] `python3 -B ed_check.py`. Tested 18 s. Expected a one-cell code with `"degeneracy": 2`, `"gap"` 0.4, `"max_port_compression"` 1e-32, second-order `Jxy` -0.28125 against the closed form, the kappa = 0 control with `Jxy` zero, residual/lambda^3 roughly constant (2.55 to 3.23), and a three-cell connected term of order 0.005 lambda^3.
- [ ] The central-charge numbers (finite-entanglement scaling c = 1.00 to 1.02 at three couplings, ring spectra M = 3 to 6 cells) come from `run_fes.py` and `ring_ed.py`, which ran on the IAS Typhon cluster (`fes_array.slurm`, `ring_array.slurm`, task lists `tasks_*.txt`). Quoted, not rerun. The M = 6 rings reach 1.55e8 states per sector.

### C2. Vacancy lift of the demonstrator (d = 11) and of the 23-letter bond at n = 1 (d = 24)

Folder `contrib/vacancy_mobile_numerics/`.

- [ ] `python3 -B brute.py`. Tested 107 s (quoted 2 min). Literal 121 x 121 bond on the physical basis against the squeezed solver. Expected agreement (`max_diff`) at most 2.1e-13 in five sectors of dimension 600 to 17640. Writes `data/brute_check.json`.
- [ ] `python3 -B d24.py check`. Tested 64 s. Literal 576 x 576 lift against the squeezed solver. Expected `max_diff` at most 2.3e-13 in sectors of dimension 780 to 3360.
- [ ] The ground-state scans, the charge-boson Casimir c = 0.998 and the FES runs (`scan.py`, `separated.py`, `modes.py`, `run_fes.py`) ran partly on Typhon (arrays listed in that README). Quoted, not rerun.

## D. What no script checks

- Assumption S of the note (Section 5.3), that is, elimination of the massive sector at fixed coupling followed by weak-coupling RG with constants independent of the gap and the wrapper ratios.
- For the d = 24 chain, inputs P (S for the frozen 23-letter chain) and M (vacancy matching), the positive dressed defect cost, and every-ring uniqueness.
- A calibrated central charge for any chain with c > d - 1, or for any small member of the 46-letter family. None exists.

These are the items for human review, together with a reading of C39 to C45 (`codex/work_C39/PROOFS.md`, `codex/work_C42/PROOF.md`, `codex/work_C44/PROOFS.md`, `codex/work_C45/PROOF.md`).
