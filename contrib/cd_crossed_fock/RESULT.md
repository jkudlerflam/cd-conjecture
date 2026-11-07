# cd_crossed_fock: the free-boson majorization of the crossed spectrum

2026-11-06. Opus 5.5 agent, route of `specs/C127_crossed_majorization.md`, with the coordinator's scope update giving priority to d = 2. Derivations are in [PROOF.md](PROOF.md); every numerical statement below is reproduced (in reduced form, about 90 s) by `run_checks.py`, which writes `checks.json`. Tags: proved, numerical, heuristic, assumed.

## Bottom line

[Proved] For XX the majorization holds with explicit constants: for every even L >= 2 and every beta, log Z~(L) <= K + 2F(L/(V' beta)) with (V'/M, K) = (4.16, 2.04) or (5.54, 1.44), from the exact crossed product. No new bound is proved for an interacting chain, in d = 2 or otherwise.

[Proved] The route reduces to one statement (H2), a trace-level domination of Z~ by the bosonic Fock space built on single-insertion levels. With the Lieb-Robinson one-particle cone of C125 and an infrared input (H3), (H2) gives A <= (d^2 - 1)[12.95 + 1.91/gamma_0] M (Theorem 5.1; 12.95 = 4e + 2.08). (H2) is at least as strong as the target.

[Proved, with one numerical inequality] The level-by-level (min-max) form of (H2) is false on the antiferromagnetic Heisenberg chain. By an SU(2) selection rule, single-site insertions reach only spin-1 crossed states. The spin-0 composite of the crossed spectrum (the dimerization primary) then lies below the lowest two-particle Fock state: eps_s/eps_t = 1.89, 1.80, 1.75 at beta = 4, 6, 8, with ratio 1 in the infrared. So no operator inequality H~ - eps_0 >= Phi^dagger H_Fock Phi with the exact single-insertion levels can exist. A proof must count composites at the trace level or degrade the one-particle energies by a factor lambda <= 1/2.

[Numerical] On eight d = 2 positive-transfer chains (XX, XXZ at Delta = 0.5, Heisenberg, XXZ at Delta = 2, critical transverse-field Ising, XY on an Ising line, XYZ, XXZ with a staggered field), (H2) holds with K_0 = 0, and the degraded level-by-level form holds globally with lambda = 0.871 or larger (sector by sector with lambda = 0.84 or larger wherever the enumeration is complete). The majorant of the spec holds with S = 2 and V'/M between 1.2 and 6.9 at K = log 2, well below the Lieb-Robinson value 4e = 10.9. The binding regime is the smallest ring, L = 2, at beta of 3 to 8. The CFT regime and the crossover are not binding.

## 1. Setting and reductions [Proved]

- log Z~(L) = S_min of the crossed Gibbs state, log Z~ <= R_L <= 2 log Z~, and a bound log Z~ <= K + kappa M beta/L gives the Casimir bound L(L e_0 - E_L) <= kappa M and A <= (6/pi) kappa M. A majorant (S, V', K) gives kappa M = S pi V'/12 and A <= S V'/2. S = 2 is one c = 1 boson (normalization check). See PROOF 1.
- Any majorant needs S pi V'/12 >= sup_L D(L), D(L) = L(L e_0 - E_L), and K >= log(top crossed multiplicity), and S V'/2 > A strictly (the zero modes cost log(beta/L) otherwise).
- The checkerboard transfer is positive for every bond with -h' = sum_a J_a A_a (x) conj(A_a) + B (x) 1 + 1 (x) conj(B), J_a >= 0. For d = 2 this covers the antiferromagnetic XYZ chains with fields that are staggered in Sx and Sz in the original frame, and the transverse-field Ising chain. A uniform field in the XXZ frame is not covered. See PROOF 2.

## 2. Task 1: smallest majorant constants, and the m dependence

[Numerical] Continuous imaginary time (m = infinity). Sources are exact ring spectra for L <= 14 and beta <= 8, with eps_0(beta) from Lanczos at m <= 11 and Richardson extrapolation (error 8e-7 on XX at beta = 8 against the exact value). The L = 2 ring at beta up to 1e4 uses the low-temperature form, and XX uses its exact product. With S = 2:

| model | M | A/M | D(2)/M | V'_min/M at K = log 2, 1, 2 | V'_min/M as K -> infinity | S = 1: V'_min/M at K = log 2 |
| :-- | --: | --: | --: | :-- | --: | --: |
| XX | 0.5000 | 2.000 | 1.454 | 4.09, 3.68, 3.07 | 2.78 | 6.49 |
| XXZ 0.5 | 0.5000 | 2.598 | 2.000 | 5.55, 5.04, 4.20 | 3.82 | 8.90 |
| Heisenberg | 0.5000 | 3.142 | 2.455 | 6.85, 6.18, 5.14 | 4.69 | 10.92 |
| XXZ 2 | 0.7500 | 0 | 2.041 | 5.81, 5.21, 4.33 | 3.90 | 9.18 |
| TFIM | 1.4142 | 0.707 | 0.399 | 1.22, 1.05, 1.00 | 0.76 | 1.85 |
| XY Ising line | 0.3953 | 0.949 | 1.082 | 3.31, 2.83, 2.29 | 2.07 | 5.00 |
| XYZ | 0.3750 | 0 | 1.719 | 5.08, 4.56, 3.64 | 3.28 | 8.03 |
| XXZ 0.5, staggered field | 0.5220 | 0 | 1.578 | 4.77, 4.17, 3.35 | 3.01 | 7.34 |

- At K = log 2 the binding point is always L = 2, at beta = 3 to 8. The K -> infinity column is the zero-temperature necessity 12 D(2)/(S pi M). Both are lattice-scale quantities. The CFT value 2A/S (the velocity a single boson needs) is smaller in every critical case (for Heisenberg, 3.14 against 6.85).
- D(L) is largest at L = 2 and decreases monotonically in L for all eight models (L <= 18), to pi A/6 for the critical ones. So the zero-temperature content of the majorant is controlled by the 2-site ring, D(2) = 4(e_0 - lambda_min(h)) <= 8M.
- With the spec's choice V' = 4eM, S = 2 needs only K <= 0.05 on the small-ring data, plus K >= log g_0 at large L (log 2 for the Neel doublet of XXZ 2). With S = 2(d^2 - 1) = 6 and V' = 4eM, K = 0 suffices on all small-ring data.

[Numerical] Finite Trotter number. Full crossed spectra, resolved by crossed momentum and charge, for m = 2 to 7 (XX to m = 8), with Trotter step beta/m <= 1. The smallest V'/M at S = 2, K = log 2 rises with m only because the accessible beta window (beta <= m) grows. For Heisenberg the values are 4.83, 6.40, 6.56, 6.66, 6.72, 6.75 for m = 2..7, and the m = infinity value is 6.85; for XX they are 2.66 to 4.04 for m = 3..8, against 4.09. At fixed beta, log Z~(L) converges like (beta/m)^2, from below at L = 2 and from above at L >= 4. For Heisenberg at beta = 4, L = 4 the values are 0.829, 0.756, 0.714, 0.688 for m = 4..7 against 0.614 at m = infinity (m^2 times the difference is 3.4 to 3.6). No constant drifts with m beyond this Trotter error. Larger runs (m = 8 for all U(1) chains, m = 9 for XX and Heisenberg) are prepared as packed SLURM jobs (`typhon_job1.slurm`, `typhon_job3.slurm`); the Typhon queue gave an estimated start of 2026-11-10, so they were cancelled unrun. Heisenberg at m = 8 was run locally for the insertion analysis.

## 3. Single insertions [Proved structure, numerical content]

- [Proved] Insertion states are |Y_{O,j}> = K_e^{1/2} Ohat_j K_o K_e^{1/2} |psi_0>, which reproduce the insertion correlators exactly at finite L with no inverse of K_e. The single-insertion span per crossed momentum has dimension at most 2(d^2 - 1) at finite m (two inequivalent insertion points per Trotter step) and at most d^2 - 1 in continuous time. For U(1) chains the span lies in Q = 0, +1, -1.
- [Proved] Selection rule: for an on-site symmetry G of the positive-frame bond, the column operator (x)_j (conj u on even slices, u on odd slices) commutes with K_e and K_o (verified numerically to 4e-16, check I3), and single insertions reach only the isotypic components occurring in the adjoint representation. For Heisenberg that is spin 1 only.
- [Numerical] For the interacting chains, single insertions are nearly sharp one-particle states. In each (Q, n) sector the dim W eigenstates of largest overlap have overlap 0.92 to 0.996 (Heisenberg, m = 7, beta = 2..8; 0.933 and 0.904 at m = 8, beta = 6 and 8, slowly decreasing in m). For the free-fermion chains they are not sharp (XX, TFIM, XY: 0.34 to 0.6), because a local insertion creates a pair continuum.
- [Numerical] Share of log Z~ carried by the single-insertion sector, log(1 + Z_1)/log Z~ with Z_1 = sum_j f_j exp(-L eps_j): 0.49 to 0.99 at L = 2, 0.61 to 1.0 at L = 4 and 0.72 to 1.0 at L = 8 over all models, m = 4..7 and beta <= 8. It decreases with beta at fixed L. In the Cardy regime L << v beta the single-insertion part grows only logarithmically while log Z~ grows like beta/L, so there the multi-particle (Fock) structure carries essentially all of Z~ [heuristic, from the scaling].

## 4. Task 2: are multi-insertion states additive? [Numerical]

No. Two-insertion states, with the vacuum and the single-insertion span removed, at m = 7 and beta = 6:

- Like charges or like signs of frequency give energies well above the sum: Heisenberg Q = 2 from two zero-momentum triplet levels at 0.259 has its dominant level at 0.927, 1.8 times the sum. XX Q = 2 sits at 0.969 against 0.525.
- Opposite charges or frequencies give most of the weight far below the sum. For Heisenberg (Q, n) = (1, 1) with (-1, -1), 67 percent of the weight lies below 0.9 times the sum, dominantly on the spin-0 composite at 0.467 against the sum 1.188. For XX, two n = 1 bosons into n = 2 put 99 percent of the weight at 1.629 against 1.938, from the concavity of asinh.

The Fock count nevertheless dominates Z~ in every case, because each true low state is reached by many insertion labels. Additivity is not what makes the bound true.

## 5. Task 3: proof attempt

- [Proved] XX: Theorem 4.1 of PROOF (exact product, concavity of asinh below a frequency X, log(2 nu) above it, and the Jacobi triple product for the fermion zero modes). This is the XX limit and one-boson normalization check of the route.
- [Proved] Conditional theorem 5.1: (H1) + (H2) + (H3) imply log Z~ <= K + kappa M beta/L with kappa = (d^2 - 1)[2 pi e/3 + 9.28/(pi e) + 1/gamma_0]. (H1) is C125's Lieb-Robinson cone for detected single-insertion levels.
- [Proved obstruction, Heisenberg] The min-max version of (H2) fails (Proposition 6.1). The failing state is a crossed composite of two spin-1 insertions, invisible to every single insertion by symmetry. Its energy ratio to the one-particle level tends to 1 in the infrared (Affleck, Gepner, Schulz and Ziman, assumed), so the degradation factor needed is lambda <= 1/2 asymptotically. The numbers are lambda_crit = 0.945, 0.902, 0.877 at beta = 4, 6, 8 (m = 7) and 0.871 at beta = 8, m = 8, decreasing toward 1/2.
- [Heuristic] Why this is the general obstruction. A many-boson crossed theory realized with few local operators would hide its extra species as composites (currents in the product of a vertex operator with its conjugate) at exactly the cone energies allowed by (H1). The single-insertion count bounds one-particle states, not species. This is C125's lesson 1 inside the lattice framework.
- [Heuristic] Minimal extra input. A crossed binding-energy bound (H2_lambda): the k-th crossed level is at least lambda(d) times the k-th Fock level built on single insertions, with lambda(d) in (0, 1/2]. It implies the target with V replaced by V/lambda. It holds numerically on all eight d = 2 chains with lambda = 0.871 at the accessible temperatures (beta <= 8). Whether fixed-d locality implies it is as open as the target. One lattice fact a proof could use: two insertions at coinciding times on one site collapse to a single insertion of the product operator, so genuine composites come from finite time separations of order 1/M.

## 6. What is not done

[Assumed or open] No proof of (H2), (H2_lambda) or (H3) for any interacting chain, including the d = 2 class. No m = infinity crossed spectra for interacting chains, only m <= 8 spectra, m = infinity partition functions for L <= 14 and beta <= 8, and the zero-temperature limit. The two-slice insertion span at finite m has no exact continuous-time counterpart (in continuous time one insertion time per frequency remains). The Typhon jobs for m = 8, 9 were cancelled before starting (queue estimate four days); m = 12, 13 Lanczos runs for eps_0 were not needed after the m <= 11 extrapolation reached 1e-6.

## Files

`cflib.py` (transfer, sectors, direct traces), `models.py` (d = 2 registry), `majorant.py` (F, XX exact product), `xx_theorem.py`, `xx_exact_majorant.py`, `run_spectra.py`, `run_spectra_model.py`, `analysis_majorant.py`, `compact_majorant.py`, `zeroT.py`, `eps0_lanczos.py`, `infinite_m.py`, `minimal_constants.py`, `insertion.py`, `analysis_insertion.py`, `two_insertion.py`, `make_tables.py`, `run_checks.py`, SLURM scripts `typhon_job*.slurm`. Outputs are in `data/` and logs in `logs/`.
