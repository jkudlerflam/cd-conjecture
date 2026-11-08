This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C158 result: infinite-size central charge near the XZ point

Derivations and the reading of the numbers are in PROOF.md, the log in PROGRESS.md. Scripts in this folder: vumps.py (C153, unchanged), fesx.py (models, MPO for an n-site cell, damped-VUMPS driver, Hadamard x reflection symmetry diagnostic), fesfit.py (estimators), chain.sh (local chi-by-chi restarts), check_alg.py, dmrg.py, dmrg_run.py, xzlib.py (C144, unchanged). Local data in local/ (chi <= 56, stopped at the coordinator's request), Typhon inputs in typhon/.

## Bottom line

No point near XZ shows c clearly above 1. Where the bond dimension gets past the crossover, every estimate is 1 within about 0.07. The naive finite-entanglement readings of 2 to 6 at small field are an artifact, explained below. hp <= 0.2 and the tuned hyperplane are not settled locally and wait on Typhon job 17595611.

## Method validation (damped VUMPS, S = (c/6) log xi; estimators: least squares over the last 3 / 4 chi, 1/delta of Rams-Czarnik-Cincio, Pollmann-Moore kappa) [numerical]

| point | exact c | c (xi, last 3 / 4) | c (1/delta) | c (kappa) | max chi |
|---|---|---|---|---|---|
| critical Ising | 1/2 | 0.498 / 0.503 | 0.487 / 0.497 | 0.56 | 24 |
| XX | 1 | 1.009 / 1.018 | 0.93 / 0.97 | 0.98 | 32 |
| Heisenberg | 1 | 0.995 / 1.019 | 1.02 / 0.76 | 0.95 | 32 |
| XZ + ZX (two Majoranas) | 1 | 1.009 / 1.018 | 0.93 / 0.97 | 0.98 | 32 |
| XZ regulated (ZX 0.2 on even bonds) | 1/2 exactly | 0.487 / 0.486 | 0.49 | 0.48 | 24 |
| XZ + 0.3 ZX (velocities 1 : 0.3) | 1 | 1.058 / 1.022 | 0.69 / 0.80 | 1.13 | 40 |
| XZ + 0.1 ZX (velocities 1 : 0.1) | 1 | 0.90 / 0.99 | fails | 0.78 | 56 |

- [proved] At XZ, a_j = X_j Z_{j+1} and b_j = Z_j X_{j+1} generate two commuting Majorana-bilinear algebras. So XZ + eps ZX is two free Majoranas with velocity ratio eps (c = 1), and XZ + eps ZX on even bonds only is A critical with B in a stabilizer state (c = 1/2 exactly). The degeneracy at XZ is handled by this regulator. A (-1)^j ZX regulator is a gauge copy of XZ + eps ZX, not a regulator.
- The two-velocity cases calibrate the method for the problem at hand. Single pairwise slopes scatter by about +-0.25 (0.83 to 1.41 at eps = 0.1) because the slow sector's share of chi grows in integer steps, while least-squares fits over 4 values land within 0.01 to 0.1 of 1. The 1/delta estimator fails when two sectors' levels interleave.

## (a) Bisector-field line h = XZ + B (x + z)/sqrt2, |B| = 2 hp (hp in K units of arXiv:2209.06221)

| hp | c | basis |
|---|---|---|
| 0.45 | 1.01 +- 0.03 | c_pair 0.95, 1.01, 1.02, 1.00 (chi 16 to 48); least squares 1.008; kappa 1.333 -> 1.014; 1/delta 0.90 to 1.06 |
| 0.3 | 1.07 +- 0.05 (pair scatter), plus the +-0.25 single-pair systematic | c_pair 1.065, 1.026, 1.131 (chi 24 to 48) after the crossover at chi = 24; kappa -> 1.15 |
| 0.2 | undetermined locally | crossover not passed at chi = 40 |
| 0.4 | not run locally | Typhon |
| 0.6 (control) | gapped, symmetry broken | S = 0.11, xi = 1.8, sym_HR = 0.616 |

- [numerical] Structure. A gapped k = pi sector (xi_pi of about 9 sites at hp = 0.3, about 20 at hp = 0.2, saturating with chi) plus one critical incommensurate sector, with wavevector 0.596 pi (hp = 0.3) and 0.332 pi (hp = 0.45) per two-site cell. The Hadamard x reflection symmetry is unbroken on the critical line (sym_HR = 1.000000 at hp = 0.45) and broken at hp = 0.6.
- [numerical] Scaling dimensions from transfer-matrix ratios at chi = 48, hp = 0.45: x_k/x_1 = 1 (+-Q), 2.05 (k = 0), 4.28 (+-Q), 5.19 (0), 5.54 (+-2Q). At hp = 0.3: 1 (+-Q), 2.09 (0), then the gapped k = pi pair at 2.36, still rising. [heuristic] The 1 : 2 : 4 pattern fits one compact boson near the free-Dirac radius, i.e. A and B Majoranas combined into one Dirac fermion (C156's guess). It does not fit an Ising + boson c = 3/2 theory, which would need an extra low level with x = 1/8.
- [numerical, the adversarial point] For chi <= 24 (hp = 0.3) and chi <= 40 (hp = 0.2) the largest xi belongs to the saturating k = pi mode while S keeps growing, so 6 dS/dlog xi_max reads 2 to 6 (single pairs up to 14, with sign changes), and 1/delta fails. This is a crossover artifact, not c > 1. Only bond dimensions where the growing incommensurate mode leads, with xi_pi saturated, are used above. Open chains at L <= 64 sit inside the same crossover (C156: S(L/2) flat between L = 32 and 48 at hp = 0.3).
- [derived] X_j and Z_j are products of one Ising string operator of each Majorana sector, so the bisector field couples A and B through dimension-1/8 operators of the critical A chain. This explains why C144's second-order B-hopping picture does not apply here.

## (b) YX and ZY at XZ (six-site cell)

- [numerical] YX eps = 0.3: c = 1.02 +- 0.05 (pairs 0.97 to 1.08, kappa -> 1.08). YX eps = 0.2: c = 1.00 +- 0.05 (pairs 0.95 to 1.07, least squares 1.005 / 0.981, kappa -> 0.984). S is identical on all six bonds, the leading mode is at k = 0, and the spectrum ratios are 1, 1, 1.10, 2.5 (SU(2)_1-like, split by the marginal coupling). This agrees with C154's mapping to a period-3 XYZ chain with SU(2)_1. ZY is the mirror image of YX (C154). ZY 0.2 and YX 0.05, 0.1 (C154's dimerization window) are in the Typhon job.

## (c) Ashkin-Teller plane and tuned hyperplane

- [numerical] At beta = 0.3 and alpha = 0.17 to 0.23 the plane is gapped at every sampled point with chi <= 48. The two gapped phases, with and without a log 2 bond-entropy alternation, meet near alpha = 0.197 to 0.2, where C144 put its c = 1/2 line. The line itself was not resolved, and metastable branches appear (energies not concave in alpha). No sign of c near or above 1.
- The tuned hyperplane T (alpha = (2/pi) beta, beta = 0.3) with YX 0.1, ZY 0.1, XX + ZZ 0.1, B_y 0.1 and bisector field hp = 0.15 admixed: Typhon only, not yet run.

## Task 3 (characterisation of c > 1 points)

Not triggered: no point has c clearly above 1. The spectra recorded for the c = 1 sectors are above.

## Pending Typhon job

- 17595611 (c158_fes): partition typhon, one 48-core node, 23:55 limit, 43 single-thread tasks. Status PENDING (Priority); squeue estimated start 2026-11-12 13:59. Directory /scratch/lustre/jkudlerflam/c158/.
- Bisector: hp = 0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5 (two-site cell; hp = 0.2, 0.3, 0.4 to chi = 384, the others to 256); hp = 0.2, 0.3, 0.4, 0.45 with a four-site cell; hp = 0.3 seed 1 and one-site cell; hp = 0.6 control.
- Also: validation (tfi, xx, heis, xzzx to chi = 128; xzdim 0.2, 0.05; xzeps 0.3, 0.1, 0.03 to chi = 256), YX 0.05 to 0.5 and ZY 0.2 (six-site cell, chi to 160), AT alpha = 0.09, 0.18, 0.198, 0.3, the five tuned-hyperplane points, and open-chain DMRG of the bisector chain at hp = 0.2, 0.4 (L = 256) and 0.1 (L = 192), chi to 160.
- Outputs: fes/<tag>.json (written after every chi, with sym_HR, 12 transfer eigenvalues, Schmidt values and magnetizations), state/<tag>.npz, dmrg/bis_h*_L*.jsonl, logs/<tag>.log, logs/status.txt, slurm_17595611.out. Analysis: python3 fesfit.py fes/*.json. Read the bisector files from the bond dimension at which xi_ic exceeds xi_pi (fesfit prints both).
- Related jobs, not mine: C157 17595396 (bisector hp = 0.3, 0.45 at L = 192, 256; pending), C156 17594698 (bisector L = 48 to 128; running), C144 17580192 and 17580432 (YX/ZY L = 128 to 256; running).
