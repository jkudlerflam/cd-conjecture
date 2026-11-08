# C158 progress log

This run was done by a Claude Opus agent (2026-11-07), not by Codex.

## Setup

- Read specs/C158_xz_idmrg.md, specs/_qubit_c_context.md, specs/C157_xz_counterexample.md, codex/work_C144/RESULT.md, codex/work_C156/RESULT.md, codex/work_C153/PROGRESS.md and vumps.py/fes.py, codex/work_C157/PROGRESS.md, and the YX/ZY part of codex/work_C154/PROGRESS.md (XZ + b YX + g ZY is a period-3 twisted XYZ chain; L = 24 ED gives SU(2)_1 for eps >= 0.2 and a dimerization tendency at eps = 0.05).
- Copied codex/work_C153/vumps.py (unchanged) and codex/work_C144/dmrg.py (unchanged) and cfit.py (as cfit_c144.py).
- New scripts: fesx.py (model builder with site-dependent bonds in an n-site cell, MPO of bond dimension 5, driver around vumps.run in damped mode, records S per bond, 12 transfer eigenvalues, Schmidt values, magnetizations, and optionally the MPS in an npz for restarts), fesfit.py (c from 6 dS/dlog xi pairwise and by least squares over the last k bond dimensions, the same with xi -> 1/delta, delta = eps_2 - eps_1 (Rams et al.), and from the Pollmann-Moore exponent).
- Conventions. Pauli bond h = sum J_ab s^a s^b + B.s (one field per site). Bisector line: J = e_x e_z^T, B = -2 hp (1, 0, 1)/sqrt2 with hp the field of arXiv:2209.06221 in units of K (C156 equivalence; the sign of B is a pi rotation about y, which fixes XZ). h_c1 = 0.511 corresponds to |B| = 1.022.
- Machine: local load average near 290 on 10 cores (other agents), so local runs are kept to chi <= 32.

## Exact facts used for validation [proved]

- At XZ the operators a_j = X_j Z_{j+1} and b_j = Z_j X_{j+1} commute with each other for all j, k; each family anticommutes only with its nearest neighbours. So XZ + eps ZX is two decoupled Majorana chains (A with hopping 1, B with hopping eps), c = 1 for eps != 0 with velocity ratio eps, e0 = -(2/pi)(1 + eps).
- XZ + eps (-1)^j ZX is the same as XZ + eps ZX: a sign alternation of Majorana hoppings is removed by gamma_j -> s_j gamma_j with s_j s_{j+1} = sign(t_j). My first "regulator" was therefore not one (c = 1); kept as model xzalt.
- XZ + eps ZX on even bonds only: the b_{2j} commute, so B sits in a stabilizer state with gap 2 eps and A is the critical Majorana chain. c = 1/2 exactly, e0 = -2/pi - eps/2. This is the regulated XZ point (model xzdim).

## Local validation (chi <= 32, damped VUMPS, two-site cell)

- tfi: e converges to -4/pi (de = 3e-8 at chi = 24); c_pair = 0.509, 0.511, 0.490; least squares over the last 3 and 4 chi: 0.498, 0.503 (xi) and 0.487, 0.497 (1/delta). Pollmann-Moore inversion gives 0.56 (kappa less reliable at small chi). [numerical]

## Typhon job 17595611 (c158_fes), submitted 2026-11-07 about 21:45

- One packed node (partition typhon, 48 cores, 23:55 limit), 43 single-thread tasks run by xargs -P 48 (typhon/tasks.txt, typhon/run.sh, typhon/c158.slurm). Directory /scratch/lustre/jkudlerflam/c158/. Outputs fes/<tag>.json (written after every chi, so partial results are usable), state/<tag>.npz (last MPS), dmrg/bis_h*_L*.jsonl, logs/<tag>.log, logs/status.txt, slurm_17595611.out. Analysis: python3 fesfit.py fes/*.json.
- Contents. Validation: tfi, xx, heis, xzzx (chi to 128), xzdim eps = 0.2, 0.05, xzalt 0.1 (chi to 192), xzeps eps = 0.3, 0.1, 0.03 (chi to 256). Bisector line hp = 0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5 (two-site cell, chi to 256), hp = 0.2, 0.3, 0.45 with a four-site cell (chi to 192), hp = 0.3, 0.45 with plain VUMPS, hp = 0.3 with a one-site cell, hp = 0.6 (soliton phase control). YX eps = 0.05, 0.1, 0.2, 0.3, 0.5 and ZY 0.2 with a six-site cell (chi to 160), YX 0.2 with a three-site cell. Ashkin-Teller plane beta = 0.3, alpha = 0.09, 0.18, 0.198, 0.3. Tuned hyperplane T (beta = 0.3, alpha = (2/pi) beta = 0.191) plus YX 0.1, ZY 0.1, (XX + ZZ) 0.1, B_y 0.1, bisector field hp = 0.15. Open-chain DMRG (C144 dmrg_run.py, chi to 160) of the bisector chain at hp = 0.2 and 0.4 (L = 256) and hp = 0.1 (L = 192), chosen not to duplicate C156 (L <= 128) or C157 job 17595396 (hp = 0.3, 0.45 at L = 192, 256).
- Status at submission: PENDING (Priority); the cluster is full.

## Local results, first round (chi <= 40, files now in local/old/)

- check_alg.py (check_alg.out) [proved, checked on an L = 8 ring]: a_j and b_k commute for all j, k; each family anticommutes only with nearest neighbours; X_j anticommutes with exactly a_{j-1} and b_j, Z_j with exactly a_j and b_{j-1}, Y_j with a_{j-1}, a_j, b_{j-1}, b_j. So Y_j is a local Majorana bilinear mixing A and B (hybridization), while X_j and Z_j are products of one string (order or disorder) operator of each Majorana chain. The bisector field therefore couples the sectors through operators of dimension 1/8 in the critical A chain.
- xzdim eps = 0.2 (regulated XZ, c = 1/2 exactly): e within 6e-8 of exact at chi = 24, bond entropies differ by log 2 (cut through a B stabilizer pair), c_pair 0.489, 0.464, 0.504, least squares 0.487 (xi), 0.492 (1/delta), Pollmann-Moore 0.481. [numerical]
- xzeps eps = 0.3: c_pair 0.91 to 1.13, least squares 1.058 (last 3), 1.022 (last 4); kappa inversion 1.13. xzeps eps = 0.1: c_pair oscillates 0.87, 1.41, 0.84, 1.02 (the chi = chi_A chi_B staircase of two decoupled sectors), least squares 0.90 to 0.99. The 1/delta estimator fails for decoupled sectors (two independent towers in the transfer spectrum). [numerical]
- Plain VUMPS (mode vumps) does not converge on the bisector chain (err 0.2 after 600 iterations at chi = 8); the damped mode does. The two plain-VUMPS Typhon tasks were replaced in tasks.txt before the job started (bis0.3_n2s1: seed 1; bis0.4_n4: four-site cell).
- Energy cross-check: the bulk energy from C156's open-chain DMRG, (E_64 - E_48)/16 at hp = 0.45, is -0.189810 K per site = -0.75924 in Pauli units; VUMPS gives -0.759275 at chi = 40. Normalizations agree. [numerical]
- Bisector hp = 0.45: c_pair 1.02, 1.35, 0.95, 1.005, 1.018; least squares 1.010 (last 3), 0.983 (last 4), 1/delta 0.84 to 1.06, kappa 1.336 -> c = 1.01. The leading transfer eigenvalues are a complex pair at arg = +-0.332 pi per two-site cell (incommensurate, wavevector about 0.166 pi per site), then a k = 0 mode with xi ratio 2.02, then another +-0.335 pi pair (ratio 4.3). [numerical]
- Bisector hp = 0.3 and 0.2: crossover. The largest xi is pinned by a k = pi mode (about 8 at hp = 0.3 and 20 at hp = 0.2) while S keeps growing, so 6 dS/dlog xi_max reads 2 to 6 (hp = 0.3: c_ls 1.8 to 2.3; hp = 0.2: 4 to 6, with negative pairs). This is the artifact the spec warns about. At hp = 0.3 an incommensurate pair (arg 0.6 pi per cell) overtakes the k = pi mode at chi = 24, after which c_pair = 1.07; at hp = 0.2 the incommensurate pair (arg 0.75 pi) grows (xi_ic = 3.4 to 7.9) but is still below the k = pi mode at chi = 40. Using xi_ic instead gives 1.18 (hp = 0.3) and 1.42 (hp = 0.2), both still contaminated by the growing k = pi sector. Not converged. [numerical]
- YX eps = 0.3, six-site cell: S uniform over the cell, k = 0 leading mode, c_pair 0.98, 0.98, 1.08, 0.97; least squares 1.024 (3), 1.015 (4); kappa 1.284 -> 1.08. [numerical] Consistent with C154's SU(2)_1.
- Typhon job 17595611: squeue --start estimates 2026-11-12 13:59 (Priority). So the decisive large-chi data will come late. Local chaining started: chain.sh runs fesx.py one chi at a time (each call capped at 290 s, single thread), restarting from local/<tag>.npz. Tags bis0.45, bis0.3, bis0.2, xzeps0.1 (chi to 96), yx0.2_n6 (to 56), at_a0.198 (AT line beta = 0.3, alpha = 0.198, C144's c = 1/2 line).

## Coordinator message (mid-run) and response

- Requests: bisector line first (hp = 0.2, 0.3, 0.4, control 0.6); Rams-Czarnik-Cincio scaling against log(1/delta) as well as xi; push xi well beyond 1e3; report transfer-spectrum ratios; use the exact Hadamard x reflection symmetry; entanglement only (not Casimir) for c; production on Typhon only.
- Done: local chains stopped (at chi = 40 to 56). The 1/delta estimator was already in fesfit.py (c_d). fesx.py now records sym_HR, the spectral radius per cell of the mixed transfer matrix between |psi> and (prod H) R |psi> (H Hadamard, R reflection, with and without a one-site shift); 1 means the state keeps the symmetry. Test: hp = 0.45 gives 1.000000 (unshifted), hp = 0.6 gives 0.616 (symmetry spontaneously broken, the soliton phase), xzdim gives 0.988 (its regulator breaks the symmetry). Typhon tasks.txt updated before the job started: bisector hp = 0.2, 0.3, 0.4 (two-site cell) now run to chi = 384, bisector tasks listed first; fesx.py and fesfit.py re-uploaded. Job 17595611 still PENDING.

## Local results, second round (chained runs, chi to 48 or 56; local/*.json)

- bis hp = 0.45: c_pair (chi = 24 to 48) 0.951, 1.005, 1.018, 0.997; least squares 1.009 (3), 1.008 (4); 1/delta 1.057, 0.901; xi_ic 1.009; kappa 1.333 -> 1.014. Transfer spectrum at chi = 48 (eps_k/eps_1 @ arg/pi per cell): 1.000 @ +-0.332, 2.047 @ 0, 4.278 @ +-0.334, 5.189 @ 0, 5.538 @ +-0.666. [numerical]
- bis hp = 0.3: incommensurate pair leads from chi = 24; c_pair 1.065, 1.026, 1.131 (chi 24 -> 48); least squares 1.066; kappa 1.236 -> 1.145; 1/delta fails (the k = pi levels cross the incommensurate ones). xi_pi saturates (8.4, 8.7, 8.9, 9.0) while xi_ic grows (9.6 to 21.3): a gapped k = pi sector plus a critical incommensurate sector. Spectrum at chi = 48: 1.000 @ +-0.596, 2.086 @ 0, 2.36 @ pi (two levels, gapped sector), 3.44 @ +-0.400. [numerical]
- bis hp = 0.2: unchanged from round one (chi = 40), crossover not passed.
- yx 0.2 (six-site cell): c_pair 1.058, 0.946, 1.069, 0.950, 0.979, 1.043; least squares 1.005, 0.981; kappa 1.357 -> 0.984; spectrum 1, 1, 1.10, 2.5 (an SU(2)_1-like near-degenerate quartet split by the marginal coupling). [numerical]
- xzeps 0.1 to chi = 56: c_pair 0.87, 1.41, 0.84, 1.02, 1.26, 0.97, 0.83; least squares 0.90 (3), 0.99 (4); kappa 1.568 -> 0.78. This calibrates the scatter of pairwise estimates in a two-velocity system at these chi: about +-0.25. [numerical]
- Validation at chi <= 32 (local/xx, xzzx, heis): XX 1.009 / 1.018 (least squares 3 / 4), 1/delta 0.93 / 0.97, kappa 0.98; XZ + ZX identical to XX to all printed digits; Heisenberg 0.995 / 1.019, kappa 0.95 (first pair 0.23, dimerization tendency of the two-site cell at chi = 8). [numerical]
- AT plane beta = 0.3, alpha = 0.17 to 0.23 (local/at/, chi <= 48): gapped everywhere sampled (S saturated, xi <= 32 at alpha = 0.196). Energies are not concave in alpha (e(0.198) lies below the chord of e(0.196) and e(0.2)), so some runs are trapped in metastable branches; the c = 1/2 line of C144 was not located here. [numerical]

## Close (2026-11-07)

- RESULT.md and PROOF.md written. No local processes left running. Typhon job 17595611 PENDING (estimated start 2026-11-12 13:59); outputs will appear in /scratch/lustre/jkudlerflam/c158/ as listed in RESULT.md.
