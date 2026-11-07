# C132 result

[Proved] The periodization construction extends to every CPTP map R with R(rho_(l-1)) = rho_l. Optimizing its bond energy is an SDP, including an exact reduction to the supports of those marginals. Claims 1 to 3 survive, with the L = 2 convention and the vacuum-length quantifier qualified in PROOF.md. The general Lemma P prime remains open, including at d = 2.

[Proved] Two additional conclusions are available at d = 2. For the full ferromagnetic Heisenberg ground-space mixture, Petz is optimal among all admissible channels and

    delta_opt(l) = delta_Petz(l) = 1/[2l(l+1)],
    min_l L^2 delta/M = L/(L-1),  L >= 3.

Thus P and P prime hold on this class with K = 2 including the product convention at L = 2. This extends by positive scaling, gauge changes and on-site rotations, and to the Delta = -1 XXZ mixture on even rings by a staggered rotation. The actual energy inequality on this frustration-free class has K = 0.

[Proved] For all nearest-neighbor quadratic fermion chains, the source's momentum-box construction and ||h_cross|| <= M give

    E_L >= e_0 L - pi^2 M/(6L).

The smallest nontrivial case is d = 2, including XX and TFIM. The argument extends to m on-site modes, d = 2^m, including C106 at d = 4, with both physical Jordan-Wigner parity sectors included. It is a uniform finite-L result for variant (a), not a proof of P prime for those chains. Classical diagonal ground-space mixtures have zero recovery defect at l = 2 for L >= 3, at every d. No new result for all interacting qubit bonds is claimed.

## Audit of the proposed obstruction

[Proved] The source's claim that trace norm cannot prove a dimension-only bound does not follow from its argument. Its continuity inequality leaves log L, so an arbitrarily long crossover can accommodate an arbitrarily large c. Also, the dual-M norm on differences with vanishing gauge pairings equals the trace norm. Consequently the marginal reformulation itself is a two-site trace-norm statement. Small error on the full recovered window is a different, stronger requirement. Energy weighting is a useful strategy, not a necessity proved by the cited obstruction. Details are in PROOF.md sections 1 and 5.

## Reproduced checks

[Numerical] An unmodified local copy of contrib/cd_free_route/run_checks.py passed every quick check. See reproduction.txt. At l = L/2 the reproduced scaled defects Q = L^2 delta/M are

| Chain | L = 6 | L = 8 | L = 10 | L = 12 |
| :--- | ---: | ---: | ---: | ---: |
| XX | 2.533747 | 2.250437 | 2.1666 | 2.1379 |
| Heisenberg | 4.053601 | 3.585604 | 3.4408 | 3.3900 |
| XXZ, Delta = 0.5 | 3.319857 | 2.939222 | 2.8245 | 2.7845 |
| XXZ, Delta = 2 | 3.285582 | 2.746558 | 2.5345 | 2.3403 |

[Numerical] SU(3) gives 9.0000 and 7.4943 at L = 3,6. TFIM at g = 1 gives 0.7710, 0.7592, 0.7532 at L = 8,10,12; g = 1.6 gives 0.1632, 0.0829, 0.0404. C106 at epsilon = 0.03 gives 3.312591 and 2.530338 at L = 4,6. The decoupled-stack identity, gauge cancellation, frozen-marginal identity, XX box smoothing and the source's full-multiplet checks pass.

[Assumed continuation, not executed] The original --long run, including L = 20 and SU(3) L = 12, was not run locally. long_checks.slurm reproduces it on Typhon. The source's larger-size numbers have not been independently reproduced here. No Typhon job was submitted because this sandbox cannot access it.

## Optimal recovery

[Numerical] The following entries use all allowed l at L <= 6 where feasible, with central-block solutions refined using Clarabel. The normalization is Q = L^2 delta/M. Values are numerical optima, not exact arithmetic certificates.

| Chain | L | Petz Q | Optimal Q |
| :--- | ---: | ---: | ---: |
| XX | 4 | 3.313708 | 3.313708 |
| XX | 6 | 2.533747 | 2.533747 |
| Heisenberg | 4 | 5.333333 | 5.333333 |
| Heisenberg | 6 | 4.053601 | 4.053601 |
| XXZ, Delta = 0.5 | 6 | 3.319857 | 3.319599 |
| XXZ, Delta = 2 | 6 | 3.285582 | 3.284746 |
| TFIM, g = 1 | 4 | 1.013635 | 0.658778 |
| TFIM, g = 1 | 6 | 0.804623 | 0.448159 |
| TFIM, g = 1.6 | 4 | 0.720166 | 0.498809 |
| TFIM, g = 1.6 | 6 | 0.306063 | 0.208667 |
| SU(3) swap | 3 | 9.000000 | 9.000000 |
| SU(3) swap | 6 | 7.494320 | 7.4943 |
| C106, epsilon = 0.03 | 4 | 3.312591 | 3.312591 |

[Numerical limitation] At L = 6, C106's central recovery SDP has Choi dimension 1024. Only its l = 2,5 SDPs were solved, both giving 6.109498, which is worse than the central Petz value 2.530338. Thus 6.109498 is not the optimum over l, and no L = 6 optimal C106 value is reported. The central problem is included in sdp_large.slurm.

[Numerical validation] For the refined L = 6 qubit solves, equality residuals and negative Choi eigenvalues are at most a few times 10^-9, and corrected dual gaps are below 10^-8 in energy units. TFIM g = 1.6 nevertheless returned optimal_inaccurate, so its displayed value should be treated accordingly. The SU(3) l = 4 solve has a roughly 10^-6 cone violation; its apparent last-digit improvement is not evidence against Petz optimality. Its l = 3 solve is much cleaner. The data files retain solver status and residuals rather than silently treating all outputs as feasible channels.

[Numerical] At L = 8, limited-iteration central-block solves for XX, Heisenberg and XXZ did not reach adequate primal feasibility. Their raw objectives are excluded as optima. The corrected dual values and the feasible Petz channels instead give the following numerical brackets for l = 4 alone.

| Chain | Dual lower Q | Petz upper Q |
| :--- | ---: | ---: |
| XX | 2.23975 | 2.25044 |
| Heisenberg | 3.56632 | 3.58560 |
| XXZ, Delta = 0.5 | 2.92565 | 2.93922 |
| XXZ, Delta = 2 | 2.72845 | 2.74656 |

[Numerical limitation] These are floating-point dual checks, not rigorous intervals, and they do not bound the minimum over uncomputed l = 5. The TFIM L = 8 central solves also remain insufficiently converged; their objectives 0.41857 and 0.10864 are exploratory only and are not used as feasible upper bounds.

[Heuristic conclusion] Optimization does not generically halve the defect on the sizes tested. XX and Heisenberg show essentially no gain; XXZ gains are small; TFIM gains are substantial and consistent with a different asymptotic coefficient. The large-L statement Q_opt -> pi A/(6M) has not been established or refuted for these fixed critical chains. In particular, finite-L agreement with Petz does not prove equality of asymptotic coefficients. The known targets are pi/3 for XX, pi^2/6 for Heisenberg and pi/(6 sqrt(2)) for critical TFIM. A claim that the optimal channel removes the factor 2 for all critical chains is not supported by these data.

## Structure of the optimizing maps

[Proved] Symmetry averaging preserves feasibility and objective, so an optimizer can be chosen covariant under the common on-site symmetries of h and rho. The ferromagnetic optimizer is explicitly the universal symmetric-spin extension R(X) = l P_sym,l (X tensor 1) P_sym,l/(l+1). Its Choi support is the total-spin-1/2 subspace, which proves its optimality.

[Numerical] At L = 6, l = 3, the XX, Heisenberg and XXZ 0.5 optimizers have two nonzero Choi eigenvalues, both approximately 2. For Heisenberg the Choi matrix agrees with Petz to about 2 times 10^-8 in Frobenius norm. For XXZ 0.5 the Frobenius difference is 0.01418, with a small but resolved energy improvement. The critical TFIM optimizer has two dominant Choi eigenvalues 3.931924 and 0.068075, and Frobenius distance 2.96854 from Petz. The g = 1.6 solution has four resolved Choi eigenvalues. Thus the improvement need not come from increasing the Kraus rank. These observations suggest looking for a covariant channel with adjustable Kraus weights. They do not establish spatial locality of its action. The matrices are expressed on reduced-state supports, not on individual physical sites.

[Proved] Rotated Petz maps and the universal twirled recovery of Junge, Renner, Sutter, Wilde and Winter are admissible. PROOF.md section 4 gives their formula and checks the marginal constraint directly, with the source https://arxiv.org/html/1509.07127v3, equations (18)-(21).

## Falsification search

[Numerical] The scan covered Delta = -0.9, -0.99, -0.9999 and -1, L = 4,6,8,10,12; TFIM g = 1 + c/L for c = -2,-0.5,0.5,2,5 on the same lengths; and 16 seeded generic Hermitian qubit bonds projected to the gauge-quotient M = 1 sphere, L = 4,6,8,10. The latter are radial projections of Gaussian coefficients in the nine two-site Pauli terms and three equally split fields, with no claim of uniform sampling on the M sphere. Small recovery SDPs were also solved at L = 4,6. Petz windows in the initial scan are l = floor(L/2) and floor(L/2)+1, clipped to the allowed range.

| Chain | L = 4 | L = 6 | L = 8 | L = 10 | L = 12 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| XXZ -0.9 | 1.88210 | 1.17713 | 0.79840 | 0.65196 | 0.59896 |
| XXZ -0.99 | 1.78771 | 1.10070 | 0.67287 | 0.46212 | 0.36048 |
| XXZ -0.9999 | 1.77788 | 1.09392 | 0.66138 | 0.44242 | 0.33266 |
| TFIM c = -2 | 0.47341 | 0.41484 | 0.40124 | 0.40124 | 0.40307 |
| TFIM c = -0.5 | 0.97327 | 0.78351 | 0.74589 | 0.73274 | 0.72616 |
| TFIM c = 0.5 | 0.99422 | 0.77084 | 0.74116 | 0.73074 | 0.72542 |
| TFIM c = 2 | 0.78238 | 0.52874 | 0.49555 | 0.47984 | 0.47020 |
| TFIM c = 5 | 0.41991 | 0.18976 | 0.15063 | 0.12942 | 0.11605 |
| random bond 5 | 3.24484 | 1.97891 | 1.88201 | 1.67707 | 1.27440 |
| random bond 13 | 2.56640 | 3.12115 | 2.67147 | 2.40149 | 2.08828 |

[Numerical] Bonds 5 and 13 were followed through L = 12, with all Petz windows tested through L = 8 and the two central windows at L = 10,12. Full small-l optimization gives Q_opt = 3.111376, 1.746836 for bond 5 at L = 4,6, and 1.914274, 2.735891 for bond 13. Thus both Petz and optimal defects can increase with L at small sizes. This is nonmonotonicity, not an unbounded-growth obstruction to P or P prime.

[Numerical reproducibility] Bond 13 has the representative h = sum_ab J_ab sigma_a tensor sigma_b + sum_a b_a(sigma_a tensor 1 + 1 tensor sigma_a)/2, with M = 1 to solver accuracy, where

    J = [[-0.034374340, -0.198108255,  0.330871226],
         [-0.507190933,  0.206203148, -0.047973979],
         [-0.216822016, -0.096748752, -0.206495638]],
    b = [0.030058455, -0.342363859, 0.122453127].

The seed and full precision matrices are in falsify.py and random_bonds.npy. Rounding this displayed representative changes M slightly.

[Proved and numerical] At Delta = -1 the full ground-space mixture is essential. Its optimal-window Q is L/(L-1), giving 1.33333,1.20000,1.14286,1.11111,1.09091 at these lengths. The two central windows instead give 1.33333,1.80000,2.13333,2.38095,2.57143, which increase toward the finite value 4. An initial sparse diagonalization missed ground-space vectors at L = 10,12. Those entries were excluded and replaced by the explicit staggered Dicke basis, with dimensions 11 and 13. The excluded records are preserved separately for audit. The analytic result, not incomplete eigensolver output, determines this endpoint.

[Heuristic assessment] No tested route constant grows without bound. The ferromagnetic crossover is compatible with lesson 3, and the g = 1+c/L scan probes lesson 2 without violating a bounded defect. The modest upturn for c = -2 is not a counterexample, and that scan changes the Hamiltonian with L. Random bond 13 is a single fixed Hamiltonian and is nonmonotone, then decreases through L = 12. None of these observations proves boundedness at arbitrarily large L.

## Remaining proof and computation

[Proved analysis] Ground-state stationarity supplies a lower bound on total channel energy cost, not the needed upper bound on one frozen bond. The modular formula is a product of noncommuting, potentially singular powers, with no lattice estimate replacing it by a local modular mismatch. A replacement channel built from another ground state loses the old boundary correlations; another ring's Petz map generally fails the prescribed-marginal constraint. The relevant toolbox results do not remove these gaps. PROOF.md section 7 gives the identities and precise failures of these particular arguments. It does not claim that every possible stationarity or modular argument is impossible.

[Assumed, single remaining sufficient lemma] For the normalized projector onto the L-ring ground space, the minimum energy defect p_L over the admissible channels and block lengths satisfies p_L <= K(2) M/L^2 with K(2) independent of the bond and L. Proving this at d = 2 alone would settle the requested qubit variant. Replacing 2 by d gives the general version. This is the remaining Lemma P prime, not a proved estimate.

[Numerical workflow] optimal_recovery.py contains the support-reduced primal, dual diagnostics and gauge-scale SDP. scan_sdp.py, refine_small.py, refine_other_windows.py, targeted_sdp.py and finish_special.py record the recovery comparisons. falsify.py and worst_cases.py record the searches. verify_c132.py passes the analytic ferromagnetic formula, a genuinely non-Petz periodization test, and an independent complex Choi contraction. All generated Python and outputs are inside this folder.

[Assumed continuation, not executed] From a copied work_C132 directory with numpy, scipy, cvxpy, SCS and Clarabel available, submit long_checks.slurm for the original long reproduction, sdp_large.slurm for larger recovery SDPs, and falsify_large.slurm for L = 14 through 20 continuation. These are prepared jobs, not submitted jobs. No network shell command, external-file edit or cluster computation was performed. STATE.md here records the handoff; the parent project's STATE.md was not edited because this task explicitly forbids writes outside the folder.
