# C117 checks

2026-10-29.

[Proved] All 15 local lower-bound certificates passed independent integer reconstruction and residual verification. The displayed decimals below are rounded downward, so they remain lower bounds. Exact rational endpoints are in certificates/.

| Bond | k=2 marginal | k=3 marginal | k=3 ground constraints |
| --- | ---: | ---: | ---: |
| Heisenberg | -0.7500002 | -0.5000002 | -0.5000002 |
| XXZ_half | -0.6250002 | -0.4215354 | -0.4215354 |
| SU3_swap | -1.0000002 | -1.0000003 | -1.0000003 |
| random_real_d2 | -0.8657721 | -0.7815390 | -0.7815401 |
| random_real_d3 | -1.9497288 | -1.6404815 | -1.6404750 |

[Numerical limitation] The random d=3 k=3 searches, with and without ground constraints, hit the 500-iteration cap at tau=0.006, 0.0012, and 0.00024. The exact certificates remain valid. No claim of optimality or significant improvement from the ground constraints is made. Other solver termination messages and exact endpoints are in numerics.json.

[Proved] All 10 full-ring intervals passed integer verification. These are the k=L reference SDPs. Decimal endpoints below are rounded outward.

| Bond | L | Lower E_0(L) | Upper E_0(L) |
| --- | ---: | ---: | ---: |
| Heisenberg | 4 | -2.0000003 | -2.0000000 |
| Heisenberg | 6 | -2.8027764 | -2.8027756 |
| XXZ_half | 4 | -1.6861410 | -1.6861406 |
| XXZ_half | 6 | -2.3680346 | -2.3680339 |
| SU3_swap | 3 | -3.0000004 | -3.0000000 |
| SU3_swap | 4 | -2.0000005 | -1.9999999 |
| random_real_d2 | 4 | -3.1058004 | -3.1057998 |
| random_real_d2 | 6 | -4.6243579 | -4.6243561 |
| random_real_d3 | 3 | -4.8637567 | -4.8637554 |
| random_real_d3 | 4 | -6.4780815 | -6.4780776 |

[Proved] Exact interval widths are below 3.83e-6. Their endpoints are rational and checked through a residual row bound and a rational Rayleigh quotient. Full-ring certificates are in ring_certificates/. This certifies finite energies, not the asymptotic Casimir limit. The SU(3) L=4 result is outside the singlet vacuum sequence and was not used as a Casimir estimator.

[Numerical] Independent floating ground-state diagnostics give maximum overlap-consistency residual 5.56e-17, maximum stationarity residual 1.45e-15, and minimum one-site ground Gram eigenvalue -3.04e-16. The last negative value is at floating roundoff. Every local lower certificate lies below the corresponding four-site ring energy. These diagnostics are separate from the integer certificate proofs.

[Numerical] The massive Ising q^4 check gives a_4/(g J^2) equal to 0.2499414146, 0.2497657572, and 0.2490645144 at J/g=0.025, 0.05, and 0.1, approaching the derived value 1/4. The direct double commutator and its degree-two cosine polynomial agree below 3e-18. See moment_check.json and PROOF.md Section 6 for the small finite-ring a_2 qualification.

[Proved reproducibility] Local commands used were python3 bootstrap.py --small, python3 verify.py, python3 ring_checks.py, and python3 moment_check.py. No cvxpy or other package installation was used. The extended bootstrap driver requires SLURM_JOB_ID. No Typhon job ran because the SSH connection timed out during banner exchange.
