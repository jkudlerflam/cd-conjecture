# C105 result

[Assumed, unresolved] I did not prove A <= C M f(d) for general interacting chains. The smooth or integrated formulation of that bound remains open. The literal unsmoothed lattice statement needs correction before it can be used as an equivalent formulation.

[Proved] A lower-energy certificate with logarithmic loss exists for XX. For nonnegative site amplitudes a_i of finite support and w_i = a_i a_(i+1),

    E0(K_w) >= -(1/pi) sum_i a_i^2,
    e0 sum_i w_i - E0(K_w)
        <= (1/(2 pi)) sum_i (a_(i+1)-a_i)^2.

The proof sandwiches the weighted one-particle hopping matrix between plus and minus D_a |T| D_a, where T is the uniform hopping operator. An explicit gamma-function factorization of w_x = x(N-x)/N gives

    Delta_N <= (1/(4 pi)) log N + O(1).

This is weaker than the existing quasi-free result A <= 2M. It provides a concrete nonlocal positive certificate for the weighted problem. The naive replacement of the XX coefficient by the centered interacting ground energy is disproved by an explicit four-site Heisenberg singlet.

[Proved] Allowing every positive, overlapping uniform-interval decomposition does not cure the boundary loss. For XX, each open interval has defect at least

    g_* = 1/sqrt(2) - 2/pi = 0.07048700881897... .

The sum of decomposition coefficients is at least max w, so every allowance constructed this way is at least g_* max w. For quadratic weights this grows linearly with N. This is a restriction on that proof method, not on the actual defect.

[Proved] A second obstruction applies to fixed-range local-marginal SDP certificates. The XX vacuum has full-rank marginals on every finite interval. Moving a small distance from its r-site marginal away from the maximally mixed matrix preserves positivity and every local overlap constraint while lowering the bond energy below e0. The resulting feasible point forces every such fixed-r dual certificate to lose a positive constant times W, which grows as N^2. The obstruction does not apply to ranges growing with N, additional ground-state constraints, or nonlocal covariance certificates.

[Proved] For the half-line XX chain the raw local deficit has the exact spectral formula

    delta_x = 2 integral_0^infinity
       P_(x-1)(lambda) P_x(lambda)
       sech(pi lambda)/(1+exp(2 pi lambda)) dlambda,

where (n+1)P_(n+1) = 2lambda P_n - nP_(n-1). In particular delta_1 = 4G/pi^2 - 1/pi = 0.05291698653, above the sharp proposed all-bond allowance 1/48.

[Numerical] The large-distance raw and filtered profiles are very different. In the exact infinite-half-line calculation,

| bond x | x^2 delta_x | x^2 (delta_(x-1)+2delta_x+delta_(x+1))/4 |
| --- | --- | --- |
| 33 | 0.84788490 | 0.01280006 |
| 257 | 4.95491062 | 0.01321859 |
| 1025 | 16.93641809 | 0.01325355 |
| 8191 | 111.42982189 | 0.01326197 |

The smooth CFT value is 1/(24 pi) = 0.0132629119243. Quadrature orders 400 and 800 agree to 2.08e-15 in the unscaled deficit through x = 8193. This scalar recurrence took 0.12 seconds, without a finite-chain diagonalization.

[Heuristic] The spectral expansion predicts a staggered term (-1)^(x-1)/(2 pi x log x). I derived its coefficient at the level of asymptotic spectral integration, but did not complete the uniform remainder estimate. Therefore divergence of the raw scaled profile is not claimed as a proved theorem here. The raw observable cannot be replaced silently by the filtered one in a proof.

[Numerical] The following small finite-interval checks were run for quadratic weights. Each slope uses the two largest sizes in the same residue class. The interval allowance is optimized over all positive interval decompositions.

| model | largest N | actual defect | optimized interval allowance | last log slope | A/(12 pi) |
| --- | --- | --- | --- | --- | --- |
| XX | 12 | 0.13122785 | 0.33718102 | 0.02396377 | 0.02652582 |
| Heisenberg | 12 | 0.22135420 | 0.58774499 | 0.03964901 | 0.04166667 |
| XXZ, anisotropy 1/2 | 12 | 0.17912500 | 0.47395777 | 0.03189360 | 0.03445806 |
| SU(3) swap | 9 | 0.50656751 | 0.83178586 | 0.10340385 | 0.11111111 |

The corresponding ramp checks also satisfy the interval inequalities. The largest eigenvector residual across these computations was 6.67e-12. These are finite-size numerical checks in balanced sectors, not asymptotic proofs or certified SDP lower bounds.

[Proved, infrastructure observation] Typhon was unreachable on two short SSH attempts, both ending with "Connection timed out during banner exchange". No job was submitted. The large finite-chain Gaussian runs are prepared in checks.py and cluster_job.sbatch but remain unrun. No heavy numerical work ran locally or on a login node.

[Assumed, next mathematical steps] A general proof still needs either (i) an interacting replacement for the nonlocal positive completion in the XX certificate, (ii) an SDP hierarchy with a controlled growing range, or (iii) an upper bound on the inverse-energy gradient response that does not introduce a gap or crossover length. The present calculation supplies none of those three inputs.

[Proved, file map] Full derivations and primary-source links are in PROOF.md. Reproduction details are in NUMERICS.md. The session memory is STATE.md and the work record is PROGRESS.md. No files outside work_C105 were modified.
