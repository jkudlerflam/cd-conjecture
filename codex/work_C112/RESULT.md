# C112 result

[Proved within CFT] The pure stress-tensor coefficient is

    kappa = 1/(24 pi).

It multiplies A = sum_a c_a v_a, including neutral sectors and unequal velocities. The derivation uses no on-site symmetry. Single-cosine perturbation theory, the full circle diffeomorphism calculation, and the exact thermodynamic XX particle-hole response agree. A Gaussian bump has the same coefficient.

[Proved] The smooth-modulation lemma as written is false. Its right side uses gauge-invariant M, while its deficit depends on the energy-density representative. There are exact counterexamples at every d >= 2 with unbounded normalized deficit and A = 0. There is also a fixed-bond counterexample at the XXZ ferromagnetic endpoint. These statements do not refute a bound on A or a revised lemma that specifies an appropriate minimizing representative.

[Proved, infrastructure observation] The requested moderate-size unrestricted numerical survey was not run. Two Typhon SSH attempts timed out during banner exchange. No job was submitted, and no heavy computation was moved to this machine. The larger survey is prepared in cluster_job.sbatch. Completed local checks took about four seconds in total and are described below.

## CFT normalization and the ring correction

[Proved within CFT] Put W = integral w, ell_w = integral dx/w, and G = integral w'^2/w. With D_L defined relative to the uniform ring energy,

    D_L(w) = A G/(24 pi) + (pi A/6)(1/ell_w-W/L^2).

For w = 1+epsilon cos(2pi n x/L),

    D_L(w) = (pi A/(6L))(n^2-1)[1-sqrt(1-epsilon^2)],
    D_L/G = (A/(24 pi))(1-1/n^2).

In particular, n = 1 gives zero CFT adaptation deficit. Increasing L at fixed n does not remove this factor. The bulk-reference deficit additionally contains the uniform Casimir contribution. For a Gaussian of width sigma on the line, the second-order deficit is A epsilon^2/(48 sqrt(pi) sigma).

[Numerical] Exact free-fermion diagonalization confirms the finite-ring formula. The table uses epsilon = 0.02 and reports D_lattice/D_CFT, with the Casimir correction included. The Gaussian has sigma = L/10.

| L | XX, cosine n = 2 | XX, Gaussian |
| --- | --- | --- |
| 32 | 1.00807963 | 1.01258964 |
| 64 | 1.00201099 | 1.00311337 |
| 128 | 1.00050219 | 1.00077630 |

[Numerical] The fundamental XX cosine has deficit at most 7.11e-15 in absolute value, consistent with zero. Independently, the exact thermodynamic susceptibility gives effective kappa/kappa_CFT = 1.01078489, 1.00267399, 1.00066712, 1.00016670 at q = 0.4,0.2,0.1,0.05. The closed-form susceptibility and direct scalar quadrature agree to 8.6e-18.

[Numerical] Balanced-sector ED on the standard spin bonds gives the following D_lattice/D_CFT for n = 2, epsilon = 0.02. These are small-size checks, not extrapolated CFT determinations.

| L | XX | XXZ, Delta = 1/2 | Heisenberg |
| --- | --- | --- | --- |
| 8 | 1.14137811 | 1.18444523 | 1.19849529 |
| 10 | 1.08748252 | 1.11316813 | 1.12348321 |
| 12 | 1.05959812 | 1.07696881 | 1.08548837 |

[Numerical] Standard SU(3) and SU(4) swap chains were also checked in balanced color sectors, through L = 9 and L = 8. Their n = 2 ratios D_lattice/D_CFT at those largest sizes are 1.21561975 and 1.41525056. These wavelengths are too short for a sharp continuum test. Gaussian checks and the smaller sizes are retained in the numerical output.

## The counterexample and the largest ratios

[Proved] Let a = diag(-delta,delta,...,delta), X = |0><1|+|1><0|, and

    h = (a tensor I+I tensor a)/2 + t(X tensor I-I tensor X).

On a uniform ring H = sum_i a_i, so the gap is 2delta, M = delta, and A = 0. Set bar_w_i = (w_i+w_(i-1))/2 and g_i = w_i-w_(i-1). The modulated chain is a sum of independent one-site operators, giving the exact formula

    D_L = sum_i [sqrt(delta^2 bar_w_i^2+t^2 g_i^2)-delta bar_w_i].

For a slowly stretched profile,

    D_L/[M sum_i (w_(i+1)-w_i)^2/w_i] -> t^2/(2delta^2).

This disproves any finite all-bonds coefficient at fixed d. It is not a finite-size Casimir artifact. Removing the divergence gives the minimizing representative, whose deficit is zero.

[Numerical and proved] The largest explicitly evaluated value in this exact family uses t = 1, delta = 0.001, L = 4096 and a Gaussian. Only a scalar sum is needed, not a 4096-site diagonalization. The cosine gives the same limit.

| d | largest evaluated ratio | A | supremum over the stated bond class |
| --- | --- | --- | --- |
| 2 | 499999.999980 | 0 | infinite, proved |
| 3 | 499999.999980 | 0 | infinite, proved |
| 4 | 499999.999980 | 0 | infinite, proved |

[Proved] The standard XXZ bond has exactly zero deficit for every positive w at Delta <= -1. Add the divergence t(S_i^z-S_(i+1)^z) at Delta = -1. The uniform chain and M = 1/2 are unchanged, but for a fundamental cosine

    lim_(epsilon -> 0) D_L/(M G_discrete)
      = t^2 L/[4(L-1)sin(pi/L)^2]
      ~ t^2 L^2/(4pi^2).

[Numerical] At t = 1, epsilon = 0.001 and L = 6,8,10,12 the measured ratios are 1.19999907, 1.95097646, 2.90892016, 4.07131480, agreeing with the second-order formula. The linear-velocity coefficient approaches zero from the critical XXZ side. The endpoint itself has quadratic magnons and uniform energy -L/4, with no 1/L Casimir term.

[Proved] The endpoint failure also persists at fixed epsilon > 0. A smooth spin-coherent product state has a field-energy gain tending to a positive constant while its exchange cost is O(1/L). Since G_discrete is O(1/L), the ratio diverges at least linearly. This gives a fixed-bond, fixed-amplitude counterexample, without exchanging small-amplitude and long-wavelength limits. Section 5.3 of PROOF.md gives the bound. Scalar evaluation of that trial state at epsilon = 0.1 yields a lower bound on the ratio of 94.93, 399.74, 1618.98, 6495.97 at L = 128,512,2048,8192. These are variational lower bounds, not computed ground-state ratios.

## Random bonds and limits of the completed survey

[Numerical] Nine generic complex Hermitian bonds were drawn, three at each d. Every weighted Hamiltonian used the full Hilbert space. Each bond was evaluated both as drawn and after gauge minimization, with cosine and Gaussian profiles. A primal representative and a projected dual witness bracket M. The table reports the largest ratios within these small samples.

| d | L | largest ratio, raw representative | largest ratio, optimized representative |
| --- | --- | --- | --- |
| 2 | 8 | 0.0773776745 | 0.0446512756 |
| 3 | 5 | 0.2207099063 | 0.1987820769 |
| 4 | 4 | 0.0766413 to 0.0766465 | 0.0533731 to 0.0533768 |

[Numerical limitation] These are implementation checks. Three draws per dimension and these short rings do not establish a useful envelope, and their Gaussian widths are not in a continuum regime. A finite-ring gap does not classify a random bond's thermodynamic phase. The independently controlled critical examples are XX, XXZ, and the swap chains; the exactly solved gapped example is the onsite family above. No claim about uniformity of the optimized-representative lemma follows from the random table.

[Numerical] The largest many-body eigenvector residual was 4.40e-14. Gauge optimization changed the uniform ground energy by at most 1.16e-14. The analytic onsite spectrum agrees with independent full-space ED to 1.2e-16. Tiny negative deficits of order 1e-15, including the four-site SU(4) singlet, are retained as floating-point subtraction error.

[Assumed, unresolved] The remaining meaningful question is a revised lemma for a specified representative, together with a justified stress-tensor matching statement. Its validity is not settled here. The all-representatives version requested in the specification is disproved, and the failure does not make A large.

[Proved, file map] PROOF.md contains the derivations and the primary CFT source. NUMERICS.md gives reproduction commands and the unrun job scope. smooth_checks.py and scalar_checks.py are the executable checks. small_results.json and scalar_results.json contain the completed output. STATE.md carries this folder's session memory. No files outside work_C112 were modified.
