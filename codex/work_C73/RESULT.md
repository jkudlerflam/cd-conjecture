# C73 result

[assumed, target unresolved] Lemma TP remains unproved for weak XXZ off half filling. No coefficient counterexample to TP was found. The results below advance the ultraviolet preparation, but do not complete the analytic integration of both chiralities or the infrared comparison.

[proved, analytic partial preparation] [PROOF Sections 2-4](PROOF.md#2-all-long-chiral-cycles-have-uniform-absolute-moments) integrate one high chirality of the C65 trace-defined reference, with all its closed cycles and open chains retained. The expansion converges on a coupling disk independent of the sea cutoff. It includes the nonlocal six-field and higher kernels. In every fixed polynomial moment norm, its thermal replacement obeys

    ||partial_A^j(F_+^T-Per_beta F_+^0)||_(H,m)
                                     <= C |g|^2 T^P,
    j<=2,   P>1.

The bubble uses C71's prepared distribution. Every longer cycle is controlled by generalized Young inequalities with exponent r/(r-1)<2. Nonzero cycle windings give the thermal gain with only a polynomial cost in r. Open chains periodize exactly. The moments and field radius are fixed once. These are thermodynamic kernel estimates, not a uniform finite-spatial-volume theorem.

[proved, complete fixed-order ultraviolet statement] [PROOF Section 6](PROOF.md#6-exact-cutoff-stabilization-at-each-fixed-order) treats every zero-temperature coefficient of the full reference preparation, including vacuum coefficients. In a fixed slow-field representative,

    U_(n,Lambda)^0 = U_(n,R_n)^0 for Lambda>=R_n,
    R_n=max(K,K_e)+(n+1)kappa+1.

Open chains are anchored to bounded external momenta. A closed chiral loop vanishes after frequency integration unless its line energies straddle the Fermi point. Bounded interaction transfers then confine all contributing momenta to an interval of size O(n kappa). Mixed high/low cycles are bounded in momentum by their low insertion. The finite-cutoff representative has finite fixed moments. Thus there is no finite-order zero-temperature ultraviolet coefficient divergence in this reference construction. The constants can still grow with n, so this does not prove analyticity of the full prepared action.

[proved, Ward correction retained] The high-covariance identity is

    D(q)a(k)a(k+q)=a(k)-a(k+q)
                           +f(k)a(k+q)-f(k+q)a(k).

The low-cutoff terms survive sea-cutoff removal. Odd high cycles cancel by reversal, but a four-cycle is generally nonzero. PROOF Section 5 gives an explicit positive example, and approximates its cutoff by allowed compact plateau cutoffs. This is the first surviving symmetrized high cycle beyond the bubble. It is finite and satisfies the thermal bound, so it is not the requested uncancelled obstruction to TP.

[assumed, exact remaining work] Prove an exponential-in-order bound for the coupled two-chirality preparation, with its full vacuum scalar and thermal derivatives. The one-chirality coefficient norm does not justify the remaining singular covariance integration. Then establish the complete lattice/reference defect and mixed thermal histories, common parameter domains, and vacuum derivative/Hessian estimates in C59 TP. C62's scale sum and the collaborator's fixed thermal-mark allocation apply only after that constructive bound is supplied. No interacting XXZ thermal theorem or retarded-action extension follows yet.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. Controls cover the contour identity, cutoff plateau, Ward remainder, nonzero four-cycle and cycle winding decomposition. The largest final winding residual is 1.08e-16. The initial insufficient image truncation and its correction are recorded in PROOF Section 9. No interacting flow was simulated.

[proved work scope] PROOF.md contains the derivations, PROGRESS.md the work record, SOURCE_AUDIT.md the attribution boundary, and STATE.md the local handoff. All writes remain inside work_C73. No independent or human review is claimed.
