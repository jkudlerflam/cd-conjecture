# C135 result

[Proved] The LTI route closes for every nearest-neighbor quadratic qubit bond, including pairing, transverse fields and DM terms in a Jordan-Wigner frame:

    e0(h)-E_LTI(n) <= M(h)[1-cos(pi/(n+1))]
                    <= pi² M(h)/(2(n+1)²).

Thus (L-LTI) holds on this class with K=pi²/2, uniformly in gaps and crossover scales. The proof smooths the finite LTI covariance with a normalized sine window, preserving covariance positivity and the particle-hole constraint. It does not assume ring saturation or a CFT limit. See [PROOF.md section 5](PROOF.md).

[Proved] For XX there is an explicit telescoping certificate with K=pi. Set v_j=sqrt(2/(n+1))sin(pi j/(n+1)), s=cos(pi/(n+1)), a_j=v_jv_(j+1)/s, and

    X_n=sum_(j=1)^(n-2) [sum_(k=1)^j a_k-1] h_(j,j+1).

Then

    h_12+X_n tensor I-I tensor X_n
      >= -I/[pi cos(pi/(n+1))]
      >= [e0-pi M/n²]I.

The counterterm contains only weighted translated bonds. Positivity follows from an integral of positive quadratic operators built from extended fermionic modes. This applies C105's weighted XX inequality to a sine window and converts it to the requested LTI certificate. It is not a claim of historical priority.

[Proved] Ring saturation is false, even for a real reflection-symmetric free bond, and allowing all on-site symmetry twists does not repair it. For

    h=(XX+YY)/4+(ZI+IZ)/8,  M=1/2,

the exact three-site optimum is E_LTI(3)=-5/12. The minimum over e0 and every untwisted ring of length at least three is -3/8. The gap is 1/24, or 1/12 in units of M. At n=4 an explicit rational Toeplitz covariance gives energy -7553/20000=-0.37765, whereas every symmetry-twisted ring with L>=4 has density at least -3/8. This gives an exact feasible gap of 0.00265, or 0.0053 in units of M. Positive integer principal minors certify the covariance, and an analytic bound handles the entire ring-length tail. These are proofs, not finite-ring extrapolations.

[Numerical] The optimized n=4 value is -0.37837433475, so the observed gap from the best twisted ring is about 0.00337433475. The rational feasible point proves the counterexample independently of that decimal optimum. The LTI hierarchy imposes a stronger finite-level requirement than the ring envelope already within the free class, even though both classes of bounds have a uniform n^(-2) rate there.

[Proved] F1 to F4 hold, with support completion for singular marginals. Twisted rings require a qualification. A twisted-translation invariant ground state, in a fixed linear ordering avoiding the seam, has LTI windows. Arbitrary ground vectors and ordinary cyclic windows crossing the seam do not automatically have this property. An explicit three-site twisted XX ground state refutes the unqualified statement. Reversal and complex-conjugation symmetrizations preserve LTI under the stated bond symmetries.

## Dual structure and periodization

[Numerical] Optimal full-spin dual matrices were extracted for XX, Heisenberg and critical TFIM at every n=4,...,8. Their matrices are in matrices/, and their Pauli decompositions and bond fits are in dual_structure.json. The following entries give the Frobenius residual after fitting X by translated copies of the original bond, divided by ||X||_F. Scalar identity components are removed first.

| Model | n=4 | n=6 | n=8 |
| :--- | ---: | ---: | ---: |
| XX | 0.08723 | 0.03075 | 0.03826 |
| Heisenberg | 0.000218 | 0.04019 | 0.04773 |
| Critical TFIM | 0.08462 | 0.13038 | 0.08220 |

[Numerical] These solver-selected duals contain Pauli strings acting on more than two sites. At n=8 their squared-coefficient fractions are about 0.00120, 0.00219 and 0.00527, respectively. Optimal duals need not be unique. A second optimization restricting X to weighted translated bonds reproduces XX's unrestricted optimum through n=8, and Heisenberg's through n=6. At Heisenberg n=8 the restricted value is about 1.8e-6 lower. For TFIM the restricted value at n=8 is -1.2832824, below the unrestricted -1.2814577 by about 0.001825. This does not show that multi-site terms are necessary after allowing independently weighted on-site and interaction terms. The analytic XX certificate already proves they are unnecessary for a Casimir-rate bound.

[Proved limitation] Complementary slackness gives S rho=0 for S=h+D(X)-E_LTI(n)>=0. For the n-site marginal eta of a final TI periodized state, the energy defect is Tr(S eta). This is nonnegative but has no upper bound from slackness alone. It depends on the recovered state's weight outside the kernel of S. The replacement channel supplies an order-one defect despite exact slackness, so admissibility plus slackness is insufficient. No new Petz or optimal-channel estimate is proved for interacting qubit bonds. Classical conditional extension has zero defect; the known symmetric ferromagnetic mixture has defect 1/[2n(n+1)]. Covariance smoothing proves the free LTI rate without proving that a channel with the prescribed marginal constraint attains it.

## Falsification computations actually run

[Numerical] The local scan includes all requested bond families and 30 random bonds, 62 fixed bonds in total. [TABLES.md](TABLES.md) reports every requested L=3,...,20 and n=2,...,8 entry, marking uncomputed entries NA. Machine-readable records are in tables.json, results/ and results_distance/. [NUMERICS.md](NUMERICS.md) gives the formulas, repairs, solver diagnostics and reproduction instructions.

[Numerical] Completed coverage is 738 ring-energy entries, 345 LTI entries at n=2,...,8, another 48 free-model entries at n=9,10, 372 trace-distance hull upper bounds and 161 trace-distance lower bounds. The shared dictionary has 952 points, including six uniform product MPS. Free ring energies reach L=20 through their analytic spectra. Interacting and generic rings reach L=10, with Heisenberg through L=12. Interacting symmetry-reduced LTI solves reach n=8, real random solves n=6 and generic complex solves n=4, plus one n=5 solve. Larger missing cases are implemented in the SLURM continuations and were not run.

[Numerical] Selected bounds on the LTI rate are below. Entries are [n²(dual_lower(n+2)-primal_upper(n))/M, n²(e0_up-dual_lower(n))/M]. NA means the n+2 solve is missing. The upper density comes from the shared variational dictionary, so it includes a finite dictionary error.

| Model | n=2 | n=4 | n=6 | n=8 |
| :--- | :--- | :--- | :--- | :--- |
| Heisenberg | [2.0000, 2.5026] | [1.0519, 2.0103] | [0.7733, 2.1565] | [NA, 2.4587] |
| XX | [1.1716, 1.4836] | [0.6470, 1.2480] | [0.4819, 1.3522] | [0.3883, 1.5473] |
| XX with b=1/2 | [0.9730, 1.1571] | [0.2037, 0.7363] | [0.5083, 1.1983] | [0.2757, 1.2267] |
| Critical TFIM | [0.3045, 0.4094] | [0.2111, 0.4195] | [0.1640, 0.4688] | [0.1343, 0.5418] |
| XXZ, Delta=-0.9999 | [0.0001333, 0.0004000] | [0.0001066, 0.0010667] | [0.0000902, 0.0021600] | [NA, 0.0036573] |

[Numerical] Across the computed entries, the largest directional lower endpoint is about 1.39563, attained at L=4 for the chiral Heisenberg bond, tied with Heisenberg within numerical error. The largest two-level LTI lower endpoint is 2 at n=2 for Heisenberg and its chiral counterpart. The largest dictionary upper endpoint for the LTI constant is about 3.25154, for the chiral bond at n=8. No unbounded family is established or refuted by these finite data.

[Numerical] The largest computed L² d_up is 3.94536 at random_29, L=6. For Heisenberg, the [L² d_low,L² d_up] brackets at L=4,6,8 are approximately [1.05186,2.01029], [0.77347,2.15649], [NA,2.45865]. The dictionary upper bounds can increase because the dictionary has a fixed approximation error. In particular the XX directional upper bound reaches 4.05638 at L=20, while the analytic free-chain answer remains bounded. That increase is not a lower witness against (K).

[Numerical, exploratory] Four-start eight-angle searches gave the best Petz upper-endpoint values

| L | Best value | Location |
| :--- | ---: | :--- |
| 6 | 4.053601 | Heisenberg |
| 8 | 3.648031 | DM-and-field point below |
| 10 | 3.440815 | Heisenberg |

At L=8 the normalized coordinates (Jx,Jy,Jz,Dx,Dy,Dz,bx,by,bz), with spin operators S=sigma/2, are

    (0.16150235,-0.17347004,-0.17534896,
     -1.98308643,0.07663464,-0.13286185,
     1.83501071,-0.04930471,0.08610785).

This point is near the U(1) family with a strong DM term and field along the same axis, but it does not establish a systematic drift. The n=2 to n=4 LTI-difference search has best value 2 at Heisenberg. Every search reached its evaluation limit, so none is reported as a converged global maximum. The largest-level all-angle search, n=6 to n=8, is prepared for SLURM and unrun. Positive directional lower witnesses at L=8 and 10 require hierarchy levels above those lengths, which the generic local scan does not supply.

[Numerical hygiene] Dual eigenvalues and repaired primal states, rather than raw objective subtraction, determine the brackets. The largest retained LTI equality residual is about 2.1e-8 and the largest repaired energy primal-dual gap about 3.4e-5. Eleven retained LTI entries at n<=8 and 28 distance lower solves have status optimal_inaccurate and are exploratory. The worst distance equality residual is about 9.2e-8. The largest ring eigenvector residual is 4.6e-11 and the largest Petz trace residual is 5.1e-11. Cross-checks found no reversed brackets or corrected LTI lower bound above a compatible computed ring energy. These are floating-point validations, not interval certificates or a proof that every sparse eigensolve found the ground eigenvalue.

## What remains

[Assumed, open] The single sufficient lemma for general interacting qubit bonds is a uniform telescoping certificate

    h tensor I+D_n(X)>=[e0(h)-K M(h)/n²]I

with K independent of h and n. It is proved here for the quadratic class, not for all d=2 bonds. Exact ring saturation cannot serve as an intermediate lemma because both proposed forms are false. The recovery version still needs an upper estimate on the energy carried outside the optimal dual slack's kernel.

[Assumed continuation, not executed] energy_scan.slurm, deeper_lti.slurm, distance_scan.slurm and search.slurm cover the missing numerical work. No cluster job was submitted. The parent project's STATE.md was not edited because the task forbids writes outside this folder; [STATE.md](STATE.md) records the complete handoff here.
