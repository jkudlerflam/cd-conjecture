# C63 progress

2026-10-28

[proved, scope] The proposed exclusion of two critical Ising sectors is false. The nearest-neighbor XY chain can have two Majoranas and total c=1. The useful target is exclusion of total c>1, including a third Majorana or an independent neutral sector above a Luttinger liquid.

Read the project memory first, then the d=2 STATUS entries, C6 parametrization and symmetry/free proofs, C15 weak/strong-coupling results, C16 continuum and covering limitations, C17 anomaly analysis, C34 symmetry partition and exclusions, and the four requested numerical contribution summaries. The gap-reduction notes distinguish the staggered-charge family from the much larger generic residual.

[assumed] All external exact-model, gap-theorem, and CFT inputs will be inherited from these local sources. No network access is used. No files outside this folder will be modified, including the parent STATE.md.

Work underway

1. Write a disjoint stratified covering with explicit unresolved rows, rather than claim that perturbative neighborhoods exhaust parameter space.
2. Audit the weak-coupling mode count and anomaly constraints, including the allowed Ising times boson theory.
3. Give a catalog of candidate CFT classes and distinguish actual range-one qubit realizations from symmetry compatibility. No exhaustive c>1 classification is available.
4. Test the staggered-field family using small charge-resolved ED, interval entropy, finite-size energy, and neutral/charged gaps, calibrated against known one-boson points. Focus on the BKT region and the quadratic factorization boundary where false c>1 values are most plausible.

No full proof or violating candidate has been established.

[numerical] The main ED scan is complete at L=8,10,12,14,16,18 for nine points, with four selected points extended to L=20. Controls include XX and XXZ. A separate Ising plus XX control tests sensitivity to c=3/2 without changing the target family. Every charge sector is checked at L=18. No supported c>1 candidate appeared. Near u=sqrt(1+b^2/4), chord and size-entropy estimates disagree with the velocity/Casimir extrapolation, as expected near the nonconformal boundary.

[proved, new explicit diagnostic] The boundary filtered-Dicke state has a Schmidt spectrum obtained from products of generating-polynomial coefficients. Its entropy is (1/2) log[L x(1-x)] plus a constant at fixed b and block fraction x. Therefore the periodic size-entropy estimator tends to 3/2 even though the exact defect branch has z=2. The large-L chord fit instead tends to about 1.241 in the chosen quarter-to-half window. This is a concrete false-positive mechanism, not a violating CFT.

[proved conditional on standard CFT input] A nontrivial U(1) current in a unitary compact CFT contributes Sugawara c=1. Unitarity of its commuting Virasoro coset implies that any nontrivial neutral coset has c>=1/2. Thus within 1<c<=3/2, the charged staggered family would require c=3/2 and an Ising coset, possibly with a non-product extension. This does not apply if the microscopic charge acts trivially on every gapless field.

Completed

- RESULT.md now contains the full stratified covering, the free/integrable weak-coupling audit, the anomaly and operator analysis, a qualified CFT catalog, the finite staggered-family test, and all steps of the filtered-state Schmidt derivation. Unresolved rows retain explicit status tags.
- The main 54 ED records and four L=20 extensions are complete. The 28-record pilot is retained separately. All charge sectors at L=18, the Ising/product controls, and the exact-boundary coefficient calculations completed. No job remains running.
- The strongest tested overshoot, (u,b)=(1.4,2), fails consistency between c_chord=1.11421, c_size=1.24686, and the conditional one-speed energy estimate 0.75143. The clean point (0.9,0.5) gives 1.02458, 1.02540, and 1.00391. These numbers do not establish a global exclusion.
- Exact boundary-state residuals and direct entropy comparisons pass at L=8. Largest scan eigenvector residual is 8.5e-10. All Python files parse, all expected records and outputs are present, and the prose restriction check passes. validation.json records versions and the convention-check residuals.
- STATE.md contains the local working memory. The project-level STATE.md, sibling directories, and repository files were not modified. No network or cluster commands were used.

Final status: the requested global proof is not established, and no c>1 candidate is supported. The new concrete outputs are the calibrated staggered-family test and the exact c_size=3/2 false-positive mechanism. The corrected sector lemma still needs a dynamical exclusion, and a non-product-CFT loophole remains even beyond it.
