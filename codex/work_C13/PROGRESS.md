# C13 progress

Status: source audit in progress. No earlier PROGRESS.md existed. All writes stay in work_C13. The project STATE.md was consulted but will not be modified under the task restriction.

Read STATUS.md, joint_fluid.md, helix_equalizer.md, length_fluid.md, and the ir_control README and operator analysis. Source claims remain source claims until independently derived or checked here.

Key issues to resolve
1. The physical length parameter is K = 1/(b+nu)^2, so treating the extra mode as an unconstrained ordinary fermion can miss low-dimension harmonics.
2. Distinct fillings alone do not prohibit zero-momentum pairing. The parity-only target needs shifted Fermi-sea centers, as already required by ir_control.
3. The joint cyclic constraint survives interactions as a Hilbert-space constraint, but separate gap and payload momenta need not survive.
4. Generic chemical potentials can avoid exact finite-ring charge crossings, but no uniform gap from those crossings follows.

Next: finish the source audit, derive the physical zero-mode lattice and resonance conditions, then run small interacting tests with calibrated entropy and a spectral/thermodynamic estimator.

## Step 1, source audit and continuum reduction

Read the requested C8 notes and RESULT sections 12 through 17, including the wrapper and the two machinery READMEs. No files outside this folder changed.

Derived the physical current twist: for payload momentum S=p.Q/(2 pi), the length current energy uses (D-S)^2/(b+nu)^2 and physical momentum is 2 pi(D-S)/(b+nu), up to an overall Fourier convention. At a translation-invariant resonance D-S=h(b+nu), its length scaling contribution is h^2. Thus the small physical K does not by itself imply infinitely many relevant uniform operators. With marker number conserved and payload parity even, the new dimension-two resonances are p.Q +/- 2 pi nu = 0 mod 2 pi for payload bilinears Q of squared norm two. Pure length locking at denominator q has dimension q^2 at the free point.

This is a bosonization derivation, not a microscopic theorem. Need retain the shifted Fermi-sea centers if flavour U(1) is absent. An all-ring unique ground state can hold generically outside countably many charge-crossing sets; it is stronger than bounded degeneracy and does not define a uniform open condition.

Numerical plan: exact joint-model calibration using entropy and Casimir; a finite local diagonal-unitary deformation supplies solvable interacting examples; small physical-graph ED with independent density and correlated-hopping interactions tests generic deviations. Report the small-size limits rather than equating toy tests with microscopic SW control.

## Step 2, first numerical tests completed

Implemented test_fluids.py using numpy/scipy only. Exact controls at M=40,64,96,128,192 compare squeezed-coordinate entropy and the shell-corrected physical Casimir coefficient. At M=192 the estimates are c=2.000301/2.000062 for n=1 and 3.000672/3.000094 for n=2. These are controls, not generic finite-coupling evidence.

Physical-graph ED for b=2, (M,r)=(6,2),(7,3),(8,3), physical L=14,17,19, compares the decoupled reference with length/payload density interactions and conditional hopping. Dimensions reach 29792. All tested fixed-charge ground states are unique. The free physical-site entropy fits overestimate c substantially (2.44 and 3.56 at L=19), so raw interacting fits cannot be quoted as accurate c estimates. The generic perturbations change the fitted entropy coefficients by less than 0.004. The complete data are graph_results.json.

Exact all-charge minimization on every L=40 through 400 finds unique interior sectors for the unperturbed model at the chosen fixed chemical potentials. This includes all primes in the interval. This test does not establish global uniqueness for the generic interacting graph.

Next: complete the mixed-margin example and contiguous-ring ED, write the detailed conditional continuum result, and separate finite-L uniqueness from a uniform untwisted CFT limit.

## Numerical verification correction

A repeat using two rather than three requested Lanczos eigenvalues exposed a missed momentum sector at n=1,L=19. The supplied free ground vector was translation invariant, so it was an inappropriate Krylov seed for searching the entire spectrum. Both scripts now use a fixed-seed random starting vector, which has support in every symmetry sector. All affected graph data will be regenerated before quoting gaps. This is a numerical-method issue, not evidence of a gap change.

## Step 3, report and zero-mode lattice

Wrote RESULT.md with the finite-order correction taxonomy, complete potentially relevant/marginal resonance list, a constructive length-mode extension of the supplied Sidon margins, Gaussian positivity, spin structures, the interacting cyclic quotient, uniqueness conditions, and the accumulation-point qualifications.

A sharper distinction emerged: irrational equilibrium charge rounding prevents convergence to one fixed untwisted vacuum sector along every integer ring length, even in an otherwise ordinary Luttinger liquid. It does not prevent each finite ring from having a unique vacuum. The report makes both conventions explicit.

Changed the generic numerical model to fixed centering constants (nu_ref=0.37,f_ref=(0.27,0.16)) on every ring, rather than using the sector's own density. This makes the generic sequence one fixed projected Hamiltonian for each n. The graph data are being regenerated with the corrected random Krylov seed.

The full accumulation point cannot be identified before the fixed-alphabet bond and its normalization exist. The bare long-cell selector does approach mu=2t with quadratic head dispersion; the equal-hopping selector has a different threshold, |mu|=4a. The report does not conflate them or assert local lattice-c bounds from the continuum c-theorem alone.

## Final status

RESULT.md is complete, with a conditional positive Gaussian-phase analysis and a negative answer to the stronger fixed-vacuum all-integer-L interpretation. The report is about 6800 words and tags each claim or paragraph by evidence status.

Final files
- test_fluids.py, controls.json, graph_results.json
- additional_checks.py, additional_results.json, contiguous_rings.json
- tower_checks.py, tower_results.json
- RESULT.md and this PROGRESS.md

Final verification
- Both ED suites reran after fixing the Krylov seed and using fixed centering constants for the interacting Hamiltonian. Repeated gaps agree within 1.43e-14.
- Twelve physical-graph entropy/ED cases and twelve contiguous-ring fixed-charge ED cases completed. Largest dimension 29792. Free graph residual at most 2.09e-15.
- All 722 grand-canonical reference checks (n=1,2, every L=40,...,400) have unique interior winning sectors; minimum intersector separation 8.80269e-5.
- Entropy/Casimir controls at M=192 give 2.0003014/2.0000625 and 3.0006723/3.0000937. Generic physical-site entropy remains biased by the same large small-size error as the free control, so no accurate generic interacting c is claimed.
- Mixed resonance enumeration at n=1,2,4,8 satisfies the conservative analytic margin bound.
- Current-lattice checks at M=40,80,160,320, for n=1,2 and independent/common spin structures, satisfy all integer gluing constraints. The current identity agrees within 8.4e-17; cosine-expansion errors have the expected bounded L^(-3) coefficients.
- All Python sources compile. No em-dash prose, network access, external writes, or sub-agents were used.

What remains unresolved
1. The actual high-channel microscopic core and its finite-coupling coefficients, including leakage control.
2. Two accurate independent estimators for a generic interacting microscopic mobile-cell chain. The present calibrated generic tests are too small.
3. All-ring uniqueness for the parity-only displaced microscopic target, beyond the proved positive-hopping subclass.
4. Whether the conjecture permits the usual size-dependent charge/flux sectors or demands one fixed vacuum sector along all integers L.
5. The normalized fixed-d accumulation Hamiltonian. The bare-selector z=2 mechanism is identified; the unfinished full family is not.

No files outside work_C13 were modified, including project STATE.md. To continue, begin with RESULT.md Sections 2, 7, 8, and 10, and the supplied core specification when one becomes available.
