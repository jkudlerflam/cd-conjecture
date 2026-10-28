# C44 result

[proved, finite-dimensional verdict] The C8 mobile wrapper applies to C39's 46-letter core for every n. It gives one repeated nearest-neighbor bond on d0=4464 states, with selected physical lengths n+3 and n+4 carrying the same complete 2^n-dimensional payload. The wrapper does not require equal vacuum amplitudes. C39's biased vacuum satisfies its actual hypotheses.

[assumed, remaining conclusion] This is a valid fixed-alphabet mobile candidate for every n, conditional on the infrared step for the claim c=n+1. That step must include volume-independent local elimination, stability of the massive auxiliary sectors at nonzero coupling, and Gaussian stability relative to the smallest payload velocity and all mixed momentum margins. Those statements remain assumptions. The audit establishes neither a microscopic central charge nor an unconditional counterexample.

The derivation is in [PROOFS.md](PROOFS.md). The main points are as follows.

1. [proved] Positivity, initialization, the complete code, the core gap, endpoint sources and exact U(1) all survive the wrapper. The occupation filter excludes every wrong occupied number. Arrow counting and the selector secant bound exclude illegal sectors, marker-free sectors and all unwanted cell lengths. These bounds tolerate arbitrarily large positive C39 express energies and filter penalties. On every ring with L>=(n+3)(n+4), the unperturbed ground space consists precisely of the selected tilings and their payload codes.

2. [proved] The compensated second-order bond has the exact source-resolvent formula in PROOFS (11), for all four ordered length pairs. Writing s for the core scale gives

       s K2_(p,q)=K_bare+E_(p,q).

   PROOFS (13) bounds the full operator norm, including its scalar part, by an explicit sum of a missing-adjacency correction and a fast-wrapper correction. The former vanishes as the square of the vacancy ratio. The latter requires s small relative to the actual charge tunneling gap, of order U rho^(2n+4). A separate nonasymptotic positive gap bound makes the parameter prescription finite and explicit.

3. [proved] Equalization is arbitrarily accurate, not exact at finite wrapper parameters. The scalar combination c_22-c_21-c_12+c_11 is a real length-density interaction. If every E_(p,q) has norm at most epsilon, this scalar interaction has magnitude at most 4epsilon, and the corresponding traceless mixed length-payload term has norm at most 8epsilon. Both retain physical prefactor lambda^2/s.

4. [proved] The mixed u_mov lambda coefficient vanishes by vacuum-versus-bus conservation. Pure marker excursions are scalar on payloads. Terms of order lambda^2 u_mov, marker-readout terms and readout leakage need not vanish; their resolvent formulas and finite-cluster bounds are included. No even-in-lambda assertion is made for odd periodic cell rings.

5. [proved] Marker transport is a positive scalar times the identity on both payloads. Its leading length Hamiltonian is hard-core hopping. The physical ring retains the joint length/payload momentum constraint, including cyclic payload relabeling when a marker crosses the origin. Finite wrapper errors perturb this reference length fluid.

6. [proved] The all-n parameter budget retains C39's visibility factor

       J0=2 g^2 q^(2n+2) D^(-2n-4)/S^2.

   Wrapper errors are chosen below J0 times the certified narrow velocity and momentum tolerances. A conservative mixed-momentum margin is pi/[4(5(4n+1)^4+1)], obtained by an explicit finite exclusion of length densities. PROOFS (19)-(23) specify the ordered choices of vacancy ratio, core scale, marker hopping and port coupling. All coefficients are fixed before L is varied. No uniform lower bound in n is asserted.

[numerical, small wrapped cells] At n=1, diagonalization of the complete legal core, vacancy and occupation-filter sectors gives the predicted two-dimensional code at both lengths. The matrix dimensions are 1728 and 4320. The next levels are 0.01475806347 and 0.00003363095421 at the stated diagnostic parameters. Code residuals are below 3.4e-16 and eigenpair residuals below 1.1e-13. Complete endpoint-helix spectra at lengths 4,5,6 also pass the kernel and gap checks. Illegal sectors are excluded analytically, not by these diagonalizations.

[numerical, source and transport checks] At n=1,2 the calculation retains every charge excitation and every helix excitation reached by the ports, for all four length pairs and four wrapper parameter choices. The largest single-cell source dimension including the helix is 7560 per bus. Direct typed-hole transport agrees with the identity on logical pair dimensions 4 and 16 to 1.4e-16; unconverted transport is zero. Every computed first-order port compression is zero. Independent source calculations and direct inverses agree to the stated precision in PROOFS Section 8.

[numerical, length fluid] Complete spectra of the projected physical marker graphs, dimensions 104,432,1024 at n=1,2, agree with the joint momentum formula to 2.4e-14. They use C39's common second-order bond, with its tiny nonscalar part rescaled for visibility. These are projected-model spectra, not finite-coupling spectra of the full 4464-state microscopic chain.

[numerical, material limitation] The diagnostic wrapper ratios do not satisfy the severe lacunary error budget. Their resolved scalar length interaction is much larger than the diagnostic payload hopping. Tiny nonscalar differences are not certified relative measurements. The all-n conclusion uses the finite norm inequalities, not numerical extrapolation. No full microscopic many-cell ED or central-charge measurement was performed.

[proved, handoff] All writes are confined to work_C44. No network, cluster computation, subagents or edits elsewhere were used. The project STATE was outside the authorized write scope; local [STATE.md](STATE.md) records the handoff. No independent human review is claimed.

[numerical, reproduction] Use Python 3, NumPy and SciPy, with one BLAS thread. Each command completes in well under a minute on the present machine.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B audit_numeric.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B audit_numeric.py cell
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B transport_ring.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B bounds_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B consistency_checks.py

[proved, files] source_bonds.npz contains the complete computed logical bonds, including scalar terms. source_results.json reports their errors and dimensions. cell_ed_results.json, transport_ring_results.json, bounds_results.json and consistency_results.json retain the other checks. core_model.py and core_literal.py are attributed C39 snapshots, with the latter modified to expose individual positive bond contributions. wrapper_tools.py is an attributed C32 snapshot. PROGRESS.md records the numerical projector correction and the scope of every check.
