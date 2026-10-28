# C24 review progress

[proved, review scope] This review is confined to work_C24. No source file outside this folder is changed and no network is used. The project memory was read first.

[proved, completed reading] Read STATUS latest continuation, research_continuation overview and infrared audit, the sparse port and circuit construction, exact rotations, scheduler and initialized-core interfaces, and the relevant C8/C13/C20 interface statements.

[proved, preliminary algebra] The Kraus normalization and parity-dressed right-channel attenuation follow from CAR. The incidence circulation argument excludes nonzero flavor transfer below six fields.

[assumed, under audit] The complete compiler implementation and the uniform local massive-sector elimination remain to be assessed independently. In particular the large-b physical-coordinate harmonics and the fixed-coupling limit need care.

[proved, next work] Run local copies of supplied scripts, add independent port, graph, clock and Gaussian checks, finish wrapper and infrared review, and write RESULT.md with severity classifications.

## Independent checks completed

[numerical] All ten requested contribution scripts and six C8 wrapper/transport scripts pass from local copies. The copied scripts write only inside reproduction/.

[numerical] independent_checks.py constructs sparse joint isometries for N=1,3,4,6, checks both ports and every occupied charge entry, and finds maximum entry error 6.7e-16. It independently generates projective planes p=2,3,5,7 and checks uniqueness of all one-edge and two-edge charge multisets.

[numerical] Explicit neutral, nonoscillating six-field cycles occur for the chosen arithmetic progression. There are 18 chiral assignments among 28 six-cycles at p=2 and 220 among 234 cycles at p=3, with one orientation per cycle. Their free dimension is three; the stated full metric bound makes them irrelevant.

[numerical] Sixty random acceptance-gap tests pass in dimensions 8,15,25. A 144-dimensional history toy has four accepted states with all-time FAIL, versus eight with terminal-only FAIL. This tests the acceptance mechanism, not the complete compiler.

[proved within the stated Gaussian dictionary] Physical detunings include division by b+nu. This can be absorbed into a fixed-p, fixed-b tolerance, but should be explicit. The reference metric satisfies K A K=B with K=T^{-1}T^{-T}.

[assumed, still under review] No complete sparse-core transition table has been independently executed here. The contribution source tree is changing during this review; current STATUS also mentions later controller prototypes, explicitly not a complete autonomous high-channel core. The review concerns the requested candidate and its displayed interfaces.

## Final result

[heuristic, verdict] Review complete. RESULT.md recommends retaining the status of conditional physics-level candidate, not an established refutation. No fatal flaw was identified. The two serious verification gaps are the complete autonomous sparse controller and the uniform local low-frequency elimination for moving cells. The absence of an optimized constant C or a simulation at a violating prime is not treated as fatal.

[proved, additional checks] The six-field resonant example has a positive Gaussian K with dimension one outside the candidate neighborhood. The full kappa > 2/3 bound, rather than positivity alone, is therefore necessary to the stated protection. The explicit half-density chemical-potential formula checks analytically. C13 already accounts for the cell-to-physical momentum conversion.

[numerical] All 16 copied contribution/C8 runs and the separate C20 algebra run passed. independent_checks.py passed after adding the positive-stiffness countercheck. Largest independent port residual is 6.7e-16 at N <= 6. No full autonomous core or finite-coupling central charge was simulated.

[proved, deliverables] RESULT.md contains the derivations, assumptions, severity table and acceptance conditions. independent_checks.py, independent_results.json and independent_checks.stdout.txt contain the new checks. reproduction/ contains local script copies and their captured results. source_manifest.json records the source hashes at the final reading phase. All work remains inside this folder.
