# C39 progress

[proved, scope] Read the project memory and requested background, candidate, C32, C31, C26, C20 and infrared files. C33 and C37 progress were consulted. C38 and C40 had no progress files at the initial check. All writes remain in work_C39. The project STATE is outside the permitted write scope.

[proved, inherited distinction] C32 supplies a positive 28-letter core and an all-n crossing theorem, but no all-n center separation. C33 has a certified n=23 witness. C37 is working on modular nonidentity certificates within the original alphabet. This run will not repeat those certificates.

[heuristic, direction] Investigate an auxiliary frustration-free ruler whose local expectations depend explicitly on position. A successful construction must distinguish a genuine band phase from a removable change of logical basis. Independent chemical potentials alone do not remove centered interflavor pair transfer. The desired completion needs both a readout mechanism and a proof that it controls the band centers after the non-diagonal hopping is diagonalized.

## Two buses and a biased vacuum

[heuristic, candidate being checked] Use two independent copies of the positive folded excited bus, with stub hopping w and w v and independent charge-exchange port phases. The proposed alphabet has 46 letters, sharing the vacuum token and data letters. After large express mass, the leading hopping is intended to be proportional to Gamma(w)+i eta Gamma(w v), where each Gamma is positive. Unlike a bandwidth hierarchy, the generalized eigenvalues of this pair affect the Fermi zeros.

[proved, auxiliary mechanism] A vacuum edge (q|j>-|j+1>)(q<j|-<j+1|) gives the unique zero vector with amplitudes q^j on the connected path. Conditional vacuum-token/data-density readout compresses to q^(2j)/sum_l q^(2l) times that data occupation. This changes the code and escapes C28's uniform-history field-rank hypothesis. It preserves exact U(1). The readout is a retained weak perturbation, not part of the zero-energy core.

[heuristic, proof route] The normalized second Gram has the form diag(v^(2j)) G(wv) diag(v^(2j)). Its generalized eigenvalues should have valuations v^(4j). A second small geometric parameter in the field separates the first derivatives of the Fermi widths. Both statements need a detailed eigenvector and remainder argument. C38's congruence warning was read and is respected: no band-center conclusion is inferred from the outer diagonal scaling alone.

## First validation and a failed ED attempt

[numerical, n=2] Direct g^2 source contraction confirms that changing the port phase from pi/2 to pi multiplies the hopping matrix by i. Same-cell quadratic entries cancel to floating-point precision. This validates the proposed independent virtual-channel phase, rather than a data basis rotation.

[numerical, n=1, two cells, Q=1 sector dimension 338] The full-in-g second-order bond plus the conditional density readout agrees with physical ED with error/lambda^4 equal to 0.03046 and 0.03037 at lambda=.1,.05. First order is zero. The full source denominator is retained.

[numerical, unsuccessful run] ARPACK on n=3, three cells, Q=1 dimension 356325 hit the 180-second budget before returning a result. It is not counted as a successful ED. A block eigensolver initialized with the known complete code in that sector is being tried.

[proved, derivation written] PROOFS.md now supplies the Gram-pencil Newton valuations, the limiting eigenvectors, the two geometric parameter choices and the exhaustive pair-sum classification. The finite-parameter selection is a terminating rational/algebraic search, not a numerical extrapolation. It remains to audit the local construction and numerical implementations against that derivation.

## Validation completed

[proved, rational certificates] certify.py certifies the limiting Gaussian reference at n=1,2,3 with exact Fraction arithmetic. All 2n roots are simple and every distinct unordered pair sum, including repeated indices, is separated modulo 2pi. Certified delta4 lower bounds are greater than 5.11e-16, 1.23e-28 and 2.36e-41. These are independent checks of the target theorem, not finite-D ED points. An initial rational angle implementation became slow during pair comparisons and was stopped. Outward rounding of already certified intervals to 250-bit dyadics reduced the completed certificate runtimes to 2.5,6.5,8.2 seconds.

[numerical, literal core n=1,2,3] The 2116-dimensional physical bond agrees with the independent legal matrices of dimensions 54,160,424 after gauging, with errors below 5.9e-16. Charge commutators vanish entrywise. The weighted vacuum is positive by its rank-one edge form; a failed diagonal-dominance assertion on that vacuum block was corrected, with diagonal dominance restricted to excited blocks.

[numerical, source ED n=1,2,3] Full-in-g two-cell bonds on all 4,16,64 code states are saved in bond_n1.npz through bond_n3.npz, with complete three-cell coefficients on 8,64,512 code states. Every source eigenstate is retained. The direct g^2 pair inverse verifies convergence to the claimed Gram coefficient at D=10,100,1000.

[numerical, finite-coupling ED] All twelve requested n=1,2,3, two/three-cell checks in charges Q=0,1 completed. The largest sector has 356325 states. Block eigensolution repaired the timed-out n=3 ARPACK run in 118 seconds. The Q=1 error/lambda^4 values at lambda=.1,.05 are (0.03046,0.03037), (0.04485,0.04478), (0.02062,0.02054), (0.03371,0.03360), (0.01399,0.01392), (0.02374,0.02364), ordered by n then cell count. The largest eigenvector residual is below 2.8e-12. Complete cluster ED of all other charges and the mobile wrapper was not performed.

[proved, coordination and remaining scope] Re-read completed C33/C37/C38/C40 progress. C37's finite original-alphabet certificates and C38's hierarchy obstruction are not duplicated. C40's exact low-frequency massive-sector remainder remains assumed here. The constructed core has 46 letters. The inherited full C8/C32 wrapper count is 4464 for frozen or mobile implementations. No unconditional thermodynamic CFT claim follows from these finite results.
