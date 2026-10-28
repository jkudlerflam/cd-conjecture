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
