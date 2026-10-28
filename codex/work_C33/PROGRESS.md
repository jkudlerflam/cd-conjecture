# C33 progress

[Proved, scope] Read project STATE, route-A background, STATUS, the spatial README and ycell/logical_chain/hermpart code, and C32 PROGRESS. C32 treats a mobile positive-core modification. This run treats frozen-cell folded-bus hopping and momentum margins. All writes stay here. No prior C33 progress file existed. Project STATE is outside the permitted write scope.

[Proved, distinction] The exact second-order interaction is not quadratic at finite conversion g. The one-particle hopping and the Hilbert-Schmidt quadratic projection used in the supplied scripts need not coincide. The coefficient of lambda^2 g^2 is quadratic and has a polynomial-size carrier-resolvent calculation. Exact finite-g one-particle extraction will be reported separately.

[In progress] Derive the scalar folded-tree filters and a quantitative finite-express positivity theorem. Build a stable structured calculation, then examine loop-induced centers and quartic momentum mismatches. A positive Gram limit alone proves neither distinct centers nor Sidon margins.

## Tree transfer identity

[Proved, derivation being written] With equal real hopping magnitudes, no loop, stub diagonal 4 after subtracting the vacuum energy, and express/marker diagonal 4+D, let C be the (n+1)-vertex express block. Its endpoint transfer is beta(z)=1/det(C+z). The exact scalar filters satisfy Laplace(f_M,i)=beta(z)^2 Laplace(f_0,i). After division by beta(0)^2 the convolution kernel is a convolution of 2(n+1) exponential densities with rates equal to the eigenvalues of C. Its mean is at most 2(n+1)/(D+2).

[Proved, consequence] A derivative bound plus a small-time Gram lower bound gives an explicit finite D(n) making the leading g^2 hopping definite for every n. A resolvent-series remainder then supplies an explicit nonzero g(n) for the exact finite-g one-particle hopping and for the quadratic projection. This addresses positivity and crossing, not distinct centers. Details and constants are being checked before the result is finalized.

## Structured calculations and an integer certificate

[Proved, algorithm] At order g^2, F_ij=2 U_i^T(A tensor I+I tensor A)^(-1)V_j, where U and V restrict a single positive source resolvent to a data row or column. The inverse has only (2n+1)^2 states. The full finite-g one-particle T uses only cell-charge pairs (0,1) and (1,0), with polynomial-size fixed-charge sectors.

[Numerical, n=1,2,3] The finite-g reduction agrees with the literal ycell calculation at D=2, shifted stub diagonal a=3, g=0.3 to relative errors 3.9e-12, 5.4e-11, 7.1e-10. At the much smaller-signal D=10, a=4, g=0.05 point the original spectral extraction lost accuracy, reaching a relative discrepancy 1.6e-4 at n=3. That failed comparison is a numerical limitation, not a successful validation. The resolved point is used for the cross-check.

[Numerical, n through 12] Full finite-g one-particle resolvents were evaluated without conversion truncation, using a convergent Neumann inverse with an exact-arithmetic truncation bound. The n=12 sector has 641875 states and took 32 seconds. The n=16 attempt hit its 160-second budget and was not completed. Fixed D=10 at a=4 becomes indefinite at n=3, already in the leading coefficient. This does not obstruct an increasing D schedule.

[Proved, computer-assisted integer certificate] For a=4, D=1000, the normalized leading Hermitian matrix is positive at n=2,4,8,12,16,23. Floor-rounded Jacobi iteration gives rational entry enclosures; fraction-free integer determinants certify positivity and an inverse-trace eigenvalue bound. At n=23 the smallest eigenvalue lies between 10^(-63.615) and 10^(-62.253). The entry-enclosure row-norm error is below 10^(-75.929). No floating-point sign test is used. The n=23 certificate took 15 seconds. The unnormalized hopping has an additional beta(0)^2 with log10 beta=-72.0416, endpoint visibility and 2g^2.

[In progress] Check larger integer certificates, write the all-n analytic schedule and finite-g remainder proof, and test the loop using a gauge in which tree matrices are real. Distinct centers and Sidon-compatible margins remain unproved for general n.

## Loop precision audit

[Proved, computer-assisted] The tree integer certificate also passes at n=32, D=1000, with 384 bits. Its smallest normalized Hermitian eigenvalue is between 10^(-88.644) and 10^(-87.139). The n=48 certificate hit its 175-second cap and is not claimed.

[Numerical, n=2 through 10] The direct imaginary same-cell formula in fast_loop.py agrees with independent source-block enumeration. The loop spectra have positive pair-sum mismatch estimates through n=10. These are leading g^2 calculations, not certified root intervals.

[Numerical, n=23, linear loop response only] Integer source arithmetic and fraction-free whitening give 46 unit-circle roots with delta2 about 0.0142 and delta4 about 2.0e-7 for a loop response of metric norm 8.5e-4. This is not yet a finite-loop result.

[Rejected numerical calculation] A full finite-loop solve at 256 fixed-point bits had an error enclosure much larger than the narrow bands. The output full_loop_256bit_rejected.json is preserved as a failed diagnostic. The cause is multiplication of tiny imaginary errors by express-sector source amplitudes proportional to beta^(-1), followed by whitening. The finite-loop solve is being repeated at 512 bits. No physical inference is drawn from the rejected roots.

## Finite n=23 band conditions certified

[Proved, computer-assisted] The 512-bit finite-loop calculation passes. A phase choice within the original 23 letters makes the scalar loop common to all three species after gauging: physical loop phases (pi,pi/2,pi/2), neutral express phase pi/(2n), port phase pi/2. At n=23 use a=4, D=1000, beta=1/det C, t_d=beta*2^(-136). Exact entry enclosures, fraction-free conditioning, rational endpoint inertia and rational pair-phase chord tests prove 46 simple roots, delta2>0.0140193, and delta4>3.1558e-7. The pair enumeration includes repeated Fermi points.

[Proved, finite conversion] A rational perturbation budget extends this certificate to the exact quadratic projection of H_2 at g=2^(-365), with common field epsilon=(2 alpha^4 beta^2 g^2)*2^(-214). The Hermitian hopping remains positive. This is not merely a g^2 asymptotic statement. T_1 and T_q are still distinguished; the full inverse was not explicitly evaluated at n=23.

[Numerical validation, n=2,3] Independent source-block enumeration agrees with the integer finite-loop g^2 matrices to relative errors below 1.1e-15. The full finite-g tree inverse at n=8, D=1000, g=.05 has positive numerical Hermitian margin. A simultaneous n=12,D=1000 run hit its 160-second cap; it is not claimed. Exact full-in-g inverses reached n=12 at D=10, and n=8 at D=1000. Leading scalar matrices were evaluated numerically through n=64, and tree positivity was certified through n=32.

[Open] The all-n distinct-center and Sidon-margin theorem is not established. No n-dependent failure of the adjustable-parameter 23-letter family was found. No alphabet enlargement was needed for the certified n=23 target. The interacting infrared conclusion and global/mobile length-sector questions remain outside this finite-band result. PROOF.md contains the derivations and certificate explanation; RESULT.md and the local STATE handoff are being finalized.
