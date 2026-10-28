# C49 local handoff

2026-10-28.

[proved, outcome] An explicit prethermal dressing and a finite-volume spectral theorem are complete. Input S and microscopic c=n+1 remain unproved. Read RESULT.md for the concise verdict and PROOF.md for the full derivations.

[assumed, imported inputs] C45's exact U0 H U0*=A+B+R on rooted ordered cells, A_j>=Delta Q_j, ||R||_a<=C_SW lambda^3, ||B||_a<=d=2lambda v0+2lambda^2G. C44 supplies the fixed-outer-marker, total-length-preserving physical interpretation. W=lambda^2 c_n, c_n=J0/s. These finite construction results were not independently re-derived.

[proved, constants] Lambda=d+2epsilon, epsilon=C_SW lambda^3. Spend h=64e Lambda/Delta in the support exponent at each step. m=floor[a Delta/(128e Lambda)] leaves exponent at least a/2. Each step has leakage at most half the previous value, so epsilon_m<=2^(-m)epsilon. The exact step remainder actually has bound o_j/7. The additional generators have summed norm <=4epsilon/Delta. The retained interaction changes by <=2epsilon.

[proved, new spectral bound] For any local all-code off-diagonal O, ||QOP||^2<=M||O||_0^2. With a block-diagonal reference threshold delta and M epsilon_m^2<=delta^2/8, low ordered eigenvalues and excitation gaps have error <=8M epsilon_m^2/delta. Ground energies alone have error <=M epsilon_m^2/delta without the size restriction. Delta-2Lambda is an admissible delta. These bounds do not give thermodynamic 1/M accuracy at fixed couplings.

[proved, physical qualifications] A generator on r cells has a local all-code factor, hence length conservation confines its physical support to at most (b+1)r+1 sites. This gives physical quasi-locality under C44's rooted representation. Spectral comparison is on the rooted tensor product. A fixed physical sector requires control of the retained-ground offset d_C, and competing marker numbers are not automatically compared.

[proved, sufficient but unverified criterion] A massive complementary compression and a weighted Schur error <=eta_M(K_M-E_M+v/M), eta_M->0, imply o(1/M) finite-size spectral errors. Section 7 proves this. The present upper bounds only give an error parameter of order M^2 epsilon_m^2/(delta v). Proving the needed state-dependent self-energy cancellation, logical RG and physical ensemble transfer remains open.

[numerical, checks] checks.py passes all controls, matrix dimensions <=256. BCH discrepancy <=1.24e-16, product-rotation formula discrepancy <=3.34e-16. No microscopic many-cell simulation and no measured central charge. checks.json records the results.

[proved, scope] All writes stayed in work_C49. No network or delegated agents. The no-network rule meant that Chen-Lucas-Yin and Yin-Lucas were available only through local reviews; no unpublished or remembered primary theorem was invoked. This file substitutes for the project STATE update, which is outside the permitted write scope.
