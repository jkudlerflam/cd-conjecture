# C37 progress

[Proved, scope] Read the project memory and the requested candidate, spatial code, C32, C26, infrared and C20 results. C33 progress and its structured leading-coefficient code were read. No prior C37 progress existed. All writes remain in this folder. Project STATE.md is outside the authorized write scope and will not be modified.

[Proved, inherited distinction] The coefficient of lambda^2 g^2 is quadratic and has polynomial-size carrier resolvents. At finite g the exact second-order bond also has non-Gaussian terms. Positivity of its Hermitian hopping part, handled by C33, is not the nonresonance assertion.

[In progress] Audit the loop-dependent leading quadratic matrix and its pair sums. Seek a structural obstruction or a certified parameter witness, without using logical spaces of size 2^n for the main calculation.

## Pair-product certificate

[Proved by finite-field computation, audit in progress] At every n=1,...,23 the leading g^2 determinant has a squarefree unordered-pair-product polynomial modulo 65521, with i mapped to 24297. The specialization is shifted stub mass 4, D=1000, t_d=10^(-2000), normalized common field 10^(-1000), total neutral express phase pi/2, neutral loop phase pi. The largest carrier-pair matrix has dimension 2209 and the pair-product polynomial degree is 1081. Files pairs_n1.json through pairs_n23.json retain the polynomial and gcd witnesses.

[Proved strategy] This excludes algebraic identities without resolving tiny Fermi-root differences. Arithmetic lifts D -> D+65521 j and exponent lifts E,L -> E,L+65520 j retain the certificate and permit arbitrarily large express energy and arbitrarily small field/loop. They can therefore be combined with the inherited positive Gram limit. An explicit discriminant bound supplies a positive, very conservative momentum margin.

[In progress] Independently audit sparse modular solves, the leading-coefficient extraction, and the symmetric-square polynomial. Write the real-parameter schedule, denominator and velocity bounds, and the conditional higher-order transfer.

## Completed finite result

[Proved, computer-assisted] Every n=1,...,23 passes the audited pair-product certificate. C33's six existing positive-hopping certificates were reused; its exact certifier filled the other 17 sizes. verify_bounds.py checks the uniform 10^(-250) normalized crossing bound by integer inequalities at all 23 sizes. All final modular solves include exact residual checks.

[Proved] RESULT.md and PROOFS.md give the concrete same-alphabet parameter choice, denominator/discriminant margins, simple-root velocities, analytic perturbation tolerance and explicit g schedule. The exact finite-g second-order quadratic projection keeps all margins. The all-n extension is open, and the higher-order thermodynamic transfer remains conditional on the prior local-elimination/dressing input.

[Numerical, audit outcome] A finite-g Richardson subtraction failed to improve the small onsite signal, with discrepancies up to 4.5e-5 of the hopping at n=3. It is retained, not used as a certificate. A replacement direct literal g^2 coefficient calculation agrees at n=1,2,3, with hopping error at most 1.4e-12 and onsite error at most 2.3e-10 in hopping units. The latter includes all 28224 excited-pair states at n=3.

[Proved, handoff] Project STATE is outside the write scope. Local STATE.md records the outcome and file map. No network or subagents were used.

[Proved, final audit] Replaced the bounded exact binary64 dot-product shortcut by NumPy int64 dot products with explicit no-overflow bounds. Recomputed all 23 modular certificates and the independent exact audits. Every polynomial and gcd agrees with the earlier result. The final n=23 run took about 33 seconds. Python syntax, certificate degrees, all gcds and the 23 crossing-bound records were checked.

[Proved/assumed distinction, final refinement] The root perturbation proof now also bounds the corrected minimum speed. Higher-order bare corrections and the retained interaction self-energy have separate assumed constants. The conversion must be reduced for a specified dressing constant before lambda is chosen. Lowering lambda alone does not suppress normalized interaction dressing. RESULT.md, PROOFS.md and local STATE.md reflect this distinction.
