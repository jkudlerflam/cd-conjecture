# C37 progress

[Proved, scope] Read the project memory and the requested candidate, spatial code, C32, C26, infrared and C20 results. C33 progress and its structured leading-coefficient code were read. No prior C37 progress existed. All writes remain in this folder. Project STATE.md is outside the authorized write scope and will not be modified.

[Proved, inherited distinction] The coefficient of lambda^2 g^2 is quadratic and has polynomial-size carrier resolvents. At finite g the exact second-order bond also has non-Gaussian terms. Positivity of its Hermitian hopping part, handled by C33, is not the nonresonance assertion.

[In progress] Audit the loop-dependent leading quadratic matrix and its pair sums. Seek a structural obstruction or a certified parameter witness, without using logical spaces of size 2^n for the main calculation.

## Pair-product certificate

[Proved by finite-field computation, audit in progress] At every n=1,...,23 the leading g^2 determinant has a squarefree unordered-pair-product polynomial modulo 65521, with i mapped to 24297. The specialization is shifted stub mass 4, D=1000, t_d=10^(-2000), normalized common field 10^(-1000), total neutral express phase pi/2, neutral loop phase pi. The largest carrier-pair matrix has dimension 2209 and the pair-product polynomial degree is 1081. Files pairs_n1.json through pairs_n23.json retain the polynomial and gcd witnesses.

[Proved strategy] This excludes algebraic identities without resolving tiny Fermi-root differences. Arithmetic lifts D -> D+65521 j and exponent lifts E,L -> E,L+65520 j retain the certificate and permit arbitrarily large express energy and arbitrarily small field/loop. They can therefore be combined with the inherited positive Gram limit. An explicit discriminant bound supplies a positive, very conservative momentum margin.

[In progress] Independently audit sparse modular solves, the leading-coefficient extraction, and the symmetric-square polynomial. Write the real-parameter schedule, denominator and velocity bounds, and the conditional higher-order transfer.
