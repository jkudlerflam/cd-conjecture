# C21 progress

2026-10-28. Initial investigation. No completed fixed-alphabet compiler.

Read the route-A background, project STATE, repository STATUS, C19 RESULT,
C18's interface construction, C8's current handoff, and the collaborator
READMEs in graded_route_a, exact_route, torus_bottleneck and
research_continuation. No prior C21 files existed. C22 was empty; C23's
progress file appeared during these reads and will be consulted before
overlapping work. All writes stay in this folder, as the task requires.

Inherited facts, not rederived: C19's graded Jacobi weights and reflection
cell are valid finite matrices; C18 already escapes the bare Bell-check
obstruction with syndrome-conditioned work energies, producing a diagonal
logical interaction. C8 already supplies the initialized 15-letter clock.
None supplies the local spectral compiler requested here. The new
first-order construction has growing alphabet and does not fill this gap.

Active work: test boundary spectral moments against nearest-neighbor
locality; construct a repeated-term excitation ladder with all poles
accessible from its endpoint; determine whether local unitary decorations
can make its source measure logical-operator-valued. Do not identify
having enough energies with having C19's operator weights.

Initial algebraic finding, proof being written: for a zero-energy code G,
H^r F G = ad_H^r(F) G. Boundary moments therefore see only a distance-r
neighborhood. C19's first Jacobi moment has nonzero coefficients for every
flavour. This excludes its literal product-register layout with distant
logical bits, but not an entangled history encoding or a change of logical
basis. A no-go statement broader than this is not justified.

## Step 1, boundary moments and a stronger product-code obstruction

[proved, being incorporated into RESULT] Disjoint boundary moment
operators commute after compression to a product subspace, even after
an arbitrary global change of logical basis. C19's left and right first
moments do not commute. For its stated target, the commutator contains
the nonzero coefficient sum_a b_a^2 t_a cos(theta_a), times logical
parity, where b_a=sum_e e X_ea. Thus its original source measures cannot
be realized on a product-subspace code in a NN cell of length >=4.
This excludes a vacuum cursor with untouched product-register code,
not entangled history codes and not route A as a whole.

## Step 2, delayed filters

[proved] Replace x(1-x)Q_j by x(1-x)^s P_j^(2s,1)(2x-1).
There are k=n+s poles, the same Cauchy orthonormality, and all source
moments through s-1 are scalar. The first nonscalar moment is s.
The inherited C19 construction then realizes the same hopping and
fields, with unchanged parity and exact abstract-code degeneracy.
The field filter d also acquires the (1-x)^s factor. Its positive
overlap is [B(2,s)-B(k+2,s)]/k. This is an algebraic repair, not a
spatial realization. For simultaneous delay s, k>=n+s is necessary.

filter_checks.py passed 20 exact rational cases, n=1,2,3,5,8 and
s=1,2,4,8. Direct positive-weight contractions at n=1,2,3 and s=1,2,3
match the intended second-order bond to <5.6e-16 absolute error.
The smallest logical scale is 3.57e-10, so the worst relative error
is 9.7e-7. Exact rational identities supply the algebraic check.
Parity commutators are zero in every numerical weight test.

## Step 3, an explicit local pole ladder and why it fails

[proved] An eight-state site, a four-state vacuum/triplet register
times a data qubit, gives a positive repeated NN hopping ladder.
A cell of m=n+1 sites has m distinct endpoint-accessible poles
Delta-2t cos(pi e/(m+1)), with sine-squared weights. All unwanted
excitation-number sectors have an analytic positive lower bound.
The complete vacuum code is exact, with 2^n states after one data pin.

[proved] Arbitrary unitary decorations of hopping on a path are removed
by a position-dependent unitary. At either endpoint the full return
source measure remains scalar. On a rooted tree the same statement
holds. This rules out a pure unitary-decorated ladder as a compiler,
even though it has enough poles and allows nontrivial data operations.

ladder_checks.py checks the decorated source measures, microscopic
bond parity and positivity, then exact finite-coupling spectra for the
undecorated control on two and three cells. These are deliberately
tests of a failed interface, not a claim of reproducing C19's bond.
The n=2 three-cell and n=3 two-cell runs were started with the smaller
cases already passing. Full results will be recorded after completion.
