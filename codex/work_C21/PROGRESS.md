# C21 progress

2026-10-28. Completed partial result. Read RESULT.md first.
No completed fixed-alphabet compiler of the intended target.

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
All ladder runs finished. The largest symmetry block has dimension
65,536 (n=2, three cells); the n=3 two-cell block has dimension 16,384.
Resolvent errors are <=1.2e-15, ED residuals <=1.3e-14, and halving-
coupling exponents are 3.0021 to 3.0038. These are scalar-bond controls.

## Step 4, a local interface that does carry all Majoranas

[proved] A different 16-state site has a vacuum, two lanes of three
cursor syndromes, a filler, and a data qubit. The cursor visits X_i
and Y_i at consecutive lane positions. Controlled-Z hopping between
sites makes the one-cursor gauge Hamiltonian a scalar path with
onsite potentials g A tensor gamma_j, where A exchanges X and Z
syndromes. Its source resolvent stays in scalar plus Clifford-vector
form under the continued-fraction recursion. A local right-port
rotation supplies the parity-dressed eta_j. No arbitrary reflection
or history controller appears in this local module.

[proved] The complete isolated-cell code has exactly 2^n states,
product parity is exact, and all nonvacuum sectors have gap at least
Delta-g-2t. The second-order bond is exactly Gaussian and has hopping
matrix J=-2 R G, where R reverses 2n indices and G is positive definite.
Strict positivity follows from distinct first nonzero Taylor orders
of the 2n filter functions. This is a genuine escape from scalar
source resolvents, with unbounded second-order hopping rank.

[proved limitation] This particular J is not C19's target. Since
J^{-1}J^T is similar to a positive-definite matrix, a zero of its
zero-field Bloch symbol can occur only at k=0 or pi. Generic choices
are gapped. The checked n=1 instance is gapped. Independent displaced
flavours and chemical potentials have NOT been realized. Changing the
source geometry is required before attempting the final selector.

[numerical] clifford_cursor_checks.py directly resolves the local
one-cursor Hamiltonian for n=1,2,3 (dimensions 12,48,144). It finds
6,12,18 accessible energies and Gaussian bond residuals below 5e-15.
The smallest n=3 hopping singular value is only 5.75e-14, so full rank
is established analytically rather than inferred from that number.
Full n=1 ED on two and three cells (dimensions 256,4096) has cubic
canonical-matrix residuals. The n=2 two-cell full ED also finished,
dimension 65536, retaining all 16 code levels. Scalar ARPACK stalled
and was interrupted; complete-parity-sector block LOBPCG succeeded.
The final script bounds each cluster calculation to 180 seconds.
Eigenpair residuals are <=1.82e-12 and halving exponents 3.0019,3.0037.
n=2 three-cell ED for this interacting cursor has not been run.

## Final handoff

RESULT.md records task-by-task completion and limitations.
LOCAL_CURSOR.md contains the entire 16-state isolated-cell rule table,
full port completions, exclusion of every nonvacuum sector, exact parity,
continued-fraction proof, growing-pole bound, full-rank proof, and the
Bloch obstruction. MOMENTS_AND_FILTERS.md contains the product-code
moment obstruction and its delayed-Jacobi algebraic repair.
SCALAR_LADDER.md records the failed eight-state ladder and tree lemma.

The remaining problem is the target's flavour structure, not merely
the existence of local nonscalar source resolvents. For the new cursor,
J^{-1}J^T has positive real spectrum, so its zero-field Fermi points
cannot leave 0 or pi. Do not promote full hopping rank to n gapless
flavours. A revised source geometry must be tested before the selector
and full ring alphabet are assembled.

Read C22 and C23's later progress during the final audit. C22 independently
found the tree-gauge lemma and is testing a distinct charged bus, with
a related wrong-band-structure obstruction. That numerical construction
was not duplicated. C23 completes a conditional assembly/IR ledger,
not the missing spectral compiler. No external files were changed.

Final audit passed. All three Python files parse, all local Markdown
links exist, and all reported numerical maxima and tested cluster sets
agree with the saved JSON records. The exact local code/gap bound,
right-port parity dressing, Clifford continued fraction, pole count,
Gaussian contraction, and Bloch limitation were checked in the write-up.
No independent referee or full n=2 three-cell cursor ED was performed.
Project STATE.md remains unchanged because the task forbids outside
writes; this file and RESULT.md carry the handoff.
