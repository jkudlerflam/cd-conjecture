# Conditional infrared source comparison and branch rank

Status: **[proved conditional sequence and matrix theorems; the required
marked GN bound is unproved for the candidate's full retarded endpoint
class]**. Neither the source-free hypotheses G1/G2 nor bounded running
couplings alone preserve the rank of infrared source vertices.

Author: the Codex GPT-6-family exact-alternative subagent, 2026-10-28.
AI-assisted derivation, not human-refereed. The source construction was
reviewed separately in [SOURCE_AUDIT.md](SOURCE_AUDIT.md).

The useful new observation is that a difference of two momentum-ray limits
removes the unknown source contact matrix. Conditional on the precise
source estimates below, its rank is exactly the number of Fermi branches.
Separate microscopic conservation of every branch is unnecessary for this
rank test. A current algebra or CFT interpretation remains additional work.

## 1. Primary-source check and the normalization issue

[read] I opened Porta, Scola and Singh,
[arXiv:2509.08665, sections 4.4 and 6.1](https://arxiv.org/html/2509.08665),
including (4.89)--(4.93), Proposition 6.3 and (6.15)--(6.26). Their source
flow distinguishes smooth diagonal density vertices from branch-changing
vertices, which can have anomalous powers. Their comparison matches source
vertices as well as action parameters. For two density insertions the
remainder is Hölder continuous, with an undetermined constant contact
term. Proposition 5.8 supplies the reference density matrix used below.
These are results for their stated microscopic class, not a theorem for
our arbitrary retarded initial action or our full set of probes.

Our N-by-N matrix below has rows indexing physical probes and columns
indexing smooth branch densities. Its off-diagonal *matrix entries* are
mixing of probe labels with density labels. They are not the anomalous
branch-changing fermion bilinears psi_b^dagger psi_c with b != c.
Distinct branch momenta remove those bilinears from smooth localization.

Write the canonically normalized branch fields as chi_b=sqrt(Z_(h,b))psi_b.
If an unnormalized source coefficient is Y_(h,ab) in front of
A_a psi_b^dagger psi_b, its normalized density vertex is

    U_(h,ab)=Y_(h,ab)/Z_(h,b).                             (1)

This identity defines our convention; a scheme which already uses
normalized fields can call its coefficient U directly. Do not divide by
Z a second time. The absolute value of a running Z need not converge and
cannot be replaced by one in this definition. The finite bare reference
normalization Zbar in section 5 is a different object.

## 2. What G1/G2 alone cannot imply

The [counterterm induction](../../direct_attack_20261028/vacancy/COUNTERTERM_INDUCTION.md)
controls source-free marginal histories and relevant counterterms. It has
no hypothesis on U. Even adding finite-scale invertibility, boundedness,
and a small bound on successive source ratios does not close this gap.

[proved counterexample to that logical inference] Set all action beta
functions to zero, so G1/G2 hold. For h<=0 set

    U_h=2^(epsilon h) I_N,    epsilon>0.

Every finite-scale matrix is invertible, ||U_h||<=1, and successive ratios
are 2^(-epsilon)I_N, as close to one as desired. Nevertheless U_h tends to
the zero matrix. A slowly rotating orthogonal U_h similarly has no limit,
although its singular values never change. These abstract recursions are
not proposed physical models; they show why a separate source cancellation
estimate is logically indispensable.

Reference Ward identities protect properly normalized smooth density
insertions. Their use must cancel the constant-history source beta term;
a bound of order epsilon per scale is insufficient. A source-dependent
comparison is needed to transfer that cancellation to the lattice action.

## 3. An explicit sufficient marked-GN hypothesis

Assume G1/G2 and their conclusions, including

    |g_h-g_*|+|nu_h| <= C epsilon 2^(theta h),  0<theta<1/4.

Let U_0 be the matched N-by-N initial source matrix. Let delta>=0 bound
the *inhomogeneous contribution to the localized one-source flow*. It is
not the raw norm of every irrelevant bare source kernel. The branch
filters have fixed momentum derivatives and higher Taylor remainders;
these need not shrink when the microscopic coupling is lowered.

Separate a finite constant C_probe controlling those bare quadratic
source remainders from delta_matching, which controls interaction-induced
source kernels after all coefficient, field-label and moment conversions.
At first order in a source, Gaussian integration of a bare quadratic
source alone leaves its two-external-fermion vertex unchanged; the only
other term is a source scalar from contraction of both fermions. Therefore
a correction to localized U involving a bare source Taylor remainder
requires at least one interaction endpoint. A marked GN proof which
retains this endpoint count may use

    delta = delta_matching + C_probe epsilon.             (1a)

This counting fact does not prove the requisite GN norm estimate. It
identifies the small quantity that estimate must use. Alternatively, a
finite-shell preprocessing step can make the source remainders small by
choosing a lower fixed initial cutoff, then choosing the interaction
last. Neither construction declares the unprocessed source remainder
small merely because lambda is small. The finite-cell matching lemma
supplies candidate bounds on delta_matching, not on C_probe. All probe
constants can be large but are fixed at the chosen finite model.

The source flow is causal and linear in the full collection of
single-source insertions. At zero external source it does not change the
action history. Write

    U_(h-1)=U_h+d_h.

The following is a sufficient *additional hypothesis*, called (S):

    ||d_h|| <= C 2^(2 theta h)
                    [epsilon max_(h<=k<=0)||U_k|| + delta]
             + C sqrt(epsilon) sum_(k=h)^0 2^[-4 theta(k-h)]
                                             ||U_k-U_h||. (S)

Norms are any fixed submultiplicative matrix norm; use the spectral norm
for the singular-value conclusions. The constants are uniform in all
cutoffs and in the source rows. The bound concerns the normalized vertices
(1), with all field normalization terms retained. An equivalent source
majorant with different positive short-memory exponents would suffice.

[proof-level reduction, not a completed extension] The expected origins
of (S) are: cancellation of the constant-history smooth-source beta by the
reference Ward identity; a short-memory gain for curvature, retarded
initial endpoints and source Taylor remainders; and a normalized marked-tree
continuity estimate for running histories. Inserting the already controlled
|g_h-g_*| and |nu_h| gives a decaying coefficient. Source-history differences
require the last term. A source-independent GN estimate or an unnormalized
vertex estimate does not establish (S).

In particular, no derivative with respect to the Fermi centers is required
by the sequence theorem below. Its constants only need to hold uniformly
at the centers selected by the existing fixed-point/compactness argument.

## 4. Conditional source-flow theorem

**[proved]** Under (S), for sufficiently small epsilon,

    sup_(h<=0) 2^(-theta h)||d_h||
                         <= C1 [epsilon ||U_0||+delta],  (2)
    U_* = lim_(h->-infinity) U_h,
    ||U_h-U_*|| <= C2 [epsilon ||U_0||+delta]2^(theta h), (3)
    ||U_*-U_0|| <= C2 [epsilon ||U_0||+delta].            (4)

Consequently, if

    C2 [epsilon ||U_0||+delta] < sigma_min(U_0),          (5)

then U_* is invertible. Arbitrary finite mixing among probe and density
labels is allowed.

**Proof.** Causality gives U_h=U_0+sum_(j=h+1)^0 d_j. Set

    D=sup_(h<=0) 2^(-theta h)||d_h||,
    S_theta=sum_(j<=0)2^(theta j)=1/(1-2^(-theta)).

For finite cutoffs, max_(h<=k<=0)||U_k||<=||U_0||+S_theta D.
Also, interchanging the finite positive sums gives

    sum_(k=h)^0 2^[-4 theta(k-h)] ||U_k-U_h||
       <= [1/(1-2^(-4theta))]
          sum_(j=h+1)^0 2^[-4theta(j-h)] ||d_j||
       <= C_theta D 2^(theta h).

Divide (S) by 2^(theta h), take the supremum, and use 2^(theta h)<=1:

    D <= C[epsilon||U_0||+delta]
             +[C epsilon S_theta+C C_theta sqrt(epsilon)]D.

For small epsilon the bracket is at most 1/2. This proves (2), uniformly
in the cutoff. Summing the absolutely convergent increments proves (3)
and (4). The smallest-singular-value perturbation inequality proves (5).
The fixed-scale causal histories and uniform bounds pass to the infinite
sequence. This proves the theorem. QED.

For the branch filters, U_0=I+E_0 after Fermi shifts and ultraviolet
matching. The explicit initial bound in BRANCH_SOURCES can make
||E_0||<1/4. Using (1a) with fixed C_probe, choosing epsilon last and making
delta_matching small then gives, for example, sigma_min(U_*)>1/2. Anomalous normalization cannot erase rank once (S)
is established; without (S), section 2 shows that it can.

## 5. A second, explicit marked-GN input for two-point comparison

The source-flow theorem is not yet a correlation-function theorem. In
addition, suppose the normalized two-source expansion matches the action
trajectory to a multichannel reference theory. Require that, after
replacing the source vertices by U_* and extracting that reference part,
the difference is a sum of shell functions E_h(p), h<=0, and a regular
ultraviolet term E_UV(p), with

    ||E_h(p)|| <= A 2^(gamma h),
    ||grad_p E_h(p)|| <= A 2^[(gamma-1)h],   0<gamma<1.     (6)

E_UV is Lipschitz near zero. The constant A can contain the fixed
bare source Taylor constants; it need not tend to zero with epsilon.
Bounds (6) are one concrete form of the
required two-source marked-GN comparison. They include source contact
vertices and cannot be inferred just by differentiating G1/G2. They follow
from an appropriate normalized short-memory/continuity theorem if that
theorem has been proved for all retained endpoints; that application is
still open here.

**[proved conditional comparison]** With intensive infinite-volume
correlators and p near zero,

    C_phys(p)=U_* C_ref(p) U_*^T+B_0+O(||p||^gamma).       (7)

**Proof.** The zero-momentum series B_0=E_UV(0)+sum_h E_h(0) converges.
The mean value theorem and (6) bound each difference by

    ||E_h(p)-E_h(0)||
       <= C A min(2^(gamma h), ||p||2^[(gamma-1)h]).

Split the sum at 2^h comparable to ||p||. Both geometric series are
bounded by C_gamma A||p||^gamma. The ultraviolet term is smaller. QED.

This is a rigorous conditional consequence of specified shell bounds,
not a claim to have proved the marked GN expansion for the candidate.
It applies to the full probe matrix; a single total-density response is
only one quadratic form of this matrix and cannot establish its rank.

## 6. A contact-free rank test

Take a reference theory with N branches, nonzero velocities v_b, finite
invertible diagonal normalization Zbar, and small forward-coupling matrix
Lambda. Set Lambda_Z=Zbar^(-1)Lambda Zbar and

    B_b(p)=1/(4 pi |v_b|) * (-i p_0+v_b p_1)/(i p_0+v_b p_1).

The reference density matrix in its scale limit has the form

    C_ref(p)=(I+B(p)Lambda_Z)^(-1) B(p) Zbar^(-2).         (8)

Here B is diagonal. Assume the displayed inverses exist; for real nonzero
Euclidean p it suffices that

    ||Lambda_Z||/(4 pi min_b |v_b|)<1.

These are reference-theory data, not an additional microscopic symmetry.
For a reference potential with nonconstant Fourier transform, (8) uses its
zero-momentum limit in the scaling covariance.

Choose rays q_0=(1,0) and q=(1,xi) with real xi != 0. Define

    D_phys(q,q_0)=lim_(rho down to 0)
                      [C_phys(rho q)-C_phys(rho q_0)].    (9)

By (7) the unknown B_0 cancels. Homogeneity of (8) gives

    D_phys=U_* D_ref U_*^T,
    D_ref=(I+B(q)Lambda_Z)^(-1)
             [B(q)-B(q_0)]
             (I+Lambda_Z B(q_0))^(-1) Zbar^(-2).        (10)

**[proved]** Every factor on the right of (10) is invertible. Hence, if
(S), (5) and (6) hold, rank(D_phys)=N.

**Proof.** For any A,B,L with the needed inverses,

    (I+AL)^(-1)A-(I+BL)^(-1)B
       =(I+AL)^(-1)(A-B)(I+LB)^(-1).

Multiplying on the left and right by I+AL and I+LB proves the identity.
Here the middle diagonal entries are

    [B(q)-B(q_0)]_b = v_b xi/[2 pi |v_b|(i+v_b xi)],      (11)

which are nonzero for every b. The other factors are invertible by the
small-coupling and Zbar hypotheses. Invertible U_* preserves rank. QED.

This argument also gives a quantitative lower bound by multiplying the
smallest singular values of the five factors in (10), including two
copies of sigma_min(U_*). It therefore excludes loss of a branch through
small source mixing or a constant contact cancellation, under the stated
hypotheses. It does not require a separate lattice continuity equation
for each filtered probe.

## 7. Relation to higher currents and the conjecture

[conditional extension] For a fixed m>=3, a marked-GN comparison with
remainder O(rho^(2-m+gamma)) transfers the reference m-point scaling
cumulant after multiplying by rho^(m-2). If the corresponding reference
Ward identity makes that cumulant vanish at nonexceptional momenta, the
same limit vanishes for the physical probes. This follows directly by
scaling the remainder. It is not an assertion of uniformity near all
coincident or exceptional momenta; distributional convergence requires
its own uniform smeared estimates.

The rank-N result is evidence of N independent smooth density directions
in the *conditional* scaling correlation theory. It is stronger than a
claim that the bare source matrix starts invertible. It does not alone
construct N chiral currents with positive levels. For the c>=N/2 argument
one still needs a local positive scaling theory and the chiral current
interpretation, or an independently established unitary conformal limit.
The full stress tensor and any extra neutral sector are not determined
by (7) or (10).

The present bottleneck is therefore precise: prove (S) and a two-source
normalized marked-GN comparison such as (6) for the retarded action, using
the reference Ward identities with the actual normalizations. G1/G2 would
supply the action trajectory needed by that proof, but do not replace it.
No c-d counterexample or complete infrared theorem is claimed here.

## 8. Conditional CFT interpretation: why two points are insufficient

The following conservative statement assumes an independently established
compact unitary CFT with a common light-cone velocity, an extremal
clustering vacuum (so its connected dimension-zero sector is removed), and
a discrete scaling spectrum of finite multiplicity below every fixed
energy. The probes are local Hermitian bosonic fields with a controlled
scaling expansion. Their connected susceptibilities are the physical
ones, with ultraviolet contact terms fixed independently of the infrared
scale. These hypotheses are not established for the candidate by the
preceding estimates.

### 8.1 Bounded response isolates dimension-one currents

[proved under these CFT hypotheses] If the susceptibility matrix stays
bounded as p=rho q tends to zero along every nonzero fixed ray, its probes
have no nonconstant dimension-less-than-one component and no scalar
primary component of dimension one.

Indeed locality gives integer spin s=h-hbar, and unitarity gives
h,hbar>=0. At total dimension Delta<1 the only allowed spin is zero.
Such a scalar's two-point Fourier transform has a positive coefficient
times |p|^(2Delta-2), which diverges. At Delta=1 a scalar with
h=hbar=1/2 contributes a positive coefficient times -log|p|. Within every
degenerate eigenspace these coefficients form a positive semidefinite
Gram matrix. Thus a nonzero scalar component is visible in at least one
Hermitian linear combination of probes, and cannot cancel against another
scalar of the same dimension. Different powers cannot cancel either.
A fixed analytic ultraviolet contact term cannot cancel these infrared
divergences. The connected identity component has already been removed.

The remaining dimension-one fields have weights (1,0) or (0,1). The
zero antiholomorphic or holomorphic weight makes the opposite derivative
a zero-norm state, so unitarity makes them chiral currents. Compactness
provides a positive gap to the next scaling dimension above one for the
retained expansion. Higher-dimensional contributions vanish in a
momentum-ray difference, apart from analytic contacts which also cancel.

Consequently the rank of the difference (9) is at most r_R+r_L, where
r_R and r_L are the dimensions of the current spaces visible in the
right- and left-chiral projections of the probes. This statement counts
visible currents; it does not yet bound c by that count.

### 8.2 The nonabelian obstruction

[proved obstruction] Current two-point rank is not the rank of an abelian
current algebra. Three SU(2) level-one currents have a nonsingular
three-by-three current two-point matrix but chiral central charge one.
Therefore the inference c_chiral >= rank(current two-point matrix) is
false. Likewise, vanishing three-point functions restricted to chosen
probes is insufficient: the span of J_x,J_y has no three-current
structure constant within that span, but [J_x,J_y] is J_z outside it.

This distinction invalidates any derivation of c>=N/2 from (10) and CFT
existence alone. No such derivation is claimed in sections 1--7.

### 8.3 Gaussian four-point functions repair that inference

Add the following hypothesis: the dimension-one probe scaling fields
have Wick Gaussian four-point functions as separated-point distributions,
and this identity passes to their holomorphic and antiholomorphic
projections. Full Gaussian scaling distributions suffice. A fixed-momentum
comparison only away from exceptional configurations has not yet supplied
this hypothesis.

[proved] The visible current spaces in each chirality are then abelian.
For Hermitian holomorphic currents J_A,J_B, write their current OPE as

    J_A(z) J_B(w) = k(A,B)/(z-w)^2
                    + i J_[A,B](w)/(z-w) + regular terms.

The positive invariant form k is the current-state Gram form. Consider
the four-point function with pairs (A,B) and (A,B), and take the two
pair-collision limits while the pair centers remain distinct. The
coefficient of (z_1-z_2)^(-1)(z_3-z_4)^(-1) is proportional to

    -k([A,B],[A,B])/(z_2-z_4)^2.

A Wick four-point function has no such double simple-pole coefficient:
the within-pair Wick contraction has double poles, while cross-pair
contractions are regular in both pair separations. Hence the bracket
has zero norm. Unitarity gives [A,B]=0 in the physical current algebra,
including components outside the initially visible span. The same
argument applies to left-chiral currents. This is why four points,
unlike three, detect the missing commutator direction.

Abelian positive-level currents admit a Sugawara stress tensor with
central charge equal to their number. The orthogonal stress sector has
nonnegative central charge in a unitary theory. Therefore

    c_R >= r_R,   c_L >= r_L,
    (c_R+c_L)/2 >= rank(D_phys)/2.                        (12)

For a nonanomalous common-velocity lattice CFT, c_R=c_L=c, so the full
rank N=46 test plus these hypotheses gives c>=23. This would suffice
against d-1=22 for the frozen d=23 candidate, *if* its CFT limit and the
required source distributional matching were established.

Different velocities and products of sectors require a separate
sector-by-sector analysis. In particular, one should not assign a single
common Euclidean spin decomposition to a multivelocity theory without
justifying it. The common-velocity theorem above is intentionally narrower
than the repository's multichannel candidate. A full Gaussian current
scaling theory with its own positive mode decomposition offers a possible
extension, but it is not proved in this note.
