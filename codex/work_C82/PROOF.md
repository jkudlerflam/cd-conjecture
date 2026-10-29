# C82 derivations

The construction and estimates below repair F5 with d=17, give a local
penalty-elimination step for F4, and prove frozen marker-sector selection
at nonzero port coupling for F3. They do not prove finite-J infrared
matching. Claim tags distinguish these statements throughout.

## 1. The grammar carried by a vacancy

[Definition] Retain C75's 15 material letters. Add two vacancies V_0,V_1.
For a letter a, let o(a) be its outgoing grammar state and I(a) its set
of acceptable incoming states.

| Letter a | o(a) | I(a) |
| :-- | :-- | :-- |
| L_0,L_1 | 0 | {0} |
| R_0,R_1 | 1 | {0,1} |
| D_(s,a),X_(s,a) | 1 | {0} |
| M_0 | 0 | {0,1} |
| M_A,M_B | 0 | {0} |
| V_q | q | {q} |

A physical pair ab is allowed if o(a) belongs to I(b). Put

    p_i = 1 - sum_(ab allowed) |ab><ab|_(i,i+1),
    N_bad = sum_i p_i,                 P = 1_(N_bad=0).

[Proved] Restricted to material letters this is C75's grammar. An L or
marker can precede every material letter. An R,D,X can precede only an
R or M_0. Both vacancies have zero internal charge and zero occupied
number. The on-site dimension is 15+2=17.

[Definition] Write h_15=K p_mat+b_15 for the entire C75 occupied bond,
including its half-site onsite terms and port. Define

    h_17 = JK p_ext + J embed(b_15)
           - t sum_(a material) sum_(q in I(a))
               (|a,V_o(a)><V_q,a| + h.c.)
           - (mu/2)(n tensor I + I tensor n).                 (1)

Here embed is zero unless both sites are occupied. The material onsite
terms are not moved onto vacancy bonds. Equation (1) is implemented in
bond17.py. Its occupied block is J h_15. The same JK coefficient penalizes
every extended mismatch, including V_0 V_1 and V_1 V_0.

[Proved] This is a one-site translation-invariant nearest-neighbor
Hermitian bond. Each hopping term moves one occupied letter without
changing its internal charge or marker status. All material terms have
C75's conserved charges. Thus occupied number, marker number, and the
internal U(1) commute with H. Vacancy color is not a new conserved number.

## 2. Why the hidden wall now costs K

[Proved] Fix a physical occupancy pattern with at least one material
site. Consider consecutive material letters a,b with r vacancies between
them. If all r+1 intervening physical bonds are allowed, each vacancy
has color o(a), and o(a) belongs to I(b). Conversely, if ab is an allowed
material pair, this unique coloring makes every intervening bond allowed.
The argument is an induction along the vacancy string. The first bond
sets its first color, each vacancy-vacancy bond copies the color, and the
last bond tests compatibility with b.

Therefore deleting vacancies maps the extended legal basis bijectively
onto legal material words with uncolored occupancy patterns. More
generally, an illegal material pair forces at least one physical mismatch
in its own intervening interval. These intervals are disjoint. As diagonal
operators in any nonempty occupied-number sector,

    N_material_bad <= N_bad.                                  (2)

This applies to arbitrarily long vacancy strings and every excited
letter. It is stronger than a test restricted to the vacuum code.

[Proved] In the referee's example R V_1 R R, the move on the final RR
creates R V_1 L R. The V_1 L bond now costs JK. Changing the vacancy to
V_0 cannot evade the penalty because R V_0 costs JK. For several
vacancies, changing color along the string creates a penalized color
boundary. The K-independent F5 mechanism is removed.

[Proved] A projected hopping event is the ordinary swap of a vacancy and
an occupied letter. For V_q a -> a V_o(a), the incoming grammar state is
q on both sides of the move and the outgoing state is o(a) on both sides.
In the reverse direction the unconstrained bond includes every q in
I(a). The predecessor's outgoing state selects exactly one of these in
P. Each projected swap consequently has amplitude -t, with no additional
multiplicity or normalization factor.

There can be off-grammar reverse hops at finite K. They have been kept in
(1). Exact invariance is claimed for PHP, not for the full finite-K bond.
The extended mismatch gap controls these processes as well as C75's
original off-grammar moves.

[Proved] The projected occupied moves also match C75's projected bond.
Every occupied move keeps the right letter's outgoing state unchanged.
The only changes in an external compatibility condition are tested by
the preceding extended bond. Onsite conversion preserves both grammar
states. Hence the coloring isometry U_col satisfies

    U_col^dagger P H P U_col
      = H_vacancy_hop + J sum_j 1_(g_j=0) (P_mat b_15 P_mat)_j
        - mu N_occ.                                         (3)

The formula is interpreted on the legal material Hilbert space. On a
finite ring it retains the physical root position and cyclic quotient,
as in the parent's squeezing construction. Coloring is equivariant under
cyclic rotation and adds no twist or extra oscillator. The all-vacancy
sector has two constant legal colorings. Equation (3) concerns nonzero
occupancy, in particular the candidate's half-density sector.

[Proved, restricted obstruction] A one-state diagonal pair penalty that
vanishes on every legal material word at every vacancy placement cannot
detect F5. Both L V and R V occur in legal embeddings, and both V L and
V R occur. Thus neither bond of R V L can be penalized while retaining
all these embeddings. The same applies to V V for longer vacancy strings.
This rules out the direct pair-penalty repair with the unchanged
one-vacancy charge space. It is not a no-go theorem for every possible
d=16 model or for a different nonperturbative charge phase.

## 3. Finite-size leakage bound and exact diagonalization

[Proved, finite volume] Let H=kappa N_bad+B and P=ker N_bad. For the
frozen chain kappa=K. For (1), kappa=JK. If kappa>2||B||, a normalized
ground vector psi with energy E satisfies E<=||B|| by a legal trial state,
and

    QHQ-E >= (kappa-2||B||)Q,
    Q psi = -(QHQ-E)^(-1) QBP psi,
    ||Q psi|| <= ||B||/(kappa-2||B||).                         (4)

Equation (2) then bounds the weight on illegal material words by the
square of (4). At fixed volume this proves inverse-square suppression.
The bound is not uniform in volume because ||B|| is extensive. Section 4
uses an interaction norm instead. In an infinite system, a small defect
density must not be confused with a large overlap with the globally
legal subspace.

[Numerical] check_vacancy.py exhausts the invariant zero-excitation
closure, with every reachable L/R orientation and vacancy color. Bits
are spectators and lambda=field=0. The diagnostic sets t=J=1. No states
are removed because of their penalty. On two n=3 cells plus one vacancy,
the uncolored and repaired dimensions are 2304 and 4608.

| K | Uncolored illegal material weight | Repaired illegal material weight |
| --: | --: | --: |
| 10 | 0.0465931 | 0.0100052 |
| 100 | 0.0267832 | 0.000107551 |
| 1000 | 0.0257054 | 0.00000108213 |
| 10000 | 0.0256079 | 0.0000000108279 |

For the repaired model K^2 times the weight approaches 1.083. A second
test, one n=2 cell and two vacancies, has dimension 480 and weight
8.60944e-10 at K=10000, with K^2 times the weight approaching 0.0861.
The projected hopping and material matrices in this test agree entry by
entry with their legal uncolored counterparts.

[Numerical] check_full.py uses the literal 17-state bond, not the reduced
vacuum rules. It checks Hermiticity and all three charges, and exhausts
all 17^4 minus 2^4 nonempty four-site words for (2). Its four-site ring
has N_occ=3, one marker, total internal charge 2, and an invariant closure
of dimension 1344. Parameters are D=4, n=3 in the phase/Delta convention,
Delta=4-2cos(pi/5), g=.07, t_d=.003, lambda=.1, mu_sel=field=mu=0,
and J=t=1. These are diagnostic parameters, not a selected-cell or
central-charge measurement. The illegal material weights at K=100,1000,
10000 are 1.48541e-5, 1.52004e-7, 1.52362e-9. All reported eigenpair
residuals across the leakage tests are below 2.2e-11.

[Numerical] check_parent.py compares the occupied block with C75's
original bond builder, after relabeling its letters. At two generic
resolved points with n=2 and n=16, the maximum entry difference is 0.0.
It reads the original source in memory and creates no external caches.

## 4. F4, local elimination at the grammar scale

[Definition] On an arbitrary finite ring write

    H = kappa N + B,    N=sum_i p_i,    B=sum_X B_X,
    nu = sup_i sum_(X contains i) ||B_X||.                    (5)

The p_i are commuting diagonal projectors on bonds. B initially has
range at most two sites. All constants below depend on a fixed order,
range, and chosen exponential weight, but not on the ring length.
The same construction applies to the original frozen 15-letter chain
and to the repaired 17-letter chain. It does not apply to a penalty on
occupied pairs alone as a test for material legality, because (2) would
then fail.

[Proved] For a term A supported on an interval X, define

    A_m = (1/(2pi)) int_0^(2pi) exp(i theta N) A exp(-i theta N)
                                      exp(-i m theta) dtheta.

Only the projectors meeting X enter this expression. All others commute
with A and with those projectors. Consequently A_m is supported on X
with at most one extra site at each end, ||A_m||<=||A||, and A_m=0
unless |m|<=|X|+1. Also

    [N,A_m]=m A_m,             A_m^dagger=(A^dagger)_(-m).    (6)

For a two-site microscopic term, at most three projectors enter, so the
first decomposition has |m|<=3 and support at most four sites. Let
B_0 be the m=0 part and B_off=sum_(m!=0) B_m. The anti-Hermitian generator

    S_1 = sum_(m!=0) B_m/(kappa m)                          (7)

satisfies [S_1,kappa N]=-B_off. Its local strength is O(nu/kappa).
There is no total-chain norm or small finite-size gap in (7).

[Proved] The Baker-Campbell-Hausdorff expansion gives, through order
nu^2/kappa,

    exp(S_1) H exp(-S_1)
      = kappa N+B_0+F_2+O_local(nu^3/kappa^2),
    F_2 = [S_1,B_0]+(1/2)[S_1,B_off].                       (8)

Take D_2=(F_2)_0 and

    S_2 = sum_(m!=0) (F_2)_m/(kappa m).

Using S=S_1+S_2 cancels the nonzero-m part of F_2. Thus

    exp(S) H exp(-S) = kappa N+B_0+D_2+R_2,
    [N,B_0+D_2]=0,
    ||D_2||_local <= C_2 nu^2/kappa,
    ||R_2||_local,exponential <= C'_2 nu^3/kappa^2.           (9)

On P the term [S_1,B_0] has zero compression. Since B_m P=0 for m<0,
the second-order effective Hamiltonian is

    P(B_0+D_2)P
       = PBP - sum_(m>0) PB_(-m) B_m P/(kappa m).            (10)

The negative sign is the usual lowering by virtual illegal states.
Although (10) is written using global P and N, its interaction before
compression is a sum of connected finite-range commutators in (8).
Disjoint terms commute. It does not introduce an extensive-range
interaction into the material theory.

[Proved] Here is the volume-independent remainder argument. Resolve all
operators into interval interactions. A nonzero nested commutator uses
a connected cluster of overlapping intervals. At order k its support
length grows at most linearly in k. Summing the norm of the next
commutator over intervals intersecting a current interval I costs at
most 2|I| times the interaction norm of the next generator. Iterating
gives a bound by k! C^k times the product of strengths, with C depending
only on the fixed generator range and exponential weight. The factor
1/k! in the exponential series cancels this factorial. The series for
the exact conjugation therefore converges in a sufficiently weak
exponential interaction norm when C||S||_local is small. Its tail
starting with three powers of B is bounded as in (9). Terms involving
kappa N have one additional commutator with S and hence the same
nu^3/kappa^2 scaling. No factor of the volume occurs in this interaction
estimate. The operator norm of R_2 on a ring can still be O(L).

[Proved] This construction iterates to any fixed order r. For an explicit
recursion, set epsilon=nu/kappa, H/kappa=N+epsilon V, and construct
T(epsilon)=sum_(j=1)^r epsilon^j T_j. At order p, calculate the known
coefficient F_p of exp(ad T_<p)(N+epsilon V), put D_p=(F_p)_0, and set
T_p=sum_(m!=0)(F_p)_m/m. The new term [T_p,N] removes the off-diagonal
coefficient. By induction all denominators are nonzero integers, all
order-p terms are supported on intervals of length O(p), and their
norms are bounded by constants depending on p alone. The preceding
connected-commutator estimate for the finite polynomial T gives the
exact finite-volume identity

    U_r H U_r^dagger = kappa N + D^(r) + R_r,
    [D^(r),N]=0,
    ||R_r||_local,exponential <= C_r nu (nu/kappa)^r.         (11)

Constants may grow rapidly with r. No convergence as r tends to infinity
and no exact all-order block diagonalization are asserted.

[Proved, scope of dressing] For a normalized phi in P, the dressed state
U_r^dagger phi has local mismatch expectation O((nu/kappa)^2). The zeroth
order is zero, and P[S,p_i]P=0 because p_i P=P p_i=0. All remaining terms
begin with two generators, with the same connected support bounds as
above. Summing gives an O(L(nu/kappa)^2) bound, not an L-independent
bound on the weight outside P. Identifying the true infrared states
with this dressed sector still requires stability under R_r.

[Numerical] check_sw.py tests (8)-(9) on the complete 96-dimensional
colored closure of one n=2 cell and one vacancy. At kappa=20,40,80,160,
kappa^2||R_2|| equals 7.41884,7.37390,7.35131,7.33998. The unitary residual
is below 4.5e-16. This checks signs and scaling, not the volume-uniform
claim, whose argument is the connected-support estimate.

[Scope] A many-body band at N=0 need not be separated from the entire
N>=1 spectrum on a long chain. Equation (11) avoids that false premise.
After F3 selects the marker pattern in the frozen model, the finite-K
cell ground space can also be used directly as an onsite code for the
next auxiliary-gap elimination. Either organization retains the K-layer
corrections. A tiny local remainder alone is not a theorem about a
gapless infrared phase.

## 5. F3, the boundary-energy argument with full cells

[Assumed, inherited proved input] C75 Section 8 establishes, on the
complete cell Hilbert space with l-1 data sites,

    E_l^infinity - (l-2)_+/(K-4) <= E_l(K) <= E_l^infinity,
    E_l^infinity = mu_sel-2cos(pi/(l+1)).                     (12)

The bound includes every illegal word and arbitrary excited-token number.
Its hypotheses include Delta>2+g and |t_d|<=1. A nonnegative data field
does not weaken the lower bound and vanishes in the empty-data trial
state for the upper bound. The marker-free ring bound is -L/(K-4).
These are imported inputs, independently checked by referee_C75. We do
not assume that individual bad-pair projectors commute with H.

[Proved] Write b for the chosen length, e_b=E_b^infinity/b, and suppose

    delta = inf_(l!=b) (E_l^infinity/l-e_b)>0,
    gamma = delta-1/(K-4)>0.                                (13)

For l!=b, (12) implies E_l(K)>=l e_b+l gamma. This is a linear-in-length
margin, so arbitrarily long wrong cells do not require a uniform walker
gap. The uniform lower cost per wrong cell gamma is a weaker consequence.

[Proved] Marker positions are exact conserved labels in the frozen
15-letter chain. At lambda=0 the full cell Hilbert spaces are tensor
factors even off grammar. Every cross-cell switch is controlled by the
identity on the preceding marker, so it belongs to the cell on its
right. All remaining coupling between cells is the port. Its norm is

    B_port = sqrt(2)|lambda|.                                (14)

For each data bit the two excited outputs and their common input form a
three-vertex star. Its nonzero singular value is sqrt(2). Disjoint stars
do not add their norms.

[Proved] Let A_m be the exact ground energy of an open chain of m full
b-cells, and let u=lim A_m/m. Cutting or inserting a single port changes
the ground energy by at most B_port. Thus

    q A_m-(q-1)B_port <= A_(qm) <= q A_m+(q-1)B_port.

The same cut bound for arbitrary concatenations proves existence of u
by almost subadditivity. Dividing the displayed inequalities by q and
taking q to infinity gives

    |A_m-m*u|<=B_port.                                       (15)

For the selected periodic-sector energy S_m,

    S_m <= m*u+2B_port,       u<=E_b(K)<=E_b^infinity.         (16)

The latter variational bound uses a product of exact finite-K vacuum
cell ground states. Every port has zero expectation because it changes
the excited-token numbers of the two incident cells.

[Proved] Take a compatible frozen ring of length L=m b with k>=1
wrong-length cells. Remove every port touching a wrong cell, at most
2k bonds. At most k selected open runs remain. Bound each run with (15),
and each wrong cell using E_l(K)>=l u/b+l gamma. Then

    E_wrong >= (u/b)L + gamma sum_wrong l - 3k B_port,
    E_wrong-S_m >= gamma sum_wrong l -3k B_port-2B_port
                  >= k(gamma-5B_port).                      (17)

The proof also covers a ring containing only wrong cells. If it has no
markers, it has no ports. Its energy is at least -L/(K-4), whereas the
selected product trial energy is L e_b. Hence it is excluded whenever
1/(K-4)<-e_b. This proves selection at nonzero lambda on every compatible
ring, including all off-grammar components in every marker sector.

[Proved] For C75's b=16 parameters, its inherited rational bounds are
delta>2^-16 and e_b<-2^-8. At K>=2^1152 one can use gamma>2^-17.
The original lambda=2^-640 satisfies gamma>5sqrt(2)|lambda| by a very
large margin. F3 is therefore closed for the original frozen candidate.

[Scope] The selected ground state is not literally in P at finite K.
There is no separate conserved illegal-word sector to exclude. Equation
(17) selects full marker sectors and allows their virtual illegal-word
admixture. This is the needed replacement for the parent's illegal-pair
cutting argument. In the mobile model the material marker separations
are conserved, but contact-dependent energies couple them to the charge
fluid. Transferring (17) to finite J requires the parent's dressed-defect
matching input; (17) is not asserted directly for that mobile problem.

## 6. A new n=16 certificate restores a violating candidate count

[Proved, inherited algebra] C75's two-species formula and factor of one
half relative to the parent hold for every n. Retaining its normalization,
the quadratic coefficient is

    T/g^2=-2 alpha^4 Re C,
    mu_quadratic/g^2=-2i alpha^4 Im A,
    E_unit=alpha^4 beta^2 g^2.

The derivation is C75 Section 3 and the referee's independent check. It
uses the legal occupied cell, which (3) reproduces. No new physical
letter is needed to increase n.

[Definition] For the new member use n=16, b=17, D=1000, and

    alpha=sqrt(2/18)sin(pi/18),       Delta=4-2cos(pi/18),
    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_17,
    phi_A=pi/32, phi_B=0, psi_A=pi, psi_B=pi/2,
    t_d=beta 2^-92,  g=2^-288,  K>=2^1280,
    0<lambda<=2^-704,       field=lambda^2 E_unit 2^-149.    (18)

Choose the midpoint marker chemical potential at b=17. The displayed K
is a lower bound for the frozen finite certificate. Section 7 imposes a
further lower bound when the fast vacancy fluid is added.

[Proved, computer-assisted] certify_n16.py makes local attributed copies
of C75's scalar-certificate workflow operative at the new n=16. It uses
256-bit integer tree and derivative calculations and a 512-bit finite
loop calculation on a 33 by 33 product graph. Approximate roots propose
rational intervals. Acceptance uses rational endpoint inertia and exact
unit-circle chord comparisons, including all 528 unordered pairs with
replacement. The output certifies

    32 simple Fermi points,
    delta_2 > 0.0145298937,
    delta_4 > 5.92937e-7,
    Hermitian hopping margin > 0.5458320,
    admissible additional normalized Bloch error = 2^-35.   (19)

The root count is exhaustive because the determinant, multiplied by
exp(16 i k), has degree at most 32. The 32 sign-changing disjoint
intervals exhaust that degree. Pair comparisons use exact rational
phases and subtract all interval radii. The raw scalar routine records
a sufficient conversion threshold g<=2^-264 in its inherited metric;
the actual choice here is the stricter g=2^-288 in budget_n16.json.

[Proved] The C75 finite-cell norm constant B_0=2^16 remains valid at
n=16. A cell has 17 sites, onsite norm below 1005, and at most twelve
unit move blocks per bond; their norm sum is below 2^16. The vacuum gap
is greater than 1/16, just as required in C75's Schur argument. Its
conservative bound ||H_2(K,g)-H_2(infinity,g)||<=2^64/K therefore carries
over with the same proof. The alpha estimate is now

    alpha^4 >=64/18^6 >2^-20.

[Proved, computer-assisted] Let W be the rational certificate congruence.
budget_n16.py evaluates with exact fractions

    err_g <=4n g^2 beta^-2 tr(W W^dagger),
    err_K <=16n tr(W W^dagger) (2^64/K)/E_unit.

At (18), these are below 10^-25.0871 and 10^-37.7304. Their sum is below
one eighth of 2^-35. Thus (19) survives finite g,K for the quadratic
projection of the exact second-order occupied-cell Hamiltonian. It does
not say that the full interacting model is Gaussian. The optional
two-cell fourth-order diagnostic in the budget file is not used as a
volume-uniform estimate.

[Proved] The all-length selector has a short general argument. Set
a_l=2[1-cos(pi/(l+1))] and delta_l=(l+1)a_l-l a_(l+1). On x>=1,

    a''(x) = 4pi sin(pi/(x+1))/(x+1)^3
             +2pi^2 cos(pi/(x+1))/(x+1)^4 >0.

Thus delta_l-delta_(l-1)=-l(a_(l-1)-2a_l+a_(l+1))<0. Since

    e_(l+1)-e_l = (2-mu_sel-delta_l)/(l(l+1)),

the midpoint 2-mu_sel=(delta_b+delta_(b-1))/2 selects b uniquely.
The smallest other density difference is attained at a neighbor.

[Proved, computer-assisted] Rational alternating-series bounds on pi and
cos give at b=17

    mu_sel = 1.912016042804466...,
    delta >1.5616829e-5 >2^-17,
    e_17 <-0.0033882037 <-2^-9,
    vacuum gap >0.09023026 >1/16.

At K>=2^1280, gamma>2^-18 suffices in (17), and lambda<=2^-704 satisfies
gamma>5B_port. F3 therefore applies also to this new member.

## 7. Order of scales and the conditional central charge

[Proved] The mobile penalty gap is kappa=JK. Its unpenalized interaction
strength nu is O(t+J B_material+|mu|), where B_material bounds the
unpenalized material terms per site, not K. The first new effective
correction from the color/grammar layer has strength

    O(nu^2/(JK)).                                           (20)

Keeping C75's old K fixed while taking J/t arbitrarily small does not
control this layer. In particular JK>>t alone is insufficient to make
(20) small relative to the smallest payload scale.

[Definition, sufficient hierarchy] Fix a desired material tolerance
tau>0, smaller than the selector margin and the retained auxiliary and
logical stability tolerances. Once g,lambda and a positive J/t have been
chosen, increase K so that

    K>=2^1280,
    JK >> nu,
    C_2 nu^2/(JK) < J tau,
    C'_2 nu^3/(JK)^2 < J tau.                               (21)

For example the third inequality requires K>C_2 nu^2/(J^2 tau).
Every quantity except K is held fixed in this choice. Increasing K
improves (12) and the finite-K certificate. All resulting couplings are
finite and independent of L. Equation (21) is a sufficient symbolic
hierarchy, not a numerical certification of a finite-J matching window.

[Assumed] There must be a positive tau for the infrared application.
This is the conventional stability/matching premise, with the velocity
of the slowest retained band included, not just the largest material
energy. At leading order (3) is the usual phase-model charge fluid
coupled through contacts to the legal material Hamiltonian. The parent
charge reference has one compact boson, K_gap=4 at half occupancy and
positive mean contact alpha_contact=1/2-2/pi^2. These reference facts are
inherited from the parent's analyzed phase model, not remeasured here.

[Assumed] Finite-J local matching while retaining the gapless charge
mode, a uniform dressed comparison for material defects, and stability
of the nonresonant multichannel Gaussian phase are still needed. Higher
order interactions, the cyclic zero-mode gluing, and the interpretation
as a collection of CFTs with different speeds belong to this infrared
input. A finite K normal form does not prove them.

[Heuristic, conditional on those inputs] With n=16 the repaired member
has sixteen payload bosons and one vacancy boson,

    d=17,       c_candidate=16+1=17 > d-1=16.                (22)

With C75's original n=15 certificate alone, the same repair instead gives
c_candidate=16=d-1. The new certificate in Section 6 is the additional
finite ingredient needed for (22). Neither count is a measured
microscopic central charge or an unconditional lattice-to-CFT result.
