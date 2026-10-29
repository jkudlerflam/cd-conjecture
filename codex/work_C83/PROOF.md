# C83. Two colored vacancies for the twelve-letter bond

[Scope] This note proves the grammar repair, its charge identities, a
new finite n=13 band certificate, frozen nonzero-port selection, and a
finite-order local penalty normal form. The central-charge conclusion
remains conditional. No microscopic central charge was measured.

[Tags] Proved means an algebraic argument or a computer-assisted rational
inequality. Numerical means a resolved finite-matrix calculation.
Heuristic means an infrared interpretation. Assumed means an input not
established here. Inherited proved inputs are identified explicitly.

## 1. The material grammar and the minimum number of colors

[Definition] Retain C80's twelve material letters

    L_0,L_1,R_0,R_1,D_A,D_B,X_A,X_B,E_0,E_1,M_A,M_B.

The markers are E_0,E_1,M_A,M_B. The two-state automaton is

| Material letter a | Outgoing state o(a) | Incoming set I(a) |
| :-- | :-- | :-- |
| L_0,L_1 | 0 | {0} |
| R_0 | 1 | {0,1} |
| R_1 | 1 | {1} |
| D_A,D_B,X_A,X_B | 1 | {0} |
| E_0 | 0 | {0,1} |
| E_1 | 0 | {1} |
| M_A,M_B | 0 | {0} |

A material pair ab is allowed if o(a) belongs to I(b). Thus an L or
marker may precede L_0,L_1,R_0,D_A,D_B,X_A,X_B,E_0,M_A,M_B. An R,D,X may
precede only R_0,R_1,E_0,E_1. This is C80's complete allowed-pair table.
The code's displaced bit is responsible for the restrictions on R_1
and E_1. Copying C82's incoming table unchanged would admit extra words.

[Definition] Add V_0,V_1 with o(V_q)=q and I(V_q)={q}. The full dimension
is d=12+2=14. Use the same allowed-pair test for all fourteen letters.
Let p_i be its diagonal mismatch projector, N_bad=sum_i p_i, and
P=1_(N_bad=0). Both vacancy letters have zero internal charge and zero
occupied number.

[Proved] Two colors suffice. Let a,b be consecutive material letters,
with r>=1 vacancies between them. The first physical bond sets the first
vacancy color to o(a). Every vacancy-vacancy bond copies it. The final
bond is allowed if and only if o(a) belongs to I(b). Thus a legal material
word at any fixed occupancy pattern with nonzero occupancy has one and
only one zero-penalty coloring. This holds for every r, every bit and
every excited letter.

Let N_matbad count bad adjacent pairs after deleting all vacancies.
Each bad material pair forces a bad physical bond in its own intervening
interval. These intervals are disjoint on the ring, including its seam.
Consequently, as diagonal operators at nonzero occupied number,

    N_matbad <= N_bad.                                      (1)

The statement includes the one-material-letter ring, where the interval
closes back onto that letter. The empty ring has two constant legal
colorings and is outside the squeezing statement used here.

[Proved, restricted minimum] One uncolored vacancy cannot implement this
repair by a nonnegative diagonal nearest-neighbor penalty while retaining
all legal material words at all vacancy placements. Legal embeddings
require both R_0 V and V L_0 to cost zero, as witnessed by legal R_0 V R_0
and L_0 V L_0 intervals. Their concatenation R_0 V L_0 then costs zero
although its squeezed pair R_0 L_0 is forbidden. Also V V must cost zero
if arbitrary gap lengths are retained. Thus one state cannot propagate
the distinction between the two outgoing states. With two states (1)
holds, proving minimality within this construction. This is not a
universal exclusion of a differently designed d=13 model.

## 2. Explicit material bond and colored lift

[Definition] Set z_A(a)=(-1)^a and z_B(a)=1. Each transition below is
accompanied by its Hermitian conjugate. Sum a=0,1, s=A,B, and all markers
m where indicated. The amplitude multiplies output ket times input bra.

| Input | Output | Amplitude |
| :-- | :-- | :-- |
| R_0 R_a | L_a R_0 | -1 |
| R_0 E_a | L_a E_0 | -1 |
| D_s R_a | L_a D_s | -z_s(a) |
| X_s R_a | L_a X_s | -exp(i phi_s) z_s(a) |
| X_s E_a | L_a M_s | -exp(i phi_s) z_s(a) |
| D_s E_a | L_a M_s | -t_d exp(i psi_s) z_s(a) |
| m D_s | m X_s | -1 |
| D_A R_1 | D_B R_0 | g |
| D_A E_1 | D_B E_0 | g |
| E_0 R_0 | M_A D_B | lambda exp(-i pi/4) |
| E_0 R_0 | M_B D_A | lambda exp(+i pi/4) |

The one-site diagonal o assigns Delta to D_s,X_s,M_s, an additional D
to X_s,M_s, mu_sel to every marker, and f a to L_a,R_a,E_a. Let w be
the transition matrix above, p_mat the material mismatch projector, and

    b_12 = w+(o tensor I+I tensor o)/2,
    h_12 = K p_mat+b_12.                                    (2)

Let embed extend b_12 by zero outside the occupied-pair block, and n_occ
be the single-site occupied projector. The complete colored bond is

    h_14 = JK p_ext + J embed(b_12)
           -t sum_(a material) sum_(q in I(a))
               (|a,V_o(a)><V_q,a| + h.c.)
           -(mu/2)(n_occ tensor I+I tensor n_occ).            (3)

Its matrix is 196 by 196. Its occupied block is J h_12. In particular,
material half-site terms remain inside occupied contacts. Moving them
to physical onsite terms would change (3). The explicit independent
implementation is bond14.py. Diagnostic parameters are ordinary floats;
the large certified hierarchy below is specified symbolically.

[Proved] The internal charge is

    q(L_a)=q(R_a)=q(E_a)=2a,
    q(D_A)=q(X_A)=q(M_A)=-1,
    q(D_B)=q(X_B)=q(M_B)=+1,
    q(V_0)=q(V_1)=0.                                        (4)

Every transition in (2) preserves total charge. For conversion the two
sums are -1+2=1+0. The port's input and outputs all have charge zero.
Hops and switches preserve the displayed sums as well. The color hop
in (3) moves the same letter a, without changing its charge or marker
status. Diagonal terms commute with all diagonal charges. Therefore the
full microscopic chain has exact conserved Q=sum q, occupied number N,
and marker number M, including all illegal words. This is an exact
U(1) statement, not a low-energy approximation. Vacancy color itself is
not conserved.

[Proved] In the frozen material chain every site's marker projector is
conserved. In (3) marker positions move with vacancies, but marker order
and the numbers of material letters between consecutive markers remain
conserved up to cyclic relabeling. The relation N=14M is a proposed
selected-sector relation, not an identity on the full Hilbert space.

[Numerical] check_lift.py reconstructs (2) independently and compares
with the local source copy of C80's bond_shift.py. At generic n=2 and
n=13 points, all matrix entries agree with error zero. The 144 material
pair permissions also agree. Hermiticity and the three charge
commutators of (3) vanish entry by entry at a generic nonzero point.
Their exact justification is the transition-by-transition argument above.

## 3. Projection, hopping and the hidden wall

[Proved] Consider V_q a -> a V_o(a). The external state entering this
pair is q and the outgoing state is o(a), on both sides. Its reverse
has one microscopic channel for every q in I(a). On P, the predecessor
selects one of them. Therefore every projected swap has amplitude -t,
with no extra multiplicity or square-root factor.

Every material transition leaves the outgoing state of its rightmost
letter unchanged. If it alters the incoming compatibility of its leftmost
letter, the preceding extended bond tests the change. The coloring map
therefore intertwines projected material moves with C80's legal material
moves. Let U_col be the coloring isometry and g_j the number of vacancies
after material letter j. On the legal material space,

    U_col^dagger P H P U_col
      = H_hop + J sum_j 1_(g_j=0) (P_mat b_12,j P_mat)
        -mu N.                                             (5)

For a ring the formula keeps the physical root and cyclic word rotation.
Coloring commutes with rotation and introduces no additional twist or
independent gapless color variable. It is a basis bijection at every
nonempty occupancy pattern, not a tensor-product enlargement of the
legal material code.

[Proved] For example, the legal interval R_0 V_1 R_0 R_0 admits the
vacuum move on its last two letters. The result is R_0 V_1 L_0 R_0,
which has a mismatch on V_1 L_0 and costs JK. Replacing V_1 by V_0
instead creates the mismatch R_0 V_0. A color change somewhere in a
longer vacancy string also costs JK. This removes the unpenalized hidden
wall of referee C75's F5.

The full finite-K Hamiltonian still mixes P and 1-P. Reverse vacancy
hops can choose a color incompatible with the predecessor. The repaired
claim is energetic suppression, not exact finite-K invariance.

[Proved, finite volume] Write H=kappa N_bad+B, with kappa=JK, and assume
that the finite invariant subspace under consideration contains a legal
trial state. Put Q_bad=1-P. If kappa>2||B||, a ground vector psi has
E<=||B||, while

    Q_bad(H-E)Q_bad >= (kappa-2||B||)Q_bad.

Projecting its eigenvalue equation gives

    Q_bad psi = -[Q_bad(H-E)Q_bad]^(-1)Q_bad B P psi,
    ||Q_bad psi|| <= ||B||/(kappa-2||B||).                   (6)

The projector onto illegal material words is bounded by Q_bad by (1).
Its weight is therefore bounded by the square of (6). This proves
K^-2 suppression at fixed volume and fixed positive J. Because ||B||
is extensive, (6) is not a volume-uniform ground-state overlap estimate.
The local statement is given in Section 7.

## 4. The n=13 quadratic certificate

[Proved, inherited C80 basis equivalence] The constrained C80 code shifts
the bit displaced by a bitless token onto the next site, with the final
bit on E_a. A legal n-data cell has (5n+3)2^n states. Its vacuum path
has n+1 clock positions and an exact 2^n-dimensional ground space at
infinite K. Its vacuum energy and endpoint amplitude are

    E_v=-2cos(pi/(n+2)),
    alpha=sqrt(2/(n+2)) sin(pi/(n+2)).

C80 PROOF Sections 2-3 prove the all-n unitary equivalence to the legal
folded carrier, including the two-site conversion and port. We use that
established equivalence rather than re-derive its constants.

[Definition] For the new member set n=13, b=n+1=14 and

    D=1000, Delta=4-2cos(pi/15),
    phi_A=pi/26, phi_B=0, psi_A=pi, psi_B=pi/2,
    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1),
    beta=1/d_14,
    t_d=beta 2^-74, g=2^-256, lambda=2^-576, K>=2^2304,
    E_unit=alpha^4 beta^2 g^2, f=lambda^2 E_unit 2^-121.     (7)

The material selector is defined in Section 6. The microscopic field
f in (2) is the physical field, with lambda^2 already included.

[Proved, inherited reduction] The scalar carrier has ordered vertices
M,X_(n-1),...,X_0,D_0,...,D_(n-1), diagonal 1004 on M,X and 4 on D,
path entries -1, and closing edge -i t_d. Write H_c for its matrix and

    R=(H_c tensor I+I tensor H_c)^(-1), x=R|M,D_0>.

U_i retains row D_i of x and V_j retains column D_j, regarded as vectors
on the product graph. With

    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j,

the g^2 coefficient, in C80's normalization, has

    T=-2 alpha^4 g^2 Re C,
    mu_quad=-2i alpha^4 g^2 Im A.                            (8)

C80's CAR calculation proves that this coefficient is quadratic. Higher
even conversion orders can contain interactions. Exact microscopic U(1)
persists, but exact Gaussianity at finite g is not claimed.

[Proved, computer-assisted] certify_n13.py regenerates the scalar data
with C80's attributed integer/rational pipeline. The product graph has
27^2=729 vertices. Tree and derivative solves use 256-bit fixed point,
and the full finite loop uses 512 bits. No n=12 or n=16 output was copied
as an n=13 result. The copied routines preserve their C75/C33 attribution.

The Jacobi contraction is at most 1/2+t_d/4. Taking
bit_length(d_14)+512+20 steps and retaining the inherited ten-unit
component radii bounds floor error and the iteration tail. The second
solve includes its source error. Products are enclosed using
10(|x_r|+|x_i|+|y_r|+|y_i|)+200 integer units. The integer LDL identity
and rational congruence transport these radii to a Bloch norm enclosure
below 9.385e-95. The exact dyadic hopping positivity margin is above
0.5608079.

Numerical roots only propose brackets. For rational t=tan(k/2), the
endpoint matrices use exact cos k=(1-t^2)/(1+t^2) and
sin k=2t/(1+t^2). A rational approximation to an eigenbasis is proved
invertible by its Gram matrix, and diagonal dominance after that
congruence proves the endpoint inertia, including entry errors.

There are 26 disjoint sign-changing brackets, with negative inertia
increasing from 0 to 13 and then decreasing to 0. The nearest-neighbor
Bloch determinant multiplied by exp(13ik) is a polynomial of degree at
most 26. The 26 distinct zeros exhaust this degree. Each is simple,
and there are no other Fermi points.

The phase midpoint of each rational bracket is exact. Subtract the
bracket chord radii from adjacent phase separations and from products
of pairs of phases. The result is

    delta_2 > 0.0149087301,
    delta_4 > 8.70283e-8,                                   (9)

for all 351 unordered pairs with replacement. The latter implies distinct
sea centers modulo pi, since equality of two center phases would give
coincident Fermi-pair sums modulo 2pi. These are finite nonresonance
inequalities, not rational independence against all integer relations.
The endpoint enclosures allow an additional normalized Hermitian Bloch
perturbation of norm 2^-35.

Exact intervals and margins are in root_certificate_exact_n13.json.
root_certificate_n13.json records their decimal summaries. Its inherited
finite_g_power=-219 is an intermediate threshold, not the schedule in
(7), which uses -256.

## 5. Finite-g and finite-K transport, including zero-port splitting

[Proved, inherited resolvent estimate] The port source has norm at most
sqrt(2) alpha^2, conversion norm at most 2g, and the two-excitation source
gap is at least three. Odd conversion orders do not return to the source
species. The even Neumann tail gives

    ||H_2(g)-H_2(0)-g^2 H_22|| <= alpha^4 g^4/4.             (10)

First order vanishes, including finite K, by excited-token number. At
second order distinct ports cannot restore both outer-cell excitation
numbers. Only onsite and neighboring-cell operators remain. The finite-K
vacuum band is not scalar. Its Hermitian second-order formula is

    (H_2)_ab=-1/2 sum_q V_aq V_qb
                [(E_q-E_a)^(-1)+(E_q-E_b)^(-1)].            (11)

Replacing all E_a by one energy would omit the C80 splitting.

[Proved] C80's cell norm bound B_0=2^16 remains valid at n=13. There are
14 sites, each onsite norm less than 1007, and fewer than 64 units of
transition norm per bond. The finite field is less than one. Two cells
therefore have unpenalized norm at most 2B_0. The legal vacuum gap is
above 1/16, certified below.

For H_K=K N_bad+B, |z|<=2B_0 and K>=8B_0, the one-cell high block has
inverse norm at most 2/K. Hence its Schur correction is at most 2B_0^2/K.
Repeating on two cells gives 2^35/K before rounding the energy bound.
C80's spectral projection and polar overlap argument then gives

    max_a |E_a-2(E_v+mu_sel)| <= 2^36/K,
    ||G_K-G_infinity|| <= 2^44/K.                           (12)

The source gap of three and the resolvent identity give a per-column
second-order change below 2^48/K, as in C80. The two-cell code now has
4^13 columns, so conversion from a column bound to an operator bound
costs at most 2^13. Allowing symmetrization and the port contraction,
a conservative enlarged bound is

    ||H_2(K,g)-H_2(infinity,g)|| <= 2^68/K.

Together with zero-port energies, field compression and subspace
alignment, the physical effective bond through second order changes
by less than 2^84/K. This bound has no lambda^2 prefactor. The larger
constants leave slack relative to C80's 2^64 and 2^80 at n=12.

[Proved, finite two-cell scope] The parity that negates the port makes
its isolated two-cell effective operator even in lambda. A contour of
radius 1/64 around its vacuum band and an analytic port disk of radius
2^-12 have uniformly convergent resolvent series. Their projection and
canonical overlap intertwiner are bounded by two. Cauchy's bound, with
slack, bounds the tail by 2^52 lambda^4 for lambda<=2^-13. This is a
finite two-cell estimate, not a many-body linked-cluster theorem.

[Proved] Let W be the rational dyadic congruence and tau_W=tr(W W^dagger).
A quadratic coefficient extracted by normalized Hilbert-Schmidt projection
is bounded by four times operator norm. Summing onsite and hopping
contributions gives the safe Bloch factor 16n. Thus the normalized errors
are bounded by

    err_g = 4n beta^-2 g^2 tau_W,
    err_K = 16n tau_W (2^84/K)/(lambda^2 E_unit),
    err_lambda,local = 16n tau_W 2^52 lambda^2/E_unit.        (13)

The physical denominator in err_K is necessary. Reducing lambda at fixed
K without rechecking this inequality is not a valid hierarchy limit.

[Proved, computer-assisted] sin(pi/15)>=2/15 gives
alpha^4>=64/15^6>2^-18. budget_n13.py uses that rational lower bound,
the exact beta recurrence and the integer congruence. Its results are

| Error | Upper bound |
| :-- | :-- |
| Finite g | 1.87e-33 |
| Finite K, including zero-port splitting | 1.12e-40 |
| Local two-cell lambda remainder | 2.60e-50 |

Their sum is below (2^-35)/8. Acceptance uses Fractions, not the decimal
summaries. Therefore the root and pair-sum certificates survive for the
quadratic projection of the finite-g, finite-K, second-order material
effective Hamiltonian. The local pair remainder also fits in the budget.
Longer-range many-cell corrections and interacting operators still
require the infrared assumptions. This calculation does not certify
them by the determinant-degree argument.

## 6. F3 transfers to complete frozen marker sectors

[Proved, inherited C80 all-cell inequality] With l-1 data sites and one
terminal marker, the full cell Hilbert space has dimension 4*8^(l-1).
It contains all illegal words and all auxiliary excitations. C80's
partial-isometry and legal-fragment argument gives

    E_l^infinity-(l-2)_+/(K-4) <= E_l(K) <= E_l^infinity,
    E_l^infinity=mu_sel-2cos(pi/(l+1)).                      (14)

It requires Delta>2+g and |t_d|<=1, which (7) satisfies. The nonnegative
field does not weaken the lower bound and vanishes on the empty-bit
trial state for the upper bound. The marker-free ring energy density is
at least -1/(K-4). These inherited all-length bounds are not inferred
from the small ED.

[Definition] Put a_l=2[1-cos(pi/(l+1))],
delta_l=(l+1)a_l-l a_(l+1), and

    mu_sel=2-(delta_14+delta_13)/2.                          (15)

For e_l=E_l^infinity/l,

    e_(l+1)-e_l=(2-mu_sel-delta_l)/(l(l+1)).

[Proved] The second derivative

    a''(x)=4pi sin(pi/(x+1))/(x+1)^3
            +2pi^2 cos(pi/(x+1))/(x+1)^4

is positive for x>=1. Therefore delta_l is strictly decreasing. Equation
(15) makes e_l decrease up to l=14 and increase thereafter. Its minimum
competing density difference is attained at l=13 or l=15. Rational
alternating-series enclosures for pi and cosine in budget_n13.py prove

    mu_sel=1.8742920463679156...,
    delta=min_(l!=14)(e_l-e_14)>3.8800592e-5>2^-15,
    e_14<-0.0058573682<-2^-8,
    vacuum gap>0.12920428>1/16.                             (16)

[Proved] C82's F3 proof needs one modification in its cell bookkeeping.
In C80 a marker followed by R_1 or E_1 is a bad pair. Nevertheless all
four markers have the same outgoing state 0, so a left boundary penalty
depends only on the first letter of the cell on its right. Assign that
term to the right cell. The switch is also the identity on its preceding
marker. All zero-port terms then act on individual complete cells. Only
the two port transitions couple cells, with norm

    B_port=sqrt(2)|lambda|.

This norm follows from the star with input E_0 R_0 and two orthogonal
outputs. No exact degeneracy of the finite-K cell code is needed.

[Proved] Let A_m be the exact ground energy of m selected full cells in
an open chain. Cutting one port changes the ground energy by at most
B_port, so

    |A_(m+r)-A_m-A_r|<=B_port.

The upper and lower almost-additive bounds give the limit u=lim A_m/m.
Apply the same inequality to q copies of an m-cell chain and divide by q
to obtain |A_m-m*u|<=B_port. The bounds used below are

    A_m >= m*u-B_port,
    S_m <= m*u+2B_port,
    u <= E_14(K) <= E_14^infinity.                          (17)

Here S_m is the periodic selected-sector energy. The last inequality
uses a product of zero-excitation finite-K cell ground states, on which
every port expectation vanishes.

Set gamma=delta-1/(K-4). For every wrong length l, (14) implies

    E_l(K)>=l e_14+l gamma >= l*u/14+l gamma.

On a compatible ring L=14m with k wrong cells, delete all ports touching
them, at most 2k. The selected cells form at most k open runs. Apply
(17) to the runs and the preceding bound to each wrong cell. This gives

    E_wrong >= (u/14)L+gamma sum_wrong l-3k B_port,
    E_wrong-S_m >= gamma sum_wrong l-3k B_port-2B_port
                 >= k(gamma-5B_port).                     (18)

It also covers a ring with only wrong cells. A ring with no markers has
no ports and is excluded by (14)'s marker-free bound and (16). At (7),
gamma>2^-16 and gamma>5sqrt(2)lambda. Thus frozen nonzero-port marker
selection holds, with illegal components retained inside each cell.

[Numerical] check_frozen.py tests two complete 1024-dimensional marker
sectors on four material sites, with lengths (2,2) and (1,3). It uses
the b=2 midpoint selector, K=100, D=4, g=.03, t_d=.001 and
lambda=0,.001,.003. No charge or excitation sector is discarded. Energy
differences are 0.58331059, 0.58331064, 0.58331102, above the respective
F3 lower bounds 0.17442881, 0.16028668, 0.13200241. Residuals are below
8.8e-12. This diagnoses the new bookkeeping, not the n=13 hierarchy.

[Scope] Equation (18) selects frozen marker sectors. It is not a mobile
finite-J defect theorem. In (3), cell terms are multiplied by fluctuating
contacts and share a gapless vacancy fluid. Transferring (18) there still
requires a uniform common-bulk-subtracted dressed-defect comparison.

## 7. F4 transfers with the actual gap JK

[Definition] On a finite ring write H=kappa N+B, N=sum_i p_i and
B=sum_X B_X, with local interaction norm

    nu=sup_i sum_(X contains i)||B_X||.

Here kappa=K in the frozen chain and kappa=JK in the lift. The p_i commute
and have integer eigenvalues, exactly as required in C82's proof.

[Proved] Resolve any local A into integer penalty transfers

    A_m=(1/(2pi)) int_0^(2pi)
          exp(i theta N) A exp(-i theta N) exp(-im theta) dtheta.

Only p_i meeting its original support X enter. Thus A_m has support in
X enlarged by one site at either end, norm at most ||A||, and vanishes
for |m|>|X|+1. Also [N,A_m]=m A_m. For a two-site term there are at most
three participating penalty projectors. The enlarged incoming sets in
Section 1 cause no change to this argument.

Write B=B_0+B_off and set

    S_1=sum_(m!=0) B_m/(kappa m).

It is anti-Hermitian and [S_1,kappa N]=-B_off. Its local norm is
O(nu/kappa), with a volume-independent constant. The BCH expansion gives

    exp(S_1)H exp(-S_1)=kappa N+B_0+F_2+O_local(nu^3/kappa^2),
    F_2=[S_1,B_0]+[S_1,B_off]/2.

Take S_2=sum_(m!=0)(F_2)_m/(kappa m). With S=S_1+S_2,

    exp(S)H exp(-S)=kappa N+B_0+(F_2)_0+R_2,
    [N,B_0+(F_2)_0]=0,
    ||(F_2)_0||_local <= C_2 nu^2/kappa,
    ||R_2||_local,exponential <= C'_2 nu^3/kappa^2.          (19)

On P=ker N, B_m P=0 for m<0 and the compression is

    P(B_0+(F_2)_0)P
      =PBP-sum_(m>0)PB_(-m)B_mP/(kappa m).                 (20)

This proves the second-order negative sign and shows why the denominators
are JK in the mobile chain, even on a vacancy bond.

[Proved] The constants in (19) are independent of chain length. Resolve
the generators into interval interactions. Disjoint terms commute, so
a nonzero nested commutator has connected overlapping support. At the
next step the sum of norms of terms meeting the current support I is
bounded by |I| times the generator's local norm, and the commutator
adds a factor two. Support grows linearly with order. Iteration gives
k! C^k times the product of local norms. BCH's 1/k! cancels this factorial;
for sufficiently small nu/kappa the remaining geometric tail converges
in a weaker exponential interaction norm. Terms from kappa N acquire
one extra generator and the same nu^3/kappa^2 bound.

At any fixed order r, remove the nonzero integer Fourier components
recursively by division by kappa m. Connected supports remain of size
O(r), giving

    U_r H U_r^dagger=kappa N+D_r+R_r,
    [D_r,N]=0,
    ||R_r||_local,exponential<=C_r nu(nu/kappa)^r.            (21)

No convergence as r tends to infinity is asserted. This is not an exact
all-order block diagonalization or a globally isolated many-body band.

For a vector phi in P, P[S,p_i]P=0. Expanding its dressed local mismatch
expectation therefore starts at second order in S, giving
O((nu/kappa)^2) per bond. The total outside-P weight need not be small
uniformly in volume. Identifying actual low-energy states with this
dressed space additionally requires stability under R_r.

[Numerical] check_sw.py uses a literal three-site, nonzero-port closure
of dimension 60, including excited letters. At K=40,80,160,320 it finds
K^2||R_2||=7.95144,7.93017,7.91948,7.91413. The second-order compression
agrees with (20) below 6.2e-21, and the unitary residual is below 4.5e-16.
This verifies the signs and finite-order scaling at resolved parameters.

## 8. Leakage checks and what the mobile budget establishes

[Numerical] check_lift.py constructs invariant closures under every
literal local transition, without truncating by mismatch number. It
exhausts all 14^4-2^4=38400 nonempty four-site basis words for (1).
Projected colored matrices agree entry by entry with the corresponding
legal uncolored matrices in all three ED tests. These tests use J=t=1,
mu=mu_sel=field=0. Excited-letter tests use n=2 in the parameter formula,
D=4, Delta=4-2cos(pi/4), g=.07, t_d=.003, lambda=.1 and the phases
in Section 2. They are diagnostic couplings, not (7).

| Test and invariant-closure dimension | K=100 | K=1000 | K=10000 |
| :-- | --: | --: | --: |
| Two n=3 zero-bit cells, one uncolored vacancy, 2304 | .0267832 | .0257054 | .0256079 |
| Same, two colors, 4608 | 1.07551e-4 | 1.08213e-6 | 1.08279e-8 |
| Nonzero-port n=2 cell, one vacancy, two colors, 268 | 2.31273e-7 | 2.67675e-9 | 2.71795e-11 |
| Nonzero-port n=2 cell, two vacancies, two colors, 1045 | 1.33164e-7 | 1.52294e-9 | 1.54444e-11 |

Entries are illegal material weights. The first pair of tests uses all
zero bits and zero port, so its invariant vacuum Hamiltonian coincides
with C82's reduced diagnostic. It was regenerated from the literal C83
bond. The excited tests have total charge Q=2 and include the displaced
bit and both conversion channels. Their uncolored counterparts tend to
2.5158e-5 and 2.6443e-5 at K=10000, respectively. Every run also records
physical mismatch expectation and its eigenpair residual. The largest
residual is below 1.55e-11. These are ground states of the stated invariant
closures, not a claim about all charge sectors of arbitrary rings.

[Proved sufficient finite-K condition] In the mobile chain
nu=O(t+J B_material+|mu|), excluding the penalty. Fix the material
parameters, a positive J/t and a positive tolerance tau in material
energy units. Choose a finite K satisfying

    K>=2^2304,
    nu/(JK) small enough for the fixed-order local expansion,
    C_2 nu^2/(JK) < J tau/2,
    C'_2 nu^3/(JK)^2 < J tau/2.                             (22)

For example the last two require

    K>2 C_2 nu^2/(J^2 tau),
    K>sqrt(2 C'_2 nu^3/(J^3 tau)).                          (23)

The constants are those of (19), fixed at a chosen local norm and order.
Every quantity on the right is independent of K and volume. A finite
choice exists for every fixed positive J,tau and improves the material
finite-K bounds. For the material quadratic tolerance alone, one can
use tau_band=lambda^2 E_unit 2^-35/(128 n tau_W), with further slack
if desired. Selector and auxiliary gaps impose additional positive
upper limits on tau.

The lower bound K=2^2304 certifies the frozen finite-dimensional
calculation. It alone is not a certification of the mobile model at an
arbitrarily small J. Merely requiring JK>>t also misses corrections
relative to the much smaller logical band scale. Equations (22)-(23)
state the additional mobile budget, without claiming an uncomputed
numerical finite-J matching window.

## 9. Conditional infrared count and remaining assumptions

[Proved, reference identification] By (5), the legal colored problem is
the usual squeezed vacancy fluid coupled through contacts to the legal
12-letter material Hamiltonian. At full occupancy it reduces to C80's
constrained bond. The color construction adds no legal degree of freedom
at a fixed occupancy pattern and material word. At nonzero occupancy it
therefore adds one vacancy-fluid mode in the reference description, not
two independent color modes.

[Assumed, inherited charge reference] At physical half occupancy the
parent charge fluid has one compact boson and mean contact per occupied
letter alpha_contact=1/2-2/pi^2. C81's cyclic-sea projection applies to
the legal material bond in (5) because coloring is rotation-equivariant.
Its physical root, cyclic quotient and all allowed charge/marker
harmonics must be retained. Half-density tuning at finite J uses positive
compressibility and remains part of the matching assumption.

[Proved conditional on C81's Gaussian dictionary] The harmonic argument
also transfers with b=14. A payload vertex of cell momentum q and
reference dimension Delta_pay has branches Q_h=q+2pi ell and

    k_phys=Q_h/28 mod 2pi,
    Delta_ref=Delta_pay+[Q_h/(28pi)]^2.

A uniform Hamiltonian term requires Q_h=56pi h. For h!=0 its charge
contribution is at least four. For h=0 the 26-root and 351-pair
certificate excludes nonforward neutral quadratic and quartic vertices.
Higher neutral fermion vertices have reference dimension at least three;
exact Q and N forbid anomalous payload pairing and number-changing charge
terms. Thus the potentially dangerous low-order vertices are removed
within the inherited dictionary, with no all-order rational-independence
assumption. For vertices below dimension three, C81's physical mismatch
bound becomes

    delta_phys >= min(delta_internal/28,pi(2-sqrt(3)))
               >3.10815e-9,
    delta_internal=min(delta_2,delta_4).

The first marker harmonic is an oscillating observable of dimension
1/196 at wavevector pi/14, not a uniform relevant perturbation. These
statements retain all integer harmonic branches. Their microscopic use
requires the dictionary and matching, as in C81.

[Assumed] A completed microscopic c=14 result still needs
(i) massive-sector matching in the required local coefficient norms,
(ii) low-frequency vacancy matching at fixed positive J,
(iii) a uniform mobile dressed-defect comparison,
(iv) multichannel stability for the interacting retarded action and
(v) sufficient current, stress or thermal control to identify the
conformal modes. The positive infrared tolerance needed for these steps
is not computed by the finite band certificate. F3 closes the frozen
sector comparison only. F4 supplies the extra finite-order grammar
layer only. Small leakage is not a substitute for these infrared steps.

C81's grading qualification also remains. exp(i pi Q/2) is logical
fermion parity on the code, but can equal +/-i on odd-Q auxiliary states.
The older microscopic Z2 grading cannot be imported unchanged before
eliminating or explicitly representing these states.

[Heuristic, conditional on those inputs] The n=13 material reference
has thirteen payload bosons and the vacancy reference contributes one.
Thus

    d=14, c_candidate=13+1=14>d-1=13.

The requested finite threshold is met, n=13=d-1. This is a certified
finite member supporting a conditional d=14, c=14 candidate, not an
unconditional counterexample. A single-color d=13 repair of this
unchanged grammar does not satisfy (1).

## 10. Reproduction and provenance

Run Python from this folder. No script needs a network call or writes
outside it. All source dependencies needed by the checks are local.

    OPENBLAS_NUM_THREADS=1 python3 -B certify_n13.py
    OPENBLAS_NUM_THREADS=1 python3 -B budget_n13.py
    OPENBLAS_NUM_THREADS=1 python3 -B check_lift.py
    OPENBLAS_NUM_THREADS=1 python3 -B check_frozen.py
    OPENBLAS_NUM_THREADS=1 python3 -B check_sw.py

The first two regenerate the integer, rational and finite-parameter
certificates. The remaining scripts regenerate lift_checks.json,
frozen_checks.json and sw_checks.json. Source copies from C80 retain
attribution to C75 and C33. bond14.py is reconstructed from C80's printed
transition table, and its comparison to bond_shift.py checks that
independent reconstruction. The local STATE.md supplies the handoff;
project-level memory was not modified because it is outside the permitted
folder. No external files, network commands, cluster jobs or subagents
were used.
