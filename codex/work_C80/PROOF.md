# C80. A twelve-letter shifted-bit encoding

[Scope] The finite statements below concern a nearest-neighbor bond, its
selected finite cells and its second-order effective Hamiltonian. The
infrared identification is conditional. No central charge is measured.
The derivations have not had an independent referee.

## 1. Explicit bond

[Definition] The twelve letters are

    L_0,L_1,R_0,R_1,D_A,D_B,X_A,X_B,E_0,E_1,M_A,M_B.

Write M={E_0,E_1,M_A,M_B} for the marker letters. Marker positions are
conserved. Define the charge by

    q(L_a)=q(R_a)=q(E_a)=2a,
    q(D_A)=q(X_A)=q(M_A)=-1,
    q(D_B)=q(X_B)=q(M_B)=+1.

The complete allowed-pair table is

| Left letter | Allowed right letters |
| :-- | :-- |
| L_0,L_1 or any marker | L_0,L_1,R_0,D_A,D_B,X_A,X_B,E_0,M_A,M_B |
| R_0,R_1,D_A,D_B,X_A,X_B | R_0,R_1,E_0,E_1 |

Every other ordered pair has penalty K. Let z_A(a)=(-1)^a and z_B(a)=1.
The following table specifies all off-diagonal terms, including their
Hermitian conjugates. Sum over a=0,1 and s=A,B where present. An amplitude
multiplies the output ket times the input bra.

| Input | Output | Amplitude |
| :-- | :-- | :-- |
| R_0 R_a | L_a R_0 | -1 |
| R_0 E_a | L_a E_0 | -1 |
| D_s R_a | L_a D_s | -z_s(a) |
| X_s R_a | L_a X_s | -exp(i phi_s) z_s(a) |
| X_s E_a | L_a M_s | -exp(i phi_s) z_s(a) |
| D_s E_a | L_a M_s | -t_d exp(i psi_s) z_s(a) |
| m D_s, for each m in M | m X_s | -1 |
| D_A R_1 | D_B R_0 | g |
| D_A E_1 | D_B E_0 | g |
| E_0 R_0 | M_A D_B | lambda exp(-i pi/4) |
| E_0 R_0 | M_B D_A | lambda exp(+i pi/4) |

The last two rows are lambda V. There are no additional port terms.
The one-site diagonal operator o assigns Delta to every D_s,X_s,M_s,
an additional D to X_s,M_s, mu_sel to every marker, and lambda^2 epsilon
times a to L_a,R_a,E_a. The repeated bond is

    h = h_table + K Pi_bad + (o tensor I + I tensor o)/2.       (1)

This is a one-site translation-invariant Hermitian bond of range one.
Its dimension is 144 by 144. `bond_shift.py`, with size=12, implements
(1). Sizes 14 and 13 implement the intermediate encodings below.

[Proved] Every row conserves the sum of the two on-site charges. For
example, D_A R_1 and D_B R_0 both have charge one, and the port output
has charge zero. Thus [h,q tensor I+I tensor q]=0 exactly, including
all illegal words. The marker projector at every site also commutes
with the full chain Hamiltonian. The excited-token number is conserved
by H_0, separately in each fixed cell. A port changes it by one in each
incident cell.

## 2. Where the letters were saved

[Definition] First use a fourteen-letter encoding

    L_a,R_a,D_v,D_A,D_B,X_A,X_B,E_a,M_v,M_A,M_B.

The tokens D_v,D_A,D_B,X_A,X_B carry no data bit. For n logical bits
(a_0,...,a_(n-1)), a token T at position j<n is encoded as

    L_(a_0)...L_(a_(j-1)) T
    R_(a_j)...R_(a_(n-2)) E_(a_(n-1)).                       (2)

At j=n the word is

    L_(a_0)...L_(a_(n-1)) M_s.                              (3)

In (2), T=D_v corresponds to the vacuum; T=D_s or X_s corresponds
to an excited species. Equation (3) has s=v,A,B. The empty marker now
stores the last bit. The token displaces the bit immediately to its
right instead of carrying it.

[Proved] A hop T R_a -> L_a T crosses precisely logical bit a_j.
The Jordan-Wigner sign is therefore z_s(a_j), as in C75. Conversion
flips that same bit using D_A R_1 <-> D_B R_0. At the final data site,
use D_A E_1 <-> D_B E_0. These are still two-site operators. Express
hops, the switch, and the port have the same scalar amplitudes as the
18-letter parent in C75. The port acts on M_v D_v and does not need
to read a data bit.

[Proved] Identify M_v with E_0, giving thirteen letters, then identify
D_v with R_0, giving twelve. These identifications are injective on the
legal words (2)-(3), although they are not injective on the on-site
Hilbert space. In a vacuum word the first non-L data letter is R_0;
it identifies the hole position. Every later R letter carries a bit.
If every data letter is L, the vacuum token is the terminal E_0.
An excited word has its unique D_s or X_s, or ends in M_s.
Thus neither the logical bits nor the clock position are lost.

[Proved] The allowed-pair table in Section 1 characterizes exactly these
words inside a cell. Before the first non-L data letter there are only
L letters. That first letter can be R_0,D_s,X_s. Every subsequent data
letter is R_0 or R_1, and the cell ends in E_a. If there is no such first
letter, the terminal letter is E_0,M_A or M_B. There are consequently

    (5n+3) 2^n

legal cell states. This is the same dimension as C75's constrained
parent. Mapping its bit-carrying token words to (2)-(3) is a permutation
of orthonormal bases. Compression of (1) to these legal words agrees
with that parent, including conversion and the two-cell port, term by
term. This proves the equivalence for every n.

[Proved limitation] The full twelve-letter H_0 does not preserve the
legal space. A vacuum move inside an R tail creates a bad pair. Its
bit shift can depend on the payload. The finite-K vacuum code is therefore
not exactly degenerate. This differs from C75 and is included in the
bounds below. A port can also act on an E_0 storing a bit instead of
a vacuum token. Such extra processes are retained in the literal bond.

[Numerical] At n=1,2, the projected H_0 and V matrices agree entry for
entry with the independent 18-letter parent builder, for all three
sizes 14,13,12. Pair dimensions are 256 and 2704. These comparisons
include the shift of the last bit into the marker and the boundary
conversion. `small_checks.json` records the results.

## 3. Constrained code and second-order reduction

[Proved] At infinite penalty the vacuum path has n+1 positions, so

    E_v=-2 cos(pi/(n+2)),
    v_j=sqrt(2/(n+2)) sin(pi(j+1)/(n+2)),
    G|a>=sum_(j=0)^n v_j |j;a>,
    alpha=v_0=v_n.                                         (4)

The marker contributes mu_sel. The vacuum gap is
2[cos(pi/(n+2))-cos(2pi/(n+2))]. On the code Q=2N_data.

Choose phi_A=pi/(2n), phi_B=0, psi_A=pi, psi_B=pi/2 and
Delta=4+E_v. Gauging the A carrier by the parity of the crossed bits,
and removing the express phases, gives the common scalar graph

    M, X_(n-1),...,X_0,D_0,...,D_(n-1).

Its diagonal is 4+D on M,X and 4 on D. Path entries are -1, and the
closing edge from D_(n-1) to M is -i t_d. Energies are above E_v.
Denote this matrix by H, and put

    R=(H tensor I+I tensor H)^(-1),  x=R|M,D_0>.

Let U_i retain row D_i of x and V_j retain column D_j, as vectors on
the product graph. Set

    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j.

[Proved] The code source has amplitude alpha^2 and species structure
|B,A>-P_left|A,B>, up to a common phase. One conversion reaches the
orthogonal AA and BB sectors. In global fermions its operator parts are

    Y_AA=sum_i U_i f_left,i^dagger-sum_j V_j f_right,j^dagger,
    Y_BB=P_left [sum_i U_i f_left,i+sum_j V_j f_right,j].

In -alpha^4 g^2(Y_AA^dagger R Y_AA+Y_BB^dagger R Y_BB), adjacent parity
factors cancel in the second term. Applying f_i f_j^dagger=
 delta_ij-f_j^dagger f_i to the first term leaves a scalar and a
number-conserving quadratic operator. There is no quartic term at g^2.
The hopping and combined onsite matrices are

    T=-2 alpha^4 g^2 Re C,
    mu=-2i alpha^4 g^2 Im A.                                (5)

These are C75's two-species matrices, without a change of scale.
The proof uses the basis equivalence, not extrapolation of the small ED.
The finite-g exact second-order operator is generally interacting.

[Proved] The pair source norm is at most sqrt(2) alpha^2, the pair
conversion norm at most 2g, and the free source-sector gap is at least
three. Odd conversion orders cannot return the source species. For
g<=1/2, the tail of the resolvent expansion gives

    ||H_2(g)-H_2(0)-g^2 H_22|| <= alpha^4 g^4/4.             (6)

For example, summing the even Neumann tail gives
(2 alpha^4/3)(2g/3)^4/[1-(2g/3)^2], bounded by the right side of (6).

[Proved] PVP=0 also at finite K, since the low vacuum subspace has zero
excited tokens. Two different neighboring ports cannot restore the
excitation numbers in their outer cells, so second order contains only
single-cell and neighboring-cell operators. On two cells, parity of
the left excited-token number negates V and commutes with H_0. The
canonical effective operator is even in lambda.

For a finite-K vacuum band with eigenvectors |a> and energies E_a,
the Hermitian second-order formula used in the checks is

    (H_2)_ab = -1/2 sum_q V_aq V_qb
               [1/(E_q-E_a)+1/(E_q-E_b)].                 (7)

All q in the full (1,1) excited sector are retained. It would be wrong
to replace all E_a by one energy at the diagnostic K=100.

## 4. The finite n=12 certificate

[Definition] The certified cell has n=12, length b=13. Set D=1000,
Delta=4-2cos(pi/14), phi_A=pi/24, phi_B=0, psi_A=pi, psi_B=pi/2.
Define integers

    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1), beta=1/d_13.

The finite schedule is

    t_d=beta 2^-66,  g=2^-216,  K=2^2048,
    lambda=2^-480,  E_unit=alpha^4 beta^2 g^2,
    epsilon=E_unit 2^-112.                                 (8)

The physical field is lambda^2 epsilon. These definitions, rather than
floating-point microscopic arrays, specify the hierarchy.

[Proved, computer-assisted] `run_certificate.py` regenerates the scalar
certificate at n=12. The scalar product graph has dimension 25^2=625.
Tree solves use 256-bit fixed-point integers; the full loop uses 512.
No C75 numerical output was reused as an n=12 result. The arithmetic
routines are local copies of C75, descended from C33.

The fixed-point Jacobi map has infinity-norm contraction at most
1/2+t_d/4. With d_13.bit_length()+512+20 iterations, floor errors and
tails are enclosed by ten integer units per real or imaginary component.
The second solve includes the first solve's source error. Product radii
use 10(|x_r|+|x_i|+|y_r|+|y_i|)+200. Fraction-free congruence transforms
these entry enclosures to a normalized row error below 4.296e-99.

The integer LDL identity is checked exactly. Its dyadic row rescaling
W gives a tree metric with diagonal entries in [1,4]. Exact rational
diagonal-dominance tests prove a Hermitian hopping margin above 0.5652
in that metric. Numerical generalized eigenvalues only propose brackets.
Acceptance uses rational t=tan(k/2), with

    cos k=(1-t^2)/(1+t^2),  sin k=2t/(1+t^2).

At each endpoint an integer approximation to a numerical eigenbasis is
proved invertible. The transformed rational matrix is strictly diagonally
dominant with certified diagonal signs, including the entry error. Its
inertia is therefore established without trusting floating-point signs.

[Proved, computer-assisted] There are 24 disjoint brackets, with inertia
changes 0 to 1 through 11 to 12, then 12 to 11 through 1 to 0. Each
contains a sign-changing root. Multiplying det h(k) by exp(12ik) gives
a polynomial of degree at most 24. These roots exhaust that degree,
so each root is simple and there are no additional roots.

The exact rational phase midpoints are sorted by orientation. Each root
phase has chord error at most its rational t-bracket width. Pair products
have error bounded by the sum of the two errors. Subtracting these errors
from neighboring max-component chord distances proves

    delta_2 > 0.01509618,
    delta_4 > 5.5894e-7.                                    (9)

All 300 unordered pairs with replacement are tested. In particular,
all sea centers are distinct modulo pi. The endpoint certificates allow
an additional normalized Hermitian Bloch perturbation of norm 2^-35.
Exact rational brackets and margins are retained in
`root_certificate_exact_n12.json`. Numerical roots and centers are in
`full_loop_bands_n12.json`.

[Scope] These are finite low-order nonresonance inequalities. They do
not establish rational independence against arbitrary integer relations.

## 5. Transport to finite K, including the zero-port splitting

[Proved] In a fixed n=12 cell write H_K=K N_bad+B. A bound B_0=2^16
covers the onsite terms and every walk and conversion, including sectors
with several excited tokens. There are 13 sites; each onsite norm is
less than 1007 and the sum of move norms per bond is less than 64.
The field in (8) is smaller than one. These loose bounds suffice.

On P=ker N_bad the compressed operator is the parent. On Q=1-P,
N_bad>=1. For |z|<=2B_0 and K>=8B_0,

    QH_KQ-z >= K/2,
    ||PBQ(QH_KQ-z)^(-1)QBP|| <= 2B_0^2/K.                  (10)

The zero-excitation parent band is a scalar E_v+mu_sel on a 2^12
space, isolated by a gap greater than 1/16. The Schur complement in
(10), followed by the spectral projection across this gap, bounds its
energy spread and its canonical subspace rotation. For two cells a
conservative set of bounds is

    max_a |E_a-2(E_v+mu_sel)| <= 2^36/K,
    ||G_K-G_infinity|| <= 2^44/K.                          (11)

To see the constants, replace B_0 by 2B_0 in (10), giving 2^35/K;
allow a factor two for the energy location. A spectral-gap denominator
of 1/16 and the high-space component 4B_0/K fit within 2^44/K. Use the
polar overlap isometry, within charge sectors, to identify the codes.

[Proved] The legal (1,1) source block has gap at least three. The pair
Schur correction and the energy shift in (11) sum to less than 2^39/K.
The resolvent identity, including the high-high block of norm at most
2/K and its mixed blocks, gives a bound 2^42/K on each source resolvent
difference. Keeping the full source includes every extra marker process.
Together with (11) and ||V||=sqrt(2), each column entering (7) changes
by at most 2^48/K. There are 4^12=2^24 columns, so passage from column
bounds to an operator norm costs at most 2^12. Symmetrization and the
source contraction remain below 2^64/K. In particular,

    ||H_2(K,g)-H_2(infinity,g)|| <= 2^64/K.                 (12)

Combining the one-cell zero-port correction, its two copies, the field
compression correction and lambda^2 times (12), the physical effective
bond through second order differs from its constrained counterpart by
less than 2^80/K for lambda<=1. This last bound includes the zero-port
splitting. It is not multiplied by lambda^2.

[Proved] A quadratic coefficient in normalized Hilbert-Schmidt projection
is bounded by four times the operator norm. Summing the two onsite and
two hopping blocks gives a conservative Bloch bound 16n times that norm.
After the congruence W, let tau=tr(W W^dagger). Then sufficient normalized
error bounds are

    error_g = 4n beta^-2 g^2 tau,
    error_K = 16n tau (2^80/K)/(lambda^2 E_unit).           (13)

The second denominator is the physical band scale. Dividing only by
E_unit would omit the principal new effect of the twelve-letter encoding.

For completeness, an isolated two-cell analytic effective operator has
remainder at most 2^50 lambda^4. One way to obtain this loose bound is
to take a contour of radius 1/64 around the parent vacuum band and an
analytic port disk of radius 2^-12. The resolvent Neumann series converges
uniformly there; the spectral projector and its canonical intertwiner
are bounded by two. The effective operator relative to its central
energy has norm below one on that disk. Evenness and Cauchy's estimate
bound its tail by 2^49 lambda^4 for lambda<=2^-13; 2^50 leaves slack
for finite K. Its normalized quadratic error on two cells is bounded by

    error_lambda,local=16n tau 2^50 lambda^2/E_unit.         (14)

Equation (14) is not a volume-uniform many-cell remainder estimate.

[Proved, computer-assisted] Since sin(pi/14)>=2/14,
alpha^4>=64/14^6>2^-17. `finite_budget.py` uses this rational lower
bound, the exact recurrence for beta, and the integer congruence to obtain

| Error in the dyadic metric | Upper bound |
| :-- | :-- |
| Finite g, equation (13) | 4.93e-18 |
| Finite K, including zero-port splitting | 1.16e-55 |
| Local two-cell lambda remainder | 3.68e-26 |

Their sum is less than one eighth of 2^-35. All acceptance comparisons
use Fractions; logarithms only summarize the output. Thus (9), the
24-root count, and distinct centers persist for the quadratic projection
of the finite-K effective Hamiltonian through second order. The quartic
and higher operators in that Hamiltonian are small, not zero.

## 6. Uniform zero-port length selection

[Definition] For any integer b>=1 put

    a_b=2[1-cos(pi/(b+1))],
    delta_b=(b+1)a_b-b a_(b+1),
    mu_sel=2-(delta_13+delta_12)/2.                         (15)

[Proved] At lambda=0 every excited hop, switch and conversion preserves
each bad-pair projector. A vacuum move only changes the projector just
to its left, when its predecessor is R,D or X. Write these terms as
-(a_j+a_j^dagger), with a_j creating that bad pair. Each a_j is a partial
isometry, including the bit it shifts, and

    a_j a_j^dagger <= Pi_bad,(j-1),
    sum_j a_j^dagger a_j <= max(k-1,0) I

in a cell with k data sites. There is no such move at the first site,
because its predecessor is a marker. The remaining H_good preserves
every bad pair. Cutting m bad pairs gives at most m+1 legal fragments.
Each fragment has at most one mobile head, with kinetic row sum at
most two. A possible conversion costs at most g per excited head and
is absorbed into Delta, since Delta>2+g. Fragments with no head have
nonnegative kinetic energy. Therefore

    H_good >= mu_sel+E_v(k)-4 N_bad.

Young's operator inequality gives

    -(a_j+a_j^dagger) >= -eta a_j a_j^dagger
                         -eta^-1 a_j^dagger a_j.

With eta=K-4 the penalty absorbs the first terms. A constrained vacuum
trial state supplies the upper bound, yielding

    mu_sel+E_v(k)-max(k-1,0)/(K-4)
       <= E_K(k) <= mu_sel+E_v(k),
    E_v(k)=-2cos(pi/(k+2)).                                (16)

This includes all data states and arbitrary illegal or excited sectors.
It does not assert degeneracy of the lowest payload levels.

[Proved] At infinite K the cell energy density is
 e_b=(mu_sel-2+a_b)/b. Direct subtraction gives

    e_(b+1)-e_b=[2-mu_sel-delta_b]/[b(b+1)].

The function a(x)=2[1-cos(pi/(x+1))] is strictly convex for x>=1,
because a''(x)=4pi sin(t)/(x+1)^3+2pi^2 cos(t)/(x+1)^4>0,
where t=pi/(x+1). Hence
 delta_b-delta_(b+1)=(b+1)(a_b-2a_(b+1)+a_(b+2))>0.
Equation (15) makes e_b strictly decrease up to b=13 and strictly
increase thereafter. The nearest competing density is a neighbor.

[Proved, computer-assisted] Rational alternating-series enclosures of
pi and cos give

    min(e_12-e_13,e_14-e_13) > 0.0000547433,
    e_13 < -0.00720578,
    vacuum gap > 0.14791808,
    mu_sel = 1.8561806413628503... .

Dividing (16) by b changes any density by at most 1/(K-4), far below
the selector gap in (8). The same chemical potential selects b=13
uniquely at finite K. On marker-free rings the fragment argument gives
energy density at least -1/(K-4), above the selected-cell trial density.
Thus the zero-port cell length is selected on compatible rings.

[Numerical] The form inequality and energy bound (16) were checked in
the full invariant two-cell spaces at n=1,2, dimensions 346 and 8224.
The lowest H_good form eigenvalues are within 4e-15 of zero.
`selector_certificate.json` contains the rational selector results and
these diagnostic checks. Exact word enumeration also checks completeness,
injectivity, total bit number and charge through n=6 in `grammar_checks.json`.

## 7. Small exact diagonalization

[Numerical] `checks.py` starts from every legal two-cell word and closes
under every literal local rule, including off-grammar processes. It uses
the entire invariant closure, not a prescribed number of violations.
The n=1,2 dimensions are 346 and 8224. The zero-excitation sectors have
16 and 196 states; retain their lowest 4 and 16 states. Equation (7)
is evaluated by sparse solves in the full source sector. Shift-invert
ED of H_0+lambda V uses a Krylov space larger than the code band.

At D=4, Delta=0.5, g=0.04, t_d=1e-4 and K=100, the errors are

| n | lambda=0.2 | lambda=0.1 | lambda=0.05 |
| :-- | :-- | :-- | :-- |
| 1, max energy error / lambda^4 | 0.00602081 | 0.00601509 | 0.00601365 |
| 2, max energy error / lambda^4 | 0.00174381 | 0.00173938 | 0.00173829 |

All code levels are retained. Eigenpair residuals are below 1.6e-13.
The constrained g^2 operator is Gaussian on the full pair Fock space
to relative residual below 2e-15.

[Numerical] At n=2 the zero-port pair-band splitting is 0.0049516990
for K=100 and 0.0004995142 for K=1000. The corresponding differences
between finite-K H_2, aligned by the polar overlap isometry, and the
constrained H_2 are 0.0037052532 and 0.0003892667. These resolve the
O(1/K) corrections and confirm that the splitting must be included.
At n=1 the zero-port code happens to remain degenerate.

[Scope] These are resolved diagnostic couplings, not floating-point ED
of the extreme hierarchy (8). The latter is handled by the finite bounds.

## 8. Alternatives and the limit of the alphabet claim

[Proved, restricted] Putting bits only on a fixed number r of tokens
cannot retain this code. For a unique clock ground state per bit string,
there are at most 2^r payload states, whereas n channels require 2^n.
Leaving the bits on the background and moving the displaced bit to the
marker avoids that obstruction. The marker is used both as the storage
site and as the right port in (1).

[Numerical] Directly merging C75's equal-charge D_(A,1),D_(B,0), or
X_(A,1),X_(B,0), was tested separately. At g=0 the merged legal-image
pair dimensions are 225 and 2304 for n=1,2. Their second-order operators
already have non-Gaussian terms. At n=2 the relative residuals are
4.63e-4 for the data merge and 4.88e-5 for the express merge. Their
non-scalar norms are 4.59e-6 and 9.79e-6 at the diagnostic parameters.
These terms survive as g tends to zero, so they cannot be bounded as
the O(g^4) remainder of the C75 expansion. `merge_c75.py` preserves the
attempt. A different nonperturbative analysis of these models is open.

[Proved, restricted; numerical values] Deleting the fold and using two
identical real endpoint carriers gives, already at n=2,
 F=[[a,b],[c,a]], with bc>a^2. Its Bloch eigenvalues are

    2a cos k +/- sqrt((b+c)^2 cos^2 k+(b-c)^2 sin^2 k).

One is positive and the other negative at every k, with disjoint ranges.
For the explicit carrier [[4,-1],[-1,4]], bc-a^2=7.08309e-9, and the
ranges are [-0.00238397,-0.00208333] and
[0.00208333,0.00238397]. A common field cannot cross both.

[Numerical] A different fold attempt replaces only the B express rail
with a fast direct path. Its scalar Hamiltonian is 4I+r times the path
Laplacian, while A retains its express fold. Speeds r=1,10,100,1000,10000
were tested at n=2,3. The onsite cancellation is lost. The onsite to
hopping norm ratio ranges from 1776 to 427033, and every common-field
band-interval overlap is negative. The smallest magnitude of the negative
overlap is 2.82e-6; the grid error bound 2||T|| pi/2048 is below 3.4e-13
there. This attempt fails at those parameters, not for every choice of
asymmetric carrier or counterterm. See `fast_carrier_checks.json`.

[Proved, restricted] A two-token replacement encoding A/B in the spatial
order of two noncrossing heads has conserved order under nearest-neighbor
hopping. More basically, AB and BA have the same total charge, whereas
the conversion coupled to a logical annihilator must change carrier
charge by two. An order-only conversion cannot conserve this U(1).
Allowing AA and BB sectors or a charged boundary state changes the
encoding and requires another analysis. A six-site, two-head matrix of
dimension 30 has zero cross-order block and zero order commutator, as
recorded in `alternative_checks.json`. This rejects that concrete
order-based fold, not all two-token cells.

[Proved, restricted; numerical checks] All twelve charge-preserving
single identifications in the twelve-letter alphabet were audited.
Merging D_s with X_s collapses 52 legal n=2 states to 44 and identifies
different rail amplitudes. Merging L_0 with R_0 collapses them to 48.
Merging L_1 with R_1 is injective on these words, but admits additional
zero-penalty paths. Its n=2 vacuum energies begin -1.80194,-1.41669,
-1.41421,-1.24698 instead of a fourfold level near -sqrt(2). The L_0/R_0
merge gives -2,-1.61803,-1.61803,-1.41421. These are order-one failures
of the scalar vacuum path, not small finite-penalty corrections.

Every remaining equal-charge identification merges a marker with a data
letter. It removes the conserved on-site marker projector. For example,
M_A=D_A makes the existing D_A R_a -> L_a D_A hop move a boundary at
unit strength. The original fixed-cell sectors then cease to be invariant.
A penalty on the old grammar does not recover that invariant. These
attempts do not supply an eleven-letter bond with the inherited proof.

[Heuristic] Twelve is the smallest plausible alphabet supported here
by a complete finite construction. It uses four background letters,
four charged rail letters and four marker letters, with the vacuum
letters already absorbed. Smaller literal quotients fail one of the
properties just checked. There is no universal lower bound a>=12.
An eleven-letter construction with a different selector, a different
charge representation, or moving markers remains possible in principle.

## 9. Infrared scope

[Assumed] Transferring the finite quadratic certificate to c=12 requires
volume-uniform elimination of the massive sectors and the usual weakly
interacting multichannel Luttinger-liquid stability input. Exact U(1),
(9), and small residual interactions provide finite ingredients. They
are not a proof of the many-body infrared theorem.

[Open] Equation (16) proves selection only at zero port. At nonzero port,
a common-bulk-subtracted dressed-defect comparison is still needed.
The local two-cell remainder bound does not supply it. Frozen cells
also retain translation-related sectors, as in C75. A vacancy lift
adds the existing vacancy matching and mobile-fluid assumptions.

[Conditional] Subject to these inputs, the frozen twelve-letter member
has c=12>11. One additional vacancy gives a thirteen-letter candidate
with c=13>12. Neither is claimed here as a completed microscopic
counterexample to the full conjecture.
